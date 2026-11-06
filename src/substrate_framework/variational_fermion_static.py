"""Finite collinear FD functional: independent scalar minimization and curvature.

For beta>0 and U>=0, spectral entropy is strictly convex on 0<gamma<I
(its Frobenius Hessian is >=4/beta); the other gamma dependence is affine
or convex quadratic. A feasible stationary point is therefore the unique
minimum, not a selected SCF saddle. Constraints fix mean N, optionally mean M;
this is not the finite canonical ensemble. No material applicability is implied.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import root
from scipy.special import expit, xlogy


def _expect(operators, gamma):
    return np.einsum("...ij,ji->...", operators, gamma).real


def _constraints(fixture, fixed_m):
    operators = [np.eye(len(fixture["h0"]))]
    targets = [fixture["N0"]]
    if fixed_m:
        operators.append(fixture["Z"])
        targets.append(fixture["M0"])
    return np.asarray(operators), np.asarray(targets)


def _ionic_hamiltonian(fixture, u):
    return (fixture["h0"] - np.einsum("a,aij->ij", u, fixture["B"])
            + .5 * np.einsum("a,b,abij->ij", u, u, fixture["W"]))


def full_functional(fixture, gamma, u):
    """Evaluate the full defining scalar, including entropy, U, W and Kion."""
    occupations = eigh(gamma, eigvals_only=True)
    if np.any(occupations <= 0) or np.any(occupations >= 1):
        raise ValueError("gamma must be strictly between zero and identity")
    entropy = np.sum(xlogy(occupations, occupations)
                     + xlogy(1 - occupations, 1 - occupations)) / fixture["beta"]
    delta = _expect(fixture["D"], gamma) - fixture["d0"]
    return float(_expect(_ionic_hamiltonian(fixture, u), gamma) + entropy
                 + .5 * delta @ fixture["U"] @ delta
                 + .5 * u @ fixture["Kion"] @ u)


def solve_minimum(fixture, u, fixed_m=False):
    """Solve FD selfconsistency plus affine constraints, without response algebra.

    Unknowns are two induced densities and number/spin multipliers. The returned
    energy is evaluated from full_functional, never from a quadratic prediction.
    Residual/domain diagnostics expose the convex-minimum justification.
    """
    u = np.asarray(u, dtype=float)
    D, U, beta = fixture["D"], fixture["U"], fixture["beta"]
    Q, targets = _constraints(fixture, fixed_m)
    external = _ionic_hamiltonian(fixture, u)
    nd = len(D)

    def state(x):
        effective = (external + np.einsum("s,sij->ij", U @ (x[:nd] - fixture["d0"]), D)
                     - np.einsum("c,cij->ij", x[nd:], Q))
        levels, vectors = eigh(effective)
        f = expit(-beta * levels)
        return (vectors * f) @ vectors.conj().T

    def residual(x):
        gamma = state(x)
        return np.r_[x[:nd] - _expect(D, gamma), _expect(Q, gamma) - targets]

    initial = np.r_[fixture["d0"], fixture["mu0"], np.zeros(len(Q) - 1)]
    solution = root(residual, initial, method="hybr", options={"xtol": 1e-11})
    gamma = state(solution.x)
    values, vectors = eigh(gamma)
    if np.any(values <= 0) or np.any(values >= 1):
        raise RuntimeError("SCF solution left the entropy domain")
    density = _expect(D, gamma)
    logit = (vectors * (np.log(values) - np.log1p(-values))) @ vectors.conj().T
    gradient = (external + np.einsum("s,sij->ij", U @ (density - fixture["d0"]), D)
                + logit / beta - np.einsum("c,cij->ij", solution.x[nd:], Q))
    residual_norm = float(np.max(np.abs(residual(solution.x))))
    stationarity = float(np.linalg.norm(gradient))
    interaction_min = float(eigh(U, eigvals_only=True)[0])
    accepted = bool(residual_norm < 1e-10 and stationarity < 1e-9
                    and interaction_min >= 0 and np.all(np.isfinite(gamma)))
    return {
        "u": u, "energy": full_functional(fixture, gamma, u), "gamma": gamma,
        "density": density, "multipliers": solution.x[nd:],
        "accepted_minimum": accepted, "solver_success": bool(solution.success),
        "solver_status": int(solution.status), "solver_message": str(solution.message),
        "solver_nfev": int(solution.nfev), "solver_residual_max": residual_norm,
        "stationarity_frobenius": stationarity,
        "constraint_residual": _expect(Q, gamma) - targets,
        "gamma_eigenvalues": values, "domain_margin": float(min(values[0], 1-values[-1])),
        "interaction_min_eigenvalue": interaction_min,
        "strict_convexity_lower_bound": 4 / beta,
    }


def joint_gibbs_curvature(fixture, fixed_m=False, grandcanonical=False):
    """Independent spectral Gibbs metric on ordered (B,D), then Q Schur projection.

    Divided differences include occupation covariance on all degenerate blocks.
    Fixed-M projects BOTH I and Z. This routine is never used by solve_minimum.
    """
    Q, _ = _constraints(fixture, fixed_m)
    physical = np.concatenate([fixture["B"], fixture["D"]])
    operators = physical if grandcanonical else np.concatenate([physical, Q])
    energies, vectors = eigh(fixture["h0"])
    f = expit(-fixture["beta"] * (energies - fixture["mu0"]))
    gap = energies[None, :] - energies[:, None]
    weight = np.divide(f[:, None] - f[None, :], gap,
                       out=np.zeros_like(gap), where=gap != 0)
    weight = np.where(gap == 0, fixture["beta"] * f[:, None] * (1-f[:, None]), weight)
    rotated = np.asarray([vectors.conj().T @ op @ vectors for op in operators])
    metric = np.einsum("aij,bji,ij->ab", rotated, rotated, weight).real
    n = len(physical)
    if not grandcanonical:
        metric = (metric[:n, :n] - metric[:n, n:]
                  @ np.linalg.solve(metric[n:, n:], metric[n:, :n]))
    return (metric + metric.T) / 2


def analytic_hessian(fixture, fixed_m=False, mutation=None):
    """Kion+<W>-C, with constrained induced-density U loop; explicit mutations.

    delta_d=(I+C_DD U)^-1 C_DB u; C=C_BB-C_BD U delta_d/du.
    wrongGC removes all mean constraints; missingW drops only <W>;
    droppedUcross zeroes only the interaction's offdiagonal spin channels.
    """
    if mutation not in (None, "wrongGC", "missingW", "droppedUcross"):
        raise ValueError("unknown curvature mutation")
    C = joint_gibbs_curvature(fixture, fixed_m, mutation == "wrongGC")
    n = len(fixture["B"])
    U = fixture["U"]
    if mutation == "droppedUcross":
        U = np.diag(np.diag(U))
    screened = C[:n, :n] - C[:n, n:] @ U @ np.linalg.solve(
        np.eye(len(U)) + C[n:, n:] @ U, C[n:, :n])
    explicit = _expect(fixture["W"], fixture["gamma0"])
    if mutation == "missingW":
        explicit = np.zeros_like(explicit)
    return fixture["Kion"] + explicit - screened

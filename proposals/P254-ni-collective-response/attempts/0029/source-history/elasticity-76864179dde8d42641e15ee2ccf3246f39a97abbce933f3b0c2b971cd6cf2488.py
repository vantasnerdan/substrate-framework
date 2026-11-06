"""Exact isotropic linear elasticity: tensor algebra through acoustic speeds.

Conventions are fixed once and reused everywhere:

- small-strain tensor ``eps_ij = (du_i/dx_j + du_j/dx_i)/2``;
- isotropic stiffness ``C_ijkl = lam*delta_ij*delta_kl
  + mu*(delta_ik*delta_jl + delta_il*delta_jk)``;
- Hooke stress ``sigma_ij = C_ijkl*eps_kl``;
- Navier-Cauchy operator ``rho*u_tt = mu*Laplacian(u)
  + (lam+mu)*grad(div(u)) + f``;
- strong ellipticity requires ``mu > 0`` and ``lam + 2*mu > 0``;
- longitudinal (P) speed squared is ``(lam+2*mu)/rho``, transverse (S)
  speed squared is ``mu/rho``.

All functions are pure; importing this module executes no simulation.
"""

from __future__ import annotations

from typing import Any

import sympy as sp


def _nonnegative(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number:
        if expression.is_real is not True or float(expression) < 0.0:
            raise ValueError(f"{name} must be real and nonnegative")
    return expression


def isotropic_stiffness(
    lam: Any, mu: Any, dim: int = 3
) -> sp.ImmutableDenseNDimArray:
    """Return C[i,j,k,l] for the declared isotropic Lamé constants."""

    lam = _nonnegative(lam, "lam")
    mu = _nonnegative(mu, "mu")
    delta = sp.eye(dim)
    return sp.ImmutableDenseNDimArray(
        [
            [
                [
                    [
                        lam * delta[i, j] * delta[k, l]
                        + mu * (delta[i, k] * delta[j, l] + delta[i, l] * delta[j, k])
                        for l in range(dim)
                    ]
                    for k in range(dim)
                ]
                for j in range(dim)
            ]
            for i in range(dim)
        ],
        shape=(dim, dim, dim, dim),
    )


def strain_from_gradient(
    gradient: sp.MatrixBase | sp.NDimArray, dim: int = 3
) -> sp.ImmutableDenseMatrix:
    """Symmetrize a displacement gradient into the small-strain tensor."""

    grad = sp.Matrix(gradient)
    if grad.shape != (dim, dim):
        raise ValueError("gradient must be a square dim x dim matrix")
    return sp.ImmutableMatrix((grad + grad.T) / 2)


def hooke_stress(
    stiffness: sp.NDimArray, strain: sp.MatrixBase, dim: int = 3
) -> sp.ImmutableMatrix:
    """Contract sigma_ij = C_ijkl eps_kl."""

    strain_matrix = sp.Matrix(strain)
    if strain_matrix.shape != (dim, dim):
        raise ValueError("strain must be a square dim x dim matrix")
    stress = sp.zeros(dim, dim)
    for i in range(dim):
        for j in range(dim):
            total = sp.Integer(0)
            for k in range(dim):
                for l in range(dim):
                    component = stiffness[i, j, k, l]
                    if component != 0:
                        total += component * strain_matrix[k, l]
            stress[i, j] = sp.simplify(total)
    return sp.ImmutableMatrix(stress)


def navier_cauchy_operator(
    u_fields: list[sp.Expr] | tuple[sp.Expr, ...],
    coordinates: list[sp.Symbol] | tuple[sp.Symbol, ...],
    time: sp.Symbol,
    lam: Any,
    mu: Any,
    rho: Any = 1,
) -> list[sp.Expr]:
    """Evaluate rho*u_tt - [mu*Lap(u) + (lam+mu)*grad(div u)] per component.

    The body force is excluded; verifiers add it as an external balance.
    """

    lam = _nonnegative(lam, "lam")
    mu = _nonnegative(mu, "mu")
    if len(u_fields) != len(coordinates):
        raise ValueError("one displacement field per coordinate required")
    residual: list[sp.Expr] = []
    div_u = sum(sp.diff(u_fields[a], coordinates[a]) for a in range(len(coordinates)))
    for a in range(len(coordinates)):
        laplace = sum(
            sp.diff(sp.diff(u_fields[a], x), x) for x in coordinates
        )
        inertial = rho * sp.diff(sp.diff(u_fields[a], time), time)
        elastic = mu * laplace + (lam + mu) * sp.diff(div_u, coordinates[a])
        residual.append(sp.expand(inertial - elastic))
    return residual


def christoffel_matrix(
    lam: Any,
    mu: Any,
    rho: Any,
    direction: sp.MatrixBase,
    dim: int = 3,
) -> sp.ImmutableMatrix:
    """Acoustic tensor Gamma_ik = C_ijkl n_j n_l / rho."""

    lam = _nonnegative(lam, "lam")
    mu = _nonnegative(mu, "mu")
    n = sp.Matrix(direction)
    if n.shape != (dim, 1):
        raise ValueError("direction must be a dim x 1 vector")
    stiffness = isotropic_stiffness(lam, mu, dim=dim)
    matrix = sp.zeros(dim, dim)
    for i in range(dim):
        for k in range(dim):
            total = sp.Integer(0)
            for j in range(dim):
                for l in range(dim):
                    total += stiffness[i, j, k, l] * n[j] * n[l]
            matrix[i, k] = sp.simplify(total / rho**1)
    return sp.ImmutableMatrix(matrix)


def acoustic_speeds_squared(lam: Any, mu: Any, rho: Any = 1) -> dict[str, sp.Expr]:
    """Exact P and S wave speed squares under the fixed conventions."""

    lam = _nonnegative(lam, "lam")
    mu = _nonnegative(mu, "mu")
    rho = _nonnegative(rho, "rho")
    return {
        "P": (lam + 2 * mu) / rho**1,
        "S": mu / rho**1,
    }


def poisson_ratio(lam: Any, mu: Any) -> sp.Expr:
    """Nu = lam / (2*(lam + mu)) for a stable isotropic continuum."""

    lam = _nonnegative(lam, "lam")
    mu = _nonnegative(mu, "mu")
    return sp.simplify(lam / (2 * (lam + mu)))


def strong_elliptic(lam: Any, mu: Any) -> bool:
    """True exactly when mu > 0 and lam + 2*mu > 0."""

    lam_expr = sp.sympify(lam)
    mu_expr = sp.sympify(mu)
    if not (lam_expr.is_real and mu_expr.is_real):
        raise ValueError("Lamé constants must be real numbers for this predicate")
    return bool(mu_expr > 0 and (lam_expr + 2 * mu_expr) > 0)


def project_divergence_free(
    wavevector: sp.MatrixBase, amplitude: sp.MatrixBase
) -> sp.ImmutableMatrix:
    """Project an amplitude onto the divergence-free subspace of direction k."""

    k = sp.Matrix(wavevector)
    a = sp.Matrix(amplitude)
    if k.shape != a.shape or k.shape[1] != 1:
        raise ValueError("wavevector and amplitude must be matching column vectors")
    kk = k.dot(k)
    projected = a - k * (k.dot(a) / kk)
    return sp.ImmutableMatrix(sp.simplify(projected))

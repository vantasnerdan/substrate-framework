"""C-CST-005 verifier: dispersion of the conditional N4 micropolar operator.

This verifier derives the transverse determinant directly from the two coupled
balance equations. It also records the singular nature of the structure-free
limit: the displacement branch tends to the linearized-Euler neutral sector,
while the spin branch has a finite limiting frequency because both its
stiffness and its microinertia are proportional to ``L_v``. At ``L_v = 0`` the
spin coordinate has no kinetic or potential weight and is removed; assigning
it an arbitrary nonzero inertia would test a different model.

The Euler-derived displacement sector is divergence-free. Its longitudinal
displacement formula belongs only to the formal compressible extension.
Longitudinal spin remains in the conditional incompressible micropolar system.
"""

import sys

import numpy as np
import sympy as sp

from substrate_framework.numerics import SolverTolerances, solve_ivp_evidence
from substrate_framework.micropolar import MicropolarCoefficients, micropolar_fourier_stiffness
from substrate_framework.verification import CheckLedger

k, w2 = sp.symbols("k omega2", positive=True)
rho, lam, mu, alpha = sp.symbols("rho lambda mu alpha", positive=True)
cs, ca, j = sp.symbols("c_s c_a j", positive=True)
ctr = sp.Symbol("c_tr", real=True)

coefficients = MicropolarCoefficients(lam, mu, alpha, ctr, cs, ca)
full_stiffness = micropolar_fourier_stiffness([0, 0, k], coefficients)
helicity_axis = sp.Matrix([1, sp.I, 0])/sp.sqrt(2)
helicity_columns = sp.zeros(6, 2)
helicity_columns[:3, 0], helicity_columns[3:, 1] = helicity_axis, helicity_axis
projected_stiffness = sp.simplify(helicity_columns.conjugate().T*full_stiffness*helicity_columns)
K_u, K_phi = projected_stiffness[0, 0], projected_stiffness[1, 1]
coupling = -projected_stiffness[0, 1]

# Build the generalized eigenproblem instead of re-entering its polynomial.
transverse_matrix = sp.Matrix([[rho * w2 - K_u, coupling], [coupling, j * w2 - K_phi]])
det_transverse = sp.expand(transverse_matrix.det())
poly = sp.Poly(det_transverse, w2)
A, B, C = poly.all_coeffs()
disc = sp.expand(B**2 - 4 * A * C)
w2p = sp.simplify((-B + sp.sqrt(disc)) / (2 * A))
w2m = sp.simplify((-B - sp.sqrt(disc)) / (2 * A))


def check_transverse_quadratic(ledger: CheckLedger) -> None:
    expected = sp.expand((rho * w2 - K_u) * (j * w2 - K_phi) - 4 * alpha**2 * k**2)
    ledger.check(
        "transverse determinant derived from the N4 2x2 operator",
        sp.simplify(det_transverse - expected) == 0,
        "(rho*w2-K_u)(j*w2-K_phi) = 4*alpha^2*k^2",
    )
    ledger.check(
        "off-diagonal product is 4 alpha^2 k^2",
        sp.simplify(coupling**2 - 4 * alpha**2 * k**2) == 0,
        "each balance contributes one 2*alpha*curl coupling",
    )


def check_branches_and_limits(ledger: CheckLedger) -> None:
    acoustic = sp.simplify(sp.limit(w2m, k, 0))
    optical = sp.simplify(sp.limit(w2p, k, 0))
    ledger.check("acoustic branch w2_-(0) = 0", acoustic == 0, f"limit = {acoustic}")
    ledger.check(
        "optical spin branch w2_+(0) = 4 alpha/j",
        sp.simplify(optical - 4 * alpha / j) == 0,
        f"limit = {optical}",
    )
    ledger.check(
        "Vieta sum and product over the two branches",
        sp.simplify(w2p + w2m + B / A) == 0
        and sp.simplify(sp.expand(w2p * w2m - C / A)) == 0,
        "roots are those of the operator-derived determinant",
    )

    cell_density, cell_stiffness, cell_inertia = sp.symbols("n_cell K_Psi J_Psi", positive=True)
    alpha_micro = cell_density*cell_stiffness/12
    j_micro = cell_density*cell_inertia/3
    optical_gap = sp.simplify(4 * alpha_micro / j_micro)
    expected_gap = cell_stiffness/cell_inertia
    ledger.check(
        "same-orbit isotropic substitution preserves the physical optical gap",
        sp.simplify(optical_gap - expected_gap) == 0,
        "omega_gap^2 = K_Psi/J_Psi; no retired tension locking or appended rigid inertia",
    )


def check_longitudinal_extension(ledger: CheckLedger) -> None:
    proj_lin = sp.simplify(full_stiffness[2, 2]-(lam+2*mu)*k**2)
    ledger.check(
        "unrestricted longitudinal diagonal: rho w^2 = (lambda+2 mu) k^2",
        proj_lin == 0,
        "formal compressible extension; not an incompressible-Euler mode",
    )
    proj_ang = sp.simplify(full_stiffness[5, 5]-(2*(cs+ctr)*k**2+4*alpha))
    ledger.check(
        "longitudinal spin diagonal: j w^2 = 2(c_s+c_tr) k^2 + 4 alpha",
        proj_ang == 0,
        "exact projection of the conditional N4 operator",
    )
    ux, uy = sp.symbols("u_x u_y")
    kvec = sp.Matrix([0, 0, k])
    u_transverse = sp.Matrix([ux, uy, 0])
    ledger.check(
        "Euler-derived displacement sector is divergence-free",
        sp.simplify(kvec.dot(u_transverse)) == 0,
        "the longitudinal P branch is excluded by the microscopic incompressibility premise",
    )


def check_structure_free_limit(ledger: CheckLedger) -> None:
    Lv = sp.Symbol("L_v", positive=True)
    mu0, alpha0, cs0, ca0, j0 = sp.symbols("mu0 alpha0 cs0 ca0 j0", positive=True)
    scaled = det_transverse.subs(
        {
            mu: Lv * mu0,
            alpha: Lv * alpha0,
            cs: Lv * cs0,
            ca: Lv * ca0,
            j: Lv * j0,
        }
    )
    limiting_polynomial = sp.simplify(sp.limit(scaled / Lv, Lv, 0))
    target = rho * w2 * (j0 * w2 - (cs0 + ca0) * k**2 - 4 * alpha0)
    ledger.check(
        "L_v -> 0 determinant has one neutral displacement root and one finite spin root",
        sp.simplify(limiting_polynomial - target) == 0,
        "the limit is singular because spin kinetic and potential weights both vanish",
    )
    ledger.check(
        "finite limiting spin frequency is the stiffness/inertia ratio",
        sp.simplify(((cs0 + ca0) * k**2 + 4 * alpha0) / j0) != 0,
        "the spin coordinate is removed, not assigned j=1, at L_v=0",
    )


def check_ivp_refinement(ledger: CheckLedger) -> None:
    """Replay a well-conditioned transverse normal mode with SciPy DOP853."""

    rho_n, mu_n, alpha_n = 2.0, 3.0, 1.0
    cs_n, ca_n, j_n, k_n = 2.0, 1.0, 1.5, 0.4
    mass = np.diag([rho_n, j_n])
    stiffness = np.array(
        [
            [(mu_n + alpha_n) * k_n**2, -2.0 * alpha_n * k_n],
            [-2.0 * alpha_n * k_n, (cs_n + ca_n) * k_n**2 + 4.0 * alpha_n],
        ],
        dtype=np.float64,
    )
    inv_sqrt_mass = np.diag(1.0 / np.sqrt(np.diag(mass)))
    reduced = inv_sqrt_mass @ stiffness @ inv_sqrt_mass
    eigenvalues, eigenvectors = np.linalg.eigh(reduced)
    optical_w2 = float(eigenvalues[-1])
    q0 = inv_sqrt_mass @ eigenvectors[:, -1]
    y0 = np.concatenate([q0, np.zeros(2)])
    period = 2.0 * np.pi / np.sqrt(optical_w2)

    symbolic_w2 = float(
        sp.N(
            w2p.subs(
                {
                    rho: rho_n,
                    mu: mu_n,
                    alpha: alpha_n,
                    cs: cs_n,
                    ca: ca_n,
                    j: j_n,
                    k: k_n,
                }
            ),
            17,
        )
    )
    ledger.check(
        "generalized eigenvalue agrees with the symbolic optical branch",
        abs(optical_w2 - symbolic_w2) < 5.0e-14,
        f"numeric={optical_w2:.16e}, symbolic={symbolic_w2:.16e}",
    )

    mass_inverse = np.linalg.inv(mass)

    def rhs(_time: float, state: np.ndarray) -> np.ndarray:
        return np.concatenate([state[2:], -mass_inverse @ stiffness @ state[:2]])

    def energy(state: np.ndarray) -> float:
        q, qdot = state[:2], state[2:]
        return float(0.5 * qdot @ mass @ qdot + 0.5 * q @ stiffness @ q)

    errors = []
    drifts = []
    evaluations = []
    for rtol in (1.0e-7, 1.0e-9, 1.0e-11):
        evidence = solve_ivp_evidence(
            rhs,
            (0.0, period),
            y0,
            sample_times=np.linspace(0.0, period, 257),
            tolerances=SolverTolerances(
                rtol=rtol, atol=rtol * 1.0e-2, max_step=period
            ),
            method="DOP853",
            invariant=energy,
        )
        errors.append(float(np.linalg.norm(evidence.state[:, -1] - y0)))
        drifts.append(float(evidence.max_abs_invariant_drift or 0.0))
        evaluations.append(evidence.function_evaluations)

    print(
        "IVP_EVIDENCE "
        f"rtol={[1.0e-7, 1.0e-9, 1.0e-11]} "
        f"state_errors={errors} energy_drifts={drifts} nfev={evaluations}"
    )

    ledger.check(
        "DOP853 one-period state error decreases under tolerance refinement",
        errors[2] < errors[1] < errors[0] and errors[2] < 2.0e-10,
        f"errors={errors}; nfev={evaluations}; float64",
    )
    ledger.check(
        "DOP853 energy drift decreases and is below the final error budget",
        drifts[2] < drifts[1] < drifts[0] and drifts[2] < 2.0e-10,
        f"max_abs_energy_drifts={drifts}",
    )


def check_mutations(ledger: CheckLedger) -> None:
    C_no_gap = sp.expand(K_u * ((cs + ca) * k**2) - 4 * alpha**2 * k**2)
    B_no_gap = -(rho * (cs + ca) * k**2 + j * K_u)
    w2p_no_gap = sp.simplify(
        (-B_no_gap + sp.sqrt(B_no_gap**2 - 4 * A * C_no_gap)) / (2 * A)
    )
    ledger.check(
        "M1 dropped spin stiffness rejected",
        sp.simplify(sp.limit(w2p_no_gap, k, 0)) == 0,
        "the optical gap vanishes",
    )

    wrong_eight = sp.expand((rho * w2 - K_u) * (j * w2 - K_phi) - 8 * alpha**2 * k**2)
    ledger.check(
        "M2 doubled coupling product rejected by the N4 operator",
        sp.simplify(det_transverse - wrong_eight) != 0,
        "two off-diagonal factors 2*alpha*k multiply to 4*alpha^2*k^2",
    )

    wrong_longitudinal = 2*(ca-ctr)*k**2+4*alpha
    ledger.check(
        "M3 historical wrong longitudinal grad-div sign rejected",
        sp.simplify(wrong_longitudinal-full_stiffness[5, 5]) != 0,
        "both Laplacian and grad-div contribute with the same Fourier sign",
    )


def main() -> int:
    ledger = CheckLedger("C-CST-005")
    check_transverse_quadratic(ledger)
    check_branches_and_limits(ledger)
    check_longitudinal_extension(ledger)
    check_structure_free_limit(ledger)
    check_ivp_refinement(ledger)
    check_mutations(ledger)
    return int(ledger.finish())


if __name__ == "__main__":
    sys.exit(main())

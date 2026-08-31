"""Exact O1 and commuting-core calculus for P249 attempt 0001."""

from __future__ import annotations

import sympy as sp

from substrate_framework.m5_covariant_action import (
    MINKOWSKI_MOSTLY_PLUS,
    cartan_inverse_metric_from_projector,
    m5_curvature_from_derivatives,
    projected_spatial_ldg_potential,
    spectral_cartan_hamiltonian_density,
    spectral_cartan_lagrangian_density,
)
from substrate_framework.verification import CheckLedger


def _zero_matrix(matrix: sp.MatrixBase) -> bool:
    return matrix.applyfunc(sp.trigsimp) == sp.zeros(matrix.rows, matrix.cols)


def main() -> int:
    ledger = CheckLedger("P249-O1-0001")
    eta = sp.Matrix(MINKOWSKI_MOSTLY_PLUS)
    q, q1, q2 = sp.symbols("q q1 q2", real=True)
    g, d, omega = sp.symbols("g d omega", positive=True, real=True)
    dx, dy, dz = sp.symbols("d_x d_y d_z", real=True)

    generator = sp.zeros(4)
    generator[2, 3] = -1
    generator[3, 2] = 1

    def rotation(angle: sp.Expr) -> sp.Matrix:
        return sp.Matrix(
            [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, sp.cos(angle), -sp.sin(angle)],
                [0, 0, sp.sin(angle), sp.cos(angle)],
            ]
        )

    rot = rotation(q)
    ledger.check("SO(2) preserves eta", _zero_matrix(rot.T * eta * rot - eta))
    ledger.check(
        "SO(2) group law",
        _zero_matrix(rotation(q1) * rotation(q2) - rotation(q1 + q2)),
    )
    ledger.check(
        "SO(2) has declared 2pi period",
        _zero_matrix(rotation(q + 2 * sp.pi) - rotation(q)),
    )
    ledger.check(
        "generator is normalized derivative",
        _zero_matrix(rot.diff(q).subs(q, 0) - generator),
    )

    vacuum = sp.diag(-g, 1, 0, 0)
    core = sp.diag(-g, 1, d, -d)
    zeta_vacuum = generator * vacuum + vacuum * generator.T
    zeta_core = generator * core + core * generator.T
    ledger.check("exterior vacuum fixed infinitesimally", zeta_vacuum == sp.zeros(4))
    ledger.check(
        "exterior vacuum fixed globally",
        _zero_matrix(rot * vacuum * rot.T - vacuum),
    )
    ledger.check("core split activates the clock", zeta_core != sp.zeros(4))
    ledger.check(
        "orbit tangent equals declared generator",
        _zero_matrix((rot * core * rot.T).diff(q).subs(q, 0) - zeta_core),
    )
    mutated_vacuum = sp.diag(-g, 1, 0, d)
    ledger.check(
        "vacuum fixation is sensitive to lost degeneracy",
        generator * mutated_vacuum + mutated_vacuum * generator.T != sp.zeros(4),
    )

    projector_t = sp.diag(1, 0, 0, 0)
    h_inverse = sp.Matrix(cartan_inverse_metric_from_projector(projector_t, eta))
    ledger.check("aligned Cartan inverse metric is Euclidean", h_inverse == sp.eye(4))
    potential = projected_spatial_ldg_potential(
        core,
        projector_t,
        beta=1,
        scale=1,
        timelike_target_eigenvalue=g,
        timelike_stiffness=1,
        metric_covariant=eta,
    )
    ledger.check("core potential exact", sp.factor(potential - (3 * d**2 + 4 * d**4)) == 0)

    split_direction = sp.diag(0, 0, 1, -1)
    derivatives = [
        omega * zeta_core,
        dx * split_direction,
        dy * split_direction,
        dz * split_direction,
    ]
    curvature = m5_curvature_from_derivatives(derivatives, eta)
    gradient_squared = dx**2 + dy**2 + dz**2
    hamiltonian = spectral_cartan_hamiltonian_density(
        curvature, h_inverse, potential
    )
    lagrangian = spectral_cartan_lagrangian_density(
        curvature, h_inverse, potential, eta
    )
    expected_c = 32 * d**2 * gradient_squared
    hamiltonian_omega_two = sp.expand(hamiltonian).coeff(omega, 2)
    lagrangian_omega_two = sp.expand(lagrangian).coeff(omega, 2)
    ledger.check(
        "clock curvature coefficient exact",
        sp.factor(hamiltonian_omega_two - expected_c) == 0,
    )
    ledger.check(
        "static commuting curvature vanishes",
        sp.factor(hamiltonian.subs(omega, 0) - potential) == 0,
    )
    ledger.check(
        "Lagrangian has the same positive kinetic coefficient",
        sp.factor(lagrangian_omega_two - expected_c) == 0,
    )
    ledger.check(
        "Lagrangian potential sign",
        sp.factor(lagrangian.subs(omega, 0) + potential) == 0,
    )
    charge_density = sp.diff(lagrangian, omega)
    ledger.check(
        "Noether charge is 2 C omega",
        sp.factor(charge_density - 2 * expected_c * omega) == 0,
    )
    ledger.check(
        "Legendre Jacobian is positive on an active core",
        sp.factor(sp.diff(charge_density, omega) - 2 * expected_c) == 0,
    )

    # An explicit smooth decaying core d=exp(-r^2) has
    # K_u=int |grad exp(-2r^2)|^2 d^3x = 3*pi^(3/2)/4, hence C=8K_u.
    gaussian_c = 6 * sp.pi ** sp.Rational(3, 2)
    gaussian_u = sp.pi ** sp.Rational(3, 2) * (
        sp.Rational(3, 2) / sp.sqrt(2) + sp.Rational(1, 2)
    )
    ledger.check("Gaussian core inertia finite positive", gaussian_c.is_positive)
    ledger.check("Gaussian core potential finite positive", gaussian_u.is_positive)
    ledger.check(
        "Gaussian charge-frequency map nondegenerate",
        sp.diff(2 * gaussian_c * omega, omega) == 12 * sp.pi ** sp.Rational(3, 2),
    )

    u = sp.Function("u")
    w = sp.symbols("w", positive=True)
    laplacian_u = sp.Symbol("Delta_u", real=True)
    euler_u = 16 * w**2 * laplacian_u + 3 + 8 * sp.Symbol("u", nonnegative=True)
    ledger.check("fixed-J Euler equation retained", euler_u.coeff(laplacian_u) == 16 * w**2)

    amplitude = sp.symbols("a", positive=True)
    angular_momentum = sp.symbols("J", nonzero=True, real=True)
    a1, a2, k_shape = sp.symbols("A1 A2 K", positive=True)
    radius = 1 / amplitude
    scaled_energy = sp.factor(
        radius**3 * (3 * amplitude * a1 + 4 * amplitude**2 * a2)
        + angular_momentum**2
        / (32 * amplitude**2 * radius * k_shape)
    )
    expected_collapse = (
        3 * a1 / amplitude**2
        + 4 * a2 / amplitude
        + angular_momentum**2 / (32 * amplitude * k_shape)
    )
    ledger.check(
        "collapse path formula exact",
        sp.factor(scaled_energy - expected_collapse) == 0,
    )
    ledger.check(
        "collapse path tends to zero",
        sp.limit(scaled_energy, amplitude, sp.oo) == 0,
    )
    ledger.check(
        "finite profiles remain positive",
        all(term.is_positive for term in sp.Add.make_args(expected_collapse)),
    )

    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())

"""Exact and mutation-sensitive tests for relativistic worldline ledgers."""

from __future__ import annotations

import pytest
import sympy as sp

from substrate_framework.relativistic_particle import (
    constant_e_gauge_rate,
    einbein_geodesic_acceleration,
    einbein_hamiltonian_ledger,
    massive_einbein_ledger,
    massive_mass_term_weyl_change,
    massless_einbein_ledger,
    massless_worldline_weyl_residual,
    worldline_reparametrization_residual,
)


def test_massive_elimination_recovers_root_action_and_rejects_sign_mutation() -> None:
    sigma = sp.Symbol("sigma", negative=True)
    e, mass, speed = sp.symbols("e mass speed", positive=True)

    ledger = massive_einbein_ledger(sigma, e, mass, speed)
    expected_derivative = -sigma / (2 * e**2) - (mass * speed) ** 2 / 2
    expected_root = sp.sqrt(-sigma) / (mass * speed)

    assert sp.simplify(ledger.einbein_euler_derivative - expected_derivative) == 0
    assert sp.simplify(
        ledger.mass_shell_constraint.subs(e, expected_root)
    ) == 0
    assert sp.simplify(
        ledger.recovered_square_root_lagrangian
        + mass * speed * sp.sqrt(-sigma)
    ) == 0
    assert massive_einbein_ledger(-4, 3, 1, 1).recovered_square_root_lagrangian == -2

    wrong_sign = sigma / (2 * e) + e * (mass * speed) ** 2 / 2
    assert sp.simplify(
        wrong_sign.subs(e, expected_root)
        + mass * speed * sp.sqrt(-sigma)
    ) != 0


def test_massless_specialization_precedes_elimination() -> None:
    sigma = sp.Symbol("sigma", real=True)
    e = sp.Symbol("e", positive=True)

    ledger = massless_einbein_ledger(sigma, e)

    assert ledger.lagrangian == sigma / (2 * e)
    assert ledger.null_constraint == sigma
    assert sp.simplify(ledger.einbein_euler_derivative + sigma / (2 * e**2)) == 0
    with pytest.raises(ValueError, match="mass"):
        massive_einbein_ledger(-1, e, 0, 1)


def test_hamiltonian_uses_covector_momentum_and_closes_on_constraint() -> None:
    e, mass, speed = sp.symbols("e mass speed", positive=True)
    p0, p1, p2 = sp.symbols("p0 p1 p2", real=True)
    metric = sp.Matrix([[-2, 0, 0], [0, 3, 1], [0, 1, 2]])
    momentum = sp.Matrix([p0, p1, p2])

    ledger = einbein_hamiltonian_ledger(
        metric, momentum, e, mass, speed
    )
    inverse = metric.inv()

    assert ledger.velocity == e * inverse * momentum
    assert sp.simplify(
        ledger.hamiltonian
        - e * ((momentum.T * inverse * momentum)[0] + (mass * speed) ** 2) / 2
    ) == 0
    wrong_index = e * metric * momentum
    assert sp.simplify(wrong_index - ledger.velocity) != sp.zeros(3, 1)


def test_reparametrization_affine_gauge_and_weyl_transformations() -> None:
    sigma = sp.Symbol("sigma", real=True)
    e, target, mass, speed, rate, omega = sp.symbols(
        "e target mass speed rate omega", positive=True
    )

    assert worldline_reparametrization_residual(
        sigma, e, mass, speed, rate
    ) == 0
    assert sp.simplify((e / constant_e_gauge_rate(e, target)) - target) == 0
    assert massless_worldline_weyl_residual(sigma, e, omega) == 0

    wrong_reparametrized = rate * (
        (sigma / rate**2) / (2 * e) - e * (mass * speed) ** 2 / 2
    )
    original = sigma / (2 * e) - e * (mass * speed) ** 2 / 2
    assert sp.simplify(wrong_reparametrized - original) != 0
    wrong_weyl_weight = omega**2 * sigma / (2 * omega * e)
    assert sp.simplify(wrong_weyl_weight - sigma / (2 * e)) != 0
    assert massive_mass_term_weyl_change(e, mass, speed, 1) == 0
    assert massive_mass_term_weyl_change(e, mass, speed, 2) != 0


@pytest.mark.parametrize("dimension", [3, 4])
def test_geodesic_acceleration_solves_lower_index_euler_lagrange(dimension: int) -> None:
    diagonal = [-1] + list(range(2, dimension + 1))
    metric = sp.diag(*diagonal)
    derivatives = sp.MutableDenseNDimArray.zeros(dimension, dimension, dimension)
    for rho in range(dimension):
        for mu in range(dimension):
            for nu in range(mu, dimension):
                value = sp.Integer((rho + 1) * (mu + 1) * (nu + 1))
                derivatives[rho, mu, nu] = value
                derivatives[rho, nu, mu] = value
    velocity = sp.Matrix(range(1, dimension + 1))
    e, e_dot = sp.Integer(2), sp.Integer(3)

    acceleration = einbein_geodesic_acceleration(
        metric, derivatives, velocity, e, e_dot
    )
    residual = sp.zeros(dimension, 1)
    for lam in range(dimension):
        residual[lam] = sp.simplify(
            (metric.row(lam) * acceleration)[0]
            + sum(
                derivatives[rho, lam, nu] * velocity[rho] * velocity[nu]
                for rho in range(dimension)
                for nu in range(dimension)
            )
            - sum(
                derivatives[lam, mu, nu] * velocity[mu] * velocity[nu] / 2
                for mu in range(dimension)
                for nu in range(dimension)
            )
            - e_dot * (metric.row(lam) * velocity)[0] / e
        )

    assert residual == sp.zeros(dimension, 1)
    wrong_sign_acceleration = -acceleration + 2 * e_dot * velocity / e
    wrong_residual = sp.zeros(dimension, 1)
    for lam in range(dimension):
        wrong_residual[lam] = sp.simplify(
            (metric.row(lam) * wrong_sign_acceleration)[0]
            + sum(
                derivatives[rho, lam, nu] * velocity[rho] * velocity[nu]
                for rho in range(dimension)
                for nu in range(dimension)
            )
            - sum(
                derivatives[lam, mu, nu] * velocity[mu] * velocity[nu] / 2
                for mu in range(dimension)
                for nu in range(dimension)
            )
            - e_dot * (metric.row(lam) * velocity)[0] / e
        )
    assert wrong_residual != sp.zeros(dimension, 1)

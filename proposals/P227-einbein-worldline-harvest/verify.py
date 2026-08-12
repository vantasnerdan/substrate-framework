#!/usr/bin/env python3
"""Exact proposal verifier for the issue #40 worldline harvest."""

from __future__ import annotations

import sympy as sp

from substrate_framework.lorentz_little_groups import (
    little_group_algebra_2plus1,
    little_group_algebra_3plus1,
)
from substrate_framework.lorentz_orbits import (
    unit_timelike_vector_orbit_metric,
)
from substrate_framework.relativistic_particle import (
    constant_e_gauge_rate,
    einbein_hamiltonian_ledger,
    massive_einbein_ledger,
    massive_mass_term_weyl_change,
    massless_einbein_ledger,
    massless_worldline_weyl_residual,
    worldline_reparametrization_residual,
)
from substrate_framework.verification import CheckLedger


def run() -> int:
    checks = CheckLedger("P227/worldline-harvest")
    sigma = sp.Symbol("sigma", negative=True)
    e, target, mass, speed, rate, omega = sp.symbols(
        "e target mass speed rate omega", positive=True
    )
    massive = massive_einbein_ledger(sigma, e, mass, speed)
    expected_root = sp.sqrt(-sigma) / (mass * speed)
    expected_action = -mass * speed * sp.sqrt(-sigma)
    checks.check(
        "positive massive root solves the auxiliary constraint",
        sp.simplify(massive.mass_shell_constraint.subs(e, expected_root)) == 0,
    )
    checks.check(
        "massive elimination exactly recovers the square-root density",
        sp.simplify(
            massive.recovered_square_root_lagrangian - expected_action
        )
        == 0,
    )
    checks.mutation_sensitive(
        "mass-term sign is load-bearing",
        lambda candidate: sp.simplify(candidate - expected_action) == 0,
        massive.recovered_square_root_lagrangian,
        (sigma / (2 * expected_root) + expected_root * (mass * speed) ** 2 / 2,),
    )

    massless = massless_einbein_ledger(sigma, e)
    checks.check(
        "massless specialization retains a regular action and null constraint",
        massless.lagrangian == sigma / (2 * e)
        and massless.null_constraint == sigma,
    )
    checks.check(
        "worldline reparametrization residual vanishes",
        worldline_reparametrization_residual(
            sigma, e, mass, speed, rate
        )
        == 0,
    )
    checks.check(
        "positive local rate attains a declared constant-e gauge",
        constant_e_gauge_rate(e, target).is_positive is True
        and sp.simplify(e / constant_e_gauge_rate(e, target) - target) == 0,
    )
    checks.check(
        "massless Weyl identity and nontrivial massive obstruction are separated",
        massless_worldline_weyl_residual(sigma, e, omega) == 0
        and massive_mass_term_weyl_change(e, mass, speed, 1) == 0
        and massive_mass_term_weyl_change(e, mass, speed, 2) != 0,
    )

    momentum = sp.Matrix(sp.symbols("p0:3", real=True))
    metric = sp.diag(-1, 2, 3)
    hamiltonian = einbein_hamiltonian_ledger(
        metric, momentum, e, mass, speed
    )
    checks.check(
        "Legendre transform closes on the covector mass-shell constraint",
        sp.simplify(
            hamiltonian.legendre_transform - hamiltonian.hamiltonian
        )
        == 0
        and hamiltonian.velocity == e * metric.inv() * momentum,
    )

    eta = sp.Symbol("eta", positive=True)
    theta = sp.Symbol("theta", positive=True)
    checks.check(
        "future unit-vector orbit metrics are exact H2 and H3 pullbacks",
        unit_timelike_vector_orbit_metric(3)
        == sp.diag(1, sp.sinh(eta) ** 2)
        and unit_timelike_vector_orbit_metric(4)
        == sp.diag(
            1,
            sp.sinh(eta) ** 2,
            sp.sinh(eta) ** 2 * sp.sin(theta) ** 2,
        ),
    )

    ledger_2plus1 = little_group_algebra_2plus1()
    ledger_3plus1 = little_group_algebra_3plus1()
    t1, t2, j3 = ledger_3plus1.massless_generators
    checks.check(
        "displayed null stabilizers fix their standard momenta",
        ledger_2plus1.massless_generators[0]
        * ledger_2plus1.standard_null_momentum
        == sp.zeros(3, 1)
        and all(
            generator * ledger_3plus1.standard_null_momentum
            == sp.zeros(4, 1)
            for generator in ledger_3plus1.massless_generators
        ),
    )
    checks.check(
        "3plus1 null stabilizers close as displayed iso2",
        t1 * t2 - t2 * t1 == sp.zeros(4)
        and j3 * t1 - t1 * j3 == -t2
        and j3 * t2 - t2 * j3 == t1,
    )
    wrong_t1 = ledger_3plus1.nonstabilizing_boosts[0] - ledger_3plus1.massive_generators[1]
    checks.check(
        "wrong little-group sign mutation fails the fixed-vector oracle",
        wrong_t1 * ledger_3plus1.standard_null_momentum != sp.zeros(4, 1),
    )
    return checks.finish()


if __name__ == "__main__":
    run()

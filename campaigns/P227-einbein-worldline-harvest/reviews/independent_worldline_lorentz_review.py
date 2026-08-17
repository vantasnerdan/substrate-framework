#!/usr/bin/env python3
"""Independent exact review of the five P227 claims.

The derivations below use direct SymPy algebra.  Canonical P227 APIs are
imported only after the independent expressions and mutations have been
constructed, so agreement with the package is a comparison rather than a
rephrasing of the implementation.
"""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def _commutator(left: sp.Matrix, right: sp.Matrix) -> sp.Matrix:
    return left * right - right * left


def _is_lorentz_generator(generator: sp.Matrix, metric: sp.Matrix) -> bool:
    return sp.simplify(metric * generator + (metric * generator).T) == sp.zeros(
        metric.rows
    )


def _pullback(
    embedding: sp.Matrix,
    coordinates: tuple[sp.Symbol, ...],
    ambient_metric: sp.Matrix,
) -> sp.Matrix:
    jacobian = embedding.jacobian(coordinates)
    return sp.simplify(jacobian.T * ambient_metric * jacobian)


def run() -> int:
    checks = CheckLedger("P227/independent-worldline-lorentz-review")

    # C-WLN-001: derive the auxiliary equation and Hamiltonian directly.
    sigma = sp.Symbol("sigma", negative=True)
    e, mass, speed = sp.symbols("e mass speed", positive=True)
    lagrangian = sigma / (2 * e) - e * (mass * speed) ** 2 / 2
    euler_derivative = sp.diff(lagrangian, e)
    constraint = sp.simplify(-2 * e**2 * euler_derivative)
    positive_root = sp.sqrt(-sigma) / (mass * speed)
    recovered = sp.simplify(lagrangian.subs(e, positive_root))
    expected_root_density = -mass * speed * sp.sqrt(-sigma)
    checks.check(
        "C-WLN-001 direct auxiliary equation and positive branch",
        sp.simplify(constraint - (sigma + e**2 * (mass * speed) ** 2)) == 0
        and sp.simplify(constraint.subs(e, positive_root)) == 0,
    )
    checks.check(
        "C-WLN-001 positive elimination recovers the root density",
        sp.simplify(recovered - expected_root_density) == 0,
    )
    wrong_mass_sign = sigma / (2 * positive_root) + positive_root * (
        mass * speed
    ) ** 2 / 2
    checks.mutation_sensitive(
        "C-WLN-001 mass-term sign is load-bearing",
        lambda candidate: sp.simplify(candidate - expected_root_density) == 0,
        recovered,
        (wrong_mass_sign,),
    )

    p0, p1, p2 = sp.symbols("p0 p1 p2", real=True)
    momentum = sp.Matrix([p0, p1, p2])
    metric = sp.Matrix([[-2, 0, 0], [0, 3, 1], [0, 1, 2]])
    inverse = metric.inv()
    velocity = e * inverse * momentum
    velocity_norm = sp.simplify((velocity.T * metric * velocity)[0])
    phase_lagrangian = velocity_norm / (2 * e) - e * (mass * speed) ** 2 / 2
    legendre = sp.simplify((momentum.T * velocity)[0] - phase_lagrangian)
    expected_hamiltonian = sp.simplify(
        e * ((momentum.T * inverse * momentum)[0] + (mass * speed) ** 2) / 2
    )
    checks.check(
        "C-WLN-001 covector Legendre transform closes on the constraint",
        sp.simplify(legendre - expected_hamiltonian) == 0,
    )
    checks.mutation_sensitive(
        "C-WLN-001 inverse metric is load-bearing",
        lambda candidate: sp.simplify(candidate - velocity) == sp.zeros(3, 1),
        velocity,
        (e * metric * momentum,),
    )

    # C-WLN-002: derive every lower-index Euler residual component on a
    # nonconstant metric.  The einbein is allowed to be nonconstant.
    x = sp.Symbol("x", real=True)
    t_dot, x_dot, e_dot = sp.symbols("t_dot x_dot e_dot", real=True)
    radial_factor = 2 + x**2
    local_metric = sp.diag(-1, radial_factor)
    local_velocity = sp.Matrix([t_dot, x_dot])
    christoffel_x_xx = x / radial_factor
    acceleration = sp.Matrix(
        [
            e_dot * t_dot / e,
            -christoffel_x_xx * x_dot**2 + e_dot * x_dot / e,
        ]
    )
    metric_derivatives = sp.MutableDenseNDimArray.zeros(2, 2, 2)
    metric_derivatives[1, 1, 1] = 2 * x

    def lower_euler_residual(candidate_acceleration: sp.Matrix) -> sp.Matrix:
        residual = sp.zeros(2, 1)
        for lam in range(2):
            residual[lam] = sp.simplify(
                (local_metric.row(lam) * candidate_acceleration)[0]
                + sum(
                    metric_derivatives[rho, lam, nu]
                    * local_velocity[rho]
                    * local_velocity[nu]
                    for rho in range(2)
                    for nu in range(2)
                )
                - sum(
                    metric_derivatives[lam, mu, nu]
                    * local_velocity[mu]
                    * local_velocity[nu]
                    / 2
                    for mu in range(2)
                    for nu in range(2)
                )
                - e_dot
                * (local_metric.row(lam) * local_velocity)[0]
                / e
            )
        return residual

    checks.check(
        "C-WLN-002 both lower-index Euler residuals vanish",
        lower_euler_residual(acceleration) == sp.zeros(2, 1),
    )
    wrong_geodesic_sign = sp.Matrix(
        [
            e_dot * t_dot / e,
            christoffel_x_xx * x_dot**2 + e_dot * x_dot / e,
        ]
    )
    missing_einbein_term = sp.Matrix([0, -christoffel_x_xx * x_dot**2])
    checks.mutation_sensitive(
        "C-WLN-002 geodesic-force sign is load-bearing",
        lambda candidate: lower_euler_residual(candidate) == sp.zeros(2, 1),
        acceleration,
        (wrong_geodesic_sign,),
    )
    checks.mutation_sensitive(
        "C-WLN-002 nonaffine einbein term is load-bearing",
        lambda candidate: lower_euler_residual(candidate) == sp.zeros(2, 1),
        acceleration,
        (missing_einbein_term,),
    )
    massless_density = sigma / (2 * e)
    checks.check(
        "C-WLN-002 massless variation gives the null constraint",
        sp.simplify(-2 * e**2 * sp.diff(massless_density, e) - sigma) == 0,
    )

    # C-WLN-003: transform density and measure directly.
    rate, target, omega = sp.symbols("rate target omega", positive=True)
    transformed_density = (sigma / rate**2) / (2 * (e / rate)) - (
        e / rate
    ) * (mass * speed) ** 2 / 2
    checks.check(
        "C-WLN-003 reparametrized density times Jacobian is invariant",
        sp.simplify(rate * transformed_density - lagrangian) == 0,
    )
    gauge_rate = e / target
    checks.check(
        "C-WLN-003 local gauge rate reaches the declared constant einbein",
        sp.simplify(e / gauge_rate - target) == 0,
    )
    checks.mutation_sensitive(
        "C-WLN-003 gauge-rate normalization is load-bearing",
        lambda candidate: sp.simplify(e / candidate - target) == 0,
        gauge_rate,
        (e / (2 * target),),
    )
    massless_weyl_residual = sp.simplify(
        omega**2 * sigma / (2 * omega**2 * e) - sigma / (2 * e)
    )
    massive_weyl_change = sp.simplify(
        -omega**2 * e * (mass * speed) ** 2 / 2
        + e * (mass * speed) ** 2 / 2
    )
    checks.check(
        "C-WLN-003 massless Weyl identity and massive obstruction separate",
        massless_weyl_residual == 0
        and massive_weyl_change.subs(omega, 1) == 0
        and massive_weyl_change.subs(omega, 2) != 0,
    )
    wrong_weyl_weight = sp.simplify(
        omega**2 * sigma / (2 * omega * e) - sigma / (2 * e)
    )
    checks.mutation_sensitive(
        "C-WLN-003 einbein Weyl weight is load-bearing",
        lambda candidate: sp.simplify(candidate) == 0,
        massless_weyl_residual,
        (wrong_weyl_weight,),
    )

    # C-LOR-001: independently pull back the two hyperboloid embeddings.
    rapidity = sp.Symbol("eta", positive=True, real=True)
    theta, phi = sp.symbols("theta phi", positive=True, real=True)
    h2_embedding = sp.Matrix(
        [
            sp.cosh(rapidity),
            sp.sinh(rapidity) * sp.cos(theta),
            sp.sinh(rapidity) * sp.sin(theta),
        ]
    )
    eta_2plus1 = sp.diag(-1, 1, 1)
    h2_metric = _pullback(h2_embedding, (rapidity, theta), eta_2plus1)
    expected_h2 = sp.diag(1, sp.sinh(rapidity) ** 2)
    checks.check(
        "C-LOR-001 H2 embedding has unit timelike norm and exact pullback",
        sp.trigsimp((h2_embedding.T * eta_2plus1 * h2_embedding)[0]) == -1
        and sp.simplify(h2_metric - expected_h2) == sp.zeros(2),
    )
    h3_embedding = sp.Matrix(
        [
            sp.cosh(rapidity),
            sp.sinh(rapidity) * sp.sin(theta) * sp.cos(phi),
            sp.sinh(rapidity) * sp.sin(theta) * sp.sin(phi),
            sp.sinh(rapidity) * sp.cos(theta),
        ]
    )
    eta_3plus1 = sp.diag(-1, 1, 1, 1)
    h3_metric = _pullback(
        h3_embedding, (rapidity, theta, phi), eta_3plus1
    )
    expected_h3 = sp.diag(
        1,
        sp.sinh(rapidity) ** 2,
        sp.sinh(rapidity) ** 2 * sp.sin(theta) ** 2,
    )
    checks.check(
        "C-LOR-001 H3 embedding has unit timelike norm and exact pullback",
        sp.trigsimp((h3_embedding.T * eta_3plus1 * h3_embedding)[0]) == -1
        and sp.simplify(h3_metric - expected_h3) == sp.zeros(3),
    )
    checks.mutation_sensitive(
        "C-LOR-001 ambient signature is load-bearing",
        lambda candidate: sp.simplify(candidate - expected_h2) == sp.zeros(2),
        h2_metric,
        (_pullback(h2_embedding, (rapidity, theta), -eta_2plus1),),
    )

    # C-LOR-002: build boosts and rotations directly in the displayed basis.
    k1_3 = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    k2_3 = sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
    j12_3 = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, -1, 0]])
    p_massive_3 = sp.Matrix([1, 0, 0])
    p_null_3 = sp.Matrix([1, 1, 0])
    null_3 = k2_3 + j12_3
    checks.check(
        "C-LOR-002 exact 2plus1 stabilizer and boost products",
        _is_lorentz_generator(j12_3, eta_2plus1)
        and _is_lorentz_generator(null_3, eta_2plus1)
        and j12_3 * p_massive_3 == sp.zeros(3, 1)
        and k1_3 * p_massive_3 == sp.Matrix([0, 1, 0])
        and k2_3 * p_massive_3 == sp.Matrix([0, 0, 1])
        and null_3 * p_null_3 == sp.zeros(3, 1),
    )
    checks.mutation_sensitive(
        "C-LOR-002 2plus1 relative boost-rotation sign is load-bearing",
        lambda candidate: candidate * p_null_3 == sp.zeros(3, 1),
        null_3,
        (k2_3 - j12_3,),
    )

    def boost4(axis: int) -> sp.Matrix:
        result = sp.zeros(4)
        result[0, axis] = 1
        result[axis, 0] = 1
        return result

    j1 = sp.Matrix(
        [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]]
    )
    j2 = sp.Matrix(
        [[0, 0, 0, 0], [0, 0, 0, -1], [0, 0, 0, 0], [0, 1, 0, 0]]
    )
    j3 = sp.Matrix(
        [[0, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 0]]
    )
    k1, k2, _k3 = (boost4(axis) for axis in (1, 2, 3))
    t1 = k1 + j2
    t2 = k2 - j1
    p_massive_4 = sp.Matrix([1, 0, 0, 0])
    p_null_4 = sp.Matrix([1, 0, 0, 1])
    checks.check(
        "C-LOR-002 3plus1 massive stabilizer closes in repository convention",
        all(_is_lorentz_generator(generator, eta_3plus1) for generator in (j1, j2, j3))
        and all(generator * p_massive_4 == sp.zeros(4, 1) for generator in (j1, j2, j3))
        and _commutator(j1, j2) == -j3
        and _commutator(j2, j3) == -j1
        and _commutator(j3, j1) == -j2,
    )
    checks.check(
        "C-LOR-002 3plus1 null stabilizer is exact iso2",
        all(_is_lorentz_generator(generator, eta_3plus1) for generator in (t1, t2, j3))
        and all(generator * p_null_4 == sp.zeros(4, 1) for generator in (t1, t2, j3))
        and _commutator(t1, t2) == sp.zeros(4)
        and _commutator(j3, t1) == -t2
        and _commutator(j3, t2) == t1,
    )
    checks.mutation_sensitive(
        "C-LOR-002 T1 relative sign is load-bearing",
        lambda candidate: candidate * p_null_4 == sp.zeros(4, 1),
        t1,
        (k1 - j2,),
    )
    checks.mutation_sensitive(
        "C-LOR-002 T2 relative sign is load-bearing",
        lambda candidate: candidate * p_null_4 == sp.zeros(4, 1),
        t2,
        (k2 + j1,),
    )
    checks.mutation_sensitive(
        "C-LOR-002 iso2 commutator sign is load-bearing",
        lambda candidate: candidate == -t2,
        _commutator(j3, t1),
        (t2,),
    )

    # Compare the independently derived objects with the canonical package.
    from substrate_framework.lorentz_little_groups import (
        little_group_algebra_2plus1,
        little_group_algebra_3plus1,
    )
    from substrate_framework.lorentz_orbits import (
        unit_timelike_vector_orbit_metric,
    )
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

    canonical_massive = massive_einbein_ledger(sigma, e, mass, speed)
    canonical_hamiltonian = einbein_hamiltonian_ledger(
        metric, momentum, e, mass, speed
    )
    checks.check(
        "C-WLN-001 independent expressions equal canonical ledgers",
        canonical_massive.mass_shell_constraint == constraint
        and canonical_massive.recovered_square_root_lagrangian == recovered
        and canonical_hamiltonian.velocity == velocity
        and canonical_hamiltonian.hamiltonian == expected_hamiltonian,
    )
    checks.check(
        "C-WLN-002 independent acceleration equals canonical helper",
        massless_einbein_ledger(sigma, e).null_constraint == sigma
        and sp.simplify(
            einbein_geodesic_acceleration(
                local_metric,
                metric_derivatives,
                local_velocity,
                e,
                e_dot,
            )
            - acceleration
        )
        == sp.zeros(2, 1),
    )
    checks.check(
        "C-WLN-003 independent transformations equal canonical residuals",
        worldline_reparametrization_residual(
            sigma, e, mass, speed, rate
        )
        == 0
        and constant_e_gauge_rate(e, target) == gauge_rate
        and massless_worldline_weyl_residual(sigma, e, omega) == 0
        and massive_mass_term_weyl_change(e, mass, speed, omega)
        == massive_weyl_change,
    )
    checks.check(
        "C-LOR-001 independent pullbacks equal canonical metrics",
        unit_timelike_vector_orbit_metric(3) == h2_metric
        and unit_timelike_vector_orbit_metric(4) == h3_metric,
    )
    canonical_3 = little_group_algebra_2plus1()
    canonical_4 = little_group_algebra_3plus1()
    checks.check(
        "C-LOR-002 independent matrices equal canonical ledgers",
        canonical_3.massive_generators == (sp.ImmutableMatrix(j12_3),)
        and canonical_3.massless_generators == (sp.ImmutableMatrix(null_3),)
        and canonical_4.massive_generators
        == tuple(sp.ImmutableMatrix(generator) for generator in (j1, j2, j3))
        and canonical_4.massless_generators
        == tuple(sp.ImmutableMatrix(generator) for generator in (t1, t2, j3)),
    )
    return checks.finish()


if __name__ == "__main__":
    run()

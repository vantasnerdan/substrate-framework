"""Independent exact checks for P253/0094."""

from __future__ import annotations

import sympy as sp

from substrate_framework.euler_action_locking import (
    action_locked_by_fixed_charge,
    carrier_complement_cohomology,
    compact_phase_charge,
    compact_phase_charge_per_action,
    electromagnetic_candidate_action_unit,
    euler_similarity_action_factor,
    leading_ring_scale_balance,
    reduced_energy_local_ledger,
    smooth_two_label_bf_ledger,
)


def check(label: str, condition) -> None:
    if condition is not True and condition != sp.true:
        raise AssertionError(f"FAIL {label}: {condition}")
    print(f"PASS {label}")


def main() -> None:
    delta, a2 = sp.symbols("delta a2", positive=True)
    energy = reduced_energy_local_ledger(-delta, a2)
    check("positive interior leading action", energy.leading_action == delta / (2 * a2))
    check("positive leading second derivative", energy.leading_second_derivative == 2 * a2)
    check("remainder and invariant reduction remain required", energy.remainder_bound_required and energy.invariant_reduction_required)

    g0, s0, number, charge = sp.symbols("g0 S0 N Q", positive=True)
    m = sp.Integer(-2)
    kappa_m = compact_phase_charge_per_action(g0, m, s0)
    phase_charge = compact_phase_charge(g0, m, number)
    check("typed compact phase moment map", sp.simplify(phase_charge - kappa_m * s0 * number) == 0)
    locked = action_locked_by_fixed_charge(charge, g0, 1, s0)
    check("fixed charge retains action-unit direction", sp.diff(locked, s0) == charge / g0)

    coefficient, epsilon, speed = sp.symbols("C epsilon c", positive=True)
    candidate = electromagnetic_candidate_action_unit(g0, epsilon, speed, coefficient)
    check("dimensionless action coefficient remains free", sp.diff(candidate, coefficient) == g0**2 / (4 * sp.pi * epsilon * speed))

    topology = carrier_complement_cohomology()
    check("point complement Gauss row", (topology.point_h1_rank, topology.point_h2_rank) == (0, 1))
    check("ring complement linking row", (topology.ring_h1_rank, topology.ring_h2_rank) == (1, 0))

    labels = smooth_two_label_bf_ledger()
    check("smooth two-label BF source vanishes", labels.exterior_derivative_vanishes_away_from_defects and not labels.supplies_smooth_bulk_electric_source)

    scale_a, scale_b = sp.symbols("A B", nonzero=True, real=True)
    check("Euler similarity action factor", euler_similarity_action_factor(scale_a, scale_b) == scale_a / scale_b**4)

    kappa, charge_q, rho, logarithm = sp.symbols(
        "kappa Q rho L", positive=True
    )
    scale = leading_ring_scale_balance(
        kappa, charge_q, rho, epsilon, logarithm
    )
    check(
        "leading fluid-Coulomb finite ring radius",
        scale.radius
        == charge_q / (2 * sp.pi * kappa * sp.sqrt(rho * epsilon)),
    )
    check(
        "leading balance impulse and strict minimum",
        scale.axial_impulse
        == charge_q**2 / (4 * sp.pi * epsilon * kappa)
        and scale.second_derivative.is_positive,
    )


if __name__ == "__main__":
    main()

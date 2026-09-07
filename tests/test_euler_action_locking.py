import pytest
import sympy as sp

from substrate_framework.euler_action_locking import (
    action_locked_by_fixed_charge,
    carrier_complement_cohomology,
    compact_phase_charge,
    compact_phase_charge_per_action,
    electromagnetic_candidate_action_unit,
    euler_similarity_action_factor,
    leading_nonzero_action,
    leading_ring_scale_balance,
    reduced_energy_local_ledger,
    smooth_two_label_bf_ledger,
)


def test_reduced_energy_leading_root_and_sign_classification():
    delta, a2 = sp.symbols("delta a2", positive=True)
    ledger = reduced_energy_local_ledger(-delta, a2)
    assert ledger.leading_action == delta / (2 * a2)
    assert ledger.leading_second_derivative == 2 * a2
    assert ledger.positive_interior_candidate
    assert not ledger.zero_is_boundary_selected
    assert ledger.remainder_bound_required
    assert ledger.invariant_reduction_required
    positive = reduced_energy_local_ledger(delta, a2)
    assert positive.zero_is_boundary_selected
    assert not positive.positive_interior_candidate
    assert leading_nonzero_action(delta, a2) == -delta / (2 * a2)


def test_compact_character_keeps_continuous_number_and_explicit_action_unit():
    g0, s0, number = sp.symbols("g0 S0 N", nonzero=True, real=True)
    m = sp.Integer(-2)
    kappa_m = compact_phase_charge_per_action(g0, m, s0)
    charge = compact_phase_charge(g0, m, number)
    action = s0 * number
    assert kappa_m == -2 * g0 / s0
    assert charge == -2 * g0 * number
    assert sp.simplify(charge - kappa_m * action) == 0
    with pytest.raises(ValueError, match="character_weight"):
        compact_phase_charge_per_action(g0, sp.Rational(1, 2), s0)


def test_fixed_charge_lock_depends_continuously_on_imported_action_unit():
    q, g0, s0 = sp.symbols("Q g0 S0", positive=True)
    locked = action_locked_by_fixed_charge(q, g0, 1, s0)
    assert locked == q * s0 / g0
    assert sp.diff(locked, s0) == q / g0


def test_electromagnetic_action_unit_retains_arbitrary_dimensionless_factor():
    coefficient, g0, epsilon, speed = sp.symbols(
        "C g0 epsilon c", positive=True
    )
    action_unit = electromagnetic_candidate_action_unit(
        g0, epsilon, speed, coefficient
    )
    assert action_unit == coefficient * g0**2 / (
        4 * sp.pi * epsilon * speed
    )
    assert sp.diff(action_unit, coefficient) == g0**2 / (
        4 * sp.pi * epsilon * speed
    )


def test_point_and_ring_complements_have_different_integral_rows():
    ledger = carrier_complement_cohomology()
    assert (ledger.point_h1_rank, ledger.point_h2_rank) == (0, 1)
    assert (ledger.ring_h1_rank, ledger.ring_h2_rank) == (1, 0)


def test_smooth_two_label_pullback_has_no_bulk_bf_source():
    ledger = smooth_two_label_bf_ledger()
    assert ledger.exterior_derivative_vanishes_away_from_defects
    assert not ledger.supplies_smooth_bulk_electric_source
    assert ledger.defect_or_extra_compact_field_required


def test_bare_euler_similarity_leaves_action_continuous():
    scale_a, scale_b = sp.symbols("A B", nonzero=True, real=True)
    assert euler_similarity_action_factor(scale_a, scale_b) == (
        scale_a / scale_b**4
    )


def test_leading_thin_ring_fluid_coulomb_balance_has_finite_radius():
    kappa, charge, rho, epsilon, logarithm = sp.symbols(
        "kappa Q rho epsilon L", positive=True
    )
    ledger = leading_ring_scale_balance(
        kappa, charge, rho, epsilon, logarithm
    )
    assert ledger.radius == charge / (
        2 * sp.pi * kappa * sp.sqrt(rho * epsilon)
    )
    assert ledger.axial_impulse == charge**2 / (
        4 * sp.pi * epsilon * kappa
    )
    assert ledger.second_derivative.is_positive

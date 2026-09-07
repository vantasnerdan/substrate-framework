import pytest
import sympy as sp

from substrate_framework.euler_gauge_action_selection import (
    charge_action_map_jacobian,
    compact_bf_route_ledger,
    compact_phase_moment_charge,
    compact_u1_topology_ledger,
    conditional_self_consistent_action,
    electromagnetic_action_scale,
    gauge_action_ratio,
    gauge_field_rescaling,
    transported_tag_charge,
)


def test_electromagnetic_action_and_field_rescaling_invariants():
    epsilon, mu, g, scale, tag = sp.symbols(
        "epsilon mu g scale tag", positive=True
    )
    ledger = gauge_field_rescaling(epsilon, mu, g, scale)
    assert ledger.maxwell_speed_squared == 1 / (epsilon * mu)
    assert ledger.charge_action_numerator == g**2 / epsilon
    c = 1 / sp.sqrt(epsilon * mu)
    q = g * tag
    q_new = ledger.coupling * tag
    assert sp.simplify(
        electromagnetic_action_scale(q_new, ledger.permittivity, c)
        - electromagnetic_action_scale(q, epsilon, c)
    ) == 0


def test_transport_tag_and_compact_phase_momentum_remain_continuous():
    g_tag, g_0, c_tag, lam, n = sp.symbols(
        "g_tag g_0 c_tag lam n", real=True
    )
    assert transported_tag_charge(g_tag, c_tag, lam) == g_tag * c_tag * lam
    assert compact_phase_moment_charge(g_0, n) == g_0 * n
    assert compact_phase_moment_charge(g_0, n, -1) == -g_0 * n
    assert sp.diff(transported_tag_charge(g_tag, c_tag, lam), lam) == (
        g_tag * c_tag
    )
    assert sp.diff(compact_phase_moment_charge(g_0, n), n) == g_0
    with pytest.raises(ValueError, match="character_weight"):
        compact_phase_moment_charge(g_0, n, sp.Rational(1, 2))


def test_compact_u1_topology_separates_whole_space_and_boundary_flux():
    ledger = compact_u1_topology_ledger()
    assert ledger.pi_3_u1 == 0
    assert ledger.h2_s3 == 0
    assert ledger.h2_s2 == "Z"
    assert ledger.whole_space_magnetic_flux == 0


def test_direct_euler_vorticity_does_not_supply_compact_bf_data():
    ledger = compact_bf_route_ledger()
    assert not ledger.euler_vorticity_is_compact_connection
    assert not ledger.euler_vorticity_periods_are_integral
    assert ledger.compact_two_form_is_added_structure
    assert not ledger.level_is_euler_selected


def test_fixed_charge_leaves_an_independent_mode_action_direction():
    g, c_tag, lam, action, epsilon, speed = sp.symbols(
        "g c_tag lam action epsilon speed", positive=True
    )
    jac = charge_action_map_jacobian(
        g, c_tag, lam, action, epsilon, speed
    )
    assert jac[:, 1] == sp.Matrix([0, 1, 0])
    assert jac.rank() == 2


def test_same_phase_identification_is_an_added_constraint_not_dynamics():
    epsilon, speed, kappa_m, action = sp.symbols(
        "epsilon speed kappa_m action", positive=True
    )
    charge = kappa_m * action
    root = conditional_self_consistent_action(kappa_m, epsilon, speed)
    equation = sp.factor(
        action - electromagnetic_action_scale(charge, epsilon, speed)
    )
    assert sp.simplify(equation.subs(action, root)) == 0
    assert sp.solve(sp.Eq(equation, 0), action) == [root]
    assert sp.simplify(gauge_action_ratio(charge, epsilon, speed, action)) == (
        kappa_m**2 * action / (4 * sp.pi * epsilon * speed)
    )


def test_conditional_root_preserves_both_coupling_signs():
    epsilon, speed, kappa_m = sp.symbols(
        "epsilon speed kappa_m", positive=True
    )
    assert conditional_self_consistent_action(-kappa_m, epsilon, speed) == (
        conditional_self_consistent_action(kappa_m, epsilon, speed)
    )
    with pytest.raises(ValueError, match="charge_per_action"):
        conditional_self_consistent_action(0, epsilon, speed)


def test_action_helpers_reject_nonpositive_domains():
    with pytest.raises(ValueError, match="permittivity"):
        electromagnetic_action_scale(1, 0, 1)
    with pytest.raises(ValueError, match="mode_action"):
        gauge_action_ratio(1, 1, 1, 0)
    with pytest.raises(ValueError, match="field_scale"):
        gauge_field_rescaling(1, 1, 1, -1)

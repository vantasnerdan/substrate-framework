from __future__ import annotations

import sympy as sp

import substrate_framework as framework
from substrate_framework.sg_gravity_chain import (
    MU_THIRD_POWER_LOWER,
    MU_THIRD_POWER_MID,
    MU_THIRD_POWER_UPPER,
    assemble_chain,
    breather_moment_third_physical_scale,
    breather_quadrupole_power,
    light_deflection_angle,
    newtonian_potential,
    schwarzschild_areal_radius_factor,
)
from substrate_framework.sg_spectral_induction import (
    newton_q_normalized,
    tower_inverse_g_normalized,
)


def test_sg_gravity_chain_api_is_exported_from_package() -> None:
    assert framework.assemble_chain is assemble_chain
    assert framework.breather_quadrupole_power is breather_quadrupole_power
    assert framework.light_deflection_angle is light_deflection_angle


def test_unit_point_recovers_the_c_med_003_scales() -> None:
    chain = assemble_chain(1, 1, 1, action_quantum=1)
    assert chain.signal_speed == 1
    assert chain.gap_frequency == 1
    assert chain.length == 1
    assert chain.energy_scale == 1
    assert chain.action_scale == 1


def test_derived_newton_constant_has_the_c_grv_001_monomial_form() -> None:
    lam, tension, mu = (sp.Symbol(s, positive=True) for s in ("lam", "T", "mu"))
    hbar = sp.Symbol("hbar", positive=True)
    chain = assemble_chain(lam, tension, mu, action_quantum=1)
    g = chain.newton_constant(hbar)
    q = newton_q_normalized(1)
    ell = sp.sqrt(tension / mu)
    c = sp.sqrt(tension / lam)
    # Exact identity: G == q * ell**2 * c**3 / hbar (the C-GRV-001 ceiling,
    # now with q supplied by the tower rather than left free).
    assert sp.simplify(g - q * ell**2 * c**3 / hbar) == 0


def test_q_is_unit_standard_invariant() -> None:
    # q must not move when the upstream coefficient ray is rescaled:
    # (lam, T, mu) -> (s*lam, s*T, s*mu) rescales every scale but not q.
    chain_a = assemble_chain(1, 1, 1, action_quantum=1)
    chain_b = assemble_chain(7, 7, 7, action_quantum=1)
    assert sp.simplify(chain_a.q_dimensionless - chain_b.q_dimensionless) == 0


def test_vacuum_exterior_factor_is_schwarzschild() -> None:
    m, x = sp.symbols("m x", positive=True)
    assert schwarzschild_areal_radius_factor(m, x) == 1 - 2 * m / x


def test_newtonian_potential_uses_the_derived_coupling() -> None:
    hbar = sp.Symbol("hbar", positive=True)
    M, r = sp.symbols("M r", positive=True)
    chain = assemble_chain(1, 1, 1, action_quantum=1)
    g = chain.newton_constant(hbar)
    phi = newtonian_potential(M, r, g, chain.signal_speed)
    s_total = tower_inverse_g_normalized(1)
    assert sp.simplify(phi - (-M / (s_total * hbar * r))) == 0


def test_deflection_angle_carries_the_derived_q() -> None:
    hbar = sp.Symbol("hbar", positive=True)
    M, b = sp.symbols("M b", positive=True)
    chain = assemble_chain(1, 1, 1, action_quantum=1)
    theta = light_deflection_angle(M, b, chain.newton_constant(hbar), chain.signal_speed)
    q = newton_q_normalized(1)
    assert sp.simplify(theta - 4 * q * M / (hbar * b)) == 0


def test_mu_third_power_midpoint_respects_the_accepted_bracket() -> None:
    assert MU_THIRD_POWER_LOWER < MU_THIRD_POWER_MID < MU_THIRD_POWER_UPPER


def test_radiative_power_is_positive_and_proportional_to_q() -> None:
    hbar = sp.Symbol("hbar", positive=True)
    chain = assemble_chain(1, 1, 1, action_quantum=1)
    p1 = breather_quadrupole_power(chain, hbar)
    assert float(p1.subs(hbar, 1)) > 0.0
    # Mutating the tower (dropping kinks) must move the power by exactly the
    # ratio of the inverse couplings: the power carries the derived G.
    chain_no_kink = assemble_chain(1, 1, 1, action_quantum=1, kink_multiplicity=0)
    p2 = breather_quadrupole_power(chain_no_kink, hbar)
    ratio = sp.N((p1 / p2).subs(hbar, 1), 12)
    expected = sp.N(
        tower_inverse_g_normalized(1, 0, 0) / tower_inverse_g_normalized(1, 0, 2), 12
    )
    assert abs(float(ratio - expected)) < 1e-10


def test_radiative_power_scaling_under_onsite_mutation() -> None:
    # Full accounting: scale**2 ~ E_scale**2 * ell**4 * omega_0**6 ~ mu**2,
    # and the derived G ~ q * ell**2 * c**3 / hbar contributes ell**2 ~ mu**-1;
    # net power scales as mu**1, so doubling the onsite coefficient doubles P.
    hbar = sp.Symbol("hbar", positive=True)
    base = breather_quadrupole_power(assemble_chain(1, 1, 1, action_quantum=1), hbar)
    doubled = breather_quadrupole_power(assemble_chain(1, 1, 2, action_quantum=1), hbar)
    ratio = sp.N((doubled / base).subs(hbar, 1), 12)
    assert abs(float(ratio) - 2.0) < 1e-10


def test_moment_scale_matches_the_declared_conversion() -> None:
    lam, tension, mu = (sp.Symbol(s, positive=True) for s in ("lam", "T", "mu"))
    chain = assemble_chain(lam, tension, mu, action_quantum=1)
    scale = breather_moment_third_physical_scale(chain, sp.Symbol("hbar", positive=True))
    c = sp.sqrt(tension / lam)
    expected = (sp.sqrt(tension * mu) / c**2) * (tension / mu) * (mu / lam) ** sp.Rational(3, 2)
    assert sp.simplify(scale - expected) == 0

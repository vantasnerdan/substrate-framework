"""Machine checks for src/substrate_framework/chromomagnetic_background.py (P229, issue #45).

Oracle strategy per the quantitative-verification skill:
- the algebraic spectrum is checked against the discretized fluctuation
  operator (never asserted);
- the heat-kernel coefficients are exact SymPy series;
- the two proper-time routes (closed hyperbolic form vs explicit Landau sum)
  must agree numerically;
- the MS-bar log coefficient and constant are extracted from the cutoff
  integral by declared subtractions and must match the analytic values;
- every check has a mutation that breaks it.
"""

import numpy as np
import sympy as sp
import pytest

from substrate_framework import chromomagnetic_background as cb
from substrate_framework.gauge_beta import GaugeFactor, product_gauge_coefficients


# ---------------------------------------------------------------------------
# Spectrum from the discretized operator
# ---------------------------------------------------------------------------


def test_orbital_landau_levels_from_operator():
    for gB in (0.7, 1.3):
        op = cb.landau_orbital_operator(gB, L=14.0, npoints=500)
        w = np.linalg.eigvalsh(op)[:5]
        exact = np.array([gB * (2 * n + 1) for n in range(5)])
        np.testing.assert_allclose(w, exact, rtol=2e-3)


def test_spin_generator_eigenvalues():
    vals = cb.spin_generator_12().eigenvals()
    assert vals == {1: 1, -1: 1, 0: 2}


def test_vector_spectrum_from_operator_includes_tachyon():
    gB = 1.0
    op = cb.vector_fluctuation_operator(gB, L=14.0, npoints=400)
    w = np.linalg.eigvalsh(op)[:8].real
    # n=0, s3=+1: E^2 = -gB (Nielsen-Olesen); then threefold E^2 = +gB
    assert w[0] == pytest.approx(-gB, rel=2e-3)
    np.testing.assert_allclose(w[1:4], gB * np.ones(3), rtol=2e-3)


def test_algebraic_spectrum_matches_encoder():
    # encoder self-consistency at the tachyon and a generic level
    assert cb.charged_vector_mass2(1.3, 0, +1) == pytest.approx(-1.3)
    assert cb.charged_vector_mass2(1.3, 2, -1, pz2=0.4) == pytest.approx(0.4 + 1.3 * 5 + 2.6)
    assert cb.charged_ghost_mass2(1.3, 0) == pytest.approx(1.3)
    gB = sp.Symbol("gB", positive=True)
    assert cb.charged_vector_mass2(gB, 0, +1) == -gB
    assert cb.charged_ghost_mass2(gB, 0) == gB


def test_mutation_spin_term_breaks_tachyon():
    """Mutation: dropping the spin term must destroy the tachyon level."""
    gB = 1.0
    op = cb.landau_orbital_operator(gB, L=14.0, npoints=400)  # no spin
    w = np.linalg.eigvalsh(op)[0]
    assert w > 0.9 * gB  # without spin the lowest level is +gB, not -gB


# ---------------------------------------------------------------------------
# Exact heat-kernel coefficients
# ---------------------------------------------------------------------------


def test_heat_kernel_b2_is_11_over_3():
    assert cb.heat_kernel_series_coefficients(4)[2] == sp.Rational(11, 3)


def test_heat_kernel_b4_exact():
    assert cb.heat_kernel_series_coefficients(4)[4] == sp.Rational(127, 180)


def test_su2_beta_background_bridge_closes_exactly():
    bridge = cb.su2_chromomagnetic_beta_bridge()

    assert bridge.adjoint_casimir == 2
    assert bridge.beta_one_loop_gauge == sp.Rational(-22, 3)
    assert bridge.background_log_coefficient == 11 / (48 * sp.pi**2)
    assert bridge.beta_implied_log_coefficient == 11 / (48 * sp.pi**2)
    assert bridge.tree_scale_derivative_coefficient == 11 / (24 * sp.pi**2)
    assert bridge.loop_scale_derivative_coefficient == -11 / (24 * sp.pi**2)
    assert bridge.rg_scale_residual == 0
    assert "mu*dg_a/dmu" in bridge.beta_convention
    assert "fixed b_field" in bridge.background_potential_convention


def test_su2_beta_background_bridge_is_mutation_sensitive():
    bridge = cb.su2_chromomagnetic_beta_bridge()

    # Independent adjoint construction from epsilon_abc, rather than the
    # fundamental traces used by the bridge implementation.
    adjoint = tuple(
        sp.Matrix(
            3,
            3,
            lambda row, column: -sp.I * sp.LeviCivita(a, row, column),
        )
        for a in range(3)
    )
    casimir = sum((generator**2 for generator in adjoint), sp.zeros(3))
    assert casimir == 2 * sp.eye(3)
    assert casimir[0, 0] == bridge.adjoint_casimir

    wrong_casimir_beta = product_gauge_coefficients(
        (GaugeFactor("mutated", sp.Integer(3)),)
    ).one_loop_gauge[0]
    assert -wrong_casimir_beta / (32 * sp.pi**2) != (
        bridge.background_log_coefficient
    )
    assert (
        -bridge.beta_one_loop_gauge / (16 * sp.pi**2)
        - bridge.background_log_coefficient
    ) != 0  # mutation: log(b_field/mu), rather than log(b_field/mu**2)
    assert (
        bridge.beta_one_loop_gauge / (16 * sp.pi**2)
        + bridge.loop_scale_derivative_coefficient
    ) != 0  # mutation: reverse the beta-function sign


def test_ghost_cancellation_leaves_pure_spin():
    s, gB = sp.symbols("s gB", positive=True)
    net = cb.net_propertime_integrand(s, gB)
    # net = gB/(4 pi sinh) * 2 cosh / (4 pi s): no additive constant term
    series = sp.series(net, gB, 0, 2).removeO()
    assert sp.simplify(series - 2 / (16 * sp.pi**2 * s**2)) == 0


# ---------------------------------------------------------------------------
# Route agreement: closed form vs explicit Landau sum
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("gB,m2", [(0.5, 4.0), (1.0, 6.0)])
def test_route_a_matches_route_b(gB, m2):
    eps = 1e-3
    va = cb.one_loop_potential_cutoff(gB, m2, eps)
    vb = cb.one_loop_potential_modesum_cutoff(gB, m2, eps)
    assert vb == pytest.approx(va, rel=2e-2)


def test_mutation_ghost_weight_changes_route_b():
    """A wrong ghost weight must move route B away from route A."""
    gB, m2, eps = 1.0, 6.0, 1e-3
    va = cb.one_loop_potential_cutoff(gB, m2, eps)
    # mutated route B: ghost weight -1 instead of -2 (leaves one s3=0 mode in)
    import substrate_framework.chromomagnetic_background as m
    n = np.arange(400)

    def mutated(sv):
        x = gB * sv
        e_up = np.exp(-sv * (gB * (2 * n + 1) - 2 * gB))
        e_dn = np.exp(-sv * (gB * (2 * n + 1) + 2 * gB))
        e_0 = np.exp(-sv * gB * (2 * n + 1))
        ssum = (gB / (2 * np.pi)) * float(np.sum(e_up + e_dn + 2 * e_0 - 1.0 * e_0))
        k0 = 2.0 / (16 * np.pi**2 * sv * sv)
        return -np.exp(-sv * m2) * (ssum / (4 * np.pi * sv) - k0) / sv

    from scipy.integrate import quad

    vb_mut, _ = quad(mutated, eps, np.inf, limit=2000)
    assert abs(vb_mut - va) > 0.05 * abs(va) + 1e-4


# ---------------------------------------------------------------------------
# MS-bar extraction: log coefficient and constant from the cutoff integral
# ---------------------------------------------------------------------------


def _extract_msbar(gB: float, ratio: float = 8.0):
    """MS-bar remainder with the IR mass tracking gB: m2 = ratio * gB.

    A fixed IR mass makes V analytic in gB at weak field; the ln(gB)
    non-analyticity only forms when m2 -> 0 with gB.  With m2 = ratio*gB the
    cutoff integral is V(eps) = C0 gB^2 (gamma_E + ln(eps*ratio*gB)) + finite,
    so rem = V - C0 gB^2 (ln eps + gamma_E) = C0 gB^2 ln(gB) + const*gB^2 + ...
    """
    c0 = float(cb.heat_kernel_series_coefficients(2)[2]) / (16 * np.pi**2)
    m2 = ratio * gB
    eps_vals = np.array([2e-3, 1e-3, 5e-4])
    rem = []
    for eps in eps_vals:
        v = cb.one_loop_potential_cutoff(gB, m2, eps)
        rem.append(v - c0 * gB**2 * (np.log(eps) + np.euler_gamma))
    return float(np.mean(rem)), float(np.std(rem))


def test_msbar_log_coefficient_from_cutoff_route():
    """Cutoff route A extraction reproduces the exact 11/(48 pi^2)."""
    gBs = np.array([0.4, 0.6, 0.8, 1.0, 1.2])
    vms = np.array([_extract_msbar(float(g))[0] for g in gBs])
    A = np.column_stack(
        [gBs**2 * np.log(gBs), gBs**2, gBs**4 * np.log(gBs), gBs**4]
    )
    coef, *_ = np.linalg.lstsq(A, vms, rcond=None)
    c_expected = float(cb.one_loop_log_coefficient())
    assert coef[0] == pytest.approx(c_expected, rel=3e-2)


def test_zeta_route_log_coefficient_exact():
    """Route C (operator zeta) reproduces 11/(48 pi^2) to high precision."""
    xs = np.array([0.4, 0.6, 0.8, 1.0, 1.3])
    re = np.array([cb.one_loop_potential_zeta(float(x)).real for x in xs])
    A = np.column_stack([xs**2 * np.log(xs), xs**2, xs**4 * np.log(xs), xs**4])
    coef, *_ = np.linalg.lstsq(A, re, rcond=None)
    assert coef[0] == pytest.approx(float(cb.one_loop_log_coefficient()), rel=1e-10)
    # the zeta-scheme constant is a distinct, documented scheme value
    assert coef[1] / coef[0] == pytest.approx(-0.9455673306104181, rel=1e-8)


def test_zeta_route_imaginary_part():
    for x in (0.5, 1.0, 1.7):
        im = cb.one_loop_potential_zeta(x).imag
        assert abs(im) == pytest.approx(x**2 / (8 * np.pi), rel=1e-8)


def test_mutation_wrong_spin_coupling_breaks_extraction():
    """A wrong spin coupling (2 cosh(gB s) instead of 2 cosh(2 gB s)) must move
    the extracted log coefficient away from 11/(48 pi^2)."""
    from scipy.integrate import quad

    def v_mutated(gB, m2, eps):
        def f(sv):
            x = gB * sv
            if x < 0.02:
                # mutated series: gB*2cosh(x)/sinh(x) = 2/s + (2/3) gB^2 s + ...
                kbmk0 = (2.0 / 3) * gB**2 * sv / (16 * np.pi**2 * sv)
            else:
                kb = 2 * gB * (1 + np.exp(-2 * x)) / (1 - np.exp(-2 * x)) / (16 * np.pi**2 * sv)  # mutated spin
                k0 = 2.0 / (16 * np.pi**2 * sv * sv)
                kbmk0 = kb - k0
            return -np.exp(-sv * m2) * kbmk0 / sv

        val, _ = quad(f, eps, np.inf, limit=2000, epsabs=1e-13, epsrel=1e-11)
        return val

    gBs = np.array([0.4, 0.6, 0.8, 1.0, 1.2])
    vms = []
    for g in gBs:
        v = v_mutated(float(g), 8.0 * float(g), 1e-3)
        c0 = float(cb.heat_kernel_series_coefficients(2)[2]) / (16 * np.pi**2)
        vms.append(v - c0 * g**2 * (np.log(1e-3) + np.euler_gamma))
    A = np.column_stack([gBs**2 * np.log(gBs), gBs**2])
    coef, *_ = np.linalg.lstsq(A, np.array(vms), rcond=None)
    assert abs(coef[0] - float(cb.one_loop_log_coefficient())) > 0.2 * float(
        cb.one_loop_log_coefficient()
    )


# ---------------------------------------------------------------------------
# Savvidy minimum
# ---------------------------------------------------------------------------


def test_savvidy_minimum_equation():
    g, mu2 = sp.symbols("g mu2", positive=True)
    ln_ratio, v_min = cb.savvidy_minimum(g, mu2)
    assert sp.simplify(ln_ratio - (-24 * sp.pi**2 / (11 * g**2))) == 0
    # V_min < 0: the chromomagnetic minimum lies below the perturbative vacuum
    assert sp.simplify(v_min).could_extract_minus_sign()


def test_minimum_at_lambda_scale_is_exactly_one():
    assert cb.minimum_at_lambda_scale(sp.symbols("g", positive=True)) == sp.Integer(1)


def test_total_potential_stationarity_numeric():
    """The symbolic minimum solves dV/dB = 0 numerically."""
    gval, mu2val = 1.0, 1.0
    B, g, mu2 = sp.symbols("B g mu2", positive=True)
    v = cb.total_potential(B, g, mu2).subs({g: gval, mu2: mu2val})
    dv = sp.diff(v, B)
    ln_ratio = float(cb.savvidy_minimum(gval, mu2val)[0])
    b_min = np.exp(ln_ratio) / gval
    assert float(dv.subs(B, b_min)) == pytest.approx(0.0, abs=1e-10)


# ---------------------------------------------------------------------------
# Literature anchors (Savvidy's open restatements, extracted verbatim in P229)
# ---------------------------------------------------------------------------


def test_savvidy_vmin_matches_literature():
    """V_min = -(11 mu^4 / 96 pi^2) exp(-48 pi^2 / (11 g^2)) (Savvidy EPJC 2020)."""
    g, mu2 = sp.symbols("g mu2", positive=True)
    _, v_min = cb.savvidy_minimum(g, mu2)
    expected = -11 * mu2**2 / (96 * sp.pi**2) * sp.exp(-48 * sp.pi**2 / (11 * g**2))
    assert sp.simplify(v_min - expected) == 0


def test_imaginary_part_is_nielsen_olesen():
    gB = sp.symbols("gB", positive=True)
    assert sp.simplify(cb.nielsen_olesen_imaginary_part(gB) - gB**2 / (8 * sp.pi)) == 0


def test_mutation_tachyon_factor_changes_imaginary_part():
    """Dropping the real-dof factor 2 must change the result."""
    gB = sp.symbols("gB", positive=True)
    pz = sp.symbols("pz", real=True)
    a = sp.sqrt(gB)
    integral = sp.integrate(sp.sqrt(a**2 - pz**2) / (2 * sp.pi), (pz, -a, a))
    mutated = sp.simplify((gB / (2 * sp.pi)) * sp.Rational(1, 2) * integral)  # no factor 2
    assert sp.simplify(mutated - gB**2 / (8 * sp.pi)) != 0

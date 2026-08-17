"""Exact-mass one-loop proper-time coefficient tests (P230, advances #76).

Oracles, per ``.agents/skills/physics-erdos-loop/references/oracles.md``:

* exact identities are verified symbolically (derivative/tail uniqueness for
  the sharp cutoff scheme, the integration-by-parts identity for its vacuum
  class, the differential recurrence plus boundary limits for the smooth
  Bessel family, and cutoff-subtraction limits for the power-subtracted
  scheme);
* the special-function closed forms are independently corroborated by
  high-precision mpmath quadrature;
* load-bearing mutations (wrong prefactor, wrong Bessel order, wrong
  exponential-integral branch, wrong xi sign, unknown regulator, defaulted
  scales) must each break a relevant check.
"""

import mpmath as mp
import pytest
import sympy as sp

import substrate_framework as framework
from substrate_framework.scalar_induced_newton import (
    SHARP_PROPER_TIME_REGULATOR,
    leading_scalar_newton_shift_coefficient,
)
from substrate_framework.scalar_one_loop_mass import (
    KNOWN_ONE_LOOP_REGULATORS,
    SMOOTH_PROPER_TIME_REGULATOR,
    ZETA_POWER_SUBTRACTED_REGULATOR,
    curvature_proper_time_integral,
    exact_mass_inverse_newton_shift,
    exact_mass_vacuum_density_shift,
    regulator_scheme_ledger,
    vacuum_proper_time_integral,
)

Lam, mu, m2 = sp.symbols("Lambda mu m2", positive=True)
tau0, t = sp.symbols("tau0 t", positive=True)
z = sp.symbols("z", positive=True)


def test_one_loop_mass_api_is_exported_from_package() -> None:
    assert framework.SMOOTH_PROPER_TIME_REGULATOR == SMOOTH_PROPER_TIME_REGULATOR
    assert framework.ZETA_POWER_SUBTRACTED_REGULATOR == ZETA_POWER_SUBTRACTED_REGULATOR
    assert framework.curvature_proper_time_integral is curvature_proper_time_integral
    assert framework.vacuum_proper_time_integral is vacuum_proper_time_integral
    assert framework.exact_mass_inverse_newton_shift is exact_mass_inverse_newton_shift
    assert framework.exact_mass_vacuum_density_shift is exact_mass_vacuum_density_shift
    assert framework.regulator_scheme_ledger is regulator_scheme_ledger


def test_sharp_curvature_closed_form_is_the_tail_integral_exactly() -> None:
    # I(tau0) = Lam^2*(exp(-z)-z*E1(z)) with z = m2*tau0, Lam^2 = 1/tau0.
    closed = (sp.exp(-m2 * tau0) - m2 * tau0 * sp.expint(1, m2 * tau0)) / tau0
    derivative_residual = sp.simplify(
        sp.diff(closed, tau0) + tau0 ** -2 * sp.exp(-m2 * tau0)
    )
    assert derivative_residual == 0
    # Tail decay pins the integration constant.  E1(x) is by definition the
    # tail integral integral_x^infty exp(-u)/u du; for u >= x > 0 the
    # integrand is bounded by exp(-u)/x, whose tail integral sympy evaluates
    # exactly, so 0 <= E1(x) <= exp(-x)/x -> 0.  The explicit exponential
    # piece vanishes directly.
    assert sp.limit(sp.exp(-m2 * tau0) / tau0, tau0, sp.oo) == 0
    x, u = sp.symbols("x u", positive=True)
    assert sp.integrate(sp.exp(-u) / x, (u, x, sp.oo)) == sp.exp(-x) / x
    mp.mp.dps = 30
    for xf in (mp.mpf(30), mp.mpf(40)):
        assert mp.e1(xf) < mp.e ** (-xf) / xf
        assert mp.e1(xf) * xf * mp.e ** xf < 1
    # The module's expression is this closed form.
    value = curvature_proper_time_integral(
        SHARP_PROPER_TIME_REGULATOR, cutoff=Lam, mass_squared=m2
    )
    assert sp.simplify(value - Lam**2 * (sp.exp(-m2 / Lam**2) - (m2 / Lam**2) * sp.expint(1, m2 / Lam**2))) == 0


def test_sharp_vacuum_closed_form_satisfies_the_by_parts_identity() -> None:
    sharp2 = curvature_proper_time_integral(
        SHARP_PROPER_TIME_REGULATOR, cutoff=Lam, mass_squared=m2
    )
    sharp3 = vacuum_proper_time_integral(
        SHARP_PROPER_TIME_REGULATOR, cutoff=Lam, mass_squared=m2
    )
    expected = sp.exp(-m2 / Lam**2) * Lam**4 / 2 - (m2 / 2) * sharp2
    assert sp.simplify(sharp3 - expected) == 0
    # Independent derivation by parts from the tail-integral form.
    by_parts = sp.integrate(sp.exp(-m2 * t) / t**3, (t, tau0, sp.oo))
    target = (sp.exp(-m2 * tau0) / (2 * tau0**2)) - (m2 / 2) * (
        (sp.exp(-m2 * tau0) - m2 * tau0 * sp.expint(1, m2 * tau0)) / tau0
    )
    # sympy may leave Ei-branch structure; compare via the derivative oracle.
    residual = sp.diff(by_parts - target, tau0).replace(
        sp.exp_polar, lambda arg: sp.exp(arg)
    )
    assert sp.simplify(residual) == 0


def test_smooth_closed_forms_satisfy_the_differential_recurrence() -> None:
    # J_p(z) = integral_0^infty t^-p exp(-z t - 1/t) dt obeys dJ_p/dz = -J_{p-1}.
    J1 = 2 * sp.besselk(0, 2 * sp.sqrt(z))
    J2 = 2 * sp.sqrt(z) * sp.besselk(1, 2 * sp.sqrt(z))
    J3 = 2 * z * sp.besselk(2, 2 * sp.sqrt(z))
    J0 = 2 * z ** sp.Rational(-1, 2) * sp.besselk(1, 2 * sp.sqrt(z))
    assert sp.simplify(sp.diff(J2, z) + J1) == 0
    assert sp.simplify(sp.diff(J3, z) + J2) == 0
    assert sp.simplify(sp.diff(J1, z) + J0) == 0
    # Boundary values pin the solutions of the recurrence.
    assert sp.limit(J2, z, 0, "+") == 1
    assert sp.limit(J3, z, 0, "+") == 1
    # The module returns these forms at z = m2/Lam^2.
    value2 = curvature_proper_time_integral(
        SMOOTH_PROPER_TIME_REGULATOR, cutoff=Lam, mass_squared=m2
    )
    value3 = vacuum_proper_time_integral(
        SMOOTH_PROPER_TIME_REGULATOR, cutoff=Lam, mass_squared=m2
    )
    assert sp.simplify(value2 - Lam**2 * J2.subs(z, m2 / Lam**2)) == 0
    assert sp.simplify(value3 - Lam**4 * J3.subs(z, m2 / Lam**2)) == 0


@pytest.mark.parametrize("zf_str", ["0.1", "1.0", "2.5"])
def test_smooth_closed_forms_match_high_precision_quadrature(zf_str) -> None:
    mp.mp.dps = 40
    zf = mp.mpf(zf_str)
    num2 = mp.quad(lambda tt: mp.e ** (-zf * tt - 1 / tt) / tt**2, [0, 1, mp.inf])
    closed2 = 2 * mp.sqrt(zf) * mp.besselk(1, 2 * mp.sqrt(zf))
    assert abs(num2 - closed2) < mp.mpf("1e-30") * abs(closed2)
    num3 = mp.quad(lambda tt: mp.e ** (-zf * tt - 1 / tt) / tt**3, [0, 1, mp.inf])
    closed3 = 2 * zf * mp.besselk(2, 2 * mp.sqrt(zf))
    assert abs(num3 - closed3) < mp.mpf("1e-30") * abs(closed3)


@pytest.mark.parametrize("zf_str", ["0.1", "1.0", "2.5"])
def test_sharp_closed_form_matches_quadrature(zf_str) -> None:
    mp.mp.dps = 40
    zf = mp.mpf(zf_str)
    num = mp.quad(lambda tt: mp.e ** (-zf * tt) / tt**2, [1.0, mp.inf])
    closed = mp.e ** (-zf) - zf * mp.e1(zf)
    assert abs(num - closed) < mp.mpf("1e-30") * abs(closed)


def test_zeta_closed_forms_are_the_power_subtraction_limits() -> None:
    sharp2 = (sp.exp(-m2 * tau0) - m2 * tau0 * sp.expint(1, m2 * tau0)) / tau0
    lim2 = sp.limit(sharp2 - 1 / tau0 - m2 * sp.log(tau0), tau0, 0, "+")
    zeta2 = curvature_proper_time_integral(
        ZETA_POWER_SUBTRACTED_REGULATOR, mass_squared=m2, renormalization_scale=mu
    )
    assert sp.simplify(lim2 - m2 * (sp.log(m2) + sp.EulerGamma - 1)) == 0
    assert sp.simplify(zeta2 - m2 * (sp.log(m2 / mu**2) + sp.EulerGamma - 1)) == 0

    sharp3 = sp.exp(-m2 * tau0) / (2 * tau0**2) - (m2 / 2) * sharp2
    lim3 = sp.limit(
        sharp3 - 1 / (2 * tau0**2) + m2 / tau0 + (m2**2 / 2) * sp.log(tau0),
        tau0,
        0,
        "+",
    )
    zeta3 = vacuum_proper_time_integral(
        ZETA_POWER_SUBTRACTED_REGULATOR, mass_squared=m2, renormalization_scale=mu
    )
    expected3 = -(m2**2 / 2) * (sp.log(m2) + sp.EulerGamma - sp.Rational(3, 2))
    assert sp.simplify(lim3 - expected3) == 0
    assert sp.simplify(zeta3 + (m2**2 / 2) * (sp.log(m2 / mu**2) + sp.EulerGamma - sp.Rational(3, 2))) == 0


def test_massless_sharp_limit_reproduces_the_accepted_leading_shift() -> None:
    N, xi = sp.Integer(3), sp.Rational(1, 12)
    shift = exact_mass_inverse_newton_shift(
        N,
        xi,
        regulator=SHARP_PROPER_TIME_REGULATOR,
        cutoff=Lam,
        mass_squared=0,
    )
    accepted = leading_scalar_newton_shift_coefficient(
        N, xi, regulator=SHARP_PROPER_TIME_REGULATOR
    )
    assert sp.simplify(shift.value - accepted.coefficient * Lam**2) == 0
    assert sp.simplify(shift.finite_mass_factor - 1) == 0

    smooth = exact_mass_inverse_newton_shift(
        N,
        xi,
        regulator=SMOOTH_PROPER_TIME_REGULATOR,
        cutoff=Lam,
        mass_squared=0,
    )
    # The smooth weight also integrates to exactly Lambda^2 when massless.
    assert sp.simplify(smooth.value - accepted.coefficient * Lam**2) == 0

    zeta = exact_mass_inverse_newton_shift(
        N,
        xi,
        regulator=ZETA_POWER_SUBTRACTED_REGULATOR,
        mass_squared=0,
        renormalization_scale=mu,
    )
    # Power divergences subtracted: no induced shift survives at m = 0.
    assert sp.simplify(zeta.value) == 0


def test_small_z_series_recovers_the_accepted_logarithmic_structure() -> None:
    F = sp.exp(-z) - z * sp.expint(1, z)
    series = sp.series(F, z, 0, 2).removeO()
    assert sp.simplify(series - (1 + z * (sp.log(z) + sp.EulerGamma - 1))) == 0


def test_scheme_ledger_reports_exact_contrasts() -> None:
    ledger = regulator_scheme_ledger(Lam, m2, mu)
    assert sp.simplify(
        ledger.sharp_over_smooth
        - (
            (sp.exp(-m2 / Lam**2) - (m2 / Lam**2) * sp.expint(1, m2 / Lam**2))
            / (2 * sp.sqrt(m2 / Lam**2) * sp.besselk(1, 2 * sp.sqrt(m2 / Lam**2)))
        )
    ) == 0
    # At z = 1 the smooth scheme induces about 1.88 times the sharp value.
    one = regulator_scheme_ledger(1, 1, 1)
    ratio = float(sp.N(one.sharp_over_smooth, 12))
    assert abs(ratio - 0.5308496410) < 1e-8
    # The PR #13 review's load-bearing factor is reproduced at z = 1.
    factor = sp.N(
        sp.exp(-1) - sp.expint(1, 1),
        10,
    )
    assert abs(factor - sp.Float("0.1484955068", 10)) < sp.Float("1e-9")


def test_exact_mass_shift_composes_the_accepted_scheme_factor() -> None:
    N, xi = sp.Integer(2), sp.Rational(1, 12)
    shift = exact_mass_inverse_newton_shift(
        N,
        xi,
        regulator=SHARP_PROPER_TIME_REGULATOR,
        cutoff=Lam,
        mass_squared=m2,
    )
    accepted = leading_scalar_newton_shift_coefficient(
        1, xi, regulator=SHARP_PROPER_TIME_REGULATOR
    )
    assert sp.simplify(shift.coefficient_per_field - accepted.coefficient_per_field) == 0
    assert sp.simplify(shift.value - N * accepted.coefficient_per_field * shift.proper_time_value) == 0


def test_vacuum_shift_composes_both_integral_classes() -> None:
    N = sp.Integer(4)
    for regulator, kwargs in (
        (SHARP_PROPER_TIME_REGULATOR, {"cutoff": Lam, "mass_squared": m2}),
        (SMOOTH_PROPER_TIME_REGULATOR, {"cutoff": Lam, "mass_squared": m2}),
        (
            ZETA_POWER_SUBTRACTED_REGULATOR,
            {"mass_squared": m2, "renormalization_scale": mu},
        ),
    ):
        vac = exact_mass_vacuum_density_shift(N, regulator=regulator, **kwargs)
        expected = -N * sp.Rational(1, 2) * (4 * sp.pi) ** -2 * (
            vac.tau_minus_three_value + m2 * vac.tau_minus_two_value
        )
        assert sp.simplify(vac.value - expected) == 0
    # Massless sharp vacuum sector: I_3 = Lambda^4/2 and m^2 I_2 = 0.
    vac0 = exact_mass_vacuum_density_shift(
        1, regulator=SHARP_PROPER_TIME_REGULATOR, cutoff=Lam, mass_squared=0
    )
    assert sp.simplify(vac0.value + Lam**4 / (4 * (4 * sp.pi) ** 2)) == 0


def test_mutations_break_the_oracles() -> None:
    # Wrong prefactor: doubling the sharp closed form fails the tail identity.
    doubled = 2 * (sp.exp(-m2 * tau0) - m2 * tau0 * sp.expint(1, m2 * tau0)) / tau0
    residual = sp.simplify(sp.diff(doubled, tau0) + tau0 ** -2 * sp.exp(-m2 * tau0))
    assert residual != 0
    # Wrong exponential-integral branch: Ei instead of E1 fails it too.
    wrong_branch = (sp.exp(-m2 * tau0) - m2 * tau0 * sp.Ei(m2 * tau0)) / tau0
    residual_branch = sp.simplify(
        sp.diff(wrong_branch, tau0) + tau0 ** -2 * sp.exp(-m2 * tau0)
    )
    assert residual_branch != 0
    # Wrong Bessel order: K_2 in the curvature class breaks the recurrence.
    wrong_order = 2 * sp.sqrt(z) * sp.besselk(2, 2 * sp.sqrt(z))
    J1 = 2 * sp.besselk(0, 2 * sp.sqrt(z))
    assert sp.simplify(sp.diff(wrong_order, z) + J1) != 0
    # Crossing the conformal value 1/6 flips the induced shift sign: the
    # minimal coupling xi = 0 has 1/6 - xi > 0 while xi = 1 flips it.
    plus = exact_mass_inverse_newton_shift(
        1, 0, regulator=SHARP_PROPER_TIME_REGULATOR, cutoff=Lam
    )
    minus = exact_mass_inverse_newton_shift(
        1, 1, regulator=SHARP_PROPER_TIME_REGULATOR, cutoff=Lam
    )
    assert plus.sign == 1 and minus.sign == -1
    assert sp.simplify(plus.value + minus.value) != 0


def test_input_contracts_reject_ambiguous_schemes() -> None:
    with pytest.raises(ValueError):
        curvature_proper_time_integral("wrong_name", cutoff=Lam)
    with pytest.raises(ValueError):
        curvature_proper_time_integral(SMOOTH_PROPER_TIME_REGULATOR)
    with pytest.raises(ValueError):
        curvature_proper_time_integral(
            SMOOTH_PROPER_TIME_REGULATOR,
            cutoff=Lam,
            renormalization_scale=mu,
        )
    with pytest.raises(ValueError):
        curvature_proper_time_integral(ZETA_POWER_SUBTRACTED_REGULATOR)
    with pytest.raises(ValueError):
        curvature_proper_time_integral(
            ZETA_POWER_SUBTRACTED_REGULATOR,
            cutoff=Lam,
            renormalization_scale=mu,
        )
    with pytest.raises(ValueError):
        curvature_proper_time_integral(
            SHARP_PROPER_TIME_REGULATOR,
            cutoff=Lam,
            mass_squared=sp.Symbol("neg", negative=True),
        )
    with pytest.raises(ValueError):
        exact_mass_inverse_newton_shift(
            sp.Rational(3, 2),
            0,
            regulator=SHARP_PROPER_TIME_REGULATOR,
            cutoff=Lam,
        )
    with pytest.raises(ValueError):
        exact_mass_inverse_newton_shift(
            1,
            sp.Symbol("xi", real=True),
            regulator=SHARP_PROPER_TIME_REGULATOR,
            cutoff=Lam,
        )


def test_module_is_target_blind() -> None:
    import inspect

    from substrate_framework import scalar_one_loop_mass

    source = inspect.getsource(scalar_one_loop_mass)
    for forbidden in ("6.674", "Planck", "M_pl", "6.708", "observed"):
        assert forbidden not in source


def test_known_regulators_are_the_three_preregistered_schemes() -> None:
    assert set(KNOWN_ONE_LOOP_REGULATORS) == {
        SHARP_PROPER_TIME_REGULATOR,
        SMOOTH_PROPER_TIME_REGULATOR,
        ZETA_POWER_SUBTRACTED_REGULATOR,
    }

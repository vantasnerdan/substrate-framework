"""Machine checks for src/substrate_framework/chromomagnetic_two_loop.py (P229, issue #46)."""

import numpy as np
import sympy as sp
import pytest

from substrate_framework import chromomagnetic_two_loop as c2
from substrate_framework.chromomagnetic_two_loop import b, g, mu2


# ---------------------------------------------------------------------------
# Level-sum identities
# ---------------------------------------------------------------------------


def test_zeta_zero_sum_cancels_divergence():
    """G(1) = 1 + zeta(0,1/2) + zeta(0,3/2) = 0: pole, ln b, ln mu all drop."""
    assert c2.zeta_zero_sum() == sp.Integer(-1)


def test_lerch_identity_numerically():
    """zeta'(0,a) = ln Gamma(a) - ln sqrt(2 pi), checked by mpmath differentiation."""
    import mpmath as mp

    with mp.workdps(40):
        for a in (mp.mpf("0.5"), mp.mpf("1.5"), mp.mpf("0.7")):
            num = mp.diff(lambda s: mp.zeta(s, a), 0)
            closed = mp.log(mp.gamma(a)) - mp.log(2 * mp.pi) / 2
            assert abs(num - closed) < mp.mpf("1e-30")


def test_log_level_sum_exact():
    assert sp.simplify(c2.log_level_sum() - (sp.log(2) - sp.I * sp.pi)) == 0


def test_log_level_sum_mpmath_independent():
    """mpmath differentiation of G(s) at s=1 vs the SymPy exact value."""
    val = c2.B2_tadpole_mpmath()
    expected = complex(-float(sp.log(2)) / (16 * np.pi**2), float(sp.pi) / (16 * np.pi**2))
    assert val == pytest.approx(expected, rel=1e-12)


def test_mutation_tachyon_branch_changes_level_sum():
    """Dropping the tachyon's -i pi must change the level sum."""
    assert sp.simplify(c2.log_level_sum() - sp.log(2)) != 0


# ---------------------------------------------------------------------------
# Two-loop term: derived vs printed
# ---------------------------------------------------------------------------


def test_re_W2_derived_exact_form():
    expected = 3 * g**2 * b**2 * (sp.log(2) ** 2 - sp.pi**2) / (512 * sp.pi**4)
    assert sp.simplify(c2.re_W2_derived() - expected) == 0


def test_re_W2_ratio_printed_over_derived_ln2_part():
    """The printed coefficient equals (4/3) x the derived ln^2 2 part, with no
    -pi^2 structure: documents finding F3 exactly."""
    ratio = sp.simplify(c2.re_W2_printed() / (3 * g**2 * b**2 * sp.log(2) ** 2 / (512 * sp.pi**4)))
    assert ratio == sp.Rational(4, 3)
    # the exact difference: printed lacks both the factor 4/3 and the -pi^2 part
    diff = sp.simplify(c2.re_W2_derived() - c2.re_W2_printed())
    assert sp.simplify(diff - (-(g**2 * b**2 * (sp.log(2) ** 2 + 3 * sp.pi**2)) / (512 * sp.pi**4))) == 0


# ---------------------------------------------------------------------------
# Minimum: F1 (the 98 vs 88 typo) and the corrected exponent
# ---------------------------------------------------------------------------


def test_printed_minimum_exponent_is_88_not_98():
    """Minimizing the printed (57) gives 3 ln^2(2) g^2/(88 pi^2): F1 confirmed,
    the exponential form printing 98 is a typo (their expansion prints 88)."""
    corr = c2.minimum_exponent(c2.W_two_loop_printed())
    expected = -(3 * sp.log(2) ** 2 * g**2) / (88 * sp.pi**2)
    assert sp.simplify(corr - expected) == 0


def test_derived_minimum_exponent():
    """With the derived W2 the exponent correction is -9(ln^2 2 - pi^2) g^2/(352 pi^2)."""
    corr = c2.minimum_exponent(c2.W_two_loop_derived())
    expected = -9 * (sp.log(2) ** 2 - sp.pi**2) * g**2 / (352 * sp.pi**2)
    assert sp.simplify(corr - expected) == 0


def test_minimum_solves_stationarity_numerically():
    """The perturbative minimum solves dW/db = 0 to the working order."""
    import scipy.optimize as so

    gval, mu2val = 0.5, 1.0
    w = sp.re(sp.expand_complex(c2.W_two_loop_derived())).subs({g: gval, mu2: mu2val})
    dwf = sp.lambdify(b, sp.diff(w, b), "numpy")
    x0 = -24 * np.pi**2 / (11 * gval**2)
    corr = float(c2.minimum_exponent(c2.W_two_loop_derived()).subs(g, gval))
    b_min = mu2val * np.exp(x0 + corr)
    wmin = float(w.subs(b, b_min))
    assert abs(dwf(b_min)) / abs(wmin) < 1e-2  # residual is O(g^4) of the potential


def test_mutation_wrong_combinatoric_factor_moves_minimum():
    """Using (g^2/2) B2^2 instead of (3 g^2/2) B2^2 must change the exponent."""
    w_mut = (
        b**2 / (2 * g**2)
        + sp.Rational(11) * b**2 / (48 * sp.pi**2) * (sp.log(b / mu2) - sp.Rational(1, 2))
        + sp.Rational(1, 2) * g**2 * c2.B2_tadpole() ** 2
    )
    corr_mut = c2.minimum_exponent(w_mut)
    corr_true = c2.minimum_exponent(c2.W_two_loop_derived())
    assert sp.simplify(corr_mut - corr_true) != 0


# ---------------------------------------------------------------------------
# RG-improved potential in Lambda units
# ---------------------------------------------------------------------------


def test_rg_improved_one_loop_minimum_at_x1():
    """With one-loop running, the Savvidy minimum sits at b = Lambda^2 exactly."""
    x = sp.symbols("x", positive=True)
    w = c2.rg_improved_potential("none")
    dw = sp.diff(w, x)
    assert sp.simplify(dw.subs(x, 1)) == 0
    assert sp.simplify(w.subs(x, 1) - (-sp.Rational(11) / (96 * sp.pi**2))) == 0


def test_rg_improved_two_loop_coefficients():
    """The x^2/ln x coefficients for the printed and derived two-loop terms."""
    x = sp.symbols("x", positive=True)
    w_p = c2.rg_improved_potential("printed")
    w_d = c2.rg_improved_potential("derived")
    coef_p = sp.simplify((w_p - c2.rg_improved_potential("none")) * sp.log(x) / x**2)
    coef_d = sp.simplify((w_d - c2.rg_improved_potential("none")) * sp.log(x) / x**2)
    assert sp.simplify(coef_p - 3 * sp.log(2) ** 2 / (176 * sp.pi**2)) == 0
    assert sp.simplify(coef_d - 9 * (sp.log(2) ** 2 - sp.pi**2) / (704 * sp.pi**2)) == 0


def test_rg_improved_two_loop_expression_is_singular_at_one_loop_minimum():
    """The running-coupling substitution cannot define a corrected x=1 minimum."""
    x = sp.symbols("x", positive=True)
    derived = c2.rg_improved_potential("derived")
    printed = c2.rg_improved_potential("printed")

    assert sp.limit(derived, x, 1, dir="-") == sp.oo
    assert sp.limit(derived, x, 1, dir="+") == -sp.oo
    assert sp.limit(printed, x, 1, dir="-") == -sp.oo
    assert sp.limit(printed, x, 1, dir="+") == sp.oo

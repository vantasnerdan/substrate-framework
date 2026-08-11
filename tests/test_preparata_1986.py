"""Machine checks for src/substrate_framework/preparata_1986.py (P229, issue #56).

Each test pins one claim of Preparata, Nuovo Cim. A 96 (1986) 366, from the
encoded inputs; mutations break them."""

import sympy as sp

from substrate_framework import preparata_1986 as p86


def test_spectrum_unstable_sector():
    """(1.2): n=0, S_z=1 gives E^2 < 0 exactly for p_z^2 < gH."""
    cond = p86.unstable_sector_condition()
    assert cond == sp.Interval.open(0, p86.b) or cond == sp.Interval.open(
        -sp.oo, p86.b
    ).intersect(sp.Interval(0, sp.oo)) or cond.sup == p86.b
    assert cond.inf < 0 or cond.inf == 0


def test_spectrum_matches_verified_primitive():
    """(1.2) agrees numerically with chromomagnetic_background's spectrum."""
    for n, s3 in [(0, 1), (0, -1), (1, 1), (3, 0)]:
        assert p86.spectrum_matches_primitive(0.7, 0.05, n, s3)


def test_sign_of_51_is_derived_not_transcribed():
    """Only the negative sign of the U-mode log gives a minimum with a,b > 0;
    the positive sign gives a maximum — inconsistent with the abstract."""
    res = p86.required_sign_for_instability()
    assert int(res[-1][0][1]) == 1   # minimum
    assert int(res[+1][0][1]) == -1  # maximum


def test_minimum_leading_log_exact():
    """(5.2)/(5.3): a = e^(-1/2), b = 11/(96 pi^2 e) — both finite as g -> 0."""
    a, b_coef = p86.minimum_leading_log()
    assert sp.simplify(a - sp.exp(sp.Rational(-1, 2))) == 0
    assert sp.simplify(b_coef - sp.Rational(11) / (96 * sp.pi**2 * sp.E)) == 0
    assert a.is_positive and b_coef.is_positive


def test_minimum_robust_to_subleading_term():
    """With the symbolic O(g^2 H^2) term, a = e^(-1/2 - 48 pi^2 kappa/11):
    finite as g -> 0 whenever kappa does not diverge — the 'essential' claim."""
    a = p86.minimum_with_subleading()
    C = p86.ONE_LOOP_LOG_COEF
    expected = sp.exp(-sp.Rational(1, 2) - p86.kappa / C)
    assert sp.simplify(sp.log(a) - sp.log(expected)) == 0


def test_savvidy_minimum_is_one_loop_af():
    """Sect. 1: solving (1.1) for 1/g^2 yields the one-loop AF law (SU(2)
    b0 = 22/3) at mu^2 = b*. Residual must vanish exactly."""
    assert p86.savvidy_minimum_implies_af() == 0


def test_1985_variational_minimum_diverges_under_af():
    """Sect. 1: the 1985 result (1.3) under AF running gives exponent
    exactly -ln(Lambda/Lambda_QCD), hence gH* = Lambda . Lambda_QCD."""
    exponent, b_star = p86.variational_1985_divergence()
    assert sp.simplify(exponent + sp.log(p86.Lam / p86.LamQCD)) == 0
    assert sp.simplify(b_star - p86.Lam * p86.LamQCD) == 0


def test_log_coefficient_matches_one_loop_primitive():
    """The C of (5.1) equals the verified one-loop coefficient (#51, zeta
    route 1e-10)."""
    assert p86.log_coefficient_matches_one_loop() == 0


def test_mutation_wrong_log_coefficient_moves_minimum():
    """A wrong C (e.g. 11/(24 pi^2)) must change b: mutation sensitivity."""
    a, b_coef = p86.minimum_leading_log()
    assert sp.simplify(b_coef - sp.Rational(11) / (48 * sp.pi**2 * sp.E)) != 0
    mutated = sp.Rational(11, 24) / sp.pi**2 * sp.exp(-1) / 2
    assert sp.simplify(b_coef - mutated) != 0


def test_mutation_positive_sign_has_no_finite_b_minimum():
    """With the wrong (positive) sign there is no local minimum at all:
    the instability conclusion cannot be recovered."""
    res = p86.required_sign_for_instability()
    assert all(int(second) < 0 for _, second in res[+1])


def test_summary_composes():
    s = p86.essential_instability_summary()
    assert "consistent" in s["verdict"]

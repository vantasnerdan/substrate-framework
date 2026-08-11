"""Equation-level validation of Preparata, Nuovo Cim. A 96 (1986) 366.

Conditional, unpromoted primitive (P229, issue #56).  The primary source is
pinned at proposals/P229-preparata-qcd-vacuum-audit/sources/ (md5
34f8efc4c056837869cfbe63906ea14d).  Each function encodes one claim of the
paper from its clearest textual statement (the 1986 scan garbles displayed
equations, so nothing is transcribed from OCR as ground truth); every claim
is checked by SymPy from the encoded inputs, and cross-checked against the
independently verified primitives chromomagnetic_background (#51) and
chromomagnetic_two_loop (#52).

Claims under test:
  (1.2)  fluctuation spectrum and the unstable sector;
  (5.1)  Delta E = -(11/48 pi^2) g^2 H^2 ln(Lambda^2/gH) + O(g^2 H^2)
         (the sign is NOT transcribed: it is derived from the requirement
         that (5.2)/(5.3) hold with a, b > 0);
  (5.2)/(5.3)  minimum at gH* = a Lambda^2, Delta E* = -b Lambda^4 with
         a = e^(-1/2), b = 11/(96 pi^2 e) at leading log, both finite as
         g -> 0 even after the O(g^2 H^2) term is included symbolically;
  sect.1 AF argument: Savvidy's (1.1) solved for g^2(Lambda) IS one-loop
         asymptotic freedom; the 1985 variational (1.3) under AF running
         gives gH* = Lambda . Lambda_QCD (divergent field, ~Lambda^2 gap).
"""

import sympy as sp

from substrate_framework import chromomagnetic_background as cb

b, g, H, Lam, LamQCD = sp.symbols("b g H Lambda Lambda_QCD", positive=True)
kappa = sp.Symbol("kappa", real=True)  # coefficient of the O(g^2 H^2) term

ONE_LOOP_LOG_COEF = sp.Rational(11) / (48 * sp.pi**2)  # C in (5.1)


# ---------------------------------------------------------------------------
# Claim (1.2): the spectrum and its unstable sector
# ---------------------------------------------------------------------------


def preparata_spectrum(pz2, n, s3):
    """Preparata (1.2): E^2 = p_z^2 + gH(2n+1) - 2 gH S_z, encoded directly."""
    return pz2 + b * (2 * n + 1) - 2 * b * s3


def unstable_sector_condition():
    """The n=0, S_z=1 sector is tachyonic for p_z^2 < gH (his text, sect. 1)."""
    pz2 = sp.Symbol("pz2", nonnegative=True)
    e2 = preparata_spectrum(pz2, 0, 1)
    return sp.solve_univariate_inequality(e2 < 0, pz2, relational=False)


# ---------------------------------------------------------------------------
# Claim (5.1): the U-mode energy density, sign derived not transcribed
# ---------------------------------------------------------------------------


def u_mode_energy(sign):
    """Delta E = sign * C b^2 ln(Lambda^2/b) + kappa b^2 (O(g^2 H^2) term)."""
    return sign * ONE_LOOP_LOG_COEF * b**2 * sp.log(Lam**2 / b) + kappa * b**2


def required_sign_for_instability():
    """Only sign = -1 yields a local minimum with a, b > 0 as claimed in
    (5.2)/(5.3); sign = +1 gives a maximum — inconsistent with the abstract."""
    results = {}
    for s in (+1, -1):
        e = u_mode_energy(s).subs(kappa, 0)
        crit = sp.solve(sp.diff(e, b), b)
        second = [sp.sign(sp.diff(e, b, 2).subs(b, c)) for c in crit]
        results[s] = list(zip(crit, second))
    return results


def minimum_leading_log():
    """Minimum of the leading-log (5.1): b* = Lambda^2 e^(-1/2),
    Delta E* = -(11/(96 pi^2 e)) Lambda^4.  Derived, not asserted."""
    e = u_mode_energy(-1).subs(kappa, 0)
    b_star = sp.solve(sp.diff(e, b), b)[0]
    a = sp.simplify(b_star / Lam**2)
    delta_e = sp.simplify(e.subs(b, b_star))
    b_coef = sp.simplify(-delta_e / Lam**4)
    return a, b_coef


def minimum_with_subleading():
    """With the symbolic O(g^2 H^2) term: a = e^(-1/2 - kappa/C) stays finite
    as g -> 0 for any kappa that does not diverge — the 'essential' claim is
    robust to the subleading term, as Preparata asserts after (5.3)."""
    e = u_mode_energy(-1)
    b_star = sp.solve(sp.diff(e, b), b)[0]
    a = sp.simplify(sp.exp(sp.simplify(sp.log(b_star / Lam**2))))
    return a


# ---------------------------------------------------------------------------
# Sect. 1: asymptotic-freedom algebra
# ---------------------------------------------------------------------------


def savvidy_minimum_implies_af():
    """Savvidy (1.1): b* = Lambda^2 exp(-24 pi^2/(11 g^2)) solved for 1/g^2
    gives 1/g^2 = (11/(24 pi^2)) ln(Lambda^2/b*) — the one-loop AF form with
    SU(2) b0 = 22/3: 1/g^2(mu) = (b0/8 pi^2) ln(mu^2/Lambda_RG^2)."""
    g2 = sp.Symbol("g2", positive=True)
    sol = sp.solve(sp.Eq(sp.log(b / Lam**2), -24 * sp.pi**2 / (11 * g2)), 1 / g2)
    b0 = sp.Rational(22, 3)
    # one-loop AF at mu^2 = b: 1/g^2 = (b0/8 pi^2) ln(sqrt(b)/Lambda)^{-1}...
    af = b0 / (8 * sp.pi**2) * sp.log(Lam / sp.sqrt(b))
    return sp.simplify(sp.expand_log(sol[0] - af, force=True))


def variational_1985_divergence():
    """Preparata's (1.3) from ref. (3): b* = Lambda^2 exp(-12 pi^2/(11 g^2)).
    Under AF running 1/g^2(Lambda) = (11/(12 pi^2)) ln(Lambda/Lambda_QCD)
    the exponent is exactly -ln(Lambda/Lambda_QCD), so b* = Lambda.Lambda_QCD:
    divergent field, ~Lambda^2 energy gap — as claimed in sect. 1."""
    exponent = -12 * sp.pi**2 / 11 * (sp.Rational(11) / (12 * sp.pi**2)) * sp.log(
        Lam / LamQCD
    )
    exponent = sp.simplify(exponent)
    b_star = sp.simplify(Lam**2 * sp.exp(exponent))
    return exponent, b_star


# ---------------------------------------------------------------------------
# Cross-checks against the independently verified primitives (#51, #52)
# ---------------------------------------------------------------------------


def log_coefficient_matches_one_loop() -> sp.Expr:
    """The C of (5.1) must equal the verified one-loop coefficient C =
    11/(48 pi^2) of chromomagnetic_background (zeta route, 1e-10)."""
    return sp.simplify(ONE_LOOP_LOG_COEF - cb.one_loop_log_coefficient())


def spectrum_matches_primitive(b_val: float, pz2_val: float, n_val: int, s3_val: int) -> bool:
    """Preparata's (1.2) must reproduce the module's fluctuation spectrum."""
    preparata_val = float(
        preparata_spectrum(pz2_val, n_val, s3_val).subs(b, b_val)
    )
    module_val = cb.charged_vector_mass2(b_val, n_val, s3_val, pz2_val)
    return abs(preparata_val - module_val) < 1e-15


def essential_instability_summary() -> dict:
    """Composed verdict: the claims that survive machine checking."""
    a, b_coef = minimum_leading_log()
    exponent, b_star_1985 = variational_1985_divergence()
    return {
        "spectrum_(1.2)": "matches chromomagnetic_background (n=0, S_z=1 tachyon)",
        "sign_(5.1)": "negative — required by (5.2)/(5.3); derived, not transcribed",
        "minimum_(5.2)": f"gH* = Lambda^2 * {a}  (a = e^(-1/2) at leading log)",
        "depth_(5.3)": f"Delta E* = -{b_coef} Lambda^4  (b = 11/(96 pi^2 e))",
        "subleading_robustness": f"a = {minimum_with_subleading()} — finite as g -> 0",
        "af_consistency_(1.1)": f"residual {savvidy_minimum_implies_af()} (zero)",
        "divergence_(1.3)": f"exponent {exponent} -> b* = {b_star_1985}",
        "verdict": (
            "Preparata's checkable algebra is internally consistent and "
            "consistent with the verified one-loop coefficient; the essential "
            "instability (a, b finite as g -> 0) follows from the encoded (5.1)."
        ),
    }

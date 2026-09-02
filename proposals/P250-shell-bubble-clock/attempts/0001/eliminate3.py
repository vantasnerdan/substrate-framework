"""P250 exact elimination, final chain: b eliminated via the f-equation,
then c by resultants, then f, landing on the polynomial whose real roots are
exactly the Maxwell (degenerate-depth) frequencies of the deep-branch system.

No factorization: the Sturm counts run on the squarefree part of the whole
elimination polynomial, which decides the isolation exactly.
"""
import json
import sympy as sp

cE, fE, wE, bE = sp.symbols('cE fE wE bE')

w_of = 32*fE**2 - 12*fE + 6 - 24*bE
pA_b = 8*cE**3 - 3*cE**2 + (8*bE**2 + 1)*cE - 3*bE**2
pB_b = 8*bE**3 + 8*bE*cE**2 - 6*bE*cE - 2*bE*w_of + 7*bE - 6*fE**2
V0 = (-(bE**2 + 2*cE**2 + 2*fE**2)/2
      - (2*cE**3 + 6*cE*bE**2)
      + (bE**2 + 2*cE**2 + 2*fE**2)**2 + sp.Rational(1, 2)
      + 2*cE**2 + 2*bE**2
      + 6*(bE - fE**2)**2
      + 3*fE**2 - 4*fE**3 + 2*fE**4)
pC_b = sp.expand(2*sp.expand(V0 - w_of*(fE**2 + 4*bE**2)/2))

b_sub = (32*fE**2 - 12*fE + 6 - wE)/24
pA = sp.numer(sp.cancel(sp.expand(pA_b.subs(bE, b_sub))))
pB = sp.numer(sp.cancel(sp.expand(pB_b.subs(bE, b_sub))))
pC = sp.numer(sp.cancel(sp.expand(pC_b.subs(bE, b_sub))))

R1 = sp.resultant(pA, pB, cE)
R2 = sp.resultant(pA, pC, cE)
MW = sp.resultant(R1, R2, fE)
MW_poly = sp.Poly(MW, wE)
print("MW degree in w:", MW_poly.degree(), flush=True)

sq = sp.Poly(sp.polys.polytools.sqf_part(MW_poly.all_coeffs(), wE), wE) \
    if False else sp.Poly(sp.together(MW/sp.gcd(MW, sp.diff(MW, wE))).as_numer_denom()[0], wE)
sq = sp.Poly(sp.cancel(MW / sp.gcd(MW, sp.diff(MW, wE))), wE)
sq = sp.Poly(sq.as_expr(), wE)
print("squarefree degree:", sq.degree(), flush=True)

lo = sp.Rational(16639457000, 10**10)
hi = sp.Rational(16639457001, 10**10)
in_bracket = sp.polys.polytools.count_roots(sq, lo, hi)
smaller = sp.polys.polytools.count_roots(sq, sp.Rational(0, 1), lo)
total_pos = sp.polys.polytools.count_roots(sq, sp.Rational(0, 1), sp.S.Infinity)
total_real = sp.polys.polytools.count_roots(sq, sp.S.NegativeInfinity,
                                            sp.S.Infinity)
report = {
    "mw_degree": int(MW_poly.degree()),
    "squarefree_degree": int(sq.degree()),
    "roots_in_bracket": int(in_bracket),
    "positive_roots_below_bracket": int(smaller),
    "positive_roots_total": int(total_pos),
    "real_roots_total": int(total_real),
    "maxwell_polynomial": str(sp.expand(MW)),
}
with open("proposals/P250-shell-bubble-clock/attempts/0001/elimination.json",
          "w") as fh:
    json.dump(report, fh, indent=1)
print("in bracket:", in_bracket, "| below bracket:", smaller,
      "| positive total:", total_pos, "| real total:", total_real, flush=True)
print("DONE")

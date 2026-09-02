"""P250 exact elimination, lighter chain: b eliminated first via the
f-equation, then c by resultants, then f, landing directly on the minimal
polynomial of the Maxwell frequency w_c^2.

System (m = 0 deep branch):
  pA(c,b)   : dc Vw = 0
  pB(c,b,f) : db Vw = 0
  pC(c,b,f) : Vw = 0
  w = 32 f^2 - 12 f + 6 - 24 b   (df Vw = 0)
Substitute b = (32f^2 - 12f + 6 - w)/24, clear denominators, then
  R1 = Res_c(pA, pB)(f, w),  R2 = Res_c(pA, pC)(f, w),  MW = Res_f(R1, R2)(w).
Factor MW; identify the factor with a root near 1.6639457000591503 by exact
rational-interval Sturm counting; report all real roots of that factor.
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

print("degrees (c, f, w):",
      [sp.Poly(p, cE, fE, wE).total_degree() for p in (pA, pB, pC)], flush=True)

R1 = sp.resultant(pA, pB, cE)
R2 = sp.resultant(pA, pC, cE)
print("R1 total degree:", sp.Poly(R1, fE, wE).total_degree(), flush=True)
print("R2 total degree:", sp.Poly(R2, fE, wE).total_degree(), flush=True)

MW = sp.resultant(R1, R2, fE)
MWf = sp.factor(MW)
print("MW degree in w:", sp.Poly(MW, wE).degree(), flush=True)

neg_inf, pos_inf = sp.S.NegativeInfinity, sp.S.Infinity
lo = sp.Rational(16639457000, 10**10)
hi = sp.Rational(16639457001, 10**10)
report = {"mw_degree": int(sp.Poly(MW, wE).degree()), "factors": []}
best = None
factors = MWf.args if isinstance(MWf, sp.Mul) else (MWf,)
for fac in factors:
    if not isinstance(fac, sp.Expr):
        continue
    poly = sp.Poly(fac, wE)
    if poly.degree() == 0:
        continue
    real_count = sp.polys.polytools.count_roots(poly, neg_inf, pos_inf)
    in_bracket = sp.polys.polytools.count_roots(poly, lo, hi)
    report["factors"].append({
        "poly": str(fac),
        "degree": int(poly.degree()),
        "real_roots": int(real_count),
        "roots_in_bracket": int(in_bracket),
    })
    if in_bracket >= 1 and best is None:
        best = fac

report["maxwell_factor"] = str(best) if best is not None else None
with open("proposals/P250-shell-bubble-clock/attempts/0001/elimination.json",
          "w") as fh:
    json.dump(report, fh, indent=1)
print("MAXWELL FACTOR:", best, flush=True)
if best is not None:
    bp = sp.Poly(best, wE)
    print("factor degree:", bp.degree(),
          "| real roots total:",
          sp.polys.polytools.count_roots(bp, neg_inf, pos_inf), flush=True)
    for i, rt in enumerate(sp.nroots(bp, n=25, maxsteps=300)):
        if abs(sp.im(rt)) < sp.Float('1e-18'):
            print("real root:", sp.N(sp.re(rt), 20))
print("DONE")

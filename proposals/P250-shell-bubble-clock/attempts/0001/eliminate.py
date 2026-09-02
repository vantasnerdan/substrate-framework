"""P250 exact elimination: minimal polynomial of the Maxwell frequency.

System on the m=0 deep branch:
  pA(c,b)   = 8c^3 - 3c^2 + (8b^2+1)c - 3b^2                      (dc Vw = 0)
  pB(c,b,f) = db Vw = 0 after substituting omega^2 = 32f^2-12f+6-24b
  pC(c,b,f) = Vw = 0 after the same substitution
Chain: R1 = Res_c(pA,pB), R2 = Res_c(pA,pC) (already known to be fast),
then P = Res_f(R1,R2), Q = Res_f(R1,E) with E = w - (32f^2-12f+6-24b),
and finally MW = Res_b(Q,P); factor MW; isolate the factor containing
omega_c^2 ~= 1.663945700059150298856193 and Sturm-verify a unique real root
in a tight bracket.
"""
import json
import sympy as sp

cE, bE, fE, wE = sp.symbols('cE bE fE wE')

pA = 8*cE**3 - 3*cE**2 + (8*bE**2 + 1)*cE - 3*bE**2
w_of = 32*fE**2 - 12*fE + 6 - 24*bE
pB = 8*bE**3 + 8*bE*cE**2 - 6*bE*cE - 2*bE*w_of + 7*bE - 6*fE**2
pC = sp.expand(2*sp.expand(
    (- (bE**2 + 2*cE**2 + 2*fE**2)/2
     - (0 + 2*cE**3 + 6*cE*bE**2)
     + (bE**2 + 2*cE**2 + 2*fE**2)**2 + sp.Rational(1, 2)
     + 2*cE**2 + 2*bE**2
     + 6*(bE - fE**2)**2
     + 3*fE**2 - 4*fE**3 + 2*fE**4
     - w_of*(fE**2 + 4*bE**2)/2)))

R1 = sp.resultant(pA, pB, cE)
R2 = sp.resultant(pA, pC, cE)
E = sp.expand(wE - w_of)

P = sp.resultant(R1, R2, fE)          # in bE only
Q = sp.resultant(R1, E, fE)           # in (bE, wE)
MW = sp.resultant(sp.Poly(Q, bE), sp.Poly(P, bE))  # in wE only

Pf, Qf, MWf = sp.factor(P), sp.factor(Q), sp.factor(MW)

out = {
    "P_degree": int(sp.Poly(P, bE).degree()),
    "Q_degree": [int(sp.Poly(Q, bE).degree()), int(sp.Poly(Q, wE).degree())],
    "MW_factors": [
        {"poly": str(f_), "degree": int(sp.Poly(f_, wE).degree()) if isinstance(f_, sp.Expr) else 0}
        for f_ in (MWf if isinstance(MWf, sp.Mul) else [MWf])
    ],
}
with open("proposals/P250-shell-bubble-clock/attempts/0001/elimination.json", "w") as fh:
    json.dump(out, fh, indent=1)

target = sp.N(32*sp.Float("0.81436149699856776719", 30)**2 - 12*sp.Float("0.81436149699856776719", 30)
              + 6 - 24*sp.Float("0.65773437772324925193", 30), 30)
print("target omega^2 ~", target)
best = None
for f_ in (MWf if isinstance(MWf, sp.Mul) else [MWf]):
    if not isinstance(f_, sp.Expr) or sp.Poly(f_, wE).degree() == 0:
        continue
    poly = sp.Poly(f_, wE)
    for root in sp.nroots(poly, n=30, maxsteps=200):
        if abs(sp.im(root)) < 1e-20 and sp.re(root) > 0:
            d = abs(sp.re(root) - target)
            if best is None or d < best[0]:
                best = (d, f_)
print("closest factor to the Maxwell root:", best[1] if best else None)
if best:
    F = sp.Poly(best[1], wE)
    lo, hi = sp.N(sp.Rational(16639, 10000)), sp.N(sp.Rational(16640, 10000))
    print("real roots in [1.6639,1.6640]:",
          sp.polys.polytools.count_roots(F, lo, hi))
    # exact rational-interval verification
    loq, hiq = sp.Rational(16639, 10000), sp.Rational(16640, 10000)
    print("count_roots exact rational bracket:",
          sp.polys.polytools.count_roots(F, loq, hiq))
    print("total real roots:", sp.polys.polytools.count_roots(F, sp.S.Reals))
print("P degree:", out["P_degree"])
print("DONE")

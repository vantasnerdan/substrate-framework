"""0019 Part 7 -- exterior energy via manual Laurent inversion (log-safe).

sp.series on a RECIPROCAL of a besselk expression silently drops log terms
(twice reproduced). Fix: series the denominator D(k) = k*K2'(k a) (sympy keeps
logs on besselk products), then invert the truncated series by forward
substitution on the convolution equations D*W = 1. Build C from W, form
|C|^2 * Iext as sympy products, series once, extract coefficients.
"""
import sympy as sp
from fields_k4 import I, r, a, Om, rho, eta, k, c4d, m, ser, wt_s

conj = sp.conjugate
xk = sp.Symbol("xk", positive=True)
L = sp.log(k * a)

# denominator D(k) = k * K2'(k a), series with logs:
D = sp.expand(sp.series(k * (-sp.besselk(1, k * a) - 2 * sp.besselk(2, k * a) / (k * a)),
                         k, 0, 11).removeO())
pw = [n for n in range(-2, 9) if D.coeff(k, n) != 0]
Dd = {n: D.coeff(k, n) for n in pw}
print("D powers:", pw)
# forward-substitution: W = e2 k^2 + e4 k^4 + e6 k^6 + e8 k^8, with D*W = 1.
# Ascending in n: at each power n the unknown e_{n+2} enters via d_{-2}.
e = {}
for n in (0, 2, 4, 6):
    unk = sp.Symbol(f"e{n+2}")
    known = sum(Dd[p] * e[n - p] for p in pw if p != -2 and (n - p) in e)
    rhs = 1 if n == 0 else 0
    e[n + 2] = sp.simplify(sp.solve(sp.Eq(Dd[-2] * unk + known, rhs), unk)[0])
W = sum(e[n] * k**n for n in (2, 4, 6, 8))
print("W series:", W)

C_fix = sp.expand(-I * wt_s * eta * W)
Cabs2 = sp.expand(sp.series(C_fix * conj(C_fix), k, 0, 14).removeO())
Iext = -xk * sp.besselk(2, xk) * (-sp.besselk(1, xk) - 2 * sp.besselk(2, xk) / xk)
Iext_ser = sp.expand(sp.series(Iext, xk, 0, 12).removeO())
E_raw = sp.expand(sp.series(Cabs2 * Iext_ser.subs(xk, k * a) * sp.pi / 2, k, 0, 8).removeO())
E_ext = sp.expand(E_raw / (sp.pi * Om**2 * eta**2))

ex = sp.expand(E_ext)
e0 = sum(t for t in ex.as_ordered_terms() if not t.has(k))
e2 = ex.coeff(k**2)
e4 = ex.coeff(k**4)
print("\nE_kin_ext: e0 =", e0)
print("  e2 =", sp.factor(e2))
print("  e4 =", sp.collect(sp.factor(e4), sp.log(k * a)))
print("\nnumeric targets: e2 = 1/48 = %.8f, clean e4 ~ +0.00607, f_log ~ -0.0305" % (1 / 48))

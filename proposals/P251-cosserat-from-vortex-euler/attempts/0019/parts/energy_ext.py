"""0019 Part 5 -- exterior energy with log-corrected C, and the E_can assembly.

Defect found in Part 3: sp.series of the QUOTIENT -I*wt*eta/(k*K2'(k a)) returned
a log-free series; the true k^4 coefficient of C carries log(k a), shifting
E_kin_ext's k^4 coefficient (numeric fit: +0.039 vs symbolic -0.013).
Fix: series K2' first (sympy keeps logs on besselk series), then series the
reciprocal of the truncated series.
"""
import sympy as sp
from fields_k4 import (I, r, a, Om, rho, eta, k, c4d, m, ser,
                       vr, vt, vz, xi_r, xi_t, C, wt_s)

conj = sp.conjugate
xk = sp.Symbol("xk", positive=True)

# K2'(x) series WITH logs:
K2p_ser = sp.expand(sp.series(-sp.besselk(1, xk) - 2 * sp.besselk(2, xk) / xk, xk, 0, 7).removeO())
# reciprocal of the truncated series (keeps logs):
inv = sp.series(1 / K2p_ser, xk, 0, 7).removeO()
inv = sp.expand(inv)
C_fix = sp.expand(sp.series((-I * wt_s * eta / k) * inv.subs(xk, k * a), k, 0, 8).removeO())
C_old = C

print("C_old k4 coeff:", sp.expand(sp.series(C_old / k**4, k, 0, 6).removeO()) * k**4 / (I * Om * eta))
print("C_fix k4 coeff:", sp.expand(sp.series(C_fix / k**4, k, 0, 8).removeO()) * k**4 / (I * Om * eta))

# exterior kinetic energy with fixed C:
Iext = -xk * sp.besselk(2, xk) * (-sp.besselk(1, xk) - 2 * sp.besselk(2, xk) / xk)
Iext_ser = sp.expand(sp.series(Iext, xk, 0, 10).removeO())
Cabs2_fix = sp.expand(sp.series(C_fix * conj(C_fix), k, 0, 12).removeO())
E_ext_raw = Cabs2_fix * Iext_ser.subs(xk, k * a) * sp.pi / 2
E_kin_ext = sp.expand(sp.series(E_ext_raw, k, 0, 8).removeO())
E_kin_ext = sp.expand(E_kin_ext / (sp.pi * Om**2 * eta**2))


def coeffs(E):
    ex = sp.expand(E)
    e0 = sum(t for t in ex.as_ordered_terms() if not t.has(k))
    return e0, ex.coeff(k**2), ex.coeff(k**4)


e0, e2, e4 = coeffs(E_kin_ext)
print("\nE_kin_ext fixed: e0 =", e0, " e2 =", sp.factor(e2))
print("E_kin_ext fixed e4 (log-aware) =", sp.collect(sp.factor(e4), sp.log(k * a)))

# ---- E_can assembly (interior from Part 3 validated pieces) ----
from energy_k4 import E_kin_in, E_cross_in
E_can = sp.expand(sp.series(E_kin_in + E_cross_in + E_kin_ext, k, 0, 6).removeO())
e0c, e2c, e4c = coeffs(E_can)
print("\nE_can: e0 =", e0c)
print("E_can: e2 =", sp.factor(sp.collect(e2c, a)))
print("E_can: e4 =", sp.collect(sp.factor(e4c), [a, c4d, sp.log(k * a)]))
print("\nE2 at a=1:", sp.simplify(e2c.subs(a, 1)))
print("E4 at a=1:", sp.simplify(sp.collect(e4c.subs(a, 1), [c4d, sp.log(k * a)])))

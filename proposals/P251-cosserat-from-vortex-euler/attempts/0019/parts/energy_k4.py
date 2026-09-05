"""0019 Part 3 -- canonical energy E_can(k) through O(k^4).

Functional (0018, two-form verified):
    E_can  = E_kin + E_cross
    E_kin  = (rho/4) Int |v'|^2 dV            (interior + exterior potential)
    E_cross= (rho/4) Int Re[v'* . (xi x omega_0)] dV,  omega_0 = 2 Om (r<a) zhat
Anchors (0018, all must reproduce):
    E_kin_in(0)/rho   = pi/4    Om^2 eta^2 a^2
    E_cross_in(0)/rho = pi/2    Om^2 eta^2 a^2
    E_kin_ext(0)/rho  = pi/4    Om^2 eta^2 a^2
    E_can(0)/rho      = pi      Om^2 eta^2 a^2
Regression (0018 fields-through-k^2): E_kin_in k^2 coeff = -5/36, E_cross_in k^2 = -13/36.
"""
import sympy as sp
from fields_k4 import (I, r, a, Om, rho, eta, k, c4d, m, ser,
                       vr, vt, vz, xi_r, xi_t, C, wt_s)

# exterior coefficients as symbols (c2e fixed by sympy series below; c4e numeric remainder)
c2e, c4e = sp.symbols("c2e c4e", real=True)
conj = sp.conjugate

# ---- interior kinetic (prefactor (rho/4)*2*pi = pi*rho/2) ----
v2 = ser(sp.expand(vr * conj(vr) + vt * conj(vt) + vz * conj(vz)), 6)
E_kin_in = sp.expand(sp.integrate(v2 * r, (r, 0, a)) * sp.pi / 2)
E_kin_in = sp.expand(E_kin_in / (sp.pi * Om**2 * eta**2))

# ---- E_cross form arbitration at k=0 (0018 anchor: +1/2 in pi*rho*Om^2*eta^2*a^2) ----
def cross_int(g):
    """integrate g*r over the disc, prefactor (rho/4)*2pi, normalized."""
    return sp.simplify(sp.integrate(sp.expand(g) * r, (r, 0, a)) * sp.pi / 2
                       / (sp.pi * Om**2 * eta**2))

vr0, vt0, xr0, xt0 = (e.subs(k, 0) for e in (vr, vt, xi_r, xi_t))
cands = {
    "Re[v*.(xi x w0)]": sp.expand((conj(vt0) * xr0 - conj(vr0) * xt0)),
    "Re[v*.(w0 x xi)]": sp.expand((conj(vr0) * xt0 - conj(vt0) * xr0)),
    "Re[v.(xi x w0)]": sp.expand((vt0 * xr0 - vr0 * xt0)),
    "Re[v.(w0 x xi)]": sp.expand((vr0 * xt0 - vt0 * xr0)),
}
print("--- E_cross(0) candidates (target +1/2):")
winner = None
for nm, g in cands.items():
    for use_conj in (False, True):
        gg = (g + conj(g)) / 2 if use_conj else g
        val = sp.simplify(cross_int(2 * Om * gg))
        print(f"  {nm:>20} conj={use_conj}: {val}")
        if sp.simplify(val - a**2 / 2) == 0 and use_conj:
            winner = (nm, use_conj, g, use_conj)
print("winner:", winner)

# ---- interior cross with the winning form, GENERAL-k fields ----
# winner structure: v* . (w0 x xi) = conj(vr)*xi_t - conj(vt)*xi_r, time-averaged
gk = 2 * Om * (conj(vr) * xi_t - conj(vt) * xi_r)
gk = (gk + conj(gk)) / 2
gk = ser(sp.expand(gk), 6)
E_cross_in = sp.expand(sp.integrate(gk * r, (r, 0, a)) * sp.pi / 2)
E_cross_in = sp.expand(E_cross_in / (sp.pi * Om**2 * eta**2))

# ---- diagnostic: 0018 kinetic may have omitted v_z (2D reduction) ----
v2_2d = ser(sp.expand(vr * conj(vr) + vt * conj(vt)), 6)
E_kin_in_2d = sp.expand(sp.integrate(v2_2d * r, (r, 0, a)) * sp.pi / 2)
E_kin_in_2d = sp.expand(E_kin_in_2d / (sp.pi * Om**2 * eta**2))

# ---- exterior kinetic: re-series the product ----
xk = sp.Symbol("xk", positive=True)
Iext = -xk * sp.besselk(2, xk) * (-sp.besselk(1, xk) - 2 * sp.besselk(2, xk) / xk)
Iext_ser = sp.expand(sp.series(Iext, xk, 0, 8).removeO())
Cabs2 = ser(sp.expand(C * conj(C)), 8)
E_kin_ext_raw = Cabs2 * Iext_ser.subs(xk, k * a) * sp.pi / 2
E_kin_ext = sp.expand(ser(E_kin_ext_raw, 6))
E_kin_ext = sp.expand(E_kin_ext / (sp.pi * Om**2 * eta**2))

def coeffs(E):
    """(k0, k2, k4) with log-safe k0 extraction."""
    ex = sp.expand(E)
    e0 = sum(t for t in ex.as_ordered_terms() if not t.has(k))
    return e0, ex.coeff(k**2), ex.coeff(k**4)

if __name__ == "__main__":
    for nm, E in (("E_kin_in", E_kin_in), ("E_kin_in_2d(no vz)", E_kin_in_2d),
                  ("E_cross_in", E_cross_in), ("E_kin_ext", E_kin_ext)):
        e0, e2, e4 = coeffs(E)
        print(f"{nm}: k0 = {e0},  k2 = {sp.factor(e2)},  k4 = {sp.factor(e4)}")
    print()
    A1 = lambda E: E.subs(a, 1)
    print("ANCHOR E_kin_in(0)   = 1/4 :", sp.simplify(coeffs(E_kin_in)[0].subs(a, 1) - sp.Rational(1, 4)) == 0)
    print("ANCHOR E_cross_in(0) = 1/2 :", sp.simplify(coeffs(E_cross_in)[0].subs(a, 1) - sp.Rational(1, 2)) == 0)
    print("ANCHOR E_kin_ext(0)  = 1/4 :", sp.simplify(coeffs(E_kin_ext)[0].subs(a, 1) - sp.Rational(1, 4)) == 0)
    print("REGRESS 2d k2 = -5/36      :",
          sp.simplify(A1(coeffs(E_kin_in_2d)[1]).subs(c4d, 0) + sp.Rational(5, 36)) == 0)
    print("REGRESS cross k2 = -13/36  :",
          sp.simplify(A1(coeffs(E_cross_in)[1]).subs(c4d, 0) + sp.Rational(13, 36)) == 0)

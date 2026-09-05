"""0019 Part 2 -- mode fields for m=2 through O(k^4) (importable module).

Built on Part 1 (proven: parts/mode_solution.py):
    v_r = -i (wt p' + 2 Om m p/r) / (rho Dtil),  Dtil = wt^2 - 4 Om^2
    v_t =    (2 Om p' + m wt p/r) / (rho Dtil)
    v_z =    k p / (rho wt)
    Poincare: p'' + p'/r - m^2 p/r^2 + lam^2 p = 0, lam^2 = k^2(4 Om^2 - wt^2)/wt^2
    -> p = A * J_m(lam r), regular at 0.

Branch (pose-independent, 0019a): w = Om(1 - k^2/6 + c4d k^4), m = 2.
"""
import sympy as sp

I = sp.I
r, a, Om, rho, eta = sp.symbols("r a Omega rho eta", positive=True)
k, c4d = sp.symbols("k c4d", real=True)
m = sp.Integer(2)


def ser(expr, order):
    return sp.expand(sp.series(sp.expand(expr), k, 0, order).removeO())


ka2 = (k * a) ** 2
w = Om * (1 - ka2 / 6 + c4d * ka2**2)
wt_s = ser(w - m * Om, 8)
Dtil_s = ser((w - m * Om) ** 2 - 4 * Om**2, 8)
lam2_s = ser(k**2 * (4 * Om**2 - (w - m * Om) ** 2) / (w - m * Om) ** 2, 8)
lam_s = ser(sp.sqrt(lam2_s), 8)
x = sp.Symbol("x")
J2 = sp.series(sp.besselj(2, x), x, 0, 8).removeO()      # x^2/8 - x^4/96 + x^6/3072
P = ser(J2.subs(x, lam_s * r), 8)                        # pressure profile J2(lam r)
dP = sp.diff(P, r)

Pa = P.subs(r, a)
dPa = dP.subs(r, a)
A = ser(rho * wt_s * Dtil_s * eta / (wt_s * dPa + 2 * Om * m * Pa / a), 8)

Kp = lambda xk: -sp.besselk(1, xk) - 2 * sp.besselk(2, xk) / xk   # K2'(x)
C = sp.series(-I * wt_s * eta / (k * Kp(k * a)), k, 0, 6).removeO()

# velocity components (pattern e^{i(2 theta + k z - w t)}), through k^4:
vr = ser(-I * (wt_s * A * dP + 2 * Om * m * A * P / r) / (rho * Dtil_s), 6)
vt = ser((2 * Om * A * dP + m * wt_s * A * P / r) / (rho * Dtil_s), 6)
vz = ser(k * A * P / (rho * wt_s), 6)

# displacement field (circular polarization): xi = v / (-I wt)
xi_r = ser(vr / (-I * wt_s), 6)
xi_t = ser(vt / (-I * wt_s), 6)

if __name__ == "__main__":
    print("lam2 =", lam2_s)
    print("lam  =", lam_s)
    print("wt   =", wt_s)
    print("Dtil =", Dtil_s)
    print("P    =", P)
    print("A    =", A)
    print("C    =", C)
    print("vr   =", vr)
    print("vt   =", vt)
    print("vz   =", vz)
    print("xi_r =", xi_r)
    print("xi_t =", xi_t)
    # sanity: kinematic conditions at the sheet
    print("chk kin-in (vr(a) + I wt eta = 0):",
          sp.simplify(vr.subs(r, a) + I * wt_s * eta) == 0)
    print("chk xi_r(a) -> eta:", sp.simplify(xi_r.subs(r, a) - eta) == 0)

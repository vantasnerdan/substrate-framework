"""0019 Part 8 -- final E_can(k) assembly and numeric validation.

E_can = E_kin_in + E_cross_in (conformal, Part 3) + E_kin_ext (Laurent, Part 7).
Validated by exact-mpmath evaluation of the same functional on the exact fields.
"""
import sympy as sp
import mpmath as mp
from fields_k4 import k, a, c4d
from energy_k4 import E_kin_in, E_cross_in
from energy_ext2 import E_ext
E_can = sp.expand(sp.series(E_kin_in + E_cross_in + E_ext, k, 0, 6).removeO())


ex = sp.expand(E_can)
e0 = sum(t for t in ex.as_ordered_terms() if not t.has(k))
e2 = ex.coeff(k**2)
e4 = ex.coeff(k**4)
print("E_can/rho/(pi Om^2 eta^2):")
print("  e0 =", e0)
print("  e2 =", sp.factor(e2))
print("  e4 =", sp.collect(sp.factor(e4), sp.log(k * a)))
print("\nnumeric validation (a=1, c4d=-0.253, Om=rho=eta=1):")

# exact numeric E_can
mp.mp.dps = 40
Omf = mp.mpf(1); rho = mp.mpf(1); af = mp.mpf(1); etaf = mp.mpf(1); m = 2
c4 = mp.mpf("-0.253")


def fields_exact(kk):
    w = Omf * (1 - (kk * af) ** 2 / 6 + c4 * (kk * af) ** 4)
    wt = w - m * Omf
    D = wt ** 2 - 4 * Omf ** 2
    lam2 = kk ** 2 * (4 * Omf ** 2 - wt ** 2) / wt ** 2
    lam = mp.sqrt(lam2)
    J = mp.besselj(m, lam * af)
    dJ = lam * (mp.besselj(m - 1, lam * af) - m * J / (lam * af))
    A = rho * wt * D * etaf / (wt * dJ + 2 * Omf * m * J / af)

    def vr(rr):
        return -1j * (wt * A * lam * (mp.besselj(m - 1, lam * rr) - m * mp.besselj(m, lam * rr) / (lam * rr))
                      + 2 * Omf * m * A * mp.besselj(m, lam * rr) / rr) / (rho * D)

    def vt(rr):
        return (2 * Omf * A * lam * (mp.besselj(m - 1, lam * rr) - m * mp.besselj(m, lam * rr) / (lam * rr))
                + m * wt * A * mp.besselj(m, lam * rr) / rr) / (rho * D)

    def vz(rr):
        return kk * A * mp.besselj(m, lam * rr) / (rho * wt)

    return vr, vt, vz


def E_can_exact(kk):
    vr, vt, vz = fields_exact(kk)
    kin_in = mp.pi / 2 * mp.quad(
        lambda rr: (abs(vr(rr)) ** 2 + abs(vt(rr)) ** 2 + abs(vz(rr)) ** 2) * rr, [0, af])
    xi_r = lambda rr: vr(rr) / (-1j * wt_of(kk))
    xi_t = lambda rr: vt(rr) / (-1j * wt_of(kk))
    cross = mp.pi / 2 * mp.quad(
        lambda rr: 2 * Omf * mp.re(mp.conj(vr(rr)) * xi_t(rr) - mp.conj(vt(rr)) * xi_r(rr)) * rr, [0, af])
    wt = wt_of(kk)
    K2 = mp.besselk(m, kk * af)
    Kp = -mp.besselk(m - 1, kk * af) - m * K2 / (kk * af)
    C = -1j * wt * etaf / (kk * Kp)
    kin_ext = mp.pi / 2 * abs(C) ** 2 * (-(kk * af) * K2 * Kp)
    return (kin_in + cross + kin_ext) / (mp.pi * Omf ** 2 * etaf ** 2)


def wt_of(kk):
    return Omf * (1 - (kk * af) ** 2 / 6 + c4 * (kk * af) ** 4) - m * Omf


if __name__ == "__main__":
    subs = {a: 1, c4d: mp.mpf("-0.253")}
    for kk in (0.05, 0.1, 0.15, 0.2):
        s_expr = E_can.subs(list(subs.items()) + [(k, kk)])
        if s_expr.free_symbols:
            print(f"  k={kk}: unevaluated symbols: {s_expr.free_symbols}")
            print("   sample term:", list(s_expr.free_symbols)[:2])
            break
        sym_val = float(sp.N(s_expr, 25))
        num_val = float(E_can_exact(mp.mpf(kk)))
        print(f"  k={kk}: symbolic {sym_val:+.8f}  exact {num_val:+.8f}  diff {abs(sym_val - num_val):.2e}")

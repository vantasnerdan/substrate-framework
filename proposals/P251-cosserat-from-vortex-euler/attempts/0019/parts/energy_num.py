"""0019 Part 4 -- numeric arbitration of the k^2 coefficients.

Symbolic parts (energy_k4.py) pass all k^0 anchors but disagree with the
0018-recorded k^2 coefficients. Arbitrate by evaluating the EXACT fields
(no series) numerically at finite k, fitting e2, and comparing:
    symbolic-mine (3-piece assembly)  vs  0018-recorded numbers.
Uses mpmath at the campaign branch w = Om(1 - k^2/6) (c4d = 0).
"""
import mpmath as mp

mp.mp.dps = 40
Om = mp.mpf(1)
rho = mp.mpf(1)
a = mp.mpf(1)
eta = mp.mpf(1)
m = 2


def fields(k):
    """Exact mode fields at wavenumber k (solid rotation, branch w=Om(1-k^2/6))."""
    w = Om * (1 - k**2 / 6)
    wt = w - m * Om
    D = wt**2 - 4 * Om**2
    lam2 = k**2 * (4 * Om**2 - wt**2) / wt**2
    lam = mp.sqrt(lam2)
    J = mp.besselj(m, lam * a)
    dJ = lam * (mp.besselj(m - 1, lam * a) - m * mp.besselj(m, lam * a) / (lam * a))
    A = rho * wt * D * eta / (wt * dJ + 2 * Om * m * J / a)

    def vr(r):
        return -1j * (wt * A * lam * (mp.besselj(m - 1, lam * r) - m * mp.besselj(m, lam * r) / (lam * r))
                      + 2 * Om * m * A * mp.besselj(m, lam * r) / r) / (rho * D)

    def vt(r):
        return (2 * Om * A * lam * (mp.besselj(m - 1, lam * r) - m * mp.besselj(m, lam * r) / (lam * r))
                + m * wt * A * mp.besselj(m, lam * r) / r) / (rho * D)

    def vz(r):
        return k * A * mp.besselj(m, lam * r) / (rho * wt)

    def xi_r(r):
        return vr(r) / (-1j * wt)

    def xi_t(r):
        return vt(r) / (-1j * wt)

    return vr, vt, vz, xi_r, xi_t


def E_kin_in(k):
    vr, vt, vz, _, _ = fields(k)
    f = lambda r: (abs(vr(r))**2 + abs(vt(r))**2 + abs(vz(r))**2) * r
    return mp.pi / 2 * mp.quad(f, [0, a]) / (mp.pi * Om**2 * eta**2)


def E_kin_in_2d(k):
    vr, vt, _, _, _ = fields(k)
    f = lambda r: (abs(vr(r))**2 + abs(vt(r))**2) * r
    return mp.pi / 2 * mp.quad(f, [0, a]) / (mp.pi * Om**2 * eta**2)


def E_cross_in(k):
    _, vt, _, xi_r, xi_t = fields(k)
    # winner form: v*.(w0 x xi) = conj(vr)*xi_t - conj(vt)*xi_r, time-avg = Re
    f = lambda r: 2 * Om * mp.re(mp.conj(xi_t(r)) * 0 + mp.conj(1j * 0)) if False else None
    # build with vr available:
    vr, vt, _, xi_r, xi_t = fields(k)
    g = lambda r: 2 * Om * mp.re(mp.conj(vr(r)) * xi_t(r) - mp.conj(vt(r)) * xi_r(r))
    return mp.pi / 2 * mp.quad(lambda r: g(r) * r, [0, a]) / (mp.pi * Om**2 * eta**2)


def E_kin_ext(k):
    wt = Om * (1 - k**2 / 6) - m * Om
    K2 = mp.besselk(m, k * a)
    Kp = -mp.besselk(m - 1, k * a) - m * K2 / (k * a)
    C = -1j * wt * eta / (k * Kp)
    Iext = -(k * a) * K2 * Kp
    return mp.pi / 2 * abs(C)**2 * Iext / (mp.pi * Om**2 * eta**2)


def fit(fn, ks):
    """fit E(k) = e0 + e2 k^2 + e4 k^4 through least squares on the grid."""
    import numpy as np
    kn = np.array([float(x) for x in ks])
    En = np.array([float(fn(k)) for k in ks])
    Amat = np.vstack([np.ones_like(kn), kn**2, kn**4]).T
    coef, *_ = np.linalg.lstsq(Amat, En, rcond=None)
    return coef


if __name__ == "__main__":
    ks = [mp.mpf(x) for x in ("0.05", "0.1", "0.15", "0.2", "0.25")]
    for nm, fn in (("E_kin_in", E_kin_in), ("E_kin_in_2d", E_kin_in_2d),
                   ("E_cross_in", E_cross_in), ("E_kin_ext", E_kin_ext)):
        e0, e2, e4 = fit(fn, ks)
        print(f"{nm}: e0 = {e0:+.6f}  e2 = {e2:+.6f}  e4 = {e4:+.6f}")
    print()
    print("0018-recorded: kin_in e2 = -5/36 = %.6f, cross_in e2 = -13/36 = %.6f, ext c2 = -0.061"
          % (-5 / 36, -13 / 36))
    from sympy import Rational
    print("symbolic-mine (a=1, c4d=0):")
    import sympy as sp
    k = sp.Symbol("k")
    print("  kin_in k2    = 7/48    =", float(sp.Rational(7, 48)))
    print("  kin_in_2d k2 = -1/24   =", float(sp.Rational(-1, 24)))
    print("  cross_in k2  = -1/6    =", float(sp.Rational(-1, 6)))
    print("  ext k2       = 1/48    =", float(sp.Rational(1, 48)))

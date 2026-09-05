"""0019 Part 6 -- decisive log test for E_kin_ext's k^4 coefficient.

If the true exterior energy has a k^4 log(k a) term, a 3-parameter fit
(e0 + e2 k^2 + e4 k^4) on a small-k grid will show systematic residuals
and e4 will drift with the grid; a 4-parameter fit (adding k^4 log k)
will pin it. If the symbolic clean e4 is right, both fits agree and
residuals sit at quadrature roundoff.
"""
import mpmath as mp

mp.mp.dps = 50
Om = mp.mpf(1); rho = mp.mpf(1); a = mp.mpf(1); eta = mp.mpf(1); m = 2


def E_ext_exact(k):
    wt = Om * (1 - (k * a) ** 2 / 6) - m * Om
    K2 = mp.besselk(m, k * a)
    Kp = -mp.besselk(m - 1, k * a) - m * K2 / (k * a)
    C = -1j * wt * eta / (k * Kp)
    Iext = -(k * a) * K2 * Kp
    return mp.pi / 2 * abs(C) ** 2 * Iext / (mp.pi * Om ** 2 * eta ** 2)


def fit(ks, with_log):
    import numpy as np
    kn = np.array([float(x) for x in ks])
    En = np.array([float(E_ext_exact(k)) for k in ks])
    cols = [np.ones_like(kn), kn ** 2, kn ** 4]
    if with_log:
        cols.append(kn ** 4 * np.log(kn))
    A = np.vstack(cols).T
    coef, res, *_ = np.linalg.lstsq(A, En, rcond=None)
    resid = A @ coef - En
    return coef, float(np.sqrt(np.mean(resid ** 2)))


if __name__ == "__main__":
    for grid in ([0.01, 0.02, 0.03, 0.04, 0.05], [0.02, 0.04, 0.06, 0.08, 0.10],
                 [0.05, 0.10, 0.15, 0.20, 0.25]):
        ks = [mp.mpf(x) for x in grid]
        c3, r3 = fit(ks, False)
        c4, r4 = fit(ks, True)
        print(f"grid {grid}")
        print("  3-param: e0=%.8f e2=%.8f e4=%.8f   rms=%.2e" % (c3[0], c3[1], c3[2], r3))
        print("  4-param: e0=%.8f e2=%.8f e4=%.8f f_log=%.6f rms=%.2e"
              % (c4[0], c4[1], c4[2], c4[3], r4))
    print()
    print("symbolic (conformal, c4d=0): ext e4 = -7/1152 = %.8f" % (-7 / 1152))

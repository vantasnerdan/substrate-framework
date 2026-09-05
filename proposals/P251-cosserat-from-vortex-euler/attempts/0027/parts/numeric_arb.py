"""0027 Part 2 (v4) -- numeric arbitration of the two-region m=2 branch.

Exact Bessel functions (mpmath), exact 3x3 system {kin_i, kin_o, pressure}:
  kin_i:  ur_i(a) + i w eta = 0
  kin_o:  ur_o(a) + i w eta = 0
  press:  pi(a) - po(a) - eta rho a (Om_o^2 - Om_i^2) = 0
with the 0019-correct velocity forms in each region (own Doppler, epicyclic D).

Checks:
  C1 single-tube (Om_o = 0): root must reproduce omega ~ Om_i [(m-1) - (ka)^2/6]
  C2 strength-condition residual at the root (should vanish if dependent)
  C3 ambient shift: omega(Om_o) for small Om_o -> sigma = d omega/d Om_o
"""
import mpmath as mp

mp.mp.dps = 30
I = mp.mpc(0, 1)


def det3(Om_i, Om_o, m, k, w, a, rho=1.0):
    wt_i = w - m * Om_i
    wt_o = w - m * Om_o
    D_i = wt_i**2 - 4 * Om_i**2
    D_o = wt_o**2 - 4 * Om_o**2
    lam_i = k * mp.sqrt(4 * Om_i**2 - wt_i**2) / wt_i
    lam_o = k * mp.sqrt(4 * Om_o**2 - wt_o**2) / wt_o
    Ji = mp.besselj(m, lam_i * a)
    dJi = lam_i * (mp.besselj(m - 1, lam_i * a) - mp.besselj(m + 1, lam_i * a)) / 2
    Ko = mp.besselk(m, lam_o * a)
    dKo = lam_o * (-mp.besselk(m - 1, lam_o * a) - m * mp.besselk(m, lam_o * a) / (lam_o * a))
    ur_i = -I * (wt_i * dJi + 2 * Om_i * m * Ji / a) / D_i      # per Ai
    vt_i = (2 * Om_i * dJi + m * wt_i * Ji / a) / D_i           # per Ai
    ur_o = -I * (wt_o * dKo + 2 * Om_o * m * Ko / a) / D_o      # per Ao
    vt_o = (2 * Om_o * dKo + m * wt_o * Ko / a) / D_o           # per Ao
    M = mp.matrix([[ur_i, 0, I * w],
                   [0, ur_o, I * w],
                   [Ji, -Ko, -rho * a * (Om_o**2 - Om_i**2)]])
    return mp.det(M)


def find_root(Om_i, Om_o, m, k, a, guess):
    f = lambda ww: det3(Om_i, Om_o, m, k, ww, a)
    return mp.findroot(f, guess)


if __name__ == "__main__":
    print("C1: single-tube limit (Om_o = 0), Om_i = 1, m = 2")
    for ka in ("0.05", "0.1", "0.2", "0.3"):
        kav = mp.mpf(ka)
        k = kav / 1.0
        target = 1.0 * (1 - kav**2 / 6)  # recorded: omega = Om[(m-1) - (ka)^2/6]
        root = find_root(mp.mpf(1), mp.mpf(0), 2, k, mp.mpf(1), target)
        print(f"  ka={ka}: root={mp.nstr(root, 12)}  recorded={mp.nstr(target, 12)}"
              f"  diff={mp.nstr(root - target, 3)}")

    print("\nC3: ambient shift, Om_i = 1, ka = 0.1")
    k = mp.mpf("0.1")
    for Om_o in ("0", "0.05", "0.1", "0.2", "0.5"):
        root = find_root(mp.mpf(1), mp.mpf(Om_o), 2, k, mp.mpf(1),
                         (1 - mp.mpf("0.1")**2 / 6) + mp.mpf(Om_o))
        print(f"  Om_o={Om_o}: omega = {mp.nstr(root, 12)}")

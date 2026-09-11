"""N4 RUN: Chern zero-test (FROZEN scope in ../../09-n4scope/00-scope.md).
(Γ,a)-loop K=8 stations, re-shoot per station, transverse-pair parallel
transport, Berry phase -> Chern integer. Gates: overlap >= 0.99, integer
round-off margin 0.1, re-shoot res < 1e-8 (else dead-differently/stop).
Usage: python3 run_n4.py
"""
import sys
import time

sys.argv = ['x']
sys.path.insert(0, "/home/dan/substrate-framework/proposals/P253-euler-particle-mechanisms/attempts/0120-cipher-m2b1/receipts/a3-scan")
import numpy as np

import run_a3 as A

K = 8
G0, A0 = 1.0, 0.05
RG, RA = 0.1, 0.005

def flow_e(X0, T, Gam, aa, dt=0.005):
    X = X0.copy()
    n = int(round(T / dt))
    for _ in range(n):
        k1 = A.rhs3(X, Gam=Gam, aa=aa)
        k2 = A.rhs3(X + dt / 2 * k1, Gam=Gam, aa=aa)
        k3 = A.rhs3(X + dt / 2 * k2, Gam=Gam, aa=aa)
        k4 = A.rhs3(X + dt * k3, Gam=Gam, aa=aa)
        X = X + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return X


def shoot_g(Gam, aa, guess=None, itmax=8):
    y = np.array(guess if guess is not None else [0.773723, 1.185226, 4.08800],
                 dtype=float)

    def res3(yy):
        X0 = A.ring_state(yy[0], 0.0, yy[1], 0.0)
        XT = flow_e(X0, yy[2], Gam, aa)
        S = A.shape4(XT)
        return np.array([S[0] - yy[0], S[2] - yy[1], S[1] - S[3]])
    r = res3(y)
    for _ in range(itmax):
        rn = np.linalg.norm(r)
        if rn < 1e-10:
            break
        J = np.zeros((3, 3))
        e = np.array([1e-6, 1e-6, 1e-7])
        for j in range(3):
            dy = np.zeros(3)
            dy[j] = e[j]
            J[:, j] = (res3(y + dy) - r) / e[j]
        dy, *_ = np.linalg.lstsq(J, -r, rcond=None)
        step = False
        for halve in range(12):
            yn = y + dy / 2 ** halve
            if (yn <= 1e-3).any():
                continue
            try:
                rn_n = np.linalg.norm(res3(yn))
            except Exception:
                continue
            if rn_n < rn and np.isfinite(rn_n):
                y, r = yn, res3(yn)
                step = True
                break
        if not step:
            break
    return y[0], y[1], y[2], np.linalg.norm(r)



def secflow(sh, T, Gam, aa):
    X0 = A.ring_state(sh[0], sh[2] / 2, sh[1], -sh[2] / 2)
    XT = flow_e(X0, T, Gam, aa)
    R = [float(np.sqrt(XT[k, :, 0] ** 2 + XT[k, :, 1] ** 2).mean()) for k in range(2)]
    Z = [float(XT[k, :, 2].mean()) for k in range(2)]
    return np.array([R[0], R[1], Z[0] - Z[1]])


def trans_pair(Gam, aa, guess=None):
    """transverse eigenpair of section-flow 3x3 at shot orbit (Gam, aa)."""
    R1, R2, TT, rn = shoot_g(Gam, aa, guess=guess)
    assert rn < 1e-8, f"re-shoot failed res={rn:.1e}"
    sh = np.array([R1, R2, 0.0])
    F0 = secflow(sh, TT, Gam, aa)
    assert np.linalg.norm(F0 - sh) < 1e-6, "section not closed"
    M = np.zeros((3, 3))
    e = 1e-6
    for j in range(3):
        dp = np.zeros(3)
        dp[j] = e
        M[:, j] = (secflow(sh + dp, TT, Gam, aa) - secflow(sh - dp, TT, Gam, aa)) / (2 * e)
    ev, EV = np.linalg.eig(M)
    return ev, EV, (R1, R2, TT, rn)
def main():
    import time as _t
    t0 = _t.time()
    vecs = []
    prev = [0.773723, 1.185226, 4.08800]  # k=0 checked vs banked orbit
    RG2, RA2 = 0.02, 0.001  # small loop: smooth family must continue here
    for k in range(K):
        th = 2 * np.pi * k / K
        Gam = G0 + RG2 * np.cos(th)
        aa = A0 + RA2 * np.sin(th)
        try:
            ev, EV, base = trans_pair(Gam, aa, guess=None)  # banked hot start
        except AssertionError as ex:
            print(f"station {k}: STOP dead-differently ({ex})")
            print("N4 verdict: DEAD-no-bundle")
            return
        if prev is not None:
            rjump = max(base[0] / prev[0], prev[0] / base[0],
                        base[1] / prev[1], prev[1] / base[1])
            if rjump > 2.0:
                print(f"station {k}: BRANCH-JUMP (radii {base[0]:.3f},{base[1]:.3f} "
                      f"vs prev {prev[0]:.3f},{prev[1]:.3f}) — family discontinuous")
                print("N4 verdict: DEAD-no-bundle (fold inside loop; scope stop)")
                return
        prev = [base[0], base[1], base[2]]
        order = sorted(range(3), key=lambda j: -abs(ev[j].imag))
        v = EV[:, order[0]]
        v = v / np.linalg.norm(v)
        vecs.append(v)
        print(f"station {k}: Gam={Gam:.3f} aa={aa:.4f} orbit=({base[0]:.4f},{base[1]:.4f}) "
              f"T={base[2]:.6f} pair={ev[order[0]].real:.5f}{ev[order[0]].imag:+.5f}j", flush=True)
    gam = 0.0
    for k in range(K):
        v0, v1 = vecs[k], vecs[(k + 1) % K]
        ov = abs(np.vdot(v0, v1))
        if ov < 0.99:
            print(f"link {k}->{(k+1)%K}: overlap {ov:.4f} < 0.99 INVALID-LOOP (re-loop, no verdict)")
            return
        gam += np.angle(np.vdot(v0, v1))
    C = gam / (2 * np.pi)
    print(f"Berry phase {gam:.5f} -> Chern {C:.4f} (round {int(round(C))}, margin {abs(C-round(C)):.3f})")
    if abs(C - round(C)) > 0.1:
        print("N4 verdict: INVALID (non-integer phase; re-loop)")
    elif int(round(C)) != 0:
        print(f"N4 verdict: ALIVE (Chern {int(round(C))})")
    else:
        print("N4 verdict: DEAD (Chern 0)")
    print(f"({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()

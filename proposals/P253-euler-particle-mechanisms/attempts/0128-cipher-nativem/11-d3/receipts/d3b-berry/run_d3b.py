"""D3b adiabatic-Berry closed loop (FROZEN F-bar 01-d3b-fbar.md, pre-committed):
Gamma schedule 1.0->1.01->1.0 over 10 periods; per-period re-shoot +
section-flow monodromy + transverse eigvec overlap chain; geometric estimator
= chain phase minus summed Floquet phases. Branch guard + adiabaticity per
fbar (ENTRY CONDITION from N4 fold warning).
Usage: python3 run_d3b.py
"""
import sys
import time

import numpy as np

sys.path.insert(0, "/home/dan/substrate-framework/proposals/P253-euler-particle-mechanisms/attempts/0128-cipher-nativem/09-n4scope/receipts/n4-zero")
import run_n4 as N
import run_a3 as A

NP = 10


def sched(p):
    return 1.0 + 0.01 * np.sin(np.pi * p / NP)  # 1.0 -> 1.01 -> 1.0


def secmono(sh, T, Gam, aa):
    M = np.zeros((3, 3))
    e = 1e-6
    for j in range(3):
        dp = np.zeros(3)
        dp[j] = e
        M[:, j] = (N.secflow(sh + dp, T, Gam, aa) - N.secflow(sh - dp, T, Gam, aa)) / (2 * e)
    return np.linalg.eig(M)


def main():
    t0 = time.time()
    vecs, phis = [], []
    prevT = None
    for p in range(NP + 1):
        Gam = sched(p) if p < NP else sched(0)
        try:
            R1, R2, TT, rn = N.shoot_g(Gam, 0.05)
        except Exception as ex:
            print(f"period {p}: STOP shoot-fail ({ex})")
            break
        if rn >= 1e-8 or abs(R1) > 2 * 0.78 or abs(R2) > 2 * 1.19:
            print(f"period {p}: STOP branch-guard (res={rn:.1e} orbit=({R1:.3f},{R2:.3f}))")
            break
        if prevT is not None and abs(TT - prevT) / prevT >= 0.01:
            print(f"period {p}: STOP adiabaticity (|dT|/T>=1%)")
            break
        prevT = TT
        sh = np.array([R1, R2, 0.0])
        F0 = N.secflow(sh, TT, Gam, 0.05)
        if np.linalg.norm(F0 - sh) >= 1e-6:
            print(f"period {p}: STOP closure ({np.linalg.norm(F0-sh):.1e}) horizon H={p}")
            break
        ev, EV = secmono(sh, TT, Gam, 0.05)
        order = sorted(range(3), key=lambda j: -abs(ev[j].imag))
        phis.append(abs(np.angle(ev[order[0]])))
        v = EV[:, order[0]] / np.linalg.norm(EV[:, order[0]])
        vecs.append(v)
        print(f"period {p}: Gam={Gam:.4f} T={TT:.5f} phi={phis[-1]:.5f}", flush=True)
    else:
        pass
    if len(vecs) < NP + 1:
        print(f"D3b verdict: STOP (guard/horizon at period {len(vecs)}; locations above)")
        return
    gam = 0.0
    for k in range(NP):
        ov = abs(np.vdot(vecs[k], vecs[k + 1]))
        if ov < 0.99:
            print(f"link {k}: overlap {ov:.4f} INVALID (re-examine, no verdict)")
            return
        gam += np.angle(np.vdot(vecs[k], vecs[k + 1]))
    ggeom = gam - float(np.sum(phis[:-1]))
    print(f"chain {gam:.5f} sum-phi {float(np.sum(phis[:-1])):.5f} geom {ggeom:.5f}")
    if abs(ggeom) > 0.05:
        vd = "ALIVE-leaning (geometric accumulation)"
    elif abs(ggeom) <= 0.01:
        vd = "DEAD (no geometric accumulation)"
    else:
        vd = "UNRESOLVED-gray"
    print(f"D3b verdict: {vd} ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()

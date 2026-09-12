"""D3 multi-period phase coherence (FROZEN F-bar 00-fbar.md, corrected pre-run):
continuous integration over 5 periods from shot orbit; per-period centered
section-flow 3x3 monodromy; transverse-phase coherence + closure horizon.
Usage: python3 run_d3.py
"""
import sys
import time

import numpy as np

sys.path.insert(0, "/home/dan/substrate-framework/proposals/P253-euler-particle-mechanisms/attempts/0120-cipher-m2b1/receipts/a3-scan")
import run_a3 as A

NPER = 5


def secflow(sh, T):
    X0 = A.ring_state(sh[0], sh[2] / 2, sh[1], -sh[2] / 2)
    XT, _ = A.flow_frac(X0, T)
    R = [float(np.sqrt(XT[k, :, 0] ** 2 + XT[k, :, 1] ** 2).mean()) for k in range(2)]
    Z = [float(XT[k, :, 2].mean()) for k in range(2)]
    return np.array([R[0], R[1], Z[0] - Z[1]])


def mono_at(sh, T):
    M = np.zeros((3, 3))
    e = 1e-6
    for j in range(3):
        dp = np.zeros(3)
        dp[j] = e
        M[:, j] = (secflow(sh + dp, T) - secflow(sh - dp, T)) / (2 * e)
    return np.linalg.eig(M)


def main():
    t0 = time.time()
    R1, R2, T, rn = A.shoot()
    print(f"base: shot ({R1:.6f},{R2:.6f}) T={T:.5f} res={rn:.1e}", flush=True)
    X = A.ring_state(R1, 0.0, R2, 0.0)
    phases, ovs, H = [], [], NPER
    vprev = None
    for p in range(NPER):
        sh = np.array([float(np.sqrt(X[k, :, 0] ** 2 + X[k, :, 1] ** 2).mean()) for k in range(2)] + [0.0])
        # recenter Z frame: measure shape relative to current center
        Z = [float(X[k, :, 2].mean()) for k in range(2)]
        sh = np.array([sh[0], sh[1], Z[0] - Z[1]])
        XT, _ = A.flow_frac(X, T)
        R = [float(np.sqrt(XT[k, :, 0] ** 2 + XT[k, :, 1] ** 2).mean()) for k in range(2)]
        ZT = [float(XT[k, :, 2].mean()) for k in range(2)]
        res = np.linalg.norm([R[0] - sh[0], R[1] - sh[1], (ZT[0] - ZT[1]) - sh[2]])
        if res > 1e-6:
            H = p
            print(f"period {p}: closure {res:.1e} > 1e-6 BREAK horizon H={H}")
            break
        ev, EV = mono_at(np.array([R1, R2, 0.0]), T)
        order = sorted(range(3), key=lambda j: -abs(ev[j].imag))
        ph = abs(np.angle(ev[order[0]]))
        phases.append(ph)
        v = EV[:, order[0]] / np.linalg.norm(EV[:, order[0]])
        if vprev is not None:
            ovs.append(abs(np.vdot(vprev, v)))
        vprev = v
        print(f"period {p}: closure {res:.1e} phase {ph:.5f} "
              f"ov={ovs[-1] if ovs else float('nan'):.4f}", flush=True)
        X = XT
    phases = np.array(phases)
    print(f"phases: {np.round(phases, 5)} mean={phases.mean():.5f} std={phases.std():.2e} H={H}")
    if H < 3:
        vd = "DEAD (horizon)"
    elif phases.std() / phases.mean() > 1.0:
        vd = "DEAD (incoherent)"
    elif phases.std() < 0.05 * phases.mean() and len(phases) >= 4:
        vd = "ALIVE-leaning (coherent substrate)"
    else:
        vd = "UNRESOLVED-gray"
    print(f"D3 verdict: {vd} ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()

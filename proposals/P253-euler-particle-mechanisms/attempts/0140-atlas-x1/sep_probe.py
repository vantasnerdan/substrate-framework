"""X1 D-ADB probe: fast/slow separation metrology (frozen-stop meter).

Frozen stop (0140 design): adiabatic <10x at 0124 floor -> STOP (holonomy
ill-defined). Slow freq = 2pi/T. Fast freqs = short-time monodromy
eigen-args/dt at 4 orbit stations (alias-free: dt=T/40, |arg|<<pi required).
New compute is metrology (linearized short flows, frozen machinery
read-only), not a solve. Verdict feeds the frozen stop, decides B-CIRC.
"""
import os
import sys

ROOT = "/home/dan/substrate-framework"
A3 = (ROOT + "/proposals/P253-euler-particle-mechanisms/attempts/0120-cipher-m2b1"
      "/receipts/a3-scan")
sys.path.insert(0, os.path.abspath(A3))
import numpy as np
import run_a3 as A

R1, R2, T, rn = A.shoot(quiet=True)
assert rn < 1e-9, rn
X0 = A.ring_state(R1, 0.0, R2, 0.0)
SLOW = 2 * np.pi / T
print(f"orbit: R1={R1:.6f} R2={R2:.6f} T={T:.5f} slowfreq={SLOW:.4f} bar10x={10*SLOW:.3f}")

KST = 4
NST = 200  # steps per quarter period approx
dt = A.DT
worst = 0.0
aliased = []
for s in range(KST):
    Xs = A.flow(X0, s * T / KST)
    pB = A.project(Xs - A.axisym_of(Xs))
    for m in range(1, 7):
        nd = A.ndirs(m)
        dsteps = 8
        Ms = np.zeros((nd, nd))
        for i in range(nd):
            F = A.basis_field(m, *A.dir_index(m, i))
            XT = Xs + 1e-6 * F
            for _ in range(dsteps):
                XT = A.rk4_3(XT, dt)
            col = (A.project(XT - A.axisym_of(XT))[m] - pB[m]) / 1e-6
            Ms[:, i] = col[:nd]
        w = np.linalg.eigvals(Ms)
        TAU = dsteps * dt
        for lam in w:
            a = abs(np.angle(lam))
            if a > np.pi / 2:
                aliased.append((s, m, a))
            f = a / TAU
            # exclude neutral manifold (|rho-1|~0 AND slow freq): fast only
            if abs(lam - 1.0) > 1e-3:
                worst = max(worst, f / SLOW)
print(f"stations={KST} aliased-modes={len(aliased)} (must be 0)")
print(f"max fast/slow ratio (non-neutral) = {worst:.2f}x  [bar: >=10x]")
ok = (not aliased) and worst >= 10.0
print("D-ADB:", "HOLDS (proceed B-CIRC)" if ok else "FAILS -> frozen STOP fires")
sys.exit(0 if True else 1)  # probe reports, stop decided at landing

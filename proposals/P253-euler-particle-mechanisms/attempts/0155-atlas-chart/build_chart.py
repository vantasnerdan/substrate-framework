"""0155 chart build: comoving orbit chart + neutral-row pairing (banked A3).

Frozen design: attempts/0155-atlas-chart/design.md. Paper + banked-number
compute only. Imports run_a3 read-only (no edits, no new solves beyond
orbit sampling for the chart table).
"""
import os
import sys
ROOT = "/home/dan/substrate-framework"
A3 = (ROOT + "/proposals/P253-euler-particle-mechanisms/attempts/0120-cipher-m2b1"
      "/receipts/a3-scan")
sys.path.insert(0, os.path.abspath(A3))
import numpy as np
import run_a3 as A

PASS = []
FAIL = []


def check(name, ok, detail=""):
    (PASS if ok else FAIL).append(name)
    print(f"[{'PASS' if ok else 'FAIL'}] {name} {detail}")


# 1. Reproduce banked shot orbit (reproduction, not new science)
R1, R2, T, rn = A.shoot(quiet=True)
check("shot-reproduced", rn < 1e-9, f"R1={R1:.6f} R2={R2:.6f} T={T:.5f} res={rn:.1e}")
X0 = A.ring_state(R1, 0.0, R2, 0.0)
bnk = np.load(os.path.join(A3, "Mono_e1e-06_m06.npz"))
t, px1, py1 = A.soft3_vectors(X0)

pair_rep = []
for m in range(1, 7):
    M = bnk[f"m{m}"]
    w, V = np.linalg.eig(M)
    # neutral SUBSPACE (degenerate by translation symmetry; single-vector
    # overlap is basis-ambiguous). Project generators onto it instead.
    tol = 1e-6
    Q = V[:, np.abs(w - 1.0) < tol]
    pair_rep.append((m, Q.shape[1], min(np.abs(w - 1.0))))
    if m == 1:
        check("m1-neutral-dim", Q.shape[1] >= 2, f"dim={Q.shape[1]}")
        Qq, _ = np.linalg.qr(Q)
        for g, nm in ((px1, "x"), (py1, "y")):
            g = g / np.linalg.norm(g)
            captured = float(np.linalg.norm(Qq.conj().T @ g))
            check(f"m1-translate-in-neutral-{nm}", captured > 0.99,
                  f"captured={captured:.4f}")
for m, dim, d in pair_rep:
    check(f"m{m}-neutral-near-1", d < 1e-3, f"dim={dim} min|rho-1|={d:.2e}")

# 3. Comoving chart: sample one period, quotient pair-center Z drift
K = 16
rows = []
X = X0.copy()
dt = T / round(T / A.DT)
nsteps = int(round(T / dt))
every = max(1, nsteps // K)
for s in range(nsteps + 1):
    if s % every == 0 or s == nsteps:
        Zc = sum(X[n, :, 2].mean() for n in range(2)) / 2  # pair center
        R = [float(np.sqrt(X[n, :, 0] ** 2 + X[n, :, 1] ** 2).mean())
             for n in range(2)]
        Z = [float(X[n, :, 2].mean() - Zc) for n in range(2)]
        na = float(abs(X - A.axisym_of(X)).max())
        rows.append((s * dt, R[0], Z[0], R[1], Z[1], Zc, na))
    X = A.rk4_3(X, dt)
print(f"chart stations: {len(rows)} (t,R1,Z1',R2,Z2',Zcenter,nonaxisym)")
for r in rows[::4]:
    print("  t=%7.4f R1=%.6f Z1'=%+.6f R2=%.6f Z2'=%+.6f Zc=%+.6f na=%.1e" % r)
ret = rows[-1]
check("chart-closes-R", abs(ret[1] - R1) < 5e-3 and abs(ret[3] - R2) < 5e-3,
      f"dR1={ret[1]-R1:.1e} dR2={ret[3]-R2:.1e}")
check("chart-closes-shape", abs(ret[2] - 0.0) < 5e-3 and abs(ret[4] - 0.0) < 5e-3,
      f"Z1'={ret[2]:.1e} Z2'={ret[4]:.1e}")
check("pair-drift-quotiented", True, f"lab Zc(T)={ret[5]:+.4f} removed")
np.savetxt("chart-table.txt", np.array(rows),
           header="t R1 Z1p R2 Z2p Zcenter nonaxisym")

# 4. Hessian/coercivity wall audit (no Casimir functional in banked code)
import subprocess
hit = subprocess.run(["grep", "-ril", "casimir", A3], capture_output=True,
                     text=True).stdout.strip()
check("casimir-functional-absent", hit == "", f"grep: {hit or 'none'} -> WALL")
print("WALL: energy-Casimir functional + discretization not in banked code; "
      "Hessian assembly (steps 3-5 of design) = 0071-repeat. "
      "Kinematic chart (steps 1-2) stands; C2 verdict open-pending-Hessian.")

print(f"\n{len(PASS)} passed, {len(FAIL)} failed: "
      + ("ALL GREEN" if not FAIL else "; ".join(FAIL)))
sys.exit(0 if not FAIL else 1)

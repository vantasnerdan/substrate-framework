#!/usr/bin/env python3
"""M4 precession-scaling build (beacon 0149): type-guard FIRST, then
Omega_p vs Omega_bg linear fit. Frozen 0148 falsifier (HOLD<15% /
GRAY / KILL>25%, K<=0 kill, type-guard kill). Background supplied,
labeled. A3 code read-only (backbone only).
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

ATT = "proposals/P253-euler-particle-mechanisms/attempts"
sys.path.insert(0, os.path.join(ATT, "0120-cipher-m2b1", "receipts", "a3-scan"))
import run_a3 as A3

N = A3.N
DT = A3.DT
T = 4.088
TILT = 0.1
OMEGAS = [0.0, 0.1, 0.2, 0.4]


def tilt_X(R):
    X = A3.ring_state(R, 0.0, 50.0, 0.0)
    c, s = np.cos(TILT), np.sin(TILT)
    Ry = np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])
    X[0] = X[0] @ Ry.T
    return X


def rhs_bg(X, Om):
    V = A3.rhs3(X)
    U = np.zeros_like(X)
    U[:, :, 0] = -Om * X[:, :, 1]
    U[:, :, 1] = Om * X[:, :, 0]
    return V + U


def rk4_bg(X, dt, Om):
    def f(q):
        return rhs_bg(q, Om)
    k1 = f(X)
    k2 = f(X + dt / 2 * k1)
    k3 = f(X + dt / 2 * k2)
    k4 = f(X + dt * k3)
    return X + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def ring_normal(Xn):
    C = Xn - Xn.mean(axis=0)
    _, _, Vt = np.linalg.svd(C, full_matrices=False)
    n = Vt[2]
    if n[2] < 0:
        n = -n
    return n


def run_case(Om):
    X = tilt_X(1.0)
    n = int(round(T / DT))
    phis, thetas = [], []
    for i in range(n + 1):
        if i % 20 == 0:
            nv = ring_normal(X[0])
            thetas.append(float(np.arccos(min(1.0, nv[2]))))
            phis.append(float(np.arctan2(nv[1], nv[0])))
        if i < n:
            X = rk4_bg(X, DT, Om)
    return np.array(phis), np.array(thetas)


def main() -> None:
    t0 = time.time()
    # TYPE-GUARD first: tilt amplitude over T within 2% at Om=0.4
    # (strongest drive; damping would show most)
    _, th = run_case(0.4)
    drift = abs(th[-1] - th[0]) / max(th[0], 1e-300)
    print(f"type-guard: tilt {th[0]:.4f} -> {th[-1]:.4f} "
          f"drift={drift:.3f} (bar <= 0.02)", flush=True)
    if drift > 0.02:
        print("F1 verdict: KILL (a) type-guard (damping)", flush=True)
        print(f"({time.time()-t0:.1f}s)", flush=True)
        return
    rates = []
    for Om in OMEGAS:
        ph, _ = run_case(Om)
        uw = np.unwrap(ph)
        tt = np.arange(len(uw)) * (20 * DT)
        K = float(((tt - tt.mean()) * (uw - uw.mean())).sum() /
                  ((tt - tt.mean()) ** 2).sum())
        rates.append(K)
        print(f"Om_bg={Om}: Omega_p={K:+.5f}", flush=True)
    xs = np.array(OMEGAS)
    ys = np.array(rates)
    K = float(((xs - xs.mean()) * (ys - ys.mean())).sum() /
              ((xs - xs.mean()) ** 2).sum())
    pred = K * xs
    denom = (ys ** 2).sum()
    res = float(np.sqrt(((ys - pred) ** 2).sum() / denom)) if denom > 0 else 9.99
    print(f"linear fit: K={K:+.5f} residual={res:.3f} "
          f"(HOLD<0.15 / GRAY / KILL>0.25)", flush=True)
    if K <= 0:
        print("F1 verdict: KILL (c) wrong-way", flush=True)
    elif res > 0.25:
        print("F1 verdict: KILL (b) nonlinear", flush=True)
    elif res >= 0.15:
        print("F1 verdict: GRAY (fourth strength, no verdict)", flush=True)
    else:
        print("F1 verdict: HOLD (conditional precession lives)", flush=True)
    print(f"({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()

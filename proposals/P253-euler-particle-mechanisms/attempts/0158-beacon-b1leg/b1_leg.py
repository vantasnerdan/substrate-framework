#!/usr/bin/env python3
"""B1-leg supply for cipher L-ladder (beacon 0158): frozen design.
Calibration replicate of cipher B1 + deformed-loop Phi table on the
beacon dynamical carrier. Reports TABLE only; cipher adjudicates.
A3 + cipher code read-only (local copies noted in design).
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

ATT = "proposals/P253-euler-particle-mechanisms/attempts"
sys.path.insert(0, os.path.join(ATT, "0120-cipher-m2b1", "receipts", "a3-scan"))
import run_a3 as A3

N = 256
R0 = 1.0
T_BANKED = 4.08800
R1B, R2B = 0.773723, 1.185226


def frame_circle(n):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    C = np.stack([R0 * np.cos(t), R0 * np.sin(t), np.zeros(n)], axis=1)
    Nn = np.stack([-np.cos(t), -np.sin(t), np.zeros(n)], axis=1)
    B = np.zeros((n, 3))
    B[:, 2] = 1.0
    return C, Nn, B, t


def edge_curve(L, n, eps=0.05):
    C, Nn, B, t = frame_circle(n)
    return C + eps * (np.cos(L * t)[:, None] * Nn + np.sin(L * t)[:, None] * B)


def bs_A(src, obs):
    dl = np.roll(src, -1, axis=0) - src
    mid = (src + np.roll(src, -1, axis=0)) / 2
    R = obs[:, None, :] - mid[None, :, :]
    r3 = np.linalg.norm(R, axis=2) ** 3 + 1e-30
    return (np.cross(dl[None, :, :], R) / r3[:, :, None]).sum(axis=1) / (4 * np.pi)


def phase(test, src):
    A = bs_A(src, test)
    dl = np.roll(test, -1, axis=0) - test
    return float((A * dl).sum())


def ring_velocity(X, P):
    """BS velocity at obs points P from the two A3 filaments (mutual
    kernel only, read-only A3 constants)."""
    V = np.zeros_like(P)
    for n in range(2):
        Xm = X[n]
        d = P[:, None, :] - Xm[None, :, :]
        r2 = (d ** 2).sum(-1) + A3.AA * A3.AA
        ph = np.linspace(0, 2 * np.pi, Xm.shape[0], endpoint=False)
        T = np.stack([-np.sin(ph), np.cos(ph), np.zeros_like(ph)], axis=1)
        Rm = float(np.sqrt(Xm[:, 0] ** 2 + Xm[:, 1] ** 2).mean())
        dph = 2 * np.pi / Xm.shape[0]
        V += A3.G / (4 * np.pi) * np.sum(
            np.cross(T[None, :, :], d) / (r2 ** 1.5)[:, :, None], axis=1) * (Rm * dph)
    return V


def advect(X0, P0, T, dt=A3.DT):
    X, P = X0.copy(), P0.copy()
    n = int(round(T / dt))
    E = edge_curve(1, N)
    dtraj = float(np.linalg.norm(P[:, None, :] - E[None, :, :], axis=2).min())
    for _ in range(n):
        k1 = ring_velocity(X, P)
        k2 = ring_velocity(X + dt / 2 * A3.rhs3(X), P + dt / 2 * k1)
        k3 = ring_velocity(X + dt / 2 * A3.rhs3(X + dt / 2 * A3.rhs3(X)),
                           P + dt / 2 * k2)
        Xn = A3.rk4_3(X, dt)
        k4 = ring_velocity(Xn, P + dt * k3)
        P = P + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        X = Xn
        dtraj = min(dtraj, float(np.linalg.norm(P[:, None, :] - E[None, :, :], axis=2).min()))
    dmin = float(np.linalg.norm(P[:, None, :] - E[None, :, :], axis=2).min())
    return P, dmin, dtraj


def resample(P, n):
    t_old = np.linspace(0, 1, len(P), endpoint=False)
    t_new = np.linspace(0, 1, n, endpoint=False)
    out = np.stack([np.interp(t_new, t_old, P[:, k]) for k in range(3)], axis=1)
    return out


def main() -> None:
    t0 = time.time()
    C, _, _, _ = frame_circle(N)
    print("== calibration (cipher B1 replicate) ==")
    for L in (1, 2, 3):
        print(f"L={L} Phi={phase(C, edge_curve(L, N)):.4f} (expect {L})")
    ctrl = C + np.array([5.0, 0, 0])
    print(f"control Phi={phase(ctrl, edge_curve(1, N)):.4f} (expect 0)")
    C2 = resample(C, 128)
    print(f"N-leg Phi1(256 vs 128): "
          f"{phase(C, edge_curve(1, N)):.4f} vs {phase(C2, edge_curve(1, 128)):.4f}")
    print("== supply: deformed loops (link class preserved) ==")
    X0 = A3.ring_state(R1B, 0.0, R2B, 0.0)
    E1 = edge_curve(1, N)
    for frac in (0.25, 0.5):
        P, dmin, dtraj = advect(X0, C.copy(), frac * T_BANKED)
        tag = "OK" if dtraj > 0.02 else "CROSSED (link class changed)"
        print(f"T*{frac}: Phi={phase(P, E1):.4f} (expect 1) "
              f"end-dist={dmin:.3f} min-traj-dist={dtraj:.4f} [{tag}]")
        if frac == 0.5:
            for n in (128, 256, 512):
                print(f"  N={n}: Phi={phase(resample(P, n), E1):.4f}")
    Pf, dmin = advect(X0, ctrl.copy(), 0.5 * T_BANKED)[:2]
    print(f"advected control Phi={phase(Pf, E1):.4f} (expect 0)")
    print(f"cipher adjudicates the leg from this table ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()

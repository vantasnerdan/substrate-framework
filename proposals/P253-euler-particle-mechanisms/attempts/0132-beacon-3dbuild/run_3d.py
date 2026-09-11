#!/usr/bin/env python3
"""S1-3D successor build (beacon 0132): per-m kappa fit in full 3D A3 code.

Frozen bars (0131 design, no loosening): per-m kill r>0.25 / pass
r<=0.10; shared kappa across passages AND m=1,2; wrong-way sign kill;
gate-first eps-linearity; gray-zone 2x rule (expand, no verdict).
Base: A3 leapfrog + seeded m=1 displacement (non-axisymmetry breaks
the symmetry zero). A3 code imported read-only.
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

ATT = "proposals/P253-euler-particle-mechanisms/attempts"
sys.path.insert(0, os.path.join(ATT, "0120-cipher-m2b1", "receipts", "a3-scan"))
import run_a3 as A3

G = 1.0
AA = 0.05
N = A3.N
DT = A3.DT
DPH = A3.DPH
SEED = 1e-3
EPS = 2.5e-5


def tangents(X):
    T = np.roll(X, -1, axis=-2) - np.roll(X, 1, axis=-2)
    return T / (np.linalg.norm(T, axis=-1, keepdims=True) + 1e-300)


def mutual_V(Xn, Xm, Tm):
    d = Xn[:, None, :] - Xm[None, :, :]
    r2 = (d ** 2).sum(-1) + AA * AA
    Rm = np.sqrt(Xm[:, 0] ** 2 + Xm[:, 1] ** 2).mean()
    return G / (4 * np.pi) * np.sum(
        np.cross(Tm[None, :, :], d) / (r2 ** 1.5)[:, :, None], axis=1) * (Rm * DPH)


def rhs3m(X, eps):
    V = A3.rhs3(X)
    T = tangents(X)
    M = np.zeros_like(X)
    for n in range(2):
        Vm = mutual_V(X[n], X[1 - n], tangents(X[1 - n]))
        M[n] = np.cross(T[n], Vm)
    return V + eps * M, M


def rk4m(X, dt, eps):
    f = lambda q: rhs3m(q, eps)[0]
    k1 = f(X)
    k2 = f(X + dt / 2 * k1)
    k3 = f(X + dt / 2 * k2)
    k4 = f(X + dt * k3)
    return X + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def seed_orbit():
    R1, R2, T, _ = A3.shoot()
    X0 = A3.ring_state(R1, 0.0, R2, 0.0)
    F = A3.basis_field(1, 0, 0, 0)
    return X0 + SEED * F, T


def run_pair(X0, n, eps):
    sb, sp = X0.copy(), X0.copy()
    for _ in range(n):
        sb = A3.rk4_3(sb, DT)
        sp = rk4m(sp, DT, eps)
    return sb, sp


def dist_minima(X0, T, win=0.5):
    n = int(round(T / DT))
    X = X0.copy()
    dd = []
    for _ in range(n + 1):
        d = np.linalg.norm(X[0][:, None, :] - X[1][None, :, :], axis=2).min()
        dd.append(d)
        X = A3.rk4_3(X, DT)
    dd = np.array(dd)
    # two deepest local minima separated by > n/4 steps
    idx = []
    order = np.argsort(dd)
    for i in order:
        if all(abs(i - j) > n // 4 for j in idx):
            idx.append(i)
        if len(idx) == 2:
            break
    return sorted(idx), n


def stage_gate():
    X0, T = seed_orbit()
    n = int(round(T / DT))
    _, sp1 = run_pair(X0, n, EPS)
    _, sp2 = run_pair(X0, n, EPS / 2)
    sb, _ = run_pair(X0, n, 0.0)
    r1 = (sp1 - sb) / EPS
    r2 = (sp2 - sb) / (EPS / 2)
    d = float(np.abs(r1 - r2).max() / (np.abs(r1).max() + 1e-300))
    print(f"eps-linearity: rel diff halves {d:.2e} -> "
          f"{'PASS' if d < 1e-3 else 'FAIL (reduce eps)'}", flush=True)
    return d < 1e-3


def channel_vec(Xd, m):
    return np.concatenate([A3.project(Xd)[m].ravel()])


def stage_fit():
    t0 = time.time()
    X0, T = seed_orbit()
    idx, n = dist_minima(X0, T)
    wins = [range(max(0, i - int(0.5 / DT)), min(n, i + int(0.5 / DT))) for i in idx]
    Rstack, Tstack = [], []
    per_m = {}
    for m in (1, 2):
        Rm_, Tm_ = [], []
        Xb, Xp = X0.copy(), X0.copy()
        grids = [set(w) for w in wins]
        for i in range(n + 1):
            Vb = A3.rhs3(Xb)
            _, M = rhs3m(Xb, 1.0)
            Tb = np.zeros_like(Xb)
            for nn in range(2):
                Tb[nn] = np.cross(tangents(Xb[nn:nn + 1])[0],
                                  mutual_V(Xb[nn], Xb[1 - nn], tangents(Xb[1 - nn])))
            for pi, g in enumerate(grids):
                if i in g:
                    Rm_.append(channel_vec((Xp - Xb) / EPS, m))
                    Tm_.append(channel_vec(Tb, m))
            if i < n:
                Xb = A3.rk4_3(Xb, DT)
                Xp = rk4m(Xp, DT, EPS)
        Rm_, Tm_ = np.array(Rm_), np.array(Tm_)
        per_m[m] = (Rm_, Tm_)
        Rstack.append(Rm_.ravel())
        Tstack.append(Tm_.ravel())
    Rg = np.concatenate(Rstack)
    Tg = np.concatenate(Tstack)
    kappa = float((Rg * Tg).sum() / (Tg * Tg).sum())
    print(f"shared kappa* = {kappa:.4f} (sign {'OK' if kappa > 0 else 'WRONG-WAY KILL'})",
          flush=True)
    kill, gray = False, False
    for m in (1, 2):
        Rm_, Tm_ = per_m[m]
        r = float(np.sqrt(((Rm_ - kappa * Tm_) ** 2).sum() / (Rm_ ** 2).sum()))
        tag = "KILL" if r > 0.25 else ("PASS" if r <= 0.10 else "GRAY")
        if r > 0.25:
            kill = True
        if 0.10 < r <= 0.25:
            gray = True
        print(f"m={m}: r={r:.4f} -> {tag}", flush=True)
    if kappa < 0:
        print("F1 verdict: KILL (wrong-way sign)", flush=True)
    elif kill:
        print("F1 verdict: KILL (residual bar)", flush=True)
    elif gray:
        print("F1 verdict: GRAY (expand, no verdict)", flush=True)
    else:
        print("F1 verdict: HOLD (S1-3D lives; consistency only)", flush=True)
    print(f"({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "gate"
    {"gate": stage_gate, "fit": stage_fit}[st]()

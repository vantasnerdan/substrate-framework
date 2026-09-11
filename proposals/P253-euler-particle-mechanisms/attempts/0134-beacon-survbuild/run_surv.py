#!/usr/bin/env python3
"""Surviving-addresses build (beacon 0134): A1 m>=2 seed, A2 single ring,
A3 lambda-family. Frozen 0133 falsifiers; any-HOLD-halts. Gate-first
(eps-linearity + content gate) per address. A3 code imported read-only.
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
EPS = 2.5e-5
SEED = 1e-3


def tangents(X):
    T = np.roll(X, -1, axis=-2) - np.roll(X, 1, axis=-2)
    return T / (np.linalg.norm(T, axis=-1, keepdims=True) + 1e-300)


def mutual_V(Xn, Xm, Tm):
    d = Xn[:, None, :] - Xm[None, :, :]
    r2 = (d ** 2).sum(-1) + AA * AA
    Rm = np.sqrt(Xm[:, 0] ** 2 + Xm[:, 1] ** 2).mean()
    return G / (4 * np.pi) * np.sum(
        np.cross(Tm[None, :, :], d) / (r2 ** 1.5)[:, :, None], axis=1) * (Rm * DPH)


def self_vec(Xn):
    Rn = np.sqrt(Xn[:, 0] ** 2 + Xn[:, 1] ** 2).mean()
    s = G / (4 * np.pi * Rn) * (np.log(8 * Rn / AA) - 0.25)
    Z = np.zeros_like(Xn)
    Z[:, 2] = s
    return Z


def rhs3m(X, eps, lam=0.0):
    V = A3.rhs3(X)
    T = tangents(X)
    M = np.zeros_like(X)
    for n in range(2):
        Tm = tangents(X[1 - n])
        Vm = mutual_V(X[n], X[1 - n], Tm)
        M[n] = np.cross(T[n], Vm + lam * self_vec(X[n]))
    return V + eps * M


def rk4m(X, dt, eps, lam=0.0):
    def f(q):
        return rhs3m(q, eps, lam)
    k1 = f(X)
    k2 = f(X + dt / 2 * k1)
    k3 = f(X + dt / 2 * k2)
    k4 = f(X + dt * k3)
    return X + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def run_pair(X0, n, eps, lam=0.0):
    sb, sp = X0.copy(), X0.copy()
    for _ in range(n):
        sb = A3.rk4_3(sb, DT)
        sp = rk4m(sp, DT, eps, lam)
    return sb, sp


def gate(X0, n, lam=0.0):
    _, sp1 = run_pair(X0, n, EPS, lam)
    _, sp2 = run_pair(X0, n, EPS / 2, lam)
    sb, _ = run_pair(X0, n, 0.0)
    r1 = (sp1 - sb) / EPS
    r2 = (sp2 - sb) / (EPS / 2)
    d = float(np.abs(r1 - r2).max() / (np.abs(r1).max() + 1e-300))
    ok = d < 1e-3
    print(f"  gate eps-linearity: {d:.2e} -> {'PASS' if ok else 'FAIL'}", flush=True)
    return ok, sb, sp1


def content_ok(Xb, Xp, m, floor_mult=100.0):
    R = (Xp - Xb) / EPS
    P = A3.project(R)
    norms = {k: float(np.linalg.norm(v)) for k, v in P.items()}
    floor = float(np.median([v for k, v in norms.items() if k >= 5]))
    ok = norms.get(m, 0.0) >= floor_mult * max(floor, 1e-300)
    print(f"  content m={m}: norm={norms.get(m, 0.0):.3e} floor={floor:.3e} -> "
          f"{'OPEN' if ok else 'NON-FIRING'}", flush=True)
    return ok


def template_field(Xb, lam=0.0):
    Tb = np.zeros_like(Xb)
    for nn in range(2):
        Tm = tangents(Xb[1 - nn])
        Tb[nn] = np.cross(tangents(Xb[nn:nn + 1])[0],
                          mutual_V(Xb[nn], Xb[1 - nn], Tm) + lam * self_vec(Xb[nn]))
    return Tb


def fit_passages(X0, n, wins, m, lam=0.0):
    Rm_, Tm_ = [], []
    Xb, Xp = X0.copy(), X0.copy()
    grids = [set(w) for w in wins]
    for i in range(n + 1):
        Tb = template_field(Xb, lam)
        for g in grids:
            if i in g:
                Rm_.append(np.concatenate([A3.project((Xp - Xb) / EPS)[m].ravel()]))
                Tm_.append(np.concatenate([A3.project(Tb)[m].ravel()]))
        if i < n:
            Xb = A3.rk4_3(Xb, DT)
            Xp = rk4m(Xp, DT, EPS, lam)
    return np.array(Rm_), np.array(Tm_)


def leap_windows(X0, T):
    n = int(round(T / DT))
    X = X0.copy()
    dd = []
    for _ in range(n + 1):
        dd.append(np.linalg.norm(X[0][:, None, :] - X[1][None, :, :], axis=2).min())
        X = A3.rk4_3(X, DT)
    dd = np.array(dd)
    idx = []
    for i in np.argsort(dd):
        if all(abs(i - j) > n // 4 for j in idx):
            idx.append(i)
        if len(idx) == 2:
            break
    idx = sorted(idx)
    return [range(max(0, i - int(0.5 / DT)), min(n, i + int(0.5 / DT))) for i in idx], n


def adjudicate(Rm_, Tm_, label):
    k = float((Rm_ * Tm_).sum() / (Tm_ * Tm_).sum())
    r = float(np.sqrt(((Rm_ - k * Tm_) ** 2).sum() / (Rm_ ** 2).sum()))
    if k < 0:
        print(f"  {label}: kappa={k:.4f} WRONG-WAY -> KILL", flush=True)
        return "KILL"
    if r > 0.25:
        print(f"  {label}: kappa={k:.4f} r={r:.4f} -> KILL", flush=True)
        return "KILL"
    if r > 0.10:
        print(f"  {label}: kappa={k:.4f} r={r:.4f} -> GRAY", flush=True)
        return "GRAY"
    print(f"  {label}: kappa={k:.4f} r={r:.4f} -> HOLD", flush=True)
    return "HOLD"


def stage_a1():
    print("A1 m>=2 seed:", flush=True)
    R1, R2, T, _ = A3.shoot()
    X0 = A3.ring_state(R1, 0.0, R2, 0.0) + SEED * A3.basis_field(2, 0, 0, 0)
    wins, n = leap_windows(X0, T)
    ok, sb, sp = gate(X0, n)
    if not ok:
        return "GATE-FAIL"
    if not content_ok(sb, sp, 2):
        return "NON-FIRING"
    Rm_, Tm_ = fit_passages(X0, n, wins, 2)
    return adjudicate(Rm_.ravel(), Tm_.ravel(), "A1 m=2 shared-kappa")


def stage_a2():
    print("A2 single ring:", flush=True)
    X0 = A3.ring_state(1.0, 0.0, 50.0, 0.0) + SEED * A3.basis_field(1, 0, 0, 0)
    T = 4.088
    n = int(round(T / DT))
    w1 = range(0, n // 2)
    w2 = range(n // 2, n)
    ok, sb, sp = gate(X0, n)
    if not ok:
        return "GATE-FAIL"
    if not content_ok(sb, sp, 1):
        return "NON-FIRING"
    Rm_, Tm_ = fit_passages(X0, n, [w1, w2], 1)
    return adjudicate(Rm_.ravel(), Tm_.ravel(), "A2 m=1 shared-kappa")


def stage_a3():
    print("A3 lambda-family:", flush=True)
    R1, R2, T, _ = A3.shoot()
    X0 = A3.ring_state(R1, 0.0, R2, 0.0) + SEED * A3.basis_field(1, 0, 0, 0)
    wins, n = leap_windows(X0, T)
    ok, sb, sp = gate(X0, n)
    if not ok:
        return "GATE-FAIL"
    best = (None, None, 1e9)
    for lam in [0.0, 0.1, 0.25, 0.5, 0.75, 1.0]:
        Rm_, Tm_ = fit_passages(X0, n, wins, 1, lam)
        Rg, Tg = Rm_.ravel(), Tm_.ravel()
        k = float((Rg * Tg).sum() / (Tg * Tg).sum())
        r = float(np.sqrt(((Rg - k * Tg) ** 2).sum() / (Rg ** 2).sum()))
        print(f"  lam={lam}: kappa={k:.4f} r={r:.4f}", flush=True)
        if r < best[2]:
            best = (lam, k, r)
    lam, k, r = best
    if lam in (0.0, 1.0):
        print(f"  A3: lam-hat={lam} ON BOUND -> KILL", flush=True)
        return "KILL"
    if k < 0:
        print(f"  A3: kappa={k:.4f} WRONG-WAY -> KILL", flush=True)
        return "KILL"
    if r > 0.25:
        print(f"  A3: lam={lam} r={r:.4f} -> KILL", flush=True)
        return "KILL"
    if r > 0.10:
        print(f"  A3: lam={lam} r={r:.4f} -> GRAY", flush=True)
        return "GRAY"
    print(f"  A3: lam={lam} kappa={k:.4f} r={r:.4f} -> HOLD", flush=True)
    return "HOLD"


if __name__ == "__main__":
    t0 = time.time()
    st = sys.argv[1] if len(sys.argv) > 1 else "all"
    stages = {"a1": stage_a1, "a2": stage_a2, "a3": stage_a3}
    if st == "all":
        for name in ("a1", "a2", "a3"):
            v = stages[name]()
            print(f"{name.upper()} verdict: {v} ({time.time()-t0:.1f}s)", flush=True)
            if v == "HOLD":
                print("PROGRAM HALT: hold reported, downstream skipped.", flush=True)
                break
    else:
        v = stages[st]()
        print(f"{st.upper()} verdict: {v} ({time.time()-t0:.1f}s)", flush=True)

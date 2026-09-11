#!/usr/bin/env python3
"""M1 interaction-energy dynamics build (beacon 0144): Neumann E_12(t)
for tagged (constitutive aa(chi)) vs untagged leapfrog pairs, three
separations. Frozen 0141 falsifier: (a) indistinct / (b) non-Coulomb /
(c) same-sign kills; HOLD = single-q Coulomb <=20% + sign correct
(conditional phenomenology only). A3 code read-only (backbone only).
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
AA0 = 0.05
QIN = 0.2
N = A3.N
DT = A3.DT
DPH = A3.DPH
SEPS = [1.0, 1.185226, 1.4]
R1FIX = 0.773723
TWIN = 4.088


def self_V(R, aa):
    return G / (4 * np.pi * R) * (np.log(8 * R / aa) - 0.25)


def rhs_aa(X, aa1, aa2):
    V = np.zeros_like(X)
    aas = [aa1, aa2]
    for n in range(2):
        Xn = X[n]
        Rn = np.sqrt(Xn[:, 0] ** 2 + Xn[:, 1] ** 2).mean()
        V[n, :, 2] = self_V(Rn, aas[n])
        Xm = X[1 - n]
        aa12 = np.sqrt(aas[n] * aas[1 - n])
        d = Xn[:, None, :] - Xm[None, :, :]
        r2 = (d ** 2).sum(-1) + aa12 * aa12
        T = np.zeros_like(Xm)
        T[:, 0] = -np.sin(A3.PH)
        T[:, 1] = np.cos(A3.PH)
        Rm = np.sqrt(Xm[:, 0] ** 2 + Xm[:, 1] ** 2).mean()
        V[n] += G / (4 * np.pi) * np.sum(
            np.cross(T[None, :, :], d) / (r2 ** 1.5)[:, :, None], axis=1) * (Rm * DPH)
    return V


def rk4_aa(X, dt, aa1, aa2):
    def f(q):
        return rhs_aa(q, aa1, aa2)
    k1 = f(X)
    k2 = f(X + dt / 2 * k1)
    k3 = f(X + dt / 2 * k2)
    k4 = f(X + dt * k3)
    return X + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def neumann_E12(X, aa12):
    dl = np.roll(X, -1, axis=1) - X
    mid = (np.roll(X, -1, axis=1) + X) / 2
    d = mid[0][:, None, :] - mid[1][None, :, :]
    r = np.sqrt((d ** 2).sum(-1) + aa12 * aa12)
    w = np.einsum("ik,jk->ij", dl[0], dl[1]) / (4 * np.pi * r)
    return float(w.sum())


def run_traj(R2, aa1, aa2):
    X = A3.ring_state(R1FIX, 0.0, R2, 0.0)
    n = int(round(TWIN / DT))
    E, dd = [], []
    for _ in range(n + 1):
        aa12 = np.sqrt(aa1 * aa2)
        E.append(neumann_E12(X, aa12))
        dd.append(np.linalg.norm(X[0][:, None, :] - X[1][None, :, :], axis=2).min())
        X = rk4_aa(X, DT, aa1, aa2)
    return np.array(E), np.array(dd)


def main() -> None:
    t0 = time.time()
    ap = AA0
    aq = AA0 * (1 + QIN)
    am = AA0 * (1 - QIN)
    rows = []
    for R2 in SEPS:
        Eu, ddu = run_traj(R2, ap, ap)
        El, _ = run_traj(R2, aq, aq)
        Eo, _ = run_traj(R2, aq, am)
        # two passage windows from untagged distance minima
        idx = []
        for i in np.argsort(ddu):
            if all(abs(i - j) > len(ddu) // 4 for j in idx):
                idx.append(i)
            if len(idx) == 2:
                break
        idx = sorted(idx)
        wins = [range(max(0, i - int(0.5 / DT)), min(len(ddu), i + int(0.5 / DT)))
                for i in idx]
        for pi, w in enumerate(wins):
            w = list(w)
            su = float(np.mean(Eu[w]))
            dm = float(np.mean(ddu[w]))
            sl = float(np.mean(El[w]) - np.mean(Eu[w]))
            so = float(np.mean(Eo[w]) - np.mean(Eu[w]))
            # passage noise scale: untagged std within window
            nz = float(np.std(Eu[w])) + 1e-300
            rows.append((dm, pi, su, sl, so, nz))
            print(f"d={dm:.3f} pass{pi}: base={su:.4f} dlike={sl:+.4e} "
                  f"dopp={so:+.4e} noise={nz:.1e}", flush=True)
    print("--- fit ---", flush=True)
    # (b) single-q Coulomb form: shifts vs measured 1/distance
    ds = np.array([1.0 / r[0] for r in rows])
    yl = np.array([r[3] for r in rows])
    yo = np.array([r[4] for r in rows])
    A = float((yl * ds).sum() / (ds * ds).sum())
    rl = float(np.sqrt(((yl - A * ds) ** 2).sum() / (yl ** 2).sum()))
    # (c) opposite-sign check: same A must fit -yo
    ro = float(np.sqrt(((yo + A * ds) ** 2).sum() / (yo ** 2).sum()))
    # (a) indistinctness: shift vs passage noise
    nzmax = max(r[5] for r in rows)
    sig = min(abs(float(np.mean([r[3] for r in rows]))),
              abs(float(np.mean([r[4] for r in rows])))) / nzmax
    print(f"single-A fit: A={A:.4e} like-res={rl:.3f} opp-res={ro:.3f} "
          f"significance={sig:.1f}x-noise", flush=True)
    if sig < 3.0:
        print("F1 verdict: KILL (a) indistinct", flush=True)
    elif rl > 0.20 or ro > 0.20:
        print("F1 verdict: KILL (b) non-Coulomb", flush=True)
    elif np.sign(np.mean([r[3] for r in rows])) == np.sign(np.mean([r[4] for r in rows])):
        print("F1 verdict: KILL (c) same-sign", flush=True)
    else:
        print("F1 verdict: HOLD (conditional phenomenology only)", flush=True)
    print(f"({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()

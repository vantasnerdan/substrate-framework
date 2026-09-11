#!/usr/bin/env python3
"""C3 shear-scan battery (beacon 0153): re-shoot under imposed shear S,
m-Floquet growth per S, available shear from trust-r3, S*-vs-available
adjudication. Frozen design.md. A3 code read-only (local copies noted).
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


def rhs_s(X, S):
    V = A3.rhs3(X)
    V[:, :, 2] += S * X[:, :, 0]
    return V


def rk4_s(X, dt, S):
    def f(q):
        return rhs_s(q, S)
    k1 = f(X)
    k2 = f(X + dt / 2 * k1)
    k3 = f(X + dt / 2 * k2)
    k4 = f(X + dt * k3)
    return X + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)


def flow_s(X0, T, S, dt=DT):
    X = X0.copy()
    n = int(round(T / dt))
    for _ in range(n):
        X = rk4_s(X, dt, S)
    frac = T - n * dt
    if frac > 1e-12:
        X = rk4_s(X, frac, S)
    return X


def shoot_s(R1g, R2g, Tg, S, quiet=True):
    y = np.array([R1g, R2g, Tg])

    def res3(yy):
        X0 = A3.ring_state(yy[0], 0.0, yy[1], 0.0)
        XT = flow_s(X0, yy[2], S)
        Sh = A3.shape4(XT)
        return np.array([Sh[0] - yy[0], Sh[2] - yy[1], Sh[1] - Sh[3]])

    r = res3(y)
    for _ in range(12):
        if np.linalg.norm(r) < 1e-10:
            break
        J = np.zeros((3, 3))
        e = np.array([1e-6, 1e-6, 1e-7])
        for j in range(3):
            dy = np.zeros(3)
            dy[j] = e[j]
            J[:, j] = (res3(y + dy) - r) / e[j]
        dy, *_ = np.linalg.lstsq(J, -r, rcond=None)
        y = y + dy
        r = res3(y)
    if not quiet:
        print(f"  shoot S={S}: R1={y[0]:.6f} R2={y[1]:.6f} T={y[2]:.5f} "
              f"|res|={np.linalg.norm(r):.2e}", flush=True)
    return y[0], y[1], y[2], np.linalg.norm(r)


def growth_s(X0, T, S, mmax=2, eps=1e-6):
    XB = flow_s(X0, T, S)
    pB = A3.project(XB - A3.axisym_of(XB))
    out = {}
    for m in range(1, mmax + 1):
        nd = A3.ndirs(m)
        M = np.zeros((nd, nd))
        for i in range(nd):
            F = A3.basis_field(m, *A3.dir_index(m, i))
            XT = flow_s(X0 + eps * F, T, S)
            col = (A3.project(XT - A3.axisym_of(XT))[m] - pB[m]) / eps
            M[:, i] = col[:nd]
        out[m] = float(np.abs(np.linalg.eigvals(M)).max())
    return out


def available_shear():
    sys.path.insert(0, os.path.join(ATT, "0117-beacon-member"))
    sys.path.insert(0, os.path.join(ATT, "0123-beacon-fitted"))
    import build_member as B
    from x18_probes import nodal_gradient
    from skfem import Basis, ElementTriP1
    r3 = np.load(os.path.join(ATT, "0117-beacon-member", "member-trust-r3.npz"))
    mesh = B.build_mesh(6.0, 3.0, 40, 20)
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0].copy()
    u = r3["u"].copy()
    g = nodal_gradient(mesh, u)
    ur, uz = g[0], g[1]
    r = np.maximum(rn, 0.15)
    Ur = -uz / r
    Uz = ur / r
    gr = nodal_gradient(mesh, Ur)
    gz = nodal_gradient(mesh, Uz)
    # axisymmetric strain: e_rr, e_zz, e_rz, e_pp=Ur/r
    err, erz1 = gr[0], gr[1]
    ezr, ezz = gz[0], gz[1]
    srate = np.sqrt(err**2 + ezz**2 + 0.5 * (erz1 + ezr)**2 + (Ur / r)**2)
    m = rn >= 0.15
    return float(srate[m].max()), float(np.abs(Ur).max() / 1.0)


def main() -> None:
    t0 = time.time()
    Sav, Uscale = available_shear()
    print(f"available shear S_av = {Sav:.4f} (U-scale {Uscale:.3f})", flush=True)
    R1g, R2g, Tg = 0.773723, 1.185226, 4.08800
    results = {}
    for S in (0.0, 0.05, 0.1, 0.2, 0.4, 0.8):
        R1, R2, T, rn_ = shoot_s(R1g, R2g, Tg, S)
        X0 = A3.ring_state(R1, 0.0, R2, 0.0)
        g = growth_s(X0, T, S)
        gm = max(g.values())
        results[S] = gm
        print(f"S={S}: " + " ".join(f"m{m}={v:.4f}" for m, v in g.items()) +
              f" max={gm:.4f} shoot-res={rn_:.1e}", flush=True)
        R1g, R2g, Tg = R1, R2, T
    Sstar = next((S for S in sorted(results) if results[S] <= 1.05), None)
    print(f"S* = {Sstar} (first S with growth <= 1.05)", flush=True)
    if Sstar is None:
        print("F1 verdict: DEAD (no stabilization in grid)", flush=True)
    elif Sstar == 0.0:
        print("F1 verdict: VACUOUS premise (stable at S=0; nothing to "
              "stabilize — C3 moot, not alive)", flush=True)
    elif Sstar > 3 * Sav:
        print("F1 verdict: DEAD (S* exceeds available)", flush=True)
    elif Sstar <= Sav:
        print("F1 verdict: ALIVE-conditional (background-supplied label "
              "permanent)", flush=True)
    else:
        print("F1 verdict: GRAY (measure better, no verdict)", flush=True)
    print(f"({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()

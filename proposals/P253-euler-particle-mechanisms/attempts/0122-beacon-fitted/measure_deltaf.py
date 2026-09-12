#!/usr/bin/env python3
"""P253/0122 R1 gate instrument (beacon): eigen-pipeline deltaF re-measure.

Drift R1: res NEVER gates (x15-60 amplification); the dF<=0.1 verdict
fires ONLY on this instrument's output. Includes R6 core/tail soft-share
decomposition. Usage: measure_deltaf.py [state.npz] (default: fitted).
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys

import numpy as np


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("state", nargs="?",
                    default="proposals/P253-euler-particle-mechanisms/"
                            "attempts/0122-beacon-fitted/member-fitted.npz")
    args = ap.parse_args()
    sys.path.insert(0, "proposals/P253-euler-particle-mechanisms/"
                       "attempts/0117-beacon-member")
    import build_member as B
    from skfem import Basis, ElementTriP1, BilinearForm, LinearForm, asm, MeshTri
    from skfem.helpers import dot, grad, inner
    att2 = "proposals/P253-euler-particle-mechanisms/attempts/0122-beacon-fitted"
    fm = np.load(f"{att2}/fitted-mesh.npz")
    d = np.load(args.state)
    mesh = MeshTri(fm["p"].T, fm["t"])
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0]
    zn = basis.doflocs[1]
    u, mu, c = d["u"], float(d["mu"]), float(d["c"])
    assert u.size == basis.N, (u.size, basis.N)

    @BilinearForm
    def stiff(a, b, w):
        return w.x[0] * dot(grad(a), grad(b))

    @BilinearForm
    def massr(a, b, w):
        return w.x[0] * inner(a, b)

    A = asm(stiff, basis).tocsr()
    Mr = asm(massr, basis)
    src = u - c * rn**2 / 2 - mu
    root = np.sqrt(src**2 + 1e-6)
    s = (src + root) / 2
    ds = 0.5 * (1 + src / root)
    f = B.EPS**-2 * s**B.P
    fi = basis.interpolator(f)

    @LinearForm
    def load(v, w):
        return (w.x[0] ** 3) * fi(w.x) * v

    R = A @ u - asm(load, basis)
    dd = B.dirichlet_dofs(basis, mesh, 6.0, 3.0)
    free = np.setdiff1d(np.arange(basis.N), dd)
    R[dd] = 0.0
    jf = 6 * B.EPS**-2 * s**5 * ds
    ji = basis.interpolator(jf)

    @BilinearForm
    def jacform(a, b, w):
        return (w.x[0] * dot(grad(a), grad(b))
                - (w.x[0] ** 3) * ji(w.x) * inner(a, b))

    Jf = asm(jacform, basis).tocsr()[free][:, free].toarray()
    w, V = np.linalg.eigh(Jf)
    cf = V.T @ R[free]
    order = np.argsort(np.abs(w))
    soft8 = order[:8]
    soft_share = float(np.sum(cf[soft8]**2) / np.sum(cf**2))
    du = np.zeros(basis.N)
    nz = np.abs(w) > 1e-12
    du[free] = V[:, nz] @ (cf[nz] / w[nz])
    h = rn * (6 * B.EPS**-2 * s**5 * ds) * du
    dF = float(np.sqrt(h @ (Mr @ h)))
    # R6: core/tail decomposition of the soft-share (tail = |z|>1.3 edge band)
    tail = np.abs(zn) > 1.0
    Vs = V[:, soft8]
    wts = (Vs**2).sum(axis=1)
    # map free-index weights back to full mesh for the mask
    wfull = np.zeros(basis.N)
    wfull[free] = wts
    tail_share = float(wfull[tail].sum() / wfull.sum())
    print(f"GATE dF = {dF:.4f} (gate: <= 0.1)", flush=True)
    print(f"soft8 share = {soft_share:.4f} "
          f"(R6 tail-mode weight = {tail_share:.4f})", flush=True)
    print(f"neg = {int((w < 0).sum())} min|.| = {np.abs(w).min():.4e}",
          flush=True)
    s_om = 5.9191
    print(f"margin: dlam <= {2 * s_om * dF:.2f} vs 11.2", flush=True)


if __name__ == "__main__":
    main()

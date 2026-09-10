#!/usr/bin/env python3
"""P253/0125 C-threshold sweep driver (beacon, banked per drift R3).

Exact nonlinear transfer |dQ| over scales of the gauged Newton step +
src=0-crossing counts -> rho* where condition C (<= 0.9) fires.
Reproduces the 0125 lemma sweep numbers (trust-r3 tensor default).
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys

import numpy as np


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", default="tensor")
    ap.add_argument("--state", default="proposals/P253-euler-particle-mechanisms/"
                    "attempts/0117-beacon-member/member-trust-r3.npz")
    ap.add_argument("--scales", default="1.0,0.5,0.25,0.125,0.0625,0.02,0.005,0.001")
    args = ap.parse_args()
    att = "proposals/P253-euler-particle-mechanisms/attempts"
    sys.path.insert(0, f"{att}/0111-beacon-ga-field")
    sys.path.insert(0, f"{att}/0117-beacon-member")
    sys.path.insert(0, f"{att}/0124-beacon-sharp")
    sys.path.insert(0, f"{att}/0123-beacon-fitted")
    import build_member as B
    from sharp_dQ import build_F, Q_of_F
    from x18_probes import nodal_gradient
    from ga_pipeline import leray, make_grid
    from skfem import Basis, ElementTriP1, BilinearForm, LinearForm, asm, MeshTri
    from skfem.helpers import dot, grad, inner

    if args.mesh == "fitted":
        att2 = f"{att}/0122-beacon-fitted"
        fm = np.load(f"{att2}/fitted-mesh.npz")
        mesh = MeshTri(fm["p"].T, fm["t"])
    else:
        mesh = B.build_mesh(6.0, 3.0, 40, 20)
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0].copy()
    zn = basis.doflocs[1].copy()
    d = np.load(args.state)
    assert d["u"].size == basis.N, (d["u"].size, basis.N)
    u, mu, c = d["u"].copy(), float(d["mu"]), float(d["c"])
    src = u - c * rn**2 / 2 - mu
    root = np.sqrt(src**2 + 1e-6)
    s = (src + root) / 2
    zeta = B.EPS**-2 * s**B.P
    n3, L3 = 64, 8.0
    F, grid = build_F(basis, rn, zn, rn * zeta, n3, L3, make_grid)
    Q, _, _ = Q_of_F(F, grid, leray)

    ds = 0.5 * (1 + src / root)
    A = asm(_stiff(), basis).tocsr()
    fi = basis.interpolator(B.EPS**-2 * s**B.P)

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
    nz = np.abs(w) > 1e-12
    du = np.zeros(basis.N)
    du[free] = V[:, nz] @ (cf[nz] / w[nz])
    gg = nodal_gradient(mesh, u)
    Qt, _ = np.linalg.qr(np.stack([gg[0][free], gg[1][free]], axis=1))
    du[free] = du[free] - Qt @ (Qt.T @ du[free])
    for sc in [float(x) for x in args.scales.split(",")]:
        dus = sc * du
        s2 = u + dus - c * rn**2 / 2 - mu
        r2 = np.sqrt(s2**2 + 1e-6)
        dz = B.EPS**-2 * (((s2 + r2) / 2)**B.P - s**B.P)
        dFe, _ = build_F(basis, rn, zn, rn * dz, n3, L3, make_grid)
        Q2, _, _ = Q_of_F(F + dFe, grid, leray)
        g = float(np.abs(np.linalg.eigvalsh(Q2 - Q)).max())
        cross = int(((src > 0) != (s2 > 0)).sum())
        print(f"scale={sc}: exact|dQ|={g:.4f} ||e||_oo={np.abs(dus).max():.4f} "
              f"xnodes={cross}", flush=True)


def _stiff():
    from skfem import BilinearForm
    from skfem.helpers import dot, grad

    @BilinearForm
    def stiff(a, b, w):
        return w.x[0] * dot(grad(a), grad(b))

    return stiff


if __name__ == "__main__":
    main()

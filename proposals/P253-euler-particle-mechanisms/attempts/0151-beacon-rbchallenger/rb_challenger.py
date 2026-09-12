#!/usr/bin/env python3
"""P253/0129 R-B mechanism experiment (beacon): semi-Lagrangian front advection.

Moves the MESH with the front (values ride nodes), reassembles Q on the
moved mesh, compares vs printed exacts. No continuity claim, no domain
linearization (G-Q3/Q4); verdict conditional on scale (G-Q1); beat-linear
trigger (G-Q5); exact-outranks-DWR (G-Q8).
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys

import numpy as np


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--scales", default="1.0,0.5,0.25")
    ap.add_argument("--ext", default="idw",
                    choices=["idw", "harm", "bump"])
    args = ap.parse_args()
    att = "proposals/P253-euler-particle-mechanisms/attempts"
    sys.path.insert(0, f"{att}/0111-beacon-ga-field")
    sys.path.insert(0, f"{att}/0117-beacon-member")
    sys.path.insert(0, f"{att}/0124-beacon-sharp")
    sys.path.insert(0, f"{att}/0123-beacon-fitted")
    sys.path.insert(0, f"{att}/0122-beacon-fitted")
    import build_member as B
    from sharp_dQ import build_F, Q_of_F
    from x18_probes import nodal_gradient
    from fitted_mesh import extract_contour
    from ga_pipeline import leray, make_grid
    from skfem import Basis, ElementTriP1, BilinearForm, LinearForm, asm, MeshTri
    from skfem.helpers import dot, grad, inner

    r3 = np.load(f"{att}/0117-beacon-member/member-trust-r3.npz")
    mesh = B.build_mesh(6.0, 3.0, 40, 20)
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0].copy()
    zn = basis.doflocs[1].copy()
    u, mu, c = r3["u"].copy(), float(r3["mu"]), float(r3["c"])
    src = u - c * rn**2 / 2 - mu
    root = np.sqrt(src**2 + 1e-6)
    s = (src + root) / 2
    ds = 0.5 * (1 + src / root)
    zeta = B.EPS**-2 * s**B.P
    n3, L3 = 64, 8.0
    F, grid = build_F(basis, rn, zn, rn * zeta, n3, L3, make_grid)
    Q, _, _ = Q_of_F(F, grid, leray)

    # gauged Newton step (0123 S1 pattern)
    @BilinearForm
    def stiff(a, b, w):
        return w.x[0] * dot(grad(a), grad(b))

    A = asm(stiff, basis).tocsr()
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
    fu = basis.interpolator(u)
    # tensor mesh covers z in [0,3]: upper-half contour only (no mirror;
    # half-field revolve convention shared by all banked Q numbers)
    cpts = np.asarray(extract_contour(fu, mu, c))
    cen = cpts.mean(axis=0)
    nxt = np.roll(cpts, -1, axis=0)
    prv = np.roll(cpts, 1, axis=0)
    tang = nxt - prv
    nrm = np.stack([tang[:, 1], -tang[:, 0]], axis=1)
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-300
    out = np.sign(((cpts - cen) * nrm).sum(axis=1)).reshape(-1, 1)
    nrm = nrm * np.where(out == 0, 1.0, out)
    # level-set normal displacement: dn = -du/|grad src| along grad direction
    gsrc = nodal_gradient(mesh, src)
    iu = basis.interpolator(du)
    igx = basis.interpolator(gsrc[0])
    igz = basis.interpolator(gsrc[1])
    duc = iu(cpts.T)
    gxv, gzv = igx(cpts.T), igz(cpts.T)
    gnorm = np.sqrt(gxv**2 + gzv**2) + 1e-300
    # project grad onto polygon normal sign for coherent orientation
    sgn = np.sign(gxv * nrm[:, 0] + gzv * nrm[:, 1]).reshape(-1, 1)
    sgn = np.where(sgn == 0, 1.0, sgn)
    dn = -(duc / gnorm) * sgn.ravel()
    Dpts = (dn.reshape(-1, 1)) * nrm
    # IDW extension, zero on outer boundary (0065-fence: smooth).
    # mesh.p is (2,N); work node-major via pn.
    p = mesh.p.copy()
    pn = p.T
    diff = pn[:, None, :] - cpts[None, :, :]
    wgt = 1.0 / (np.sum(diff**2, axis=2) + 1e-6)
    W = wgt / wgt.sum(axis=1, keepdims=True)
    edge = np.minimum.reduce([pn[:, 0], 6.0 - pn[:, 0],
                              pn[:, 1] + 3.0, 3.0 - pn[:, 1]])
    cut = np.clip(edge / 0.5, 0.0, 1.0).reshape(-1, 1)
    if args.ext == "idw":
        D = (W @ Dpts) * cut
    elif args.ext == "bump":
        dmin = np.sqrt((diff**2).sum(-1).min(axis=1)).reshape(-1, 1)
        D = (W @ Dpts) * np.exp(-((dmin / 0.3) ** 2)) * cut
    else:
        from skfem import BilinearForm, asm
        from skfem.helpers import dot, grad
        from scipy.sparse.linalg import spsolve

        @BilinearForm
        def lap(a, b, w):
            return dot(grad(a), grad(b))

        K = asm(lap, basis).tocsr()
        dix = ((diff**2).sum(-1)).argmin(axis=0)
        bnd = edge < 1e-9
        fixed = np.zeros(basis.N, dtype=bool)
        fixed[bnd] = True
        dval = np.zeros((basis.N, 2))
        dval[dix] = Dpts
        fixed[dix] = True
        free = ~fixed
        D = np.zeros((basis.N, 2))
        D[fixed] = dval[fixed]
        Kff = K[free][:, free].tocsr()
        Kfd = K[free][:, fixed]
        for c in range(2):
            D[free, c] = spsolve(Kff, -Kfd @ dval[fixed, c])
        D = D * cut
    ext_e = float(np.abs(D).max())
    print(f"extension {args.ext}: max|D|/unit-du = {ext_e:.4f}, "
          f"cpts = {len(cpts)}", flush=True)
    # native winding signs: detJ gate is RELATIVE (tensor mesh ships
    # mixed-wound 800/1600 in this convention; skfem orients internally)
    _t = mesh.t
    _v0 = p[:, _t[0]] - p[:, _t[2]]
    _v1 = p[:, _t[1]] - p[:, _t[2]]
    _s0 = np.sign(0.5 * (_v0[0] * _v1[1] - _v0[1] * _v1[0]))
    exacts = {1.0: 9.1750, 0.5: 3.4895, 0.25: 1.5141}
    lin_g1 = 5.2052
    import copy as _copy
    for sc in [float(x) for x in args.scales.split(",")]:
        pmn = pn + sc * D
        pm = pmn.T
        t = mesh.t
        v0 = pm[:, t[0]] - pm[:, t[2]]
        v1 = pm[:, t[1]] - pm[:, t[2]]
        areas = 0.5 * (v0[0] * v1[1] - v0[1] * v1[0])
        bad = int((np.sign(areas) != _s0).sum())
        if bad:
            print(f"scale={sc}: detJ GATE FAIL ({bad} flipped) — dropped",
                  flush=True)
            continue
        m2 = MeshTri(pm.copy(), mesh.t.copy())
        b2 = Basis(m2, ElementTriP1())
        # moved coordinates
        F2, _ = build_F(b2, m2.p[0].copy(), m2.p[1].copy(), rn * zeta,
                        n3, L3, make_grid)
        Q2, _, _ = Q_of_F(F2, grid, leray)
        shape_err = float(np.abs(np.linalg.eigvalsh(Q2 - Q)).max())
        ex = exacts[sc]
        lin = lin_g1 * sc
        verdict = ("HOLD-PT" if (shape_err <= 0.2 * ex and shape_err <= lin)
                   else "MISS-PT")
        print(f"scale={sc}: shape-err={shape_err:.4f} exact={ex:.4f} "
              f"linear={lin:.4f} [{verdict}]", flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""P253/0123 x18 decision probes (beacon): S1 null-projection, S2 layer-R, S3 sharpness.

Reads 0122 fitted state (+ 0117 trust r2/r3 on the tensor mesh for S3).
Prints the three probe verdicts. No solves, no mesh builds.
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys

import numpy as np


def load_fitted():
    sys.path.insert(0, "proposals/P253-euler-particle-mechanisms/"
                       "attempts/0117-beacon-member")
    import build_member as B
    from skfem import Basis, ElementTriP1, MeshTri
    att2 = "proposals/P253-euler-particle-mechanisms/attempts/0122-beacon-fitted"
    fm = np.load(f"{att2}/fitted-mesh.npz")
    d = np.load(f"{att2}/member-fitted.npz")
    mesh = MeshTri(fm["p"].T, fm["t"])
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0].copy()
    u, mu, c = d["u"].copy(), float(d["mu"]), float(d["c"])
    return B, basis, mesh, rn, u, mu, c


def residual_and_jac(B, basis, mesh, rn, u, mu, c):
    from skfem import BilinearForm, LinearForm, asm
    from skfem.helpers import dot, grad, inner

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
    return R, Jf, free, Mr, src, jf

def nodal_gradient(mesh, u):
    """Exact P1 gradient, area-weighted to nodes (numpy only)."""
    p, t = mesh.p, mesh.t
    g = np.zeros((2, p.shape[1]))
    wsum = np.zeros(p.shape[1])
    for tri in t.T:
        xy = p[:, tri]
        M = np.stack([xy[:, 1] - xy[:, 0], xy[:, 2] - xy[:, 0]])
        ge = np.linalg.solve(M, u[tri[1:]] - u[tri[0]])
        area = abs(np.linalg.det(M)) / 2
        for i in tri:
            g[:, i] += area * ge
            wsum[i] += area
    return g / wsum

def main() -> None:
    B, basis, mesh, rn, u, mu, c = load_fitted()
    R, Jf, free, Mr, src, jf = residual_and_jac(B, basis, mesh, rn, u, mu, c)
    w, V = np.linalg.eigh(Jf)
    order = np.argsort(np.abs(w))
    print(f"spectrum: min|.|={np.abs(w).min():.4e} neg={(w < 0).sum()}",
          flush=True)

    # ---- S1: translation-subspace projection ----
    # discrete translation modes: nodal gradients of u along r and z
    gg = nodal_gradient(mesh, u)
    Q, _ = np.linalg.qr(np.stack([gg[0][free], gg[1][free]], axis=1))
    v0 = V[:, order[0]]
    ov = float(np.sqrt((Q.T @ v0) @ (Q.T @ v0)))
    print(f"S1: softest-mode |overlap| with translation plane = {ov:.4f}",
          flush=True)
    # project translation subspace out of the Newton step, re-measure dF
    cf = V.T @ R[free]
    nz = np.abs(w) > 1e-12
    du = np.zeros(basis.N)
    du[free] = V[:, nz] @ (cf[nz] / w[nz])
    h = rn * jf * du
    dF = float(np.sqrt(h @ (Mr @ h)))
    du_f = du[free] - Q @ (Q.T @ du[free])
    du2 = np.zeros(basis.N)
    du2[free] = du_f
    h2 = rn * jf * du2
    dF_notrans = float(np.sqrt(h2 @ (Mr @ h2)))
    print(f"S1: dF = {dF:.4f} -> without-translation {dF_notrans:.4f}",
          flush=True)
    print(f"S1: ||du||_2 = {float(np.linalg.norm(du_f)):.4e} plain, "
          f"{float(np.sqrt(du2 @ (Mr @ du2))):.4e} r-weighted", flush=True)

    # ---- S2: layer-R decomposition ----
    delta = 2 * B.EPS
    Rfull = np.zeros(basis.N)
    Rfull[free] = R[free]
    # lumped nodal R^2 proxy via mass-lumped weights
    mdiag = np.asarray(Mr.sum(axis=1)).ravel()
    R2 = (Rfull**2) * mdiag
    layer = np.abs(src) < delta
    print(f"S2: layer |R|^2 share = {R2[layer].sum() / R2.sum():.4f} "
          f"(band |src|<{delta:.2f})", flush=True)
    print(f"S2: max|jf| = {np.abs(jf).max():.1f}, "
          f"layer nodes = {int(layer.sum())}/{basis.N}", flush=True)

    # ---- S3: bound sharpness on trust r2 -> r3 (tensor mesh) ----
    att7 = "proposals/P253-euler-particle-mechanisms/attempts/0117-beacon-member"
    from skfem import Basis as _B, ElementTriP1 as _E
    r2 = np.load(f"{att7}/member-recover-exploratory.npz")
    r3 = np.load(f"{att7}/member-trust-exploratory.npz")
    mesh_t = B.build_mesh(6.0, 3.0, 40, 20)
    bt = _B(mesh_t, _E())
    rnt = bt.doflocs[0].copy()
    du23 = r3["u"] - r2["u"]
    # F-perturbation of the r2->r3 step under r3 weights
    u3, mu3, c3 = r3["u"], float(r3["mu"]), float(r3["c"])
    src3 = u3 - c3 * rnt**2 / 2 - mu3
    root3 = np.sqrt(src3**2 + 1e-6)
    jf3 = 6 * B.EPS**-2 * ((src3 + root3) / 2)**5 * 0.5 * (1 + src3 / root3)
    from skfem import BilinearForm as _BF, asm as _asm
    from skfem.helpers import inner as _inner

    @_BF
    def _m(a, b, w):
        return w.x[0] * _inner(a, b)

    Mt = _asm(_m, bt)
    hh = rnt * jf3 * du23
    dF23 = float(np.sqrt(hh @ (Mt @ hh)))
    print(f"S3: r2->r3 dF_step = {dF23:.4f}, "
          f"bound-predicted |dlam| <= {2 * 5.9191 * dF23:.2f}", flush=True)
    print("S3: compare vs actual eigengap drift across trust rounds "
          "(see trust-report: min|.| 0.016 coarse-stable -> decision by ratio)",
          flush=True)

def s3_full() -> None:
    """S3-full: ||J(r3)-J(r2)||_2 vs the 2*s*dF bound (tensor mesh)."""
    import sys
    sys.path.insert(0, "proposals/P253-euler-particle-mechanisms/"
                       "attempts/0117-beacon-member")
    import build_member as B
    from skfem import Basis, ElementTriP1
    att7 = "proposals/P253-euler-particle-mechanisms/attempts/0117-beacon-member"
    r2 = np.load(f"{att7}/member-recover-exploratory.npz")
    r3 = np.load(f"{att7}/member-trust-exploratory.npz")
    mesh = B.build_mesh(6.0, 3.0, 40, 20)
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0].copy()
    from skfem import BilinearForm, asm
    from skfem.helpers import dot, grad, inner

    def free_jac(u, mu, c):
        src = u - c * rn**2 / 2 - mu
        root = np.sqrt(src**2 + 1e-6)
        ds = 0.5 * (1 + src / root)
        s = (src + root) / 2
        ji = basis.interpolator(6 * B.EPS**-2 * s**5 * ds)

        @BilinearForm
        def jf(a, b, w):
            return (w.x[0] * dot(grad(a), grad(b))
                    - (w.x[0] ** 3) * ji(w.x) * inner(a, b))

        J = asm(jf, basis).tocsr()
        dd = B.dirichlet_dofs(basis, mesh, 6.0, 3.0)
        free = np.setdiff1d(np.arange(basis.N), dd)
        return J[free][:, free].toarray()

    J2 = free_jac(r2["u"], float(r2["mu"]), float(r2["c"]))
    J3 = free_jac(r3["u"], float(r3["mu"]), float(r3["c"]))
    dJ2 = float(np.abs(np.linalg.eigvalsh(J3 - J2)).max())
    print(f"S3-full: ||dJ||_2 = {dJ2:.4f} vs bound 1.13 "
          f"(looseness x{1.13 / max(dJ2, 1e-12):.1f})", flush=True)

if __name__ == "__main__":
    main()
    s3_full()

#!/usr/bin/env python3
"""0138 one-assembly decider (beacon): Q at u_h + du_central (trust-r3).

Expect-dead prediction: lam_min ~= 2 (worst-case arithmetic).
Revive: lam_min >= 5. Either outcome banked. Tested pipeline only.
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys

import numpy as np


def main() -> None:
    att = "proposals/P253-euler-particle-mechanisms/attempts"
    sys.path.insert(0, f"{att}/0111-beacon-ga-field")
    sys.path.insert(0, f"{att}/0117-beacon-member")
    sys.path.insert(0, f"{att}/0124-beacon-sharp")
    sys.path.insert(0, f"{att}/0123-beacon-fitted")
    import build_member as B
    from sharp_dQ import build_F, Q_of_F
    from x18_probes import nodal_gradient
    from ga_pipeline import leray, make_grid
    from skfem import Basis, ElementTriP1, BilinearForm, LinearForm, asm
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
    n3, L3 = 64, 8.0
    F, grid = build_F(basis, rn, zn, rn * B.EPS**-2 * s**B.P, n3, L3, make_grid)
    Q, _, _ = Q_of_F(F, grid, leray)
    lam0 = float(np.linalg.eigvalsh(Q)[0])
    print(f"Q(u_h) lam_min = {lam0:.4f}", flush=True)

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

    u2 = u + du
    src2 = u2 - c * rn**2 / 2 - mu
    root2 = np.sqrt(src2**2 + 1e-6)
    s2 = (src2 + root2) / 2
    F2, _ = build_F(basis, rn, zn, rn * B.EPS**-2 * s2**B.P, n3, L3, make_grid)
    Q2, _, _ = Q_of_F(F2, grid, leray)
    lam2 = np.linalg.eigvalsh(Q2)
    g = float(lam2[0])
    print(f"Q(u_h+du) eigenvalues: {lam2}", flush=True)
    print(f"DECIDER lam_min = {g:.4f} (base {lam0:.4f}; expect-dead ~= 2; "
          f"revive >= 5.0)", flush=True)
    if g >= 5.0:
        print("VERDICT: REVIVE — alignment favorable, run amended survey",
              flush=True)
    else:
        print("VERDICT: CONFIRM — prediction converts to measurement",
              flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""P253/0122 + IDEA-05: bank soft-aware solver mode shapes (beacon).

Dense eigendecomposition of the bordered-PDE Jacobian at the trust
state (coarse 40x20); saves the 8-mode soft cluster with shape metrics
(dipole moments, core/tail weight) and overlap vs the S9 tail direction
(-z elongation). Diagnostic only: tells the (a) build WHICH modes to
deflate/target; never a QOI input (per IDEA-05 terms).
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys

import numpy as np


def main() -> None:
    sys.path.insert(0, "proposals/P253-euler-particle-mechanisms/"
                       "attempts/0117-beacon-member")
    import build_member as B
    from skfem import Basis, ElementTriP1, BilinearForm, asm
    from skfem.helpers import dot, grad, inner
    mesh = B.build_mesh(6.0, 3.0, 40, 20)
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0]
    zn = basis.doflocs[1]
    att7 = "proposals/P253-euler-particle-mechanisms/attempts/0117-beacon-member"
    att2 = "proposals/P253-euler-particle-mechanisms/attempts/0122-beacon-fitted"
    d = np.load(f"{att7}/member-trust-r3.npz")
    u, mu, c = d["u"], float(d["mu"]), float(d["c"])

    src = u - c * rn**2 / 2 - mu
    root = np.sqrt(src**2 + 1e-6)
    s = (src + root) / 2
    ds = 0.5 * (1 + src / root)
    jf = 6 * B.EPS**-2 * s**5 * ds
    ji = basis.interpolator(jf)

    @BilinearForm
    def jacform(a, b, w):
        return (w.x[0] * dot(grad(a), grad(b))
                - (w.x[0] ** 3) * ji(w.x) * inner(a, b))

    J = asm(jacform, basis).tocsr()
    dd = B.dirichlet_dofs(basis, mesh, 6.0, 3.0)
    free = np.setdiff1d(np.arange(basis.N), dd)
    w, V = np.linalg.eigh(J[free][:, free].toarray())
    order = np.argsort(np.abs(w))[:8]
    modes = np.zeros((basis.N, 8))
    modes[free] = V[:, order]
    # shape metrics: z-dipole (tail direction), core weight r<2
    dip = (modes * zn[:, None]).sum(axis=0) / np.abs(modes).sum(axis=0)
    core = (np.abs(modes[rn < 2.0]).sum(axis=0)
            / np.abs(modes).sum(axis=0))
    print("soft evals:", w[order], flush=True)
    print("z-dipoles:", dip, flush=True)
    print("core weight:", core, flush=True)
    np.savez(f"{att2}/soft-modes.npz", modes=modes,
             evals=w[order], dipoles=dip, core_weight=core,
             mu=mu, c=c)
    print("banked soft-modes.npz (diagnostic only)", flush=True)


if __name__ == "__main__":
    main()

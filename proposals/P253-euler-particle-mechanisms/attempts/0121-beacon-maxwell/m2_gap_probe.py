#!/usr/bin/env python3
"""P253/0121 M2 gap probe (beacon, ARCHIVED from session transcript).

Dense Jacobian spectrum at the trust state; ARPACK-SM failure note.
Rerun to verify: matches transcript values below.
Transcript: 1 negative (-0.506), soft cluster 0.016-0.075 (8 modes),
bulk to 43.2; ARPACK: No convergence, 8001 iters, 0/6.
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
    d = np.load("proposals/P253-euler-particle-mechanisms/attempts/"
                "0117-beacon-member/member-trust-r3.npz")
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
    w = np.linalg.eigvalsh(J[free][:, free].toarray())
    print("n =", free.size)
    print("10 smallest:", w[:10])
    print("negatives:", int((w < 0).sum()))
    print("10 largest:", w[-10:])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""R-B front/bulk decomposition challenger (beacon 0147): where does the
TRUE exact delta-zeta live? Push band-only and bulk-only parts of the
exact nonlinear field perturbation separately through the tested Q
pipeline. No mesh motion, no profile reconstruction, no winding.
Front explains <=> ||dQ_front|| ~= ||dQ_exact|| (same trigger bars).
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
    from sharp_dQ import build_F, Q_of_F
    from ga_pipeline import leray, make_grid
    import build_member as B
    from skfem import Basis, ElementTriP1

    r3 = np.load(f"{att}/0117-beacon-member/member-trust-r3.npz")
    mesh = B.build_mesh(6.0, 3.0, 40, 20)
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0].copy()
    zn = basis.doflocs[1].copy()
    # gauged du: reuse sweep numbers pattern (recompute compactly via J)
    from skfem import BilinearForm, LinearForm, asm
    from skfem.helpers import dot, grad, inner
    import sys as _s
    _s.path.insert(0, f"{att}/0123-beacon-fitted")
    from x18_probes import nodal_gradient

    u, mu, c = r3["u"].copy(), float(r3["mu"]), float(r3["c"])
    src = u - c * rn**2 / 2 - mu
    root = np.sqrt(src**2 + 1e-6)
    s = (src + root) / 2
    ds = 0.5 * (1 + src / root)
    n3, L3 = 64, 8.0
    F, grid = build_F(basis, rn, zn, rn * B.EPS**-2 * s**B.P, n3, L3, make_grid)
    Q, _, _ = Q_of_F(F, grid, leray)

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

    from sharp_dQ import build_F as _bf  # noqa (same fn, explicit)
    for sc, ex, lin in ((1.0, 9.1750, 5.2052), (0.5, 3.4895, 2.6026),
                        (0.25, 1.5141, 1.3013)):
        dus = sc * du
        s2 = u + dus - c * rn**2 / 2 - mu
        r2 = np.sqrt(s2**2 + 1e-6)
        dz = B.EPS**-2 * (((s2 + r2) / 2)**B.P - s**B.P)
        for delta in (0.1, 0.3, 1.0):
            band = np.abs(src) < delta
            eshare = float((dz[band]**2).sum() / (np.linalg.norm(dz) ** 2))
            dFe, _ = _bf(basis, rn, zn, rn * dz * band, n3, L3, make_grid)
            Qe, _, _ = Q_of_F(F + dFe, grid, leray)
            gf = float(np.abs(np.linalg.eigvalsh(Qe - Q)).max())
            print(f"scale={sc} band={delta}: dz-share={eshare:.4f} "
                  f"front|dQ|={gf:.4f} (exact {ex:.4f})", flush=True)


if __name__ == "__main__":
    main()

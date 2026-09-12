#!/usr/bin/env python3
"""R-B Eulerian front-isolation challenger (beacon 0147): NO mesh motion.

Hypothesis made exact: Q response = front POSITION change with FROZEN
cross-front profile s(d) (empirical bins) + exact original far field,
reassembled on the FIXED mesh. Same frozen trigger. If HOLD: front
motion vindicated (STOP reopens). If MISS: STOP stands fortified
(IDW/harm/profile-frozen converge; bump exposed as bulk-freeze artifact).
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys

import numpy as np


def seg_dist(px, pz, ax, az, bx, bz):
    abx, abz = bx - ax, bz - az
    denom = abx * abx + abz * abz + 1e-300
    t = np.clip(((px - ax) * abx + (pz - az) * abz) / denom, 0.0, 1.0)
    dx, dz = px - (ax + t * abx), pz - (az + t * abz)
    return np.sqrt(dx * dx + dz * dz)


def poly_dist(rn, zn, cpts):
    d = np.full_like(rn, 1e9)
    for i in range(len(cpts)):
        a, b = cpts[i], cpts[(i + 1) % len(cpts)]
        dd = seg_dist(rn, zn, a[0], a[1], b[0], b[1])
        d = np.minimum(d, dd)
    return d


def winding_inside(rn, zn, cpts):
    inside = np.zeros_like(rn, dtype=bool)
    x, y = rn, zn
    for i in range(len(cpts)):
        x1, y1 = cpts[i]
        x2, y2 = cpts[(i + 1) % len(cpts)]
        cond = ((y1 > y) != (y2 > y))
        xinters = (x2 - x1) * (y - y1) / (y2 - y1 + 1e-300) + x1
        inside ^= cond & (x < xinters)
    return inside


def main() -> None:
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
    n3, L3 = 64, 8.0
    F, grid = build_F(basis, rn, zn, rn * B.EPS**-2 * s**B.P, n3, L3, make_grid)
    Q, _, _ = Q_of_F(F, grid, leray)

    # gauged du (0123 S1 pattern, condensed)
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
    ds = 0.5 * (1 + src / root)
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

    # contour + normal displacement (upper half; tensor mesh z in [0,3])
    fu = basis.interpolator(u)
    cpts = np.asarray(extract_contour(fu, mu, c))
    cen = cpts.mean(axis=0)
    tang = np.roll(cpts, -1, axis=0) - np.roll(cpts, 1, axis=0)
    nrm = np.stack([tang[:, 1], -tang[:, 0]], axis=1)
    nrm /= np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-300
    out = np.sign(((cpts - cen) * nrm).sum(axis=1)).reshape(-1, 1)
    nrm = nrm * np.where(out == 0, 1.0, out)
    gsrc = nodal_gradient(mesh, src)
    iu = basis.interpolator(du)
    duc = iu(cpts.T)
    gxv = basis.interpolator(gsrc[0])(cpts.T)
    gzv = basis.interpolator(gsrc[1])(cpts.T)
    gnorm = np.sqrt(gxv**2 + gzv**2) + 1e-300
    sgn = np.sign(gxv * nrm[:, 0] + gzv * nrm[:, 1]).reshape(-1, 1)
    sgn = np.where(sgn == 0, 1.0, sgn)
    dn = -(duc / gnorm) * sgn.ravel()
    Dpts = dn.reshape(-1, 1) * nrm
    # empirical profile: s binned by polygon-signed distance (same
    # winding convention as advected side; bias cancels at sc -> 0)
    d0 = poly_dist(rn, zn, cpts) * np.where(winding_inside(rn, zn, cpts), 1.0, -1.0)
    bedges = np.arange(-1.0, 1.0, 0.02)
    prof = np.zeros(len(bedges) - 1)
    bi = np.clip(np.digitize(d0, bedges) - 1, 0, len(prof) - 1)
    for b in range(len(prof)):
        m = bi == b
        prof[b] = np.median(s[m]) if m.any() else 0.0

    def s_of(d):
        b = np.clip(np.digitize(d, bedges) - 1, 0, len(prof) - 1)
        return prof[b]

    exacts = {0.0: 0.0, 1.0: 9.1750, 0.5: 3.4895, 0.25: 1.5141}
    lin_g1 = 5.2052
    for sc in (0.0, 1.0, 0.5, 0.25):
        ca = cpts + sc * Dpts
        da = poly_dist(rn, zn, ca) * np.where(winding_inside(rn, zn, ca), 1.0, -1.0)
        s_new = s_of(da)
        blend = np.clip((np.abs(da) - 0.4) / 0.1, 0.0, 1.0)
        s_use = (1 - blend) * s_new + blend * s
        Fe, _ = build_F(basis, rn, zn, rn * B.EPS**-2 * s_use**B.P,
                        n3, L3, make_grid)
        Qe, _, _ = Q_of_F(Fe, grid, leray)
        err = float(np.abs(np.linalg.eigvalsh(Qe - Q)).max())
        ex, lin = exacts[sc], lin_g1 * sc
        v = "HOLD-PT" if (err <= 0.2 * ex and err <= lin) else "MISS-PT"
        print(f"scale={sc}: front-only-err={err:.4f} exact={ex:.4f} "
              f"linear={lin:.4f} [{v}]", flush=True)


if __name__ == "__main__":
    main()

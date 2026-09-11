#!/usr/bin/env python3
"""R-A execution survey (beacon 0139): min-Q over the stall tube.

Amended design (0127 amendments + drift charter): E = 0.43 tube,
8 softest J-modes + 2 param dims, 2d axis + 64 corners + 32 interior
(+32 doubling convergence), central point banked (20.33).
HOLD >= 5.0 everywhere -> halt to firewall. STOP < 5.0 -> map banked.
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys
import time

import numpy as np

E_TUBE = 0.43
E_MU = 0.004
E_C = -0.008
HOLD = 5.0
N_CORNER = 64
N_IN = 32


def main() -> None:
    att = "proposals/P253-euler-particle-mechanisms/attempts"
    sys.path.insert(0, f"{att}/0111-beacon-ga-field")
    sys.path.insert(0, f"{att}/0117-beacon-member")
    sys.path.insert(0, f"{att}/0124-beacon-sharp")
    import build_member as B
    from sharp_dQ import build_F, Q_of_F
    from ga_pipeline import leray, make_grid
    from skfem import Basis, ElementTriP1, BilinearForm, LinearForm, asm
    from skfem.helpers import dot, grad, inner

    t0 = time.time()
    r3 = np.load(f"{att}/0117-beacon-member/member-trust-r3.npz")
    mesh = B.build_mesh(6.0, 3.0, 40, 20)
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0].copy()
    zn = basis.doflocs[1].copy()
    u0, mu0, c0 = r3["u"].copy(), float(r3["mu"]), float(r3["c"])
    n3, L3 = 64, 8.0
    grid = None

    def lam_of(du, dmu, dc):
        nonlocal grid
        uu = u0 + du
        src = uu - (c0 + dc) * rn**2 / 2 - (mu0 + dmu)
        root = np.sqrt(src**2 + 1e-6)
        s = (src + root) / 2
        F, grid = build_F(basis, rn, zn, rn * B.EPS**-2 * s**B.P,
                          n3, L3, make_grid)
        Q, _, _ = Q_of_F(F, grid, leray)
        return float(np.linalg.eigvalsh(Q)[0])

    lam_h = lam_of(np.zeros_like(u0), 0.0, 0.0)
    print(f"base lam_min = {lam_h:.4f}", flush=True)

    # 8 softest J free-block modes, scattered to full mesh
    @BilinearForm
    def stiff(a, b, w):
        return w.x[0] * dot(grad(a), grad(b))

    A = asm(stiff, basis).tocsr()
    src0 = u0 - c0 * rn**2 / 2 - mu0
    root0 = np.sqrt(src0**2 + 1e-6)
    s0 = (src0 + root0) / 2
    ds0 = 0.5 * (1 + src0 / root0)
    fi = basis.interpolator(B.EPS**-2 * s0**B.P)

    @LinearForm
    def load(v, w):
        return (w.x[0] ** 3) * fi(w.x) * v

    R = A @ u0 - asm(load, basis)
    dd = B.dirichlet_dofs(basis, mesh, 6.0, 3.0)
    free = np.setdiff1d(np.arange(basis.N), dd)
    R[dd] = 0.0
    jf = 6 * B.EPS**-2 * s0**5 * ds0
    ji = basis.interpolator(jf)

    @BilinearForm
    def jacform(a, b, w):
        return (w.x[0] * dot(grad(a), grad(b))
                - (w.x[0] ** 3) * ji(w.x) * inner(a, b))

    Jf = asm(jacform, basis).tocsr()[free][:, free].toarray()
    w, V = np.linalg.eigh(Jf)
    order = np.argsort(np.abs(w))[:8]
    modes = []
    for j in order:
        v = np.zeros(basis.N)
        v[free] = V[:, j]
        nrm = np.abs(v).max()
        modes.append(v / nrm if nrm > 0 else v)

    rng = np.random.default_rng(0)
    recs = [("central-du", 20.3322)]
    # axis points: modes + param dims
    for i, v in enumerate(modes):
        for sgn in (1.0, -1.0):
            recs.append((f"axis-m{i}{'+' if sgn > 0 else '-'}",
                         lam_of(sgn * E_TUBE * v, 0.0, 0.0)))
    for sgn in (1.0, -1.0):
        recs.append((f"axis-mu{'+' if sgn > 0 else '-'}",
                     lam_of(np.zeros_like(u0), sgn * E_MU, 0.0)))
        recs.append((f"axis-c{'+' if sgn > 0 else '-'}",
                     lam_of(np.zeros_like(u0), 0.0, sgn * E_C)))
    print(f"axis done ({time.time()-t0:.1f}s) min={min(r[1] for r in recs):.4f}",
          flush=True)
    # random corners over 10 collective dims
    for k in range(N_CORNER):
        s = rng.choice([-1.0, 1.0], size=10)
        du = sum(s[i] * modes[i] for i in range(8))
        du = E_TUBE * du / max(np.abs(du).max(), 1e-300)
        recs.append((f"corner-{k}",
                     lam_of(du, s[8] * E_MU, s[9] * E_C)))
    print(f"corners done ({time.time()-t0:.1f}s) min={min(r[1] for r in recs):.4f}",
          flush=True)
    # random interior (full-space) + doubling convergence
    mins = []
    for round2, nn in ((False, N_IN), (True, N_IN)):
        vals = []
        for k in range(nn):
            du = rng.normal(size=basis.N)
            du = (E_TUBE * rng.random() * du / max(np.abs(du).max(), 1e-300))
            vals.append(lam_of(du, (rng.random() * 2 - 1) * E_MU,
                               (rng.random() * 2 - 1) * E_C))
        tag = "in2" if round2 else "in1"
        for k, v in enumerate(vals):
            recs.append((f"{tag}-{k}", v))
        mins.append(min(vals))
        print(f"{tag} done ({time.time()-t0:.1f}s) min={mins[-1]:.4f}", flush=True)
    if len(mins) == 2:
        move = abs(mins[1] - mins[0]) / max(abs(mins[0]), 1e-300)
        print(f"sampling-doubling move = {move:.3f} (bar < 0.10)", flush=True)
    lo = min(recs, key=lambda r: r[1])
    print(f"SURVEY min = {lo[1]:.4f} at {lo[0]} over {len(recs)} samples",
          flush=True)
    print(f"bottom-5: {sorted(recs, key=lambda r: r[1])[:5]}", flush=True)
    if lo[1] >= HOLD:
        print("VERDICT: HOLD -> halt to firewall before downstream use",
              flush=True)
    else:
        print("VERDICT: STOP (dead-by-measurement, min-Q map banked)",
              flush=True)


if __name__ == "__main__":
    main()

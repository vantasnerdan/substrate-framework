#!/usr/bin/env python3
"""M3 recoil-reciprocity build (beacon 0146): far-field momentum ratio
for two unit Gaussian drives on trust-r3. Frozen 0145 falsifier
(KILL>25% / GRAY 15-25% / HOLD<15%) + M2-routing. Conditional label.
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys
import time

import numpy as np


def main() -> None:
    att = "proposals/P253-euler-particle-mechanisms/attempts"
    sys.path.insert(0, f"{att}/0117-beacon-member")
    import build_member as B
    from skfem import Basis, ElementTriP1, BilinearForm, LinearForm, asm
    from skfem.helpers import dot, grad, inner

    t0 = time.time()
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

    @BilinearForm
    def farform(a, b, w):
        return np.where(w.x[0] > 3.0, (w.x[0] ** 2) * inner(a, b), 0.0)

    Mfar = asm(farform, basis)
    mdiag = np.asarray(Mfar.sum(axis=1)).ravel()

    def drive_response(rc, zc, wdt=0.25):
        df = np.exp(-((rn - rc) ** 2 + (zn - zc) ** 2) / wdt**2)
        df = df / np.linalg.norm(df)
        rhs = np.zeros(basis.N)
        rhs[free] = -df[free]
        du = np.zeros(basis.N)
        du[free] = np.linalg.solve(Jf, rhs[free])
        dP = float(mdiag @ du)
        return dP

    dPa = drive_response(1.0, 0.5)
    dPb = drive_response(2.0, 1.5)
    print(f"dP_A = {dPa:.4e} dP_B = {dPb:.4e}", flush=True)
    Ra, Rb = abs(dPa), abs(dPb)
    mean = 0.5 * (Ra + Rb)
    # (b) noise bar: assembly/solve floor ~1e-12 relative
    if min(Ra, Rb) < 1e-9:
        print("F1 verdict: KILL (b) at-noise", flush=True)
    else:
        rel = abs(Ra - Rb) / max(mean, 1e-300)
        print(f"disagreement = {rel:.3f} (KILL>0.25 / GRAY 0.15-0.25 / "
              f"HOLD<0.15)", flush=True)
        if rel > 0.25:
            print("F1 verdict: KILL (a) ratio-varies", flush=True)
        elif rel >= 0.15:
            print("F1 verdict: GRAY (third drive, no verdict)", flush=True)
        else:
            print("F1 verdict: HOLD (conditional reciprocity lives)",
                  flush=True)
    print(f"({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()

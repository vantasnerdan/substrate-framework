#!/usr/bin/env python3
"""P253/0117 feed approximate member into G-a witness pipeline (beacon).

Loads member-bordered-exploratory.npz (documents provenance below), rebuilds
zeta = f on the FEM nodes, revolves axisymmetric omega_theta = r*zeta onto a
3D Cartesian grid, and runs the Q_F/H-matrix budget (imports leray/make_grid
from the TESTED ga_pipeline; same math). B_g absent (no Maxwell stage):
H-rows only (ga-status rows 1-3). Labels: EXPLORATORY-APPROXIMATE with the
member residual as systematic floor — NOT G-a2 numbers.
"""

from __future__ import annotations

import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

sys.path.insert(0, "proposals/P253-euler-particle-mechanisms/attempts/0111-beacon-ga-field")
from ga_pipeline import leray, make_grid  # noqa: E402  (tested FFT path)

ATTEMPT = "proposals/P253-euler-particle-mechanisms/attempts/0117-beacon-member"


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default=f"{ATTEMPT}/member-bordered-exploratory.npz")
    ap.add_argument("--nr", type=int, default=80)
    ap.add_argument("--nz", type=int, default=40)
    ap.add_argument("--n3", type=int, default=48)
    ap.add_argument("--L3", type=float, default=8.0)
    args = ap.parse_args()

    sys.path.insert(0, f"{ATTEMPT}")
    import build_member as B

    dat = np.load(args.npz)
    u, mu, c = dat["u"], float(dat["mu"]), float(dat["c"])
    print(f"member provenance: kap={dat['kap']} rbar={dat['rbar']} "
          f"iz={dat['iz']} res={dat['res']} (APPROXIMATE)", flush=True)
    mesh = B.build_mesh(6.0, 3.0, args.nr, args.nz)
    basis = Basis_check(mesh, B, args)
    rn = basis.doflocs[0]
    assert u.size == basis.N, (u.size, basis.N)
    from build_member import EPS, P
    src = u - c * rn**2 / 2 - mu
    root = np.sqrt(src**2 + 1e-6)
    s = (src + root) / 2
    zeta = EPS**-2 * s ** P  # == omega_theta / r
    wth = rn * zeta
    # revolve onto 3D grid via (r,z) scatter interpolation
    from scipy.interpolate import LinearNDInterpolator
    zn = basis.doflocs[1]
    interp = LinearNDInterpolator(np.stack([rn, zn], axis=1), wth,
                                  fill_value=0.0)
    n3, L3 = args.n3, args.L3
    print(f"FEED npz={args.npz} mesh={args.nr}x{args.nz} grid={n3}^3 L3={L3}", flush=True)
    XYZ, K, dx = make_grid(n3, L3)
    X, Y, Z = XYZ
    R = np.sqrt(X**2 + Y**2)
    W = interp(np.stack([R.ravel(), Z.ravel()], axis=1)).reshape(R.shape)
    Fx = np.where(R > 1e-12, -W * Y / R, 0.0)
    Fy = np.where(R > 1e-12, W * X / R, 0.0)
    F = np.stack([Fx, Fy, np.zeros_like(Fx)], axis=0)
    dV = dx**3
    k2 = K[0]**2 + K[1]**2 + K[2]**2
    e = np.eye(3)
    S = np.zeros((3, 3, n3, n3, n3), dtype=complex)
    for j in range(3):
        Fxej = np.cross(F, e[j][:, None, None, None], axis=0)
        H = np.stack([np.fft.fftn(Fxej[a]) for a in range(3)]) * dV
        S[j] = leray(H, np.stack(K), k2)
    Sphys = np.stack([np.fft.ifftn(S[j], axes=(1, 2, 3)) / dV
                      for j in range(3)])
    Q = np.real(np.einsum("jaxyz,iaxyz->ji", np.conj(Sphys), Sphys)) * dV
    Q = (Q + Q.T) / 2
    lam = np.linalg.eigvalsh(Q)
    print(f"Q_omega eigenvalues: {lam}", flush=True)
    assert np.all(lam > 0), "Q not positive"
    print(f"lambda_omega={lam[0]:.4e} s_omega={np.sqrt(lam[-1]):.4e} "
          f"(EXPLORATORY, member-res floor {dat['res']:.1e})", flush=True)
    print("MEMBER FEED COMPLETE (H-rows only; G-rows await Maxwell stage)")


def Basis_check(mesh, B, args):
    from skfem import Basis, ElementTriP1
    return Basis(mesh, ElementTriP1())


if __name__ == "__main__":
    main()

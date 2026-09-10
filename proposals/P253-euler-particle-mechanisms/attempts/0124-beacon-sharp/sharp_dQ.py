#!/usr/bin/env python3
"""P253/0124 sharp Weyl evaluation (beacon): Q + dQ on the fitted state.

Reuses the tested ga_pipeline FFT path (same math as 0117 feed).
dQ is the EXACT 3x3 perturbation from the S1-gauge-projected
linearized step (S linear in F). Gates G1 (charter L_J<=1) + G2
(transfer remainder) per frozen README design.
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys

import numpy as np


def build_F(basis, rn, zn, wth, n3, L3, make_grid):
    from scipy.interpolate import LinearNDInterpolator
    interp = LinearNDInterpolator(np.stack([rn, zn], axis=1), wth,
                                  fill_value=0.0)
    XYZ, K, dx = make_grid(n3, L3)
    X, Y, Z = XYZ
    R = np.sqrt(X**2 + Y**2)
    W = interp(np.stack([R.ravel(), Z.ravel()], axis=1)).reshape(R.shape)
    Fx = np.where(R > 1e-12, -W * Y / R, 0.0)
    Fy = np.where(R > 1e-12, W * X / R, 0.0)
    return np.stack([Fx, Fy, np.zeros_like(Fx)], axis=0), (XYZ, K, dx)


def Q_of_F(F, grid, leray):
    XYZ, K, dx = grid
    dV = dx**3
    k2 = K[0]**2 + K[1]**2 + K[2]**2
    n3 = F.shape[1]
    e = np.eye(3)
    Sphys = np.zeros((3, 3, n3, n3, n3), dtype=complex)
    for j in range(3):
        Fxej = np.cross(F, e[j][:, None, None, None], axis=0)
        H = np.stack([np.fft.fftn(Fxej[a]) for a in range(3)]) * dV
        Sj = leray(H, np.stack(K), k2)
        Sphys[j] = np.fft.ifftn(Sj, axes=(1, 2, 3)) / dV
    Q = np.real(np.einsum("jaxyz,iaxyz->ji", np.conj(Sphys), Sphys)) * dV
    return (Q + Q.T) / 2, Sphys, dV


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--mesh", default="fitted")
    ap.add_argument("--state", default="proposals/P253-euler-particle-mechanisms/"
                    "attempts/0122-beacon-fitted/member-fitted.npz")
    args = ap.parse_args()
    sys.path.insert(0, "proposals/P253-euler-particle-mechanisms/"
                       "attempts/0111-beacon-ga-field")
    from ga_pipeline import leray, make_grid
    sys.path.insert(0, "proposals/P253-euler-particle-mechanisms/"
                       "attempts/0117-beacon-member")
    import build_member as B
    from skfem import Basis, ElementTriP1, BilinearForm, LinearForm, asm, MeshTri
    from skfem.helpers import dot, grad, inner

    if args.mesh == "fitted":
        att2 = "proposals/P253-euler-particle-mechanisms/attempts/0122-beacon-fitted"
        fm = np.load(f"{att2}/fitted-mesh.npz")
        mesh = MeshTri(fm["p"].T, fm["t"])
    else:
        mesh = B.build_mesh(6.0, 3.0, 40, 20)
    basis = Basis(mesh, ElementTriP1())
    rn = basis.doflocs[0].copy()
    zn = basis.doflocs[1].copy()
    d = np.load(args.state)
    assert d["u"].size == basis.N, (d["u"].size, basis.N)
    u, mu, c = d["u"].copy(), float(d["mu"]), float(d["c"])

    src = u - c * rn**2 / 2 - mu
    root = np.sqrt(src**2 + 1e-6)
    s = (src + root) / 2
    ds = 0.5 * (1 + src / root)
    zeta = B.EPS**-2 * s**B.P
    n3, L3 = 64, 8.0
    F, grid = build_F(basis, rn, zn, rn * zeta, n3, L3, make_grid)
    Q, Sphys, dV = Q_of_F(F, grid, leray)
    lam = np.linalg.eigvalsh(Q)
    print(f"Q_fitted eigenvalues: {lam}", flush=True)
    print(f"lam_min = {lam[0]:.4f} s_w = {np.sqrt(lam[-1]):.4f}", flush=True)

    # linearized gauge-projected step (0123 S1 machinery)
    @BilinearForm
    def stiff(a, b, w):
        return w.x[0] * dot(grad(a), grad(b))

    @BilinearForm
    def massr(a, b, w):
        return w.x[0] * inner(a, b)

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
    # S1 gauge: project discrete translation modes (nodal grads of u)
    sys.path.insert(0, "proposals/P253-euler-particle-mechanisms/"
                       "attempts/0123-beacon-fitted")
    from x18_probes import nodal_gradient
    gg = nodal_gradient(mesh, u)
    Qt, _ = np.linalg.qr(np.stack([gg[0][free], gg[1][free]], axis=1))
    du[free] = du[free] - Qt @ (Qt.T @ du[free])
    print(f"||du_gauged||_oo = {np.abs(du).max():.4e}", flush=True)
    # dF field: dzeta = (jf/6) du (P = 6)
    dzeta = (jf / 6) * du
    dF, _ = build_F(basis, rn, zn, rn * dzeta, n3, L3, make_grid)
    XYZ, K, dx = grid
    k2 = K[0]**2 + K[1]**2 + K[2]**2
    e = np.eye(3)
    dS_all = []
    # exact bilinear cross terms (vectorized cleanly)
    for j in range(3):
        dFxej = np.cross(dF, e[j][:, None, None, None], axis=0)
        dH = np.stack([np.fft.fftn(dFxej[a]) for a in range(3)]) * dV
        dS_all.append(np.fft.ifftn(leray(dH, np.stack(K), k2),
                                   axes=(1, 2, 3)) / dV)
    dS_all = np.stack(dS_all)
    dQ = (np.real(np.einsum("jaxyz,iaxyz->ji", np.conj(dS_all), Sphys))
          + np.real(np.einsum("jaxyz,iaxyz->ji", np.conj(Sphys), dS_all))) * dV
    dQ = (dQ + dQ.T) / 2
    dl = np.linalg.eigvalsh(dQ)
    g1 = float(np.abs(dl).max())
    print(f"dQ eigenvalues: {dl}", flush=True)
    print(f"G1: max|eig(dQ)| = {g1:.4f} (gate <= 0.9); "
          f"margin vs lam_min x{lam[0] / max(g1, 1e-12):.1f}", flush=True)

    # G2: quadratic remainder on tube E = 2||du||_oo
    # d2zeta = EPS^-2 P[(P-1)s^{P-2}ds^2 + s^{P-1}d2s], d2s = 1e-6/(2 root^3)
    d2s = 1e-6 / (2 * root**3)
    d2z = B.EPS**-2 * B.P * ((B.P - 1) * s**(B.P - 2) * ds**2
                             + s**(B.P - 1) * d2s)
    E = 2 * float(np.abs(du).max())
    rem_zeta = 0.5 * float(np.abs(d2z).max()) * E**2
    # remainder field scale vs linearized: compare representative scales
    lin_zeta = float(np.abs(jf / 6 * du).max())
    print(f"G2: tube E = {E:.4e}, quad remainder zeta-scale {rem_zeta:.4e} "
          f"vs linearized {lin_zeta:.4e}", flush=True)
    # G2': EXACT nonlinear transfer (no Taylor remainder theory).
    # True field perturbation at u+du with (mu,c) tube corners, pushed
    # through the same Q pipeline; sensitivity envelope 0.5x/2x.
    def exact_dQ(scale, dmu, dc):
        dus = scale * du
        s2 = u + dus - (c + dc) * rn**2 / 2 - (mu + dmu)
        r2 = np.sqrt(s2**2 + 1e-6)
        dz_exact = B.EPS**-2 * (((s2 + r2) / 2)**B.P - s**B.P)
        dFe, _ = build_F(basis, rn, zn, rn * dz_exact, n3, L3, make_grid)
        Q2, _, _ = Q_of_F(F + dFe, grid, leray)
        return float(np.abs(np.linalg.eigvalsh(Q2 - Q)).max())

    g_central = exact_dQ(1.0, 0.0, 0.0)
    g_half = exact_dQ(0.5, 0.0, 0.0)
    g_double = exact_dQ(2.0, 0.004, -0.008)
    print(f"G2': exact |dQ| central={g_central:.4f} half={g_half:.4f} "
          f"double+corner={g_double:.4f} (gate <= 0.9)", flush=True)


if __name__ == "__main__":
    main()

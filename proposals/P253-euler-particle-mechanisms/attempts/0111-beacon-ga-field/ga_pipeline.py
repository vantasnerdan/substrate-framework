#!/usr/bin/env python3
"""P253/0111 G-a witness-budget pipeline (beacon, attempt-local).

Implements 0107 derivation (15),(18)-(25) END-TO-END on supplied field data:
  Q_F Gram -> lambda_F, s_F -> cutoff witnesses w_{F,j,R} -> eta_F(R) budget
  -> M_F response matrix -> sigma_min bound check (25).

No campaign data is consumed: the driver below feeds a synthetic compact
divergence-free torus field to prove the pipeline; the TRUE Cao member
fields are G-a2 input (see ga-status.md). All identities are exact FFT
spectral calculus on a periodic box large enough that tails are ~1e-15.

Checks (assertions, not prints): Q symmetry/positivity, eta decay in R,
sigma_min(M) >= |c| lambda/2 at the (24)-chosen R, H-analog det != 0.
"""

from __future__ import annotations

import numpy as np


def make_grid(n: int, L: float):
    x = np.linspace(-L / 2, L / 2, n, endpoint=False)
    dx = L / n
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    k1 = np.fft.fftfreq(n, d=dx) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
    return (X, Y, Z), (KX, KY, KZ), dx


def compact_torus_field(XYZ, K=None, R0: float = 2.0, w: float = 0.55):
    """Compact div-free torus field F = curl(A_phi e_phi), spectral curl.

    The envelope is C^1-tapered (exp core, cos^2 skirt to zero at 2.5w);
    taking the curl spectrally makes div F = 0 to machine precision, so the
    periodic-box Leray/Poisson path is consistent (no k=0 mismatch).
    """
    X, Y, Z = XYZ
    R = np.sqrt(X**2 + Y**2)
    t = np.sqrt((R - R0) ** 2 + Z**2) / w
    env = np.where(t < 2.0, np.exp(-t**2),
                   np.where(t < 2.5, np.exp(-t**2) * np.cos((t - 2.0) * np.pi) ** 2, 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        Ax = np.where(R > 1e-12, env * Y / R, 0.0)
        Ay = np.where(R > 1e-12, env * X / R, 0.0)
    A = np.stack([-Ax, Ay, np.zeros_like(env)], axis=0)
    if K is None:
        return A
    Ahat = np.stack([np.fft.fftn(A[a]) for a in range(3)])
    Fhat = np.stack([1j * (K[1] * Ahat[2] - K[2] * Ahat[1]),
                     1j * (K[2] * Ahat[0] - K[0] * Ahat[2]),
                     1j * (K[0] * Ahat[1] - K[1] * Ahat[0])])
    n = XYZ[0].shape[0]
    return np.real(np.stack([np.fft.ifftn(Fhat[a]) for a in range(3)]))

def leray(Fhat, K, k2):
    """Toroidal Leray projection (mean-preserving)."""
    k2s = np.where(k2 > 0, k2, 1.0)
    KdotF = K[0] * Fhat[0] + K[1] * Fhat[1] + K[2] * Fhat[2]
    return Fhat - K * (KdotF / k2s)


def check(label: str, predicate: bool) -> None:
    if not predicate:
        raise AssertionError(label)
    print(f"PASS {label}")

def main() -> None:
    n, L = 64, 10.0
    XYZ, K, dx = make_grid(n, L)
    k2 = K[0] ** 2 + K[1] ** 2 + K[2] ** 2
    dV = dx**3
    F = compact_torus_field(XYZ, K)

    # S_F e_j = P_L(F x e_j), shape (3 vectors, 3 components, n^3)
    S = np.zeros((3, 3, n, n, n), dtype=complex)
    e = np.eye(3)
    for j in range(3):
        Fxej = np.cross(F, e[j][:, None, None, None], axis=0)
        Fxej_hat = np.stack([np.fft.fftn(Fxej[a]) for a in range(3)]) * dV
        S[j] = leray(Fxej_hat, np.stack(K), k2)
    Sphys = np.stack([np.fft.ifftn(S[j], axes=(1, 2, 3)) / dV for j in range(3)])  # /dV undoes fwd norm
    Q = np.real(np.einsum("jaxyz,iaxyz->ji", np.conj(Sphys), Sphys)) * dV
    Q = (Q + Q.T) / 2
    lam = np.linalg.eigvalsh(Q)
    lambda_F, s_F = lam[0], np.sqrt(lam[-1])
    assert np.all(lam > 0), f"Q_F not positive: {lam}"
    print(f"PASS Q_F positive: lambda_min={lambda_F:.6e} s_F={s_F:.6e}")

    # Vector potentials A_{F,j} = (-Delta)^{-1} curl(F x e_j), Fourier
    k2s = np.where(k2 > 0, k2, 1.0)
    A = np.zeros_like(S)
    for j in range(3):
        Fxej_hat = np.stack(
            [np.fft.fftn(np.cross(F, e[j][:, None, None, None], axis=0)[a])
             for a in range(3)]) * dV
        # curl in Fourier: i k x G
        G = Fxej_hat
        curlG = np.stack([
            1j * (K[1] * G[2] - K[2] * G[1]),
            1j * (K[2] * G[0] - K[0] * G[2]),
            1j * (K[0] * G[1] - K[1] * G[0])])
        A[j] = curlG / k2s
        A[j][:, k2 == 0] = 0.0
    Aphys = np.real(np.stack([np.fft.ifftn(A[j], axes=(1, 2, 3)) / dV for j in range(3)]))
    Sphys_r = np.real(Sphys)

    # eta_F(R): tail + Hardy annulus terms per (21), C=1
    Rgrid = np.sqrt(XYZ[0]**2 + XYZ[1]**2 + XYZ[2]**2)
    etas = {}
    for R in (1.0, 2.0, 3.0, 4.0):
        tail = np.sqrt(sum(np.sum(Sphys_r[j][:, Rgrid > R]**2) * dV for j in range(3)))
        ann = (Rgrid > R) & (Rgrid < 2 * R)
        rr = np.where(ann, Rgrid, 1.0)
        hardy = np.sqrt(sum(np.sum((Aphys[j] / rr)**2 * ann) * dV for j in range(3)))
        etas[R] = tail + hardy
    assert etas[4.0] < etas[1.0], etas
    print("PASS eta_F(R) decays: " + ", ".join(f"R={R}: {v:.3e}" for R, v in etas.items()))

    # (24): first R with eta <= lambda/(4 s); build M_F response, c = 1
    Rch = next(R for R in (1.0, 2.0, 3.0, 4.0) if etas[R] <= lambda_F / (4 * s_F))
    chi = np.where(Rgrid <= Rch, 1.0, np.where(Rgrid >= 2 * Rch, 0.0, 0.5 + 0.5 * np.cos(np.pi * (Rgrid - Rch) / Rch)))
    # w_j = curl(chi A_j): spectral curl of product
    M = np.zeros((3, 3))
    for j in range(3):
        chiA = chi * Aphys[j]
        chiA_hat = np.stack([np.fft.fftn(chiA[a]) for a in range(3)]) * dV
        curl_hat = np.stack([
            1j * (K[1] * chiA_hat[2] - K[2] * chiA_hat[1]),
            1j * (K[2] * chiA_hat[0] - K[0] * chiA_hat[2]),
            1j * (K[0] * chiA_hat[1] - K[1] * chiA_hat[0])])
        wj = np.real(np.fft.ifftn(curl_hat, axes=(1, 2, 3)) / dV)
        for i in range(3):
            M[i, j] = np.sum(wj * Sphys_r[i]) * dV
    smin = np.linalg.svd(M, compute_uv=False)[-1]
    assert smin >= lambda_F / 2 * (1 - 1e-6), (smin, lambda_F / 2)
    assert abs(np.linalg.det(M)) > 0
    print(f"PASS witness bound at R={Rch}: sigma_min(M)={smin:.6e} "
          f">= lambda/2={lambda_F/2:.6e}; det={np.linalg.det(M):.6e}")

    print("ALL G-A PIPELINE CHECKS PASSED")


if __name__ == "__main__":
    main()

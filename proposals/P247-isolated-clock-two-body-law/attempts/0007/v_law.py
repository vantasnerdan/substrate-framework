"""P247 attempt 0007: S-potential interaction law and boundary-layer check.

V(S) = -0.5 tr(S^2) - tr(S^3) + 0.5 (tr S^2)^2 + 0.5  (exact form from
the committed functional). Debox vacuum = radial projector field
P(x) = n_hat n_hat, V(P) = 0 pointwise.

Parts:
  A  r-resolved potential density of the committed roots across the
     R ladder (mu = 1 slice and mu-averaged): is the box-growth a
     boundary layer or a core effect?
  B  numeric Hessian of V at the projector diag(1,0,0) in the 5-dim
     symmetric-deviation space: channel masses (hand expansion says
     X11: +0.5, X22: +0.5, X33: +0.5, X12: -2, X13: -2, X23: +1 as
     coefficients of Xij^2 in V).
  C  validation: single-clock integral of V(S1) - V(P1) on the modal
     grid must reproduce the functional's potential component
     (17.2276 at R = 24).
  D  two-clock ladder (cylindrical grid, exact frames per clock):
     E_int(d) = int [ V(S1 + S2 - P1) - V(S1) - V(S2) + V(P1) ],
     fit the de-boxed law; box-refinement evidence.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
from numpy.polynomial.chebyshev import Chebyshev

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "0004"))
sys.path.insert(0, str(HERE.parent / "0002"))
sys.path.insert(0, str(HERE.parent / "0001"))

import c3_functional as c3  # noqa: E402

ORDER = 24
R_LADDER = [20.0, 24.0, 28.0, 30.0]
D_LADDER = [6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 21.0, 24.0, 28.0]
BOX_HALF = 40.0  # cylindrical half-extent in z around the pair center
RHO_MAX = 40.0


def potentials(s: np.ndarray) -> np.ndarray:
    """V for stacked 3x3 symmetric matrices s (..., 3, 3)."""
    tr2 = np.einsum("...ij,...ji->...", s, s)
    tr3 = np.einsum("...ij,...jk,...ki->...", s, s, s)
    return -0.5 * tr2 - tr3 + tr2**2 + 0.5


def spatial_from_root(root_coef: np.ndarray, R: float, r: np.ndarray, mu: np.ndarray):
    """Exact debox hedgehog field: returns S (..., 3, 3), eigen-mains.

    root_coef shape (3, ORDER); r and mu broadcastable arrays.
    """
    n2 = np.minimum(np.clip(r / R, 0.0, None) ** 2, 1.0)
    x = 2.0 * n2 - 1.0
    modal = [
        Chebyshev(root_coef[row], domain=[-1.0, 1.0])(x)
        for row in range(3)
    ]
    q = n2 + n2 * (1.0 - n2) * modal[0]
    tangent = (1.0 - n2) * (1.0 / 3.0 + modal[1])
    delta = n2**2 * (1.0 - n2) * modal[2] * (1.0 - mu**2)
    sin_th = np.sqrt(np.maximum(0.0, 1.0 - mu**2))
    director = np.stack([sin_th * np.ones_like(mu), np.zeros_like(np.sin(mu * 0)), mu], axis=-1)
    polar = np.stack([mu * np.ones_like(mu), np.zeros_like(mu), -sin_th], axis=-1)
    azimuthal = np.stack([np.zeros_like(mu), np.ones_like(mu), np.zeros_like(mu)], axis=-1)
    lam_n = tangent + q
    lam_p = tangent + delta
    lam_a = tangent - delta
    S = (
        lam_n[..., None, None] * _outer(director)
        + lam_p[..., None, None] * _outer(polar)
        + lam_a[..., None, None] * _outer(azimuthal)
    )
    return S, (lam_n, lam_p, lam_a)


def _outer(v: np.ndarray) -> np.ndarray:
    return v[..., :, None] * v[..., None, :]


def part_a(roots) -> dict:
    out = {}
    for R in R_LADDER:
        root = np.asarray(roots[f"R{int(R)}"], dtype=np.float64).reshape(3, ORDER)
        # mu = 1 slice (axis) and a mu = 0 ring
        r = np.linspace(0.05 * R, 0.995 * R, 160)
        for mu_val, tag in [(1.0, "axis"), (0.0, "equator")]:
            rr = np.broadcast_to(r, (160,))
            mmu = np.full(160, mu_val)
            S, _ = spatial_from_root(root, R, rr, mmu)
            v = potentials(S)
            out[f"R{int(R)}_{tag}"] = {
                "r": r.tolist(),
                "V": v.tolist(),
            }
        # integrated potential on this root's own modal grid
        flat = torch.zeros((4, ORDER), dtype=torch.float64)
        flat[:3] = torch.tensor(root)
        _total, comp = c3.extended_static(
            flat, radial_order=ORDER, radial_nodes=72, angular_nodes=36,
            radius=R, chi_angular_degrees=(0,),
        )
        out[f"R{int(R)}_functional_potential"] = float(comp["potential"])
    return out


def part_b() -> dict:
    base = np.zeros((3, 3))
    base[0, 0] = 1.0
    names = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    diag = {}
    for i, j in names:
        E = np.zeros((3, 3))
        E[i, j] = E[j, i] = 1.0
        eps = 1e-4
        vp = float(potentials(base[None] + eps * E[None])[0])
        vm = float(potentials(base[None] - eps * E[None])[0])
        v0 = float(potentials(base[None])[0])
        second = (vp + vm - 2 * v0) / eps**2
        diag[f"X{i+1}{j+1}"] = 0.5 * second
    return {"V_projector": v0, "quadratic_coefficients": diag}


def part_c_validation(roots) -> dict:
    """Integrate V(S1) on the modal grid with my own V/frame code and
    compare against the functional's potential component (R = 24)."""
    root = np.asarray(roots["R24"], dtype=np.float64).reshape(3, ORDER)
    R = 24.0
    from cpu_energy import gauss_grid  # noqa: E402  (paths set by c3 import)

    r_nodes, w_r, mu_nodes, w_mu = gauss_grid(72, 36, R)
    rr = np.broadcast_to(np.asarray(r_nodes)[:, None], (72, 36))
    mm = np.broadcast_to(np.asarray(mu_nodes)[None, :], (72, 36))
    S, _ = spatial_from_root(root, R, rr, mm)
    v = potentials(S)
    w = np.asarray(w_r)[:, None] * np.asarray(w_mu)[None, :]
    integral = float(np.sum(2 * np.pi * rr**2 * w * v))
    flat = torch.zeros((4, ORDER), dtype=torch.float64)
    flat[:3] = torch.tensor(root)
    _total, comp = c3.extended_static(
        flat, radial_order=ORDER, radial_nodes=72, angular_nodes=36,
        radius=R, chi_angular_degrees=(0,),
    )
    ref = float(comp["potential"])
    return {"my_modal_integral": integral, "functional": ref, "abs_diff": abs(integral - ref)}


def pair_field(root: np.ndarray, R: float, X: np.ndarray, center: np.ndarray):
    """S and P (projector) for one clock on Cartesian points X (N,3)."""
    rel = X - center[None, :]
    r = np.linalg.norm(rel, axis=1)
    mu = np.where(r > 1e-9, rel[:, 2] / np.maximum(r, 1e-9), 1.0)
    S, _ = spatial_from_root(root, R, r, mu)
    sin_th = np.sqrt(np.maximum(0.0, 1.0 - mu**2))
    director = np.stack([sin_th, np.zeros_like(mu), mu], axis=-1)
    P = _outer(director)
    return S, P


def part_d_ladder(roots, masses=(0.5,)) -> dict:
    root = np.asarray(roots["R24"], dtype=np.float64).reshape(3, ORDER)
    R = 24.0
    results = {}
    for d in D_LADDER:
        z0 = -d / 2.0
        z1 = d / 2.0
        n_rho, n_z = 260, 420
        rho = np.linspace(0.0, RHO_MAX, n_rho)
        z = np.linspace(z0 - BOX_HALF, z1 + BOX_HALF, n_z)
        dr = rho[1] - rho[0]
        dz = z[1] - z[0]
        RHO, Z = np.meshgrid(rho, z, indexing="ij")
        X = np.stack([RHO.ravel(), np.zeros_like(RHO).ravel(), Z.ravel()], axis=1)
        c1 = np.array([0.0, 0.0, z0])
        c2 = np.array([0.0, 0.0, z1])
        S1, P1 = pair_field(root, R, X, c1)
        S2, P2 = pair_field(root, R, X, c2)
        # superpose DEVIATIONS on the common background P1
        S_pair = P1 + (S1 - P1) + (S2 - P2)
        V_pair = potentials(S_pair)
        V_s1 = potentials(S1)
        V_s2 = potentials(S2)
        V_p1 = potentials(P1)
        integrand = V_pair - V_s1 - V_s2 + V_p1
        weights = 2 * np.pi * RHO.ravel() * dr * dz
        e_int = float(np.sum(weights * integrand))
        results[f"d{d}"] = e_int
        print(f"[D d={d}] E_int = {e_int:.6e}", flush=True)
    ds = np.array(sorted(D_LADDER))
    vals = np.array([results[f"d{d}"] for d in ds])
    # saturating fit: E(d) = Einf - A exp(-gamma d), exclude the two
    # smallest separations (core-overlap regime)
    m = ds >= 10.0
    A_mat = np.stack([np.ones(m.sum()), -np.exp(-ds[m])], axis=1)

    def resid(theta):
        Einf, logA, g = theta
        return vals[m] - (Einf - np.exp(logA) * np.exp(-g * ds[m]))

    from scipy.optimize import least_squares
    sol = least_squares(lambda th: resid(th), x0=[8.2, 1.0, 0.5])
    Einf, Aamp, gam = float(sol.x[0]), float(np.exp(sol.x[1])), float(sol.x[2])
    ss = float(np.sum(sol.fun**2))
    fit = {"E_inf": Einf, "A": Aamp, "gamma": gam, "sum_sq_resid": ss,
           "force_sign": "attractive" if vals[-1] > vals[0] else "repulsive",
           "ddependent_part_d10": float(Einf - vals[list(ds).index(10.0)])}
    print(f"[D fit] E_inf={Einf:.4f} A={Aamp:.4f} gamma={gam:.4f} ssr={ss:.2e} force={fit['force_sign']}", flush=True)
    return {"ladder": {k: v for k, v in results.items()}, "fit": fit}


def main() -> None:
    started = time.time()
    roots = np.load(HERE.parent / "0003" / "roots24.npz")
    payload = {}
    payload["a_rresolved"] = part_a(roots)
    payload["b_hessian"] = part_b()
    print("[B]", json.dumps(payload["b_hessian"]["quadratic_coefficients"], indent=1), flush=True)
    payload["c_validation"] = part_c_validation(roots)
    print("[C]", payload["c_validation"], flush=True)
    payload["d_ladder"] = part_d_ladder(roots)
    payload["minutes_total"] = round((time.time() - started) / 60.0, 2)
    (HERE / "v-law-measurement.json").write_text(json.dumps(payload, indent=2, default=float))
    print("WROTE v-law-measurement.json")


if __name__ == "__main__":
    main()

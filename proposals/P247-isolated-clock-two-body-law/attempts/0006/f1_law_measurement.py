"""P247 attempt 0006, gates F1/F2 measurement: leading two-clock law.

Inputs: committed order-24 roots (modal Chebyshev coefficients in
x = 2 (r/R)^2 - 1), the extended functional for component attribution.

  M0  components at chi=0 across the R ladder -> which sector carries
      the linear box-growth (+0.38/R).
  M1  per-row tail-mass fit on the R=24 and R=30 profiles:
      ln|phi - phi_vac| vs r on the outer region; exponential vs
      power-law discrimination.
  M2  effective source j = (-d2/dr2 - (2/r)d/dr + m^2) (phi - phi_vac)
      on the physical gauss nodes (exact polynomial derivatives).
  M3  E_int(d) = double quadrature j G_m j over the node set on a
      separation ladder; fit the decay constant gamma vs the fitted m.
  M4  sign comparison against the C-M5S-008 confined repulsive law.
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

import c3_functional as c3  # noqa: E402  (pulls in cpu_energy paths)
from cpu_energy import gauss_grid  # noqa: E402

ORDER = 24
R_LADDER = [20.0, 22.0, 24.0, 26.0, 28.0, 30.0]
D_LADDER = [6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 26.0]
N_R, N_MU = 72, 36


def base_nodes(n: int) -> np.ndarray:
    """Chebyshev-Gauss nodes in (0,1), increasing."""
    return 0.5 * (1.0 + np.cos((2 * np.arange(n) + 1) * np.pi / (2 * n)))[::-1]


def m0_channel_attribution(roots) -> dict:
    rows = {}
    for R in R_LADDER:
        flat = torch.zeros((4, ORDER), dtype=torch.float64)
        flat[:3] = torch.tensor(
            np.asarray(roots[f"R{int(R)}"], dtype=np.float64).reshape(3, ORDER)
        )
        _total, comp = c3.extended_static(
                flat, radial_order=ORDER, radial_nodes=72, angular_nodes=36,
                radius=R, chi_angular_degrees=(0,),
            )
        scalars = {}
        for k, v in comp.items():
            t = torch.as_tensor(v, dtype=torch.float64)
            if t.numel() == 1 and torch.isfinite(t):
                scalars[k] = float(t)
        rows[R] = scalars
        print(f"[M0 R={R}] " + " ".join(f"{k}={v:.4f}" for k, v in sorted(scalars.items())), flush=True)
    growth = {}
    for k in rows[R_LADDER[0]]:
        vals = np.array([rows[R][k] for R in R_LADDER])
        if np.isfinite(vals).all() and vals.std() > 1e-6:
            growth[k] = float(np.polyfit(R_LADDER, vals, 1)[0])
    return {"rows": {str(k): v for k, v in rows.items()}, "slope_vs_R": growth}


def m1_tail_fit(roots, n_tail: int = 12) -> dict:
    out = {}
    nodes01 = base_nodes(ORDER)
    for name in ["R24", "R30"]:
        root = np.asarray(roots[name], dtype=np.float64).reshape(3, ORDER)
        R = float(name[1:])
        fits = {}
        for row in range(3):
            coef = Chebyshev(root[row], domain=[-1.0, 1.0])
            r_fine = np.linspace(0.82 * R, 0.995 * R, 400)
            phi = coef(2 * (r_fine / R) ** 2 - 1.0)
            phi_vac = float(coef(2 * 0.995**2 - 1.0))
            amp = np.abs(phi - phi_vac)
            keep = amp > 1e-13 * max(amp.max(), 1e-300)
            r_k, y = r_fine[keep], np.log(amp[keep])
            n = min(n_tail * 30, len(y))
            slope, intercept = np.polyfit(r_k[-n:], y[-n:], 1)
            exp_resid = float(np.std(y[-n:] - (slope * r_k[-n:] + intercept)))
            pslope = float(np.polyfit(np.log(r_k[-n:]), y[-n:], 1)[0])
            fits[f"row{row}"] = {
                "mass_fit": float(-slope),
                "exp_fit_residual": exp_resid,
                "power_fit_exponent": float(-pslope),
                "phi_vac": phi_vac,
                "tail_amplitude_at_0.55R": float(amp[0]),
            }
            print(
                f"[M1 {name} row{row}] mass={-slope:.4f} exp_resid={exp_resid:.2e} "
                f"pow_exp={-pslope:.2f} amp={amp[0]:.3e}",
                flush=True,
            )
        out[name] = {"R": R, "rows": fits}
    return out


def m2_m3_pairing(roots, masses: dict, d_ladder=D_LADDER) -> dict:
    r_nodes, w_r, mu_nodes, w_mu = gauss_grid(N_R, N_MU, 24.0)
    r_nodes = np.asarray(r_nodes, dtype=np.float64)
    mu_nodes = np.asarray(mu_nodes, dtype=np.float64)
    w = np.asarray(w_r, dtype=np.float64)[:, None] * np.asarray(w_mu, dtype=np.float64)[None, :]
    w = w.ravel()
    sin_th = np.sqrt(np.maximum(0.0, 1.0 - mu_nodes**2))
    X = np.stack(
        [
            (r_nodes[:, None] * sin_th[None, :]).ravel(),
            np.zeros(N_R * N_MU),
            (r_nodes[:, None] * mu_nodes[None, :]).ravel(),
        ],
        axis=1,
    )

    root = np.asarray(roots["R24"], dtype=np.float64).reshape(3, ORDER)
    R = 24.0
    results = {"rows_used": [], "ladders": {}, "source_diagnostics": {}}

    for row, m in masses.items():
        coef = Chebyshev(root[row], domain=[-1.0, 1.0])
        x_of_r = 2 * (r_nodes / R) ** 2 - 1.0
        psi = coef(x_of_r) - float(coef(2 * 0.9999**2 - 1.0))
        d1 = coef.deriv(1)(x_of_r) * (4 * r_nodes / R**2)
        d2 = coef.deriv(2)(x_of_r) * (4 * r_nodes / R**2) ** 2 + coef.deriv(1)(x_of_r) * (4 / R**2)
        j = -d2 - (2.0 / np.maximum(r_nodes, 1e-6)) * d1 + m * m * psi
        # j lives on rings: broadcast over mu
        j_ring = np.repeat(j, N_MU)

        core = float(np.abs(j_ring).max())
        tail = float(np.abs(j_ring[-N_MU:]).max())
        results["source_diagnostics"][f"row{row}"] = {
            "max_abs_j": core,
            "outer_ring_abs_j": tail,
            "tail_to_core": tail / max(core, 1e-300),
            "monopole": float(np.sum(w * j_ring)),
        }
        print(f"[M2 row{row}] m={m:.4f} |j|max={core:.3e} outer/inner={tail / max(core, 1e-300):.2e}", flush=True)

        ladder = {}
        for d in d_ladder:
            diff = X[:, None, :] - (X[None, :, :] + np.array([0.0, 0.0, d]))
            dist = np.sqrt((diff * diff).sum(-1))
            G = np.exp(-m * dist) / (4.0 * np.pi * np.maximum(dist, 1e-9))
            ladder[d] = float(np.einsum("i,ij,j->", w, G * j_ring[:, None], w * j_ring))
        results["ladders"][f"row{row}"] = ladder
        ds = np.array(sorted(ladder))
        es = np.abs(np.array([ladder[d] for d in ds]))
        good = es > 1e-13 * es.max()
        gamma = float(-np.polyfit(ds[good], np.log(es[good]), 1)[0])
        # local power: fit with free power p
        A = np.stack([np.ones(good.sum()), -ds[good], -np.log(ds[good])], axis=1)
        sol, *_ = np.linalg.lstsq(A, np.log(es[good]), rcond=None)
        results[f"fit_row{row}"] = {"gamma": gamma, "gamma_with_power": float(sol[1]), "power_p": float(sol[2])}
        print(
            f"[M3 row{row}] gamma={gamma:.4f} (m={m:.4f}) gamma_p={sol[1]:.4f} p={sol[2]:.2f} "
            f"E_int(10)={ladder[10.0]:.6e} sign={'repulsive' if ladder[10.0] > 0 else 'attractive'}",
            flush=True,
        )
        results["rows_used"].append(row)
    return results


def main() -> None:
    started = time.time()
    roots = np.load(HERE.parent / "0003" / "roots24.npz")

    payload = {"m0": m0_channel_attribution(roots)}
    m1 = m1_tail_fit(roots)
    payload["m1"] = m1

    masses = {}
    for row in range(3):
        f = m1["R24"]["rows"][f"row{row}"]
        if f["exp_fit_residual"] < 0.10 and f["mass_fit"] > 1e-2:
            masses[row] = f["mass_fit"]
    payload["fitted_masses"] = {str(k): v for k, v in masses.items()}
    payload["m2_m3"] = m2_m3_pairing(roots, masses) if masses else None
    payload["minutes_total"] = round((time.time() - started) / 60.0, 2)

    (HERE / "f1-law-measurement.json").write_text(json.dumps(payload, indent=2, default=float))
    print("WROTE f1-law-measurement.json")


if __name__ == "__main__":
    main()

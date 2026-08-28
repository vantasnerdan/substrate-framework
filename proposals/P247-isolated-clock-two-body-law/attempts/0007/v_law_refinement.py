"""P247 0011: materialize the box/resolution refinement of the de-boxed law.

Runs the attempt-0007 part_d ladder at (a) 2x grid resolution and
(b) a 1.5x larger integration box, refits the saturating exponential on
each, and records the refinement next to the primary measurement.
Output: v-law-refinement.json (this directory).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import v_law  # noqa: E402


def refined_ladder(n_mult: float, box_mult: float) -> dict:
    roots = np.load(HERE.parent / "0003" / "roots24.npz")
    root = np.asarray(roots["R24"], dtype=np.float64).reshape(3, v_law.ORDER)
    R = 24.0
    results = {}
    for d in v_law.D_LADDER:
        z0 = -d / 2.0
        z1 = d / 2.0
        n_rho = int(260 * n_mult)
        n_z = int(420 * n_mult)
        rho_max = v_law.RHO_MAX * box_mult
        box_half = v_law.BOX_HALF * box_mult
        rho = np.linspace(0.0, rho_max, n_rho)
        z = np.linspace(z0 - box_half, z1 + box_half, n_z)
        dr = rho[1] - rho[0]
        dz = z[1] - z[0]
        RHO, Z = np.meshgrid(rho, z, indexing="ij")
        X = np.stack([RHO.ravel(), np.zeros_like(RHO).ravel(), Z.ravel()], axis=1)
        c1 = np.array([0.0, 0.0, z0])
        c2 = np.array([0.0, 0.0, z1])
        S1, P1 = v_law.pair_field(root, R, X, c1)
        S2, P2 = v_law.pair_field(root, R, X, c2)
        S_pair = P1 + (S1 - P1) + (S2 - P2)
        V_pair = v_law.potentials(S_pair)
        V_s1 = v_law.potentials(S1)
        V_s2 = v_law.potentials(S2)
        V_p1 = v_law.potentials(P1)
        integrand = V_pair - V_s1 - V_s2 + V_p1
        weights = 2 * np.pi * RHO.ravel() * dr * dz
        results[f"d{d}"] = float(np.sum(weights * integrand))
        print(f"[mult=({n_mult},{box_mult}) d={d}] E_int = {results[f'd{d}']:.6e}", flush=True)

    ds = np.array(sorted(v_law.D_LADDER))
    vals = np.array([results[f"d{d}"] for d in ds])
    m = ds >= 10.0
    from scipy.optimize import least_squares

    sol = least_squares(
        lambda th: vals[m] - (th[0] - np.exp(th[1]) * np.exp(-th[2] * ds[m])),
        x0=[8.2, 1.0, 0.5],
    )
    fit = {
        "E_inf": float(sol.x[0]),
        "A": float(np.exp(sol.x[1])),
        "gamma": float(sol.x[2]),
        "sum_sq_resid": float(np.sum(sol.fun**2)),
    }
    print(f"[fit mult=({n_mult},{box_mult})] {fit}", flush=True)
    return {"ladder": results, "fit": fit,
            "grid": [n_rho, n_z], "box": [float(rho_max), float(box_half)]}


def main() -> None:
    started = time.time()
    payload = {
        "primary": {"grid": [260, 420], "box": [v_law.RHO_MAX, v_law.BOX_HALF],
                    "note": "primary fit in v-law-measurement.json: gamma=0.2442, E_inf=8.1133"},
        "refine_2x_grid": refined_ladder(2.0, 1.0),
        "refine_15x_box": refined_ladder(1.0, 1.5),
        "minutes_total": round((time.time() - started) / 60.0, 2),
    }
    (HERE / "v-law-refinement.json").write_text(json.dumps(payload, indent=2, default=float))
    print("WROTE v-law-refinement.json")


if __name__ == "__main__":
    main()

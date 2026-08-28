"""P247 attempt 0002, milestone 1: curvature-class continuation R = 20..30.

Continues the clean order-16 root chain from the attempt-0001 R=18 root under
the campaign aliasing gate, extracts the radius-free x-space components
C[c], V[c] on every rung R=12..30 (uniform evaluator across the whole ladder),
and applies the interpretation bands frozen in the proposal and in this
attempt's manifest BEFORE execution:

- SATURATION: C[c(R)] approaches a finite limit within the fitted band and
  E(R) turns over or plateaus.
- PERSISTENT_GROWTH: C[c(R)] keeps growing without saturation through R=30.
- OTHER: inconclusive with the named mechanism.

Decision procedure (fixed here, before execution): fit C(R) on the full
R=12..30 set with (a) a linear model and (b) a saturating model
A - B*exp(-R/xi) (bounded fits). Decisive saturation requires the saturating
model to win by dAICc >= 10 with a resolved finite asymptote (xi <= half the
R-span and 0 < A < max(C)) and a vanishing tail derivative; decisive
persistent growth requires a positive slope at R=30 with the linear model
competitive (dAICc > -10) or no resolved asymptote. Otherwise OTHER.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
PREV = HERE.parent / "0001"
sys.path.insert(0, str(PREV))

import debox_common as base  # noqa: E402  (installs P240 0041/0042 paths, pins threads)

import solve_radial_1d  # noqa: E402
import xspace_energy  # noqa: E402

ORDER = base.ORDER
SOLVE_NODES = base.SOLVE_NODES  # (48, 24) committed solve nodes at order 16
ALIAS_NODES = (96, 48)  # committed aliasing gate nodes (clean_ladder.py)
REL_GRAD_GATE = 1.0e-10
ALIAS_GATE = 1.0e-6
COMPONENT_GATE = 1.0e-6
NEW_RADII = [20.0, 22.0, 24.0, 26.0, 28.0, 30.0]


def solve_rung(radius: float, seed: np.ndarray):
    settings = dict(
        radial_order=ORDER,
        radial_nodes=SOLVE_NODES[0],
        angular_nodes=SOLVE_NODES[1],
        radius=radius,
    )
    return solve_radial_1d.solve_order(ORDER, seed, settings)


def alias_energy(values: np.ndarray, radius: float) -> float:
    oracle = solve_radial_1d.Oracle(
        dict(
            radial_order=ORDER,
            radial_nodes=ALIAS_NODES[0],
            angular_nodes=ALIAS_NODES[1],
            radius=radius,
        )
    )
    return oracle.evaluate(np.asarray(values, dtype=np.float64))[0]


def components_at(values, nodes: tuple[int, int]):
    import torch

    flat = torch.tensor(np.asarray(values, dtype=np.float64), dtype=torch.float64)
    curvature, potential, inertia = xspace_energy.xspace_components(
        flat, ORDER, nodes[0], nodes[1]
    )
    return float(curvature), float(potential), float(inertia)


def fit_models(radii: np.ndarray, curvature: np.ndarray) -> dict:
    span = float(radii[-1] - radii[0])

    def linear(r, c0, c1):
        return c0 + c1 * r

    def saturating(r, amp, drop, xi):
        return amp - drop * np.exp(-r / xi)

    lin_p, _ = curve_fit(linear, radii, curvature, p0=(float(curvature[0]), 1.0))
    lin_res = float(np.sum((curvature - linear(radii, *lin_p)) ** 2))
    best_sat = None
    amp_hi = float(np.max(curvature)) * 50.0
    for xi0 in (2.0, 5.0, 10.0, 20.0):
        try:
            p, _ = curve_fit(
                saturating,
                radii,
                curvature,
                p0=(
                    float(np.max(curvature)),
                    float(np.max(curvature)) - float(curvature[0]),
                    xi0,
                ),
                bounds=([0.0, 0.0, 0.5], [amp_hi, amp_hi, span]),
                maxfev=20000,
            )
        except Exception:
            continue
        res = float(np.sum((curvature - saturating(radii, *p)) ** 2))
        if best_sat is None or res < best_sat[1]:
            best_sat = (p, res)
    n = radii.size

    def aicc(res: float, k: int) -> float:
        return n * np.log(res / n) + 2 * k + (2 * k * (k + 1)) / max(1, n - k - 1)

    out = {
        "linear_params": [float(v) for v in lin_p],
        "linear_rss": lin_res,
        "linear_aicc": float(aicc(lin_res, 2)),
    }
    if best_sat is not None:
        (amp, drop, xi), res = best_sat
        out.update(
            saturating_params=[float(amp), float(drop), float(xi)],
            saturating_rss=res,
            saturating_aicc=float(aicc(res, 3)),
            saturating_asymptote=float(amp),
            saturating_xi_over_span=float(xi / span),
            daicc_sat_minus_lin=float(aicc(res, 3) - aicc(lin_res, 2)),
        )
    else:
        out["saturating_fit_failed"] = True
    return out


def decide(fits: dict, radii: np.ndarray, curvature: np.ndarray) -> tuple[str, str]:
    """Returns (band, mechanism). Bands frozen; arithmetic fixed pre-execution."""
    if "saturating_fit_failed" in fits:
        if fits["linear_params"][1] > 0:
            return "PERSISTENT_GROWTH", "saturating fit failed to converge; linear slope positive"
        return "OTHER", "both fit models inconclusive"
    daicc = fits["daicc_sat_minus_lin"]
    xi_rel = fits["saturating_xi_over_span"]
    amp = fits["saturating_asymptote"]
    cmax = float(np.max(curvature))
    asymptote_resolved = xi_rel <= 0.5 and 0.0 < amp < cmax and daicc <= -10.0
    _, drop_param, xi_param = fits["saturating_params"]
    tail_derivative = drop_param * np.exp(-radii[-1] / xi_param) / xi_param
    if asymptote_resolved and tail_derivative <= max(1e-3 * abs(amp), 1e-12):
        return "SATURATION", "resolved finite asymptote wins by dAICc >= 10; tail derivative ~ 0"
    if fits["linear_params"][1] > 0 and (daicc > -10.0 or xi_rel > 0.5):
        return (
            "PERSISTENT_GROWTH",
            "positive slope at R=30; saturating model not decisively better or asymptote unresolved",
        )
    return "OTHER", "fit evidence mixed; no isolation verdict drawable"


def main() -> None:
    started = time.time()
    with open(PREV / "clean-ladder.json") as handle:
        previous = json.load(handle)
    prior = {row["radius"]: row for row in previous["rungs"] if row.get("accepted")}
    last_values = np.asarray(prior[18.0]["values"], dtype=np.float64)

    rows: list[dict] = []
    for radius in NEW_RADII:
        rung_started = time.time()
        solution = solve_rung(radius, last_values)
        accepted = solution["relative_gradient"] < REL_GRAD_GATE
        row: dict = {
            "radius": radius,
            "success": bool(solution["success"]),
            "relative_gradient": solution["relative_gradient"],
            "accepted_root": accepted,
        }
        if accepted:
            values = np.asarray(solution["values"], dtype=np.float64)
            energy = float(solution["energy"])
            alias = alias_energy(values, radius)
            alias_relative = abs(alias - energy) / abs(energy)
            curvature_a, _, _ = components_at(values, (48, 16))
            curvature_b, potential_b, inertia_b = components_at(values, (72, 24))
            component_gap = abs(curvature_a - curvature_b) / max(1e-300, abs(curvature_a))
            phi = 1.0 / (4.0 * inertia_b)
            reconstructed = radius**3 * potential_b + (curvature_b + phi) / radius
            decomposition_gap = abs(reconstructed - energy) / abs(energy)
            a_matrix, _ = base.component_hessians(values, 72, 24)
            lam_a, vector_a = base.lambda_min(a_matrix)
            fractions, _ = solve_radial_1d.analyze_mode(vector_a)
            row.update(
                {
                    "energy_total": energy,
                    "alias_relative_gap": alias_relative,
                    "alias_pass": bool(alias_relative <= ALIAS_GATE),
                    "curvature_C_nodes_48x16": curvature_a,
                    "curvature_C_nodes_72x24": curvature_b,
                    "component_quadrature_gap": component_gap,
                    "component_pass": bool(component_gap <= COMPONENT_GATE),
                    "potential_V": potential_b,
                    "inertia": inertia_b,
                    "decomposition_reconstruction_gap": decomposition_gap,
                    "lambda_min_A": lam_a,
                    "morse_index_A": base.committed_morse_gate(
                        np.linalg.eigvalsh((a_matrix + a_matrix.T) / 2.0)
                    ),
                    "split_fraction": float(fractions[2]),
                    "values": values.tolist(),
                }
            )
            row["accepted"] = bool(
                row["alias_pass"] and row["component_pass"] and decomposition_gap < 1e-8
            )
            if row["accepted"]:
                last_values = values
        else:
            row["accepted"] = False
        row["minutes"] = round((time.time() - rung_started) / 60.0, 2)
        rows.append(row)
        print(
            f"[rung R={radius}] accepted={row.get('accepted')} "
            f"relgrad={row['relative_gradient']:.3e} "
            f"alias={row.get('alias_relative_gap', float('nan')):.3e} "
            f"C={row.get('curvature_C_nodes_72x24', float('nan')):.6e} "
            f"lamA={row.get('lambda_min_A', float('nan')):.6e} "
            f"E={row.get('energy_total', float('nan')):.8f} "
            f"({row['minutes']} min)",
            flush=True,
        )

    # Uniform C, V extraction across the FULL ladder R=12..30 from stored roots.
    potential_map: dict[float, float] = {}
    curvature_list: list[float] = []
    energy_list: list[float] = []
    radii_list: list[float] = []
    for radius in sorted(prior):
        c_b, v_b, _ = components_at(prior[radius]["values"], (72, 24))
        c_a, _, _ = components_at(prior[radius]["values"], (48, 16))
        gap = abs(c_a - c_b) / max(1e-300, abs(c_a))
        if gap > COMPONENT_GATE:
            raise RuntimeError(f"prior rung R={radius} fails component cross-check ({gap:.2e})")
        radii_list.append(float(radius))
        curvature_list.append(c_b)
        potential_map[float(radius)] = v_b
        energy_list.append(float(prior[radius]["energy_total"]))
    for row in rows:
        if row.get("accepted"):
            radii_list.append(float(row["radius"]))
            curvature_list.append(float(row["curvature_C_nodes_72x24"]))
            potential_map[float(row["radius"])] = float(row["potential_V"])
            energy_list.append(float(row["energy_total"]))
    order = np.argsort(radii_list)
    radii = np.asarray(radii_list, dtype=np.float64)[order]
    curvature = np.asarray(curvature_list, dtype=np.float64)[order]
    energies = np.asarray(energy_list, dtype=np.float64)[order]

    fits = fit_models(radii, curvature)
    band, mechanism = decide(fits, radii, curvature)
    slope_fit = np.polyfit(radii, energies, 1)
    payload = {
        "ladder": {
            "radii": radii.tolist(),
            "curvature_C": curvature.tolist(),
            "potential_V": [potential_map[r] for r in radii.tolist()],
            "energies": energies.tolist(),
        },
        "C_fits": fits,
        "E_linear_fit_slope": float(slope_fit[0]),
        "E_linear_fit_intercept": float(slope_fit[1]),
        "band": band,
        "mechanism": mechanism,
        "new_rungs": [{k: v for k, v in row.items() if k != "values"} for row in rows],
        "minutes_total": round((time.time() - started) / 60.0, 2),
    }
    (HERE / "m1-continuation.json").write_text(json.dumps(payload, indent=2))
    print(f"C(R) fits: {json.dumps(fits)}")
    print(f"E(R) linear slope on R=12..30: {payload['E_linear_fit_slope']:.6f}")
    print(f"MILESTONE-1 BAND: {band} ({mechanism})")
    print(f"WROTE m1-continuation.json runtime={payload['minutes_total']} min")


if __name__ == "__main__":
    main()

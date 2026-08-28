"""P247 attempt 0001, gate G2: clean de-boxing ladder R = 12..18.

Reproduces the reported 0044 ladder (PR #153 / issue #151 comments, transfer
references) from committed merged machinery:

- clean roots at order 16 on the committed biaxial-hedgehog functional,
  continued in R from the committed largeR-roots.json backgrounds,
  accepted only at relative gradient < 1e-10;
- the aliasing gate reported for 0044 and adopted as a campaign-wide rule:
  every accepted root reproduces its energy at doubled quadrature to a
  relative 1e-6;
- the potential-curvature operator A = d^2 V[c] on each clean background,
  quadrature-checked between two node sets;
- lambda_min(A), the committed Morse-index convention, and the committed
  channel-fraction decomposition of the lowest eigenvector.

Frozen gates (memory contract): slope in [0.25, 0.45] per unit R with
R^2 >= 0.999; lambda_min(A) < 0 at every rung; |lambda_min(A)| at the
outermost rung in [1e-7, 2e-6]; Morse index exactly 1; split-channel
fraction > 0.9.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from debox_common import (
    CERT_NODES_A,
    CERT_NODES_B,
    ORDER,
    REPORTED,
    SOLVE_NODES,
    component_hessians,
    committed_morse_gate,
    lambda_min,
    load_committed_roots,
)

HERE = Path(__file__).resolve().parent

RADII = [12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0]
ALIAS_NODES = (96, 48)
REL_GRAD_GATE = 1.0e-10
ALIAS_GATE = 1.0e-6
ACURACY_GATE = 1.0e-6


def solve_rung(radius: float, seed: np.ndarray):
    import solve_radial_1d

    settings = dict(
        radial_order=ORDER, radial_nodes=SOLVE_NODES[0], angular_nodes=SOLVE_NODES[1], radius=radius
    )
    solution = solve_radial_1d.solve_order(ORDER, seed, settings)
    return solution


def alias_energy(values: np.ndarray, radius: float, nodes: tuple[int, int]) -> float:
    import solve_radial_1d

    oracle = solve_radial_1d.Oracle(
        dict(radial_order=ORDER, radial_nodes=nodes[0], angular_nodes=nodes[1], radius=radius)
    )
    return oracle.evaluate(np.asarray(values, dtype=np.float64))[0]


def main() -> None:
    started = time.time()
    roots = load_committed_roots()
    seeds = {
        12.0: np.asarray(roots["R12"]["values"], dtype=np.float64),
        14.0: np.asarray(roots["R14"]["values"], dtype=np.float64),
    }
    rows: list[dict] = []
    backgrounds: dict[str, dict[str, np.ndarray]] = {}
    last_values: np.ndarray | None = None

    for radius in RADII:
        rung_started = time.time()
        seed = seeds.get(radius, last_values)
        if seed is None:
            raise RuntimeError(f"no seed for radius {radius}")
        solution = solve_rung(radius, seed)
        # Acceptance follows the committed phase1 gate (relative gradient
        # < 1e-10); hybr's success flag is recorded but not gating, because at
        # xtol=1e-14 it can report no-progress at the rounding floor of a
        # fully converged root (observed: relgrad 1.17e-14, success=False).
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
            row["values"] = values.tolist()
            alias = alias_energy(values, radius, ALIAS_NODES)
            alias_relative = abs(alias - energy) / abs(energy)
            a_a, _ = component_hessians(values, *CERT_NODES_A)
            a_b, d_b = component_hessians(values, *CERT_NODES_B)
            lam_a, vector_a = lambda_min(a_b)
            lam_a_coarse, _ = lambda_min(a_a)
            accuracy = abs(lam_a - lam_a_coarse) / max(1.0e-300, abs(lam_a), abs(lam_a_coarse))
            import solve_radial_1d

            eigenvalues_full = np.linalg.eigvalsh((a_b + a_b.T) / 2.0)
            morse_index = committed_morse_gate(eigenvalues_full)
            fractions, nodes = solve_radial_1d.analyze_mode(vector_a)
            row.update(
                {
                    "energy_total": energy,
                    "static_energy": solution["components"]["static"],
                    "inertia": solution["components"]["inertia"],
                    "omega": solution["components"]["frequency"],
                    "alias_energy_doubled_quad": alias,
                    "alias_relative_gap": alias_relative,
                    "alias_pass": bool(alias_relative <= ALIAS_GATE),
                    "lambda_min_A_nodes_48x16": lam_a_coarse,
                    "lambda_min_A_nodes_72x24": lam_a,
                    "lambda_min_A_quadrature_gap": accuracy,
                    "lambda_min_A_pass": bool(accuracy <= ACURACY_GATE),
                    "morse_index_A": morse_index,
                    "mode_fractions_q_tangent_split": [float(f) for f in fractions],
                    "split_fraction": float(fractions[2]),
                    "split_radial_nodes": int(nodes),
                    "lambda_min_total_at_root": float(solution["lambda_min"]),
                    "morse_index_total_at_root": int(solution["morse_index"]),
                }
            )
            row["accepted"] = bool(
                row["alias_pass"] and row["lambda_min_A_pass"]
            )
            if row["accepted"]:
                backgrounds[f"R{int(radius)}"] = {"A": a_b, "D": d_b}
                last_values = values
        else:
            row["accepted"] = False
        row["minutes"] = round((time.time() - rung_started) / 60.0, 2)
        rows.append(row)
        print(f"[rung R={radius}] accepted={row.get('accepted')} "
              f"relgrad={row['relative_gradient']:.3e} "
              f"alias={row.get('alias_relative_gap', float('nan')):.3e} "
              f"lamA={row.get('lambda_min_A_nodes_72x24', float('nan')):.6e} "
              f"index={row.get('morse_index_A')} split={row.get('split_fraction', float('nan')):.6f} "
              f"E={row.get('energy_total', float('nan')):.8f} ({row['minutes']} min)", flush=True)

    accepted_rows = [row for row in rows if row.get("accepted")]
    radii = np.array([row["radius"] for row in accepted_rows])
    energies = np.array([row["energy_total"] for row in accepted_rows])
    slope, intercept, r_squared = float("nan"), None, float("nan")
    if radii.size >= 3:
        fit_slope, fit_intercept = np.polyfit(radii, energies, 1)
        fitted = fit_slope * radii + fit_intercept
        ss_res = float(np.sum((energies - fitted) ** 2))
        ss_tot = float(np.sum((energies - energies.mean()) ** 2))
        slope, intercept, r_squared = float(fit_slope), float(fit_intercept), 1.0 - ss_res / ss_tot

    outermost = accepted_rows[-1] if accepted_rows else None
    gates = {
        "slope_band": {
            "slope": float(slope),
            "band": REPORTED["slope_band"],
            "pass": bool(REPORTED["slope_band"][0] <= slope <= REPORTED["slope_band"][1]),
        },
        "r_squared": {
            "r_squared": float(r_squared),
            "min": REPORTED["r2_min"],
            "pass": bool(r_squared >= REPORTED["r2_min"]),
        },
        "lambda_min_negative_all_rungs": {
            "values": [row["lambda_min_A_nodes_72x24"] for row in accepted_rows],
            "pass": bool(all(row["lambda_min_A_nodes_72x24"] < 0.0 for row in accepted_rows)),
        },
        "flattening_band_outermost": {
            "value": outermost["lambda_min_A_nodes_72x24"] if outermost else None,
            "band": REPORTED["lambda_min_A_flattening_band"],
            "pass": bool(
                outermost is not None
                and REPORTED["lambda_min_A_flattening_band"][0]
                <= abs(outermost["lambda_min_A_nodes_72x24"])
                <= REPORTED["lambda_min_A_flattening_band"][1]
            ),
        },
        "morse_index_one_all_rungs": {
            "values": [row["morse_index_A"] for row in accepted_rows],
            "pass": bool(all(row["morse_index_A"] == REPORTED["morse_index"] for row in accepted_rows)),
        },
        "split_fraction_all_rungs": {
            "values": [row["split_fraction"] for row in accepted_rows],
            "min": REPORTED["split_fraction_min"],
            "pass": bool(all(row["split_fraction"] > REPORTED["split_fraction_min"] for row in accepted_rows)),
        },
    }
    gates_all_pass = all(entry["pass"] for entry in gates.values())
    payload = {
        "gate": "G2_ladder",
        "order": ORDER,
        "solve_nodes": SOLVE_NODES,
        "alias_nodes": ALIAS_NODES,
        "certification_nodes": [CERT_NODES_A, CERT_NODES_B],
        "rel_grad_gate": REL_GRAD_GATE,
        "alias_gate_relative": ALIAS_GATE,
        "committed_seed_provenance": "largeR-roots.json R12/R14; continuation by nearest accepted root",
        "rungs": rows,
        "accepted_rung_count": len(accepted_rows),
        "linear_fit": {"slope": float(slope), "intercept": float(intercept) if radii.size >= 3 else None,
                        "r_squared": float(r_squared)},
        "reported_slope_transfer_reference": REPORTED["slope_per_unit_R"],
        "gates": gates,
        "gates_all_pass": gates_all_pass,
        "verdict": "ESTABLISHED" if gates_all_pass else "MIXED",
        "runtime_seconds": round(time.time() - started, 1),
    }
    (HERE / "clean-ladder.json").write_text(json.dumps(payload, indent=2) + "\n")
    np.savez(HERE / "backgrounds.npz", **{
        f"{tag}_{name}": matrix for tag, pair in backgrounds.items() for name, matrix in pair.items()
    })
    print(json.dumps({"gates": gates, "linear_fit": payload["linear_fit"],
                      "verdict": payload["verdict"]}, indent=2))
    print(f"G2 VERDICT: {payload['verdict']}")


if __name__ == "__main__":
    main()

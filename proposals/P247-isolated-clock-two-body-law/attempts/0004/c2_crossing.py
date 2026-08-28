"""P247 attempt 0004, gate C2: order-24 bisection of the static crossing.

Bracket from attempt-0003 (g9-w1w3-order24.json): lambda_min(A) < 0 at
R = 20 (-2.078540e-09, above budget, q-dominated) and > 0 at R = 22
(+4.297582e-10, above budget, split-channel). Bisection at order 24 with
the attempt-0003 protocol: R-continuation seeding from the negative-side
bracket end's order-24 root (roots24.npz), root-continuity gate 0.25,
aliasing gate 1e-6, relgrad 1e-10, lambda_min(A) at (108,36) with doubled
(216,72) quadrature gauge and jitter sign tests. Stop after 3 midpoints
(final bracket width 0.25). Continuity reference for intermediate radii:
linear interpolation of the bracket ends' order-24 energies (no order-16
root exists at intermediate radii; the bracket ends were continuity-
accepted at order 24 in attempt 0003).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "0002"))
sys.path.insert(0, str(HERE.parent / "0001"))

import debox_common as base  # noqa: E402

import order_study as osd  # noqa: E402
import solve_radial_1d  # noqa: E402

REL_GRAD_GATE = 1.0e-10
ALIAS_GATE = 1.0e-6
CONTINUITY_GATE = 0.25
ORDER = 24
PRIMARY = (108, 36)
DOUBLED = (216, 72)
MIDPOINTS = 3


def channel_split(matrix: np.ndarray) -> dict:
    symmetric = (np.asarray(matrix, dtype=np.float64) + np.asarray(matrix).T) / 2.0
    _eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    fractions, nodes = solve_radial_1d.analyze_mode(eigenvectors[:, 0])
    return {
        "fractions_q_tangent_split": [float(f) for f in fractions],
        "split_radial_nodes": int(nodes),
    }


def main() -> None:
    started = time.time()
    roots = np.load(HERE.parent / "0003" / "roots24.npz")
    with open(HERE.parent / "0003" / "g9-w1w3-order24.json") as handle:
        prior = json.load(handle)
    e24 = {
        float(r): row["energy"]
        for r, row in prior["roots24"].items()
    }
    lo, hi = 20.0, 22.0
    lo_sign = "negative"
    lo_values = roots["R20"]
    values_by_radius = {20.0: lo_values, 22.0: roots["R22"]}

    records = []
    for _level in range(MIDPOINTS):
        mid = 0.5 * (lo + hi)
        reference = e24[lo] + (e24[hi] - e24[lo]) * (mid - lo) / (hi - lo)
        solution = osd.solve_at(ORDER, mid, lo_values)
        rel_grad = solution["relative_gradient"]
        if not rel_grad < REL_GRAD_GATE:
            records.append({"radius": mid, "converged": False, "relative_gradient": rel_grad})
            print(f"[C2 R={mid:.3f}] solve failed relgrad={rel_grad:.2e}", flush=True)
            break
        values = np.asarray(solution["values"], dtype=np.float64)
        energy = float(solution["energy"])
        alias_relative = abs(osd.alias_energy(values, ORDER, mid) - energy) / abs(energy)
        if not alias_relative < ALIAS_GATE:
            records.append({"radius": mid, "converged": False, "alias_gap": alias_relative})
            print(f"[C2 R={mid:.3f}] alias failed gap={alias_relative:.2e}", flush=True)
            break
        continuity = abs(energy - reference) / abs(reference)
        if not continuity < CONTINUITY_GATE:
            records.append({"radius": mid, "converged": False, "continuity": continuity})
            print(f"[C2 R={mid:.3f}] continuity failed {continuity:.2e}", flush=True)
            break
        a_primary, _ = base.component_hessians(values, PRIMARY[0], PRIMARY[1], order=ORDER)
        metrics = osd.soft_metrics(a_primary)
        a_doubled, _ = base.component_hessians(values, DOUBLED[0], DOUBLED[1], order=ORDER)
        lam_doubled = base.lambda_min(a_doubled)[0]
        quad_gauge = abs(lam_doubled - metrics["lambda_min"])
        jit = osd.jitter_spread(values, ORDER, PRIMARY)
        budget = quad_gauge + jit["spread"]
        sign_stable = (jit["min"] > 0.0 and metrics["lambda_min"] > 0.0) or (
            jit["max"] < 0.0 and metrics["lambda_min"] < 0.0
        )
        resolved = bool(abs(metrics["lambda_min"]) > 10.0 * budget and sign_stable)
        sign = "unresolved"
        if resolved:
            sign = "negative" if metrics["lambda_min"] < 0.0 else "positive"
        record = {
            "radius": mid,
            "converged": True,
            "energy": energy,
            "continuity": continuity,
            "alias_relative": alias_relative,
            **metrics,
            "lambda_min_doubled": lam_doubled,
            "budget": budget,
            "sign": sign,
            "softest_channel": channel_split(a_primary),
        }
        records.append(record)
        values_by_radius[mid] = values
        print(
            f"[C2 R={mid:.3f}] lamA={metrics['lambda_min']:+.6e} sign={sign} "
            f"budget={budget:.1e} cont={continuity:.2e}",
            flush=True,
        )
        if sign == "unresolved":
            break
        if sign == lo_sign:
            lo = mid
            lo_values = values
            e24[mid] = energy
        else:
            hi = mid
            e24[mid] = energy

    payload = {
        "final_bracket": [lo, hi],
        "records": [
            {k: v for k, v in row.items() if k != "values"} for row in records
        ],
        "minutes_total": round((time.time() - started) / 60.0, 2),
    }
    (HERE / "c2-crossing.json").write_text(json.dumps(payload, indent=2))
    print(f"C2 final bracket: [{lo}, {hi}], width {hi - lo:.3f}")
    print("WROTE c2-crossing.json")


if __name__ == "__main__":
    main()

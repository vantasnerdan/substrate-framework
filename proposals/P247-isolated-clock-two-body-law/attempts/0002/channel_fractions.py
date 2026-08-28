"""P247 attempt 0002 addendum: channel attribution of the refined negative mode.

Re-solves the deterministic order-24 roots at R=14 and R=18 (identical seed
and node protocol to order_study.py) and reports the channel fractions
(q, tangent, split) and radial nodal counts of the three lowest eigenvectors
of A. Preregistered adjudication target: the reported 0044 claim that the
negative direction is 'purely split-channel' with Morse index 1.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PREV = HERE.parent / "0001"
sys.path.insert(0, str(PREV))

import debox_common as base  # noqa: E402

import solve_radial_1d  # noqa: E402

ORDER = 24
NODES = (72, 36)
RADII = [14.0, 18.0]


def main() -> None:
    started = time.time()
    with open(PREV / "clean-ladder.json") as handle:
        previous = json.load(handle)
    prior = {row["radius"]: row for row in previous["rungs"] if row.get("accepted")}

    # reuse the exact reprojection + solver protocol from order_study
    from order_study import reproject  # noqa: E402  (same directory)

    records: dict = {}
    for radius in RADII:
        root16 = np.asarray(prior[radius]["values"], dtype=np.float64)
        seed = reproject(root16, ORDER)
        settings = dict(
            radial_order=ORDER, radial_nodes=NODES[0], angular_nodes=NODES[1], radius=radius
        )
        solution = solve_radial_1d.solve_order(ORDER, seed, settings)
        rel_grad = solution["relative_gradient"]
        if not rel_grad < 1.0e-10:
            print(f"[R={radius}] solve FAILED relgrad={rel_grad:.3e}", flush=True)
            continue
        values = np.asarray(solution["values"], dtype=np.float64)
        a_matrix, _ = base.component_hessians(values, 108, 36, order=ORDER)
        symmetric = (a_matrix + a_matrix.T) / 2.0
        eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
        rows = []
        for index in range(3):
            vector = eigenvectors[:, index]
            fractions, nodes = solve_radial_1d.analyze_mode(vector)
            rows.append(
                {
                    "eigenvalue": float(eigenvalues[index]),
                    "fractions_q_tangent_split": [float(f) for f in fractions],
                    "split_radial_nodes": int(nodes),
                }
            )
            print(
                f"[R={radius} mode {index}] lam={eigenvalues[index]:+.6e} "
                f"fractions={np.round(fractions, 6).tolist()} nodes={nodes}",
                flush=True,
            )
        records[str(radius)] = {
            "relative_gradient": rel_grad,
            "morse_gate": base.committed_morse_gate(eigenvalues),
            "modes": rows,
            "values": values.tolist(),
        }

    payload = {
        "records": records,
        "order": ORDER,
        "nodes": list(NODES),
        "minutes_total": round((time.time() - started) / 60.0, 2),
    }
    (HERE / "g8-channel-fractions.json").write_text(json.dumps(payload, indent=2))
    print("WROTE g8-channel-fractions.json")


if __name__ == "__main__":
    main()

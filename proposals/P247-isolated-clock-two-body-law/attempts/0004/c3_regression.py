"""P247 attempt 0004, gate C3-0: chi=0 regression of the extended functional.

The extended 4x4 functional must reproduce the committed 3x3 functional
exactly at chi = 0 on a committed root: curvature, potential, inertia
(relative 1e-9) and the derived fixed-J and omega.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "0002"))
sys.path.insert(0, str(HERE.parent / "0001"))

import c3_functional as c3  # noqa: E402
import solve_radial_1d  # noqa: E402

ORDER = 16
NODES = (48, 24)


def components_committed(values: np.ndarray, nodes, radius: float):
    oracle = solve_radial_1d.Oracle(
        dict(
            radial_order=ORDER,
            radial_nodes=nodes[0],
            angular_nodes=nodes[1],
            radius=radius,
        )
    )
    result = oracle.evaluate(values)
    total, components = result[0], result[3]
    return float(total), {k: float(v) for k, v in components.items()}


def main() -> None:
    with open(HERE.parent / "0001" / "clean-ladder.json") as handle:
        previous = json.load(handle)
    row18 = [
        row
        for row in previous["rungs"]
        if row.get("accepted") and row["radius"] == 18.0
    ][0]
    root16 = np.asarray(row18["values"], dtype=np.float64)

    radius = 18.0
    total_c, comp_c = components_committed(root16, NODES, radius)

    flat = np.zeros((4, ORDER), dtype=np.float64)
    flat[:3] = root16.reshape(3, ORDER)
    flat_tensor = torch.tensor(flat, dtype=torch.float64, requires_grad=False)
    total_e, comp_e = c3.extended_static(
        flat_tensor,
        radial_order=ORDER,
        radial_nodes=NODES[0],
        angular_nodes=NODES[1],
        radius=radius,
    )

    record = {
        "committed": {"total": total_c, **comp_c},
        "extended_chi0": {k: float(v) for k, v in comp_e.items()},
        "relative_diffs": {
            key: abs(float(comp_e[key]) - comp_c[key])
            / max(abs(comp_c[key]), 1e-300)
            for key in ("curvature", "potential", "inertia")
            if key in comp_c
        },
    }
    record["gates_pass"] = bool(
        all(v < 1.0e-9 for v in record["relative_diffs"].values())
    )
    (HERE / "c3-regression.json").write_text(json.dumps(record, indent=2))
    print(json.dumps(record, indent=2))
    print("C3-0 REGRESSION:", "PASS" if record["gates_pass"] else "FAIL")


if __name__ == "__main__":
    main()

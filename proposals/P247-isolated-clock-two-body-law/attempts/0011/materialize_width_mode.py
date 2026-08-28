"""P247 0011 correction pass: materialize repaired width-mode evidence.

Writes attempts/0009/width-mode.json (the repaired-class bottom-mode
characterization at the R=24 root) plus the family width ladder over
the actual minimum-branch roots R=20..30.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ATT9 = HERE.parent / "0009"
ATT10 = HERE.parent / "0010"
sys.path.insert(0, str(HERE.parent / "0004"))
sys.path.insert(0, str(ATT10))

from narrow_branch_final import REPAIR_MARKER, REPAIRED  # noqa: E402

ORDER = 24


def main() -> None:
    src = (HERE.parent / "0004" / "c3_functional.py").read_text()
    assert src.count(REPAIR_MARKER) == 1
    src = src.replace(REPAIR_MARKER, REPAIRED)
    src = src.replace("def extended_static(", "def ext_rep(", 1)
    src = src.replace("__file__",
                      repr(str((HERE.parent / "0004" / "c3_functional.py").resolve())))
    ns: dict = {"__name__": "rep"}
    exec(compile(src, "rep", "exec"), ns)

    coef24 = np.asarray(np.load(HERE.parent / "0003" / "roots24.npz")["R24"],
                        dtype=np.float64).reshape(3, ORDER)
    flat = torch.zeros((4, ORDER), dtype=torch.float64)
    flat[:3] = torch.tensor(coef24)
    flat.requires_grad_(True)
    total, comp = ns["ext_rep"](flat, radial_order=ORDER, radial_nodes=48,
                                angular_nodes=24, radius=24.0,
                                chi_angular_degrees=(0,))
    total = total + comp["fixed_j"]
    g = torch.autograd.grad(total, flat, create_graph=True)[0].reshape(-1)
    rows = [torch.autograd.grad(g[i], flat, retain_graph=True)[0].reshape(-1)
            for i in range(g.numel())]
    H = torch.stack(rows).detach().numpy()
    H = (H + H.T) / 2.0
    ev, evec = np.linalg.eigh(H)

    coef30 = np.asarray(np.load(HERE.parent / "0003" / "roots24.npz")["R30"],
                        dtype=np.float64).reshape(3, ORDER)
    dc = np.concatenate([(coef30 - coef24).reshape(-1), np.zeros(ORDER)])
    dc_hat = dc / np.linalg.norm(dc)

    from fixed_width_probe import width_proxy

    roots = np.load(HERE.parent / "0003" / "roots24.npz")
    widths = {}
    for name in roots.files:
        c = np.asarray(roots[name], dtype=np.float64).reshape(3, ORDER)
        R = float(name[1:])
        widths[name] = width_proxy(c, R)

    payload = {
        "class": "sigma-repaired extended functional",
        "root": "R24 (roots24.npz), quadrature (48, 24)",
        "spectrum_bottom6": [float(v) for v in ev[:6]],
        "lambda_min": float(ev[0]),
        "num_negative": int((ev < 0).sum()),
        "bottom_mode": {
            "row_fractions": {
                "m0": float(np.linalg.norm(evec[:ORDER, 0])),
                "m1": float(np.linalg.norm(evec[ORDER:2 * ORDER, 0])),
                "m2": float(np.linalg.norm(evec[2 * ORDER:3 * ORDER, 0])),
                "chi": float(np.linalg.norm(evec[3 * ORDER:, 0])),
            },
            "split_row_m2_fraction": float(
                np.linalg.norm(evec[2 * ORDER:3 * ORDER, 0])
                / np.linalg.norm(evec[:3 * ORDER, 0])),
            "family_drift_overlap_R24_to_R30": float(
                np.dot(evec[:, 0], dc_hat)),
            "second_eigenvalue": float(ev[1]),
        },
        "family_core_width_ladder": widths,
        "width_proxy_definition": "r where the summed |lambda_i - lambda_i(boundary)| first falls below half its maximum, measured from the peak",
    }
    (ATT9 / "width-mode.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload["bottom_mode"], indent=1))
    print("family widths:", {k: round(v, 3) for k, v in widths.items()})
    print("WROTE", ATT9 / "width-mode.json")


if __name__ == "__main__":
    main()

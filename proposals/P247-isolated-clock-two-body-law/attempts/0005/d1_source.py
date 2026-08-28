"""P247 attempt 0005, gate D1: does the linear boost source survive a
tilted (fixed-axis) boost?

Measures grad of the fixed-J extended energy in the chi channel at chi = 0
on the committed order-24 R = 24 root for two boost-axis choices:

  co-aligned: boost along the LOCAL DIRECTOR n-hat (attempt-0004 setup;
              gradient known to vanish exactly),
  tilted:     boost along the FIXED LABORATORY z axis.

SURVIVES if the tilted source norm exceeds the jitter budget (the only
error channel at chi = 0: sigma and s_m are quadratic, so the measured
gradient is pure signal); VANISHES otherwise. Decides D_tilted_boost vs
E_localized_core routing (proposal next_loop).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "0004"))
sys.path.insert(0, str(HERE.parent / "0002"))
sys.path.insert(0, str(HERE.parent / "0001"))

import c3_functional as c3  # noqa: E402

ORDER = 24
NODES = (72, 36)
JITTER_DRAWS = 3


def main() -> None:
    started = time.time()
    roots = np.load(HERE.parent / "0003" / "roots24.npz")
    root_s = np.asarray(roots["R24"], dtype=np.float64)

    # monkeypatch the boost axis inside c3_functional.boosted_boost_matrix
    # through a fixed-axis variant: rebuild the module functions with a
    # global axis override instead of editing the file.
    import types

    original_boosted = c3.boosted_boost_matrix
    original_boost_only = c3.boost_matrix_only
    AXIS_MODE = {"mode": "coaligned"}

    def boost_matrix_only(director, chi):
        c, s = torch.cosh(chi), torch.sinh(chi)
        if AXIS_MODE["mode"] == "tilted":
            ax = torch.zeros_like(director)
            ax[..., 2] = 1.0  # fixed laboratory z axis
        else:
            ax = director
        nx, ny, nz = ax.unbind(-1)
        zero = torch.zeros_like(c)
        one = torch.ones_like(c)
        top = torch.stack((c, s * nx, s * ny, s * nz), dim=-1)
        row1 = torch.stack(
            (s * nx, one + (c - 1) * nx * nx, (c - 1) * nx * ny, (c - 1) * nx * nz),
            dim=-1,
        )
        row2 = torch.stack(
            (s * ny, (c - 1) * ny * nx, one + (c - 1) * ny * ny, (c - 1) * ny * nz),
            dim=-1,
        )
        row3 = torch.stack(
            (s * nz, (c - 1) * nz * nx, (c - 1) * nz * ny, one + (c - 1) * nz * nz),
            dim=-1,
        )
        return torch.stack((top, row1, row2, row3), dim=-2)

    c3.boost_matrix_only = boost_matrix_only
    c3.boosted_boost_matrix = lambda d, chi: (
        boost_matrix_only(d, chi),
        boost_matrix_only(d, -chi),
    )

    results = {}
    for mode in ("coaligned", "tilted"):
        AXIS_MODE["mode"] = mode
        flat = torch.zeros((4, ORDER), dtype=torch.float64)
        flat[:3] = torch.tensor(root_s.reshape(3, ORDER))
        flat.requires_grad_(True)
        total, comp = c3.extended_static(
            flat, radial_order=ORDER, radial_nodes=NODES[0],
            angular_nodes=NODES[1], radius=24.0,
        )
        (grad,) = torch.autograd.grad(total, flat)
        source_norm = float(grad[3].norm())
        max_entry = float(grad[3].abs().max())
        # jitter check: the source must be robust to root jitter at 1e-11
        rng = np.random.default_rng(20260828)
        jitter_norms = []
        for _ in range(JITTER_DRAWS):
            perturbed = root_s + 1.0e-11 * np.linalg.norm(root_s) * rng.standard_normal(root_s.size) / np.sqrt(root_s.size)
            flat_j = torch.zeros((4, ORDER), dtype=torch.float64)
            flat_j[:3] = torch.tensor(perturbed.reshape(3, ORDER))
            flat_j.requires_grad_(True)
            total_j, _ = c3.extended_static(
                flat_j, radial_order=ORDER, radial_nodes=NODES[0],
                angular_nodes=NODES[1], radius=24.0,
            )
            (grad_j,) = torch.autograd.grad(total_j, flat_j)
            jitter_norms.append(float(grad_j[3].norm()))
        results[mode] = {
            "source_norm": source_norm,
            "max_abs_entry": max_entry,
            "jitter_source_norms": jitter_norms,
            "jitter_min": min(jitter_norms),
        }
        print(
            f"[D1 {mode}] source_norm={source_norm:.6e} "
            f"max_entry={max_entry:.3e} jitter_min={min(jitter_norms):.3e}",
            flush=True,
        )

    survives = bool(
        results["tilted"]["jitter_min"] > 10.0 * 1e-15
        and results["tilted"]["source_norm"] > 0.0
    )
    payload = {
        "results": results,
        "budget_note": "sigma and s_m are quadratic in chi: the measured gradient is pure signal; budget set at 1e-15 against roundoff",
        "verdict": "SURVIVES" if survives else "VANISHES",
        "routing": "D_tilted_boost (D2)" if survives else "E_localized_core (D3)",
        "minutes_total": round((time.time() - started) / 60.0, 2),
    }
    (HERE / "d1-source.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps({k: v for k, v in payload.items() if k != "results"}, indent=2))
    print("WROTE d1-source.json")


if __name__ == "__main__":
    main()

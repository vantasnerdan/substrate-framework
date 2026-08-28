"""P247 attempt 0010: finalize narrow-branch diagnostics from saved roots."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "0004"))
sys.path.insert(0, str(HERE.parent / "0001"))
sys.path.insert(0, str(HERE.parent / "0002"))

import debox_common as base  # noqa: F401,E402
import solve_radial_1d  # noqa: E402
from fixed_width_probe import ORDER, width_proxy  # noqa: E402

REPAIR_MARKER = """    sigma_density = KAPPA * (
        -torch.sum(du_r * (du_r @ ETA), dim=-1)
        - torch.sum(du_theta * (du_theta @ ETA), dim=-1)"""
REPAIRED = """    sigma_density = KAPPA * (
        torch.sum(du_r * (du_r @ ETA), dim=-1)
        + torch.sum(du_theta * (du_theta @ ETA), dim=-1)"""


def repaired_hessian_min(coef: np.ndarray, radius: float) -> dict:
    src = (HERE.parent / "0004" / "c3_functional.py").read_text()
    assert src.count(REPAIR_MARKER) == 1
    src = src.replace(REPAIR_MARKER, REPAIRED)
    src = src.replace("def extended_static(", "def ext_rep(", 1)
    src = src.replace("__file__", repr(str(HERE.parent / "0004" / "c3_functional.py")))
    ns: dict = {"__name__": "rep"}
    exec(compile(src, "rep", "exec"), ns)
    flat = torch.zeros((4, ORDER), dtype=torch.float64)
    flat[:3] = torch.tensor(coef.reshape(3, ORDER))
    flat.requires_grad_(True)
    total, comp = ns["ext_rep"](flat, radial_order=ORDER, radial_nodes=48,
                                angular_nodes=24, radius=radius,
                                chi_angular_degrees=(0,))
    total = total + comp["fixed_j"]
    g = torch.autograd.grad(total, flat, create_graph=True)[0].reshape(-1)
    rows = [torch.autograd.grad(g[i], flat, retain_graph=True)[0].reshape(-1)
            for i in range(g.numel())]
    H = torch.stack(rows).detach().numpy()
    H = (H + H.T) / 2.0
    ev = np.linalg.eigvalsh(H)
    return {"E_repaired": float(total.detach()), "lambda_min": float(ev[0]),
            "num_negative": int((ev < 0).sum()),
            "bottom5": [float(v) for v in ev[:5]]}


def main() -> None:
    started = time.time()
    roots = np.load(HERE.parent / "0003" / "roots24.npz")
    payload: dict = {}

    coef_prev = None
    for R in (26, 28):
        path = HERE / f"narrow{R}.npy"
        if path.exists():
            coef = np.load(path)
        else:
            from fixed_width_probe import embed_and_reproject
            seed = embed_and_reproject(coef_prev, float(R - 2), float(R))
            print(f"[R{R}] solving from embedded narrow R{R-2} seed...", flush=True)
            sol = solve_radial_1d.solve_order(ORDER, seed.reshape(-1),
                                              dict(radial_nodes=48, angular_nodes=24,
                                                   radius=float(R)))
            coef = np.asarray(sol["values"], dtype=np.float64).reshape(3, ORDER)
            np.save(path, coef)
            print(f"[R{R}] solved: E={sol['energy']:.6f} "
                  f"relgrad={sol['relative_gradient']:.3e}", flush=True)
        coef_prev = coef
        oracle = solve_radial_1d.Oracle(dict(radial_order=ORDER, radial_nodes=48,
                                             angular_nodes=24, radius=float(R)))
        total, grad, _hess, _comp = oracle.evaluate(coef.reshape(-1))
        rel = float(np.max(np.abs(grad)) / max(1.0, abs(float(total))))
        rec = {
            "E_solver_functional": float(total),
            "relative_gradient": rel,
            "width": width_proxy(coef, float(R)),
            "repaired_stability": repaired_hessian_min(coef, float(R)),
        }
        payload[f"R{R}"] = rec
        print(f"[R{R}] E={rec['E_solver_functional']:.6f} relgrad={rel:.3e} "
              f"width={rec['width']:.3f}", flush=True)
        print(f"[R{R}] repaired:", rec["repaired_stability"], flush=True)
        payload[f"R{R}"] = rec

    oracle26 = solve_radial_1d.Oracle(dict(radial_order=ORDER, radial_nodes=48,
                                           angular_nodes=24, radius=26.0))
    payload["family_R26_energy"] = float(oracle26.evaluate(
        np.asarray(roots["R26"], dtype=np.float64).reshape(-1))[0])
    print("family R26 energy:", payload["family_R26_energy"], flush=True)

    payload["minutes_total"] = round((time.time() - started) / 60.0, 2)
    (HERE / "narrow-branch.json").write_text(json.dumps(payload, indent=2, default=float))
    print("WROTE narrow-branch.json")


if __name__ == "__main__":
    main()

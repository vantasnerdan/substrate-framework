"""P247 attempt 0010: narrow-branch tracing + repaired stability.

Follows the fixed-width probe's R=26 narrow convergence:
  1. Re-solve at R=26 with full diagnostics (relative_gradient,
     solver Morse index, lambda_min).
  2. Repaired-Hessian Morse index at the narrow R=26 solution
     (the attempt-0009 sigma-repair machinery).
  3. Seed R=28 from the narrow R=26 solution; solve; diagnose.
Artifacts: narrow26.npy, narrow28.npy, narrow-branch.json
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
sys.path.insert(0, str(HERE.parent / "0001"))
sys.path.insert(0, str(HERE.parent / "0002"))

import debox_common as base  # noqa: F401,E402
import solve_radial_1d  # noqa: E402
from fixed_width_probe import ORDER, embed_and_reproject, width_proxy  # noqa: E402


def repaired_hessian_min(coef: np.ndarray, radius: float):
    """Full repaired-class Hessian (sigma sign fix) at a 3-row root."""
    import c3_functional as c3  # noqa: F401  (ensures module context)

    src = (HERE.parent / "0004" / "c3_functional.py").read_text()
    marker = """    sigma_density = KAPPA * (
        -torch.sum(du_r * (du_r @ ETA), dim=-1)
        - torch.sum(du_theta * (du_theta @ ETA), dim=-1)"""
    assert src.count(marker) == 1
    src = src.replace(marker, """    sigma_density = KAPPA * (
        torch.sum(du_r * (du_r @ ETA), dim=-1)
        + torch.sum(du_theta * (du_theta @ ETA), dim=-1)""")
    src = src.replace("def extended_static(", "def ext_rep(", 1)
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


def solve_diagnose(coef_seed: np.ndarray, radius: float, tag: str):
    sol = solve_radial_1d.solve_order(ORDER, coef_seed.reshape(-1),
                                      dict(radial_nodes=48, angular_nodes=24,
                                           radius=radius))
    values = np.asarray(sol["values"], dtype=np.float64).reshape(3, ORDER)
    rec = {
        "success": sol["success"],
        "E_converged": sol["energy"],
        "relative_gradient": sol["relative_gradient"],
        "solver_morse_index_unrepaired": sol["morse_index"],
        "solver_lambda_min_unrepaired": sol["lambda_min"],
        "width_converged": width_proxy(values, radius),
        "coef_l2_from_seed": float(np.linalg.norm(values - coef_seed)),
    }
    np.save(HERE / f"narrow{int(radius)}.npy", values)
    print(f"[{tag}] relgrad={rec['relative_gradient']:.3e} E={rec['E_converged']:.6f} "
          f"width={rec['width_converged']:.3f} morse(unsolved-fn)={rec['solver_morse_index_unrepaired']} "
          f"moved={rec['coef_l2_from_seed']:.3e}", flush=True)
    return values, rec


def main() -> None:
    started = time.time()
    payload: dict = {}

    coef24 = np.asarray(np.load(HERE.parent / "0003" / "roots24.npz")["R24"],
                        dtype=np.float64).reshape(3, ORDER)
    seed26 = embed_and_reproject(coef24, 24.0, 26.0)
    narrow26, rec26 = solve_diagnose(seed26, 26.0, "R26")
    rec26["repaired_stability"] = repaired_hessian_min(narrow26, 26.0)
    print("[R26] repaired:", rec26["repaired_stability"], flush=True)
    payload["R26"] = rec26

    seed28 = embed_and_reproject(narrow26, 26.0, 28.0)
    e_seed28 = width_proxy(seed28, 28.0)
    narrow28, rec28 = solve_diagnose(seed28, 28.0, "R28")
    rec28["width_seed"] = e_seed28
    rec28["repaired_stability"] = repaired_hessian_min(narrow28, 28.0)
    print("[R28] repaired:", rec28["repaired_stability"], flush=True)
    payload["R28"] = rec28

    # family comparison
    roots = np.load(HERE.parent / "0003" / "roots24.npz")
    oracle = solve_radial_1d.Oracle(dict(radial_order=ORDER, radial_nodes=48,
                                         angular_nodes=24, radius=26.0))
    payload["family_R26_energy"] = float(oracle.evaluate(
        np.asarray(roots["R26"], dtype=np.float64).reshape(-1))[0])
    print("family R26 energy:", payload["family_R26_energy"], flush=True)

    payload["minutes_total"] = round((time.time() - started) / 60.0, 2)
    (HERE / "narrow-branch.json").write_text(json.dumps(payload, indent=2, default=float))
    print("WROTE narrow-branch.json")


if __name__ == "__main__":
    main()

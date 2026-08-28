"""P247 attempt 0008: width-mode diagnosis of the spreading branch.

The W3 ladder found lambda_min ~ +5e-10 at R = 24..30: a quasi-zero
mode. Question: is that eigenvector the FAMILY-DRIFT (width) direction?
If yes, the spreading branch slides along a marginally-stiff valley and
the isolation obstruction is the near-flat width stiffness itself.

  M1  full Hessian of the extended fixed-J energy (4x24 = 96 params) at
      the R=24 and R=30 roots; bottom-6 eigenvalues.
  M2  family-drift direction dc = c(R+2) - c(R) (same modal layout);
      overlap |<v_bot, dc_hat>| for the bottom eigenvectors.
  M3  physical profile of the bottom eigenvector (q/tangent/split
      lambda profiles) to see the width character directly.

Cross-check: lambda_min at R=24 should reproduce the W3 value
+4.805453e-10 (same functional, same nodes as the w1w3 run used).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
from numpy.polynomial.chebyshev import Chebyshev

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "0004"))
sys.path.insert(0, str(HERE.parent / "0002"))
sys.path.insert(0, str(HERE.parent / "0001"))

import c3_functional as c3  # noqa: E402

ORDER = 24
NODES = (72, 36)


def total_energy(flat: torch.Tensor, radius: float):
    total, comp = c3.extended_static(
        flat, radial_order=ORDER, radial_nodes=NODES[0],
        angular_nodes=NODES[1], radius=radius, chi_angular_degrees=(0,),
    )
    return total + comp["fixed_j"]


def hessian_at(radius: float, root: np.ndarray) -> np.ndarray:
    flat = torch.zeros((4, ORDER), dtype=torch.float64)
    flat[:3] = torch.tensor(root.reshape(3, ORDER))
    flat.requires_grad_(True)
    total = total_energy(flat, radius)
    grad = torch.autograd.grad(total, flat, create_graph=True)[0].reshape(-1)
    rows = []
    for i in range(grad.numel()):
        rows.append(torch.autograd.grad(grad[i], flat, retain_graph=True)[0].reshape(-1))
    H = torch.stack(rows).detach().numpy()
    return (H + H.T) / 2.0


def lambda_profiles(coef_row: np.ndarray, R: float, r: np.ndarray, mu: float):
    x = np.clip(2.0 * (r / R) ** 2 - 1.0, -1.0, 1.0)
    return Chebyshev(coef_row, domain=[-1.0, 1.0])(x)


def main() -> None:
    started = time.time()
    roots = np.load(HERE.parent / "0003" / "roots24.npz")
    payload = {}

    for Rname, Rnext in [("R24", 26.0), ("R30", None)]:
        R = float(Rname[1:])
        root = np.asarray(roots[Rname], dtype=np.float64)
        H = hessian_at(R, root)
        evals, evecs = np.linalg.eigh(H)
        bottom = evals[:6]
        print(f"[M1 {Rname}] bottom eigenvalues: {bottom}", flush=True)
        entry = {"bottom_eigenvalues": bottom.tolist()}

        if Rnext is not None:
            root_next = np.asarray(roots[f"R{int(Rnext)}"], dtype=np.float64)
            dc = np.concatenate([root_next.reshape(-1), np.zeros(ORDER)])
            dc_hat = dc / np.linalg.norm(dc)
            overlaps = []
            for k in range(6):
                v = evecs[:, k]
                sign = np.sign(np.dot(v, dc_hat)) or 1.0
                overlaps.append(float(np.dot(v * sign, dc_hat)))
            entry["drift_overlaps_bottom6"] = overlaps
            print(f"[M2 {Rname}] |<v_k, dc_hat>| = {np.round(np.abs(overlaps), 4)}", flush=True)

            # physical profile of the bottom eigenvector (mu=1 slice)
            v = evecs[:, 0] * (np.sign(np.dot(evecs[:, 0], dc_hat)) or 1.0)
            v_s = v.reshape(4, ORDER)[:3]
            r = np.linspace(0.3 * R, 0.99 * R, 120)
            prof = {row: lambda_profiles(v_s[row], R, r, 1.0).tolist() for row in range(3)}
            entry["bottom_mode_profiles"] = {
                "r": r.tolist(),
                "rows": prof,
                "max_abs": float(max(np.abs(p).max() for p in prof.values())),
            }
        payload[Rname] = entry

    payload["minutes_total"] = round((time.time() - started) / 60.0, 2)
    (HERE / "width-mode.json").write_text(json.dumps(payload, indent=2, default=float))
    print("WROTE width-mode.json")


if __name__ == "__main__":
    main()

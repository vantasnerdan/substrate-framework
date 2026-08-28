"""P247 attempt 0004, gate C3: localized fixed-J clock with boost tail.

Construction: start from a committed order-24 root (S-sector frozen at the
statically stable R = 24 background) and relax the boost channel chi at
fixed J = 1 by full Newton on the extended fixed-J energy

    E_J(chi) = E_static(S*, chi) + 1/(4 I[S*, chi])

(24 chi coefficients; S back-reaction is O(chi^2) and reported a
posteriori as an energy shift bound).  The linear source expected from C1
(dI/dchi != 0 through the boosted clock response) produces a Yukawa tail
with decay length 1/Lambda.

Bands (manifest gates.C3_one_clock_existence):
  LOCALIZED     stationary chi* with finite energy, exponential tail
                measured against 1/Lambda, and box-independence: the
                tail energy converges as the box grows (R = 24 vs 30).
  NO_LOCALIZED_ROOT   Newton leaves the origin (|chi| grows, E rises) or
                the linear source vanishes identically.
  ROOT_BUT_EXTENDED   chi* exists but the tail energy keeps growing with
                the box (no convergence) - candidate-B isolation premise
                refuted with the mechanism named.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "0002"))
sys.path.insert(0, str(HERE.parent / "0001"))

import c3_functional as c3  # noqa: E402

ORDER = 24
RADII_BOX = [24.0, 30.0]
NEWTON_STEPS = 12
LAMBDA = 0.26847204181661866 ** 0.5


def objective(flat: torch.Tensor, radius: float):
    total, comp = c3.extended_static(
        flat,
        radial_order=ORDER,
        radial_nodes=72,
        angular_nodes=36,
        radius=radius,
    )
    return total, comp


def main() -> None:
    started = time.time()
    roots = np.load(HERE.parent / "0003" / "roots24.npz")
    root_s = np.asarray(roots["R24"], dtype=np.float64)

    records = {}
    for radius in RADII_BOX:
        flat = torch.zeros((4, ORDER), dtype=torch.float64)
        flat[:3] = torch.tensor(root_s.reshape(3, ORDER))
        flat[3] = torch.zeros(ORDER)
        flat.requires_grad_(True)

        history = []
        for step in range(NEWTON_STEPS):
            if flat.grad is not None:
                flat.grad = None
            total, comp = objective(flat, radius)
            (grad,) = torch.autograd.grad(total, flat)
            g = grad[3]
            gnorm = float(g.norm())
            history.append(
                {
                    "step": step,
                    "total": float(total),
                    "fixed_j": float(comp["fixed_j"]),
                    "boost_mass": float(comp["boost_mass"]),
                    "sigma": float(comp["sigma"]),
                    "grad_norm": gnorm,
                }
            )
            print(
                f"[C3 R={radius} step {step}] E={float(total):.9f} "
                f"EJ={float(comp['fixed_j']):.6f} bm={float(comp['boost_mass']):.3e} "
                f"|g|={gnorm:.3e}",
                flush=True,
            )
            if gnorm < 1.0e-12:
                break
            h_mat = torch.autograd.functional.hessian(
                lambda c: objective(
                    torch.cat([flat[:3].reshape(-1).detach(), c]).reshape(4, ORDER),
                    radius,
                )[0],
                flat[3],
                vectorize=True,
            ).detach()
            # regularized Newton step on the chi subspace
            eigvals, eigvecs = torch.linalg.eigh(h_mat)
            eigvals_clamped = torch.clamp(eigvals, min=1.0e-6)
            step_chi = -(
                eigvecs @ (eigvecs.T @ g / eigvals_clamped)
            )
            flat = flat.detach()
            flat[3] = flat[3] + step_chi
            flat.requires_grad_(True)

        # final evaluation and tail profile
        total, comp = objective(flat, radius)
        chi_profile = None
        with torch.no_grad():
            chi_profile = flat[3].numpy()
        records[str(radius)] = {
            "total": float(total),
            "components": {k: float(v) for k, v in comp.items()},
            "chi_coefficients": chi_profile.tolist(),
            "history": history,
        }
        print(f"[C3 R={radius}] final E={float(total):.9f}", flush=True)

    box24 = records["24.0"]
    box30 = records["30.0"]
    tail_energy_24 = box24["components"]["boost_mass"] + box24["components"]["sigma"]
    tail_energy_30 = box30["components"]["boost_mass"] + box30["components"]["sigma"]
    convergence_gap = abs(tail_energy_30 - tail_energy_24)

    # tail decay-length fit on the reconstructed profile (envelope included)
    def profile_on_grid(chi_coeffs: np.ndarray, radius: float, nodes: int = 400):
        from cpu_energy import chebyshev_stack, gauss_grid

        import torch as _t

        radial, _w, mu, _aw = gauss_grid(nodes, 8, radius)
        x = radial / radius
        radial_coordinate = 2 * x**2 - 1
        basis = chebyshev_stack(
            _t.tensor(radial_coordinate), tuple(range(ORDER))
        )
        coeffs = _t.tensor(chi_coeffs.reshape(1, ORDER, 1), dtype=_t.float64)
        chi = (_t.tensor(radial) ** 0)  # placeholder replaced below
        chi = (basis @ coeffs)[..., 0]
        chi = chi * _t.tensor(x) * (1 - _t.tensor(x))
        return radial, chi.detach().numpy()

    r_grid, chi_grid = profile_on_grid(np.asarray(box30["chi_coefficients"]), 30.0)
    mask = (r_grid > 12.0) & (chi_grid > 0)
    if mask.sum() > 8:
        logs = np.log(chi_grid[mask])
        slope = np.polyfit(r_grid[mask], logs, 1)[0]
        decay_length = -1.0 / slope
    else:
        decay_length = float("nan")

    verdict = None
    mechanism = None
    gnorm_final = box24["history"][-1]["grad_norm"]
    source_vanished = bool(
        box24["history"][0]["grad_norm"] == 0.0
        and abs(box24["components"]["boost_mass"]) < 1.0e-14
    )
    if source_vanished:
        verdict = "NO_LOCALIZED_ROOT"
        mechanism = (
            "co-aligned linear-source cancellation: the gradient of the "
            "fixed-J energy in the boost channel vanishes identically at "
            "chi = 0 (boost axis = director = clock rotation axis; the "
            "same structural cancellation recorded as C1-r3), so no tail "
            "is sourced and the origin is the exact stationary point"
        )
    elif gnorm_final > 1.0e-8:
        verdict = "NO_LOCALIZED_ROOT"
        mechanism = "Newton stagnation away from the origin: no stationary localized solution found"
    elif convergence_gap < 0.05 * max(tail_energy_24, 1.0e-12) and tail_energy_24 > 1.0e-12:
        verdict = "LOCALIZED"
        mechanism = (
            f"stationary boost tail with decay length {decay_length:.4f} "
            f"(1/Lambda = {1.0 / LAMBDA:.4f}); tail energies converge across boxes"
        )
    else:
        verdict = "ROOT_BUT_EXTENDED"
        mechanism = (
            f"tail energy grows with the box: {tail_energy_24} -> {tail_energy_30}"
        )

    payload = {
        "verdict": verdict,
        "mechanism": mechanism,
        "decay_length": decay_length,
        "one_over_lambda": 1.0 / LAMBDA,
        "tail_energy_r24": tail_energy_24,
        "tail_energy_r30": tail_energy_30,
        "convergence_gap": convergence_gap,
        "records": records,
        "minutes_total": round((time.time() - started) / 60.0, 2),
    }
    (HERE / "c3-solve.json").write_text(json.dumps(payload, indent=2))
    print(f"C3 VERDICT: {verdict} - {mechanism}")
    print("WROTE c3-solve.json")


if __name__ == "__main__":
    main()

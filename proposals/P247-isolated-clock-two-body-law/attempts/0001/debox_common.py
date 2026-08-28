"""Shared loaders and x-space component Hessians for P247 attempt 0001.

Reuses committed, merged attempt infrastructure as source-audited inputs
(P244 precedent, P247 proposal allowed_imports):

- proposals/P240-m5-kinetic-axis/attempts/0041/cpu_energy.py  (committed functional helpers)
- proposals/P240-m5-kinetic-axis/attempts/0041/solve_radial_1d.py  (Oracle, solve_order, analyze_mode)
- proposals/P240-m5-kinetic-axis/attempts/0042/xspace_energy.py  (fixed-node x-space component evaluator)
- proposals/P240-m5-kinetic-axis/attempts/0042/largeR-roots.json  (committed window-family roots R10/R12/R14)

Numerics follow the small-ratio-numerics prescriptions: float64 end to end,
pinned BLAS/threads, compensated reductions via math.fsum where sums are
load-bearing, and declared gates.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
ATTEMPTS_0041 = HERE / ".." / ".." / ".." / "P240-m5-kinetic-axis" / "attempts" / "0041"
ATTEMPTS_0042 = HERE / ".." / ".." / ".." / "P240-m5-kinetic-axis" / "attempts" / "0042"

for path in (str(ATTEMPTS_0041), str(ATTEMPTS_0042)):
    if path not in sys.path:
        sys.path.insert(0, path)

torch.set_num_threads(1)  # small-ratio-numerics: pin threads for reproducible reductions

ORDER = 16
SOLVE_NODES = (48, 24)  # >= 3x basis order radial; aliasing gate re-evaluates at (96, 48)
CERT_NODES_A = (48, 16)
CERT_NODES_B = (72, 24)

# Frozen reproduction gates (memory contract milestone_gates_attempt_0001).
REPORTED = {
    "slope_per_unit_R": 0.33,
    "slope_band": (0.25, 0.45),
    "r2_min": 0.999,
    "lambda_min_A_flattening_band": (1.0e-7, 2.0e-6),
    "morse_index": 1,
    "split_fraction_min": 0.9,
    "window": (8.0, 34.0),
    "window_relative_halfwidth": 0.20,
    "pinned_R8_energy": 60.17,
    "pinned_R8_energy_rel": 0.02,
    "free_R8_energy": 50.44,
    "free_R8_energy_rel": 0.05,
    "wall_value_band": (0.98, 1.01),
    "inertia_pair": (0.67638, 0.67662),
    "inertia_rel": 1.0e-2,
    "omega_pair": (0.73923, 0.73929),
    "omega_rel": 5.0e-3,
}


def load_committed_roots() -> dict:
    with open(ATTEMPTS_0042 / "largeR-roots.json") as handle:
        return json.load(handle)


def component_hessians(
    flat: np.ndarray, radial_nodes: int, angular_nodes: int, order: int = ORDER
) -> tuple[np.ndarray, np.ndarray]:
    """Exact Hessians (A, D) of the x-space components V and C + Phi.

    With E(R)[c] = R^3 V[c] + (C[c] + Phi[c]) / R and Phi = 1/(4 I)
    (committed xspace_energy.py decomposition), the total Hessian is
    H(R) = R^3 A + R^{-1} D with A = d^2 V, D = d^2 (C + Phi).
    """
    import torch as _torch

    import xspace_energy

    flat_tensor = _torch.tensor(np.asarray(flat, dtype=np.float64), dtype=_torch.float64, requires_grad=True)
    curvature, potential, inertia = xspace_energy.xspace_components(
        flat_tensor, order, radial_nodes, angular_nodes
    )
    phi = 1.0 / (4.0 * inertia)

    def hessian(expression: _torch.Tensor) -> np.ndarray:
        gradient = _torch.autograd.grad(expression, flat_tensor, create_graph=True)[0]
        rows = [
            _torch.autograd.grad(gradient[index], flat_tensor, retain_graph=True)[0]
            for index in range(flat_tensor.numel())
        ]
        return torch.stack(rows).detach().numpy().astype(np.float64)

    a_matrix = hessian(potential)
    d_matrix = hessian(curvature + phi)
    return a_matrix, d_matrix


def lambda_min(matrix: np.ndarray) -> tuple[float, np.ndarray]:
    symmetric = (matrix + matrix.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    return float(eigenvalues[0]), eigenvectors[:, 0]


def committed_morse_gate(eigenvalues: np.ndarray) -> int:
    """Committed index convention (solve_radial_1d.py:231)."""
    return int(np.sum(eigenvalues < -1.0e-8 * max(1.0, float(np.max(np.abs(eigenvalues))))))

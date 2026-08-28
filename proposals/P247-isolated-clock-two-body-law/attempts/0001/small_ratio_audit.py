"""Attempt-0001 addendum: small-ratio-numerics audit of the ladder eigenvalues.

Implements the prescriptions from .agents/skills/small-ratio-numerics/SKILL.md
for every floor-adjacent quantity in the G2/G3 adjudication:

- pre-symmetrization asymmetry ||H - H^T|| recorded before symmetrizing;
- softest eigenpair residual ||H v - lambda v|| against the UNSYMMETRIZED
  operator, relative to |lambda|;
- lambda_min / lambda_2 reported with every verdict;
- in-situ quadrature gauge: lambda_min(A) recomputed at doubled quadrature
  (144 x 48 vs 72 x 24); the shift is a measured discretization-plus-quadrature
  error bar at this exact basis order;
- roundoff-jitter sign test: coefficient perturbations at 1e-13 and 1e-11
  relative scale (fixed seed); the lambda_min spread they induce is compared
  against |lambda_min|;
- an error-budget verdict per rung: a lambda_min is RESOLVED only when it
  exceeds the sum of quadrature gauge, jitter spread, and residual scale by
  the customary order of magnitude.

This script records measurements; it does not re-adjudicate attempt 0001,
whose floor-scoped wording (flatness-to-floor, not positivity) stands.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from debox_common import component_hessians

HERE = Path(__file__).resolve().parent

JITTER_SCALES = (1.0e-13, 1.0e-11)

JITTER_DRAWS = 3


def soft_spectrum(matrix: np.ndarray, count: int = 2) -> tuple[np.ndarray, list[float]]:
    """Eigenvalues, and residuals of the softest pairs vs the raw matrix."""
    symmetric = (matrix + matrix.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    residuals = [
        float(np.linalg.norm(matrix @ eigenvectors[:, i] - eigenvalues[i] * eigenvectors[:, i]))
        for i in range(count)
    ]
    return eigenvalues[:count], residuals


def audit_rung(radius: float, values: np.ndarray) -> dict:
    record: dict = {"radius": radius}

    a_coarse, d_coarse = component_hessians(values, 72, 24)
    a_fine, d_fine = component_hessians(values, 144, 48)

    for label, a_matrix, d_matrix in (("coarse_72x24", a_coarse, d_coarse), ("fine_144x48", a_fine, d_fine)):
        pencil = radius**3 * a_matrix + (1.0 / radius) * d_matrix
        for name, matrix in (("A", a_matrix), ("pencil", pencil)):
            scale = float(np.max(np.abs(matrix)))
            record[f"asym_max_{name}_{label}"] = float(np.max(np.abs(matrix - matrix.T))) / scale
        eigenvalues, residuals = soft_spectrum(a_matrix)
        lam_min, lam_2 = float(eigenvalues[0]), float(eigenvalues[1])
        record[f"lambda_min_A_{label}"] = lam_min
        record[f"lambda_2_A_{label}"] = lam_2
        record[f"ratio_lambda_min_over_lambda_2_{label}"] = lam_min / lam_2
        record[f"residual_softest_A_{label}"] = residuals[0] / max(abs(lam_min), 1e-300)
        pencil_eigenvalues, pencil_residuals = soft_spectrum(pencil)
        record[f"lambda_min_pencil_{label}"] = float(pencil_eigenvalues[0])
        record[f"residual_softest_pencil_{label}"] = (
            pencil_residuals[0] / max(abs(float(pencil_eigenvalues[0])), 1e-300)
        )
        record[f"morse_gate_A_{label}"] = int(np.sum(eigenvalues[0] < -1e-8 * max(1.0, lam_2)))

    record["quadrature_gauge_lambda_min_A"] = abs(
        record["lambda_min_A_fine_144x48"] - record["lambda_min_A_coarse_72x24"]
    )

    # Roundoff-jitter sign test at two perturbation scales, fixed seed.
    rng = np.random.default_rng(20260828)
    norm = float(np.linalg.norm(values))
    for jitter_scale in JITTER_SCALES:
        minima = []
        for _ in range(JITTER_DRAWS):
            perturbed = values + jitter_scale * norm * rng.standard_normal(values.shape) / np.sqrt(
                values.size
            )
            a_jitter, _ = component_hessians(perturbed, 72, 24)
            eigenvalues, _ = soft_spectrum(a_jitter)
            minima.append(float(eigenvalues[0]))
        record[f"jitter_min_A_scale_{jitter_scale:g}"] = float(np.min(minima))
        record[f"jitter_max_A_scale_{jitter_scale:g}"] = float(np.max(minima))
        record[f"jitter_spread_A_scale_{jitter_scale:g}"] = float(np.max(minima) - np.min(minima))

    budget = record["quadrature_gauge_lambda_min_A"] + record["jitter_spread_A_scale_1e-11"]
    softest_residual = max(
        record["residual_softest_A_coarse_72x24"], record["residual_softest_pencil_coarse_72x24"]
    )
    record["error_budget_lambda_min_A"] = budget
    record["resolved_above_budget"] = abs(record["lambda_min_A_coarse_72x24"]) > 10.0 * budget
    record["softest_relative_residual"] = softest_residual
    return record


def main() -> None:
    started = time.time()
    with open(HERE / "clean-ladder.json") as handle:
        rungs = json.load(handle)["rungs"]
    records = []
    for row in rungs:
        if not row.get("accepted"):
            continue
        values = np.asarray(row["values"], dtype=np.float64)
        record = audit_rung(row["radius"], values)
        records.append(record)
        print(
            f"[R={row['radius']:.0f}] lamA={record['lambda_min_A_coarse_72x24']:+.4e} "
            f"lam2={record['lambda_2_A_coarse_72x24']:+.4e} "
            f"ratio={record['ratio_lambda_min_over_lambda_2_coarse_72x24']:+.3e} "
            f"quad_gauge={record['quadrature_gauge_lambda_min_A']:.3e} "
            f"jit_spread={record['jitter_spread_A_scale_1e-11']:.3e} "
            f"res={record['softest_relative_residual']:.2e} "
            f"asym={record['asym_max_A_coarse_72x24']:.1e} "
            f"lam_pencil={record['lambda_min_pencil_coarse_72x24']:+.4e} "
            f"resolved={record['resolved_above_budget']}"
        )
    payload = {
        "records": records,
        "jitter_scales": list(JITTER_SCALES),
        "jitter_draws": JITTER_DRAWS,
        "runtime_s": time.time() - started,
    }
    with open(HERE / "g5-small-ratio-audit.json", "w") as handle:
        json.dump(payload, handle, indent=1)
    print(f"WROTE g5-small-ratio-audit.json runtime={time.time() - started:.1f}s")


if __name__ == "__main__":
    main()

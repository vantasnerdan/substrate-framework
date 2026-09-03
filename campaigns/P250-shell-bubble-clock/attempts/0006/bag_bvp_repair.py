#!/usr/bin/env python
"""Repair and independently replay the P250 spherical bag family.

The moving-coordinate/phase-fixed formulation removes the nearly singular
translation mode that drove attempt 0003 into ``solve_bvp`` node ceilings.
Pressure, charge, physical energy, and the radial virial are evaluated with
their canonical definitions; see ``receipt.md`` in this directory.
"""

from __future__ import annotations

import importlib.util
import contextlib
import io
import json
import math
import os
from pathlib import Path
import sys

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import CubicSpline
from scipy.linalg import sqrtm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "src"))

from substrate_framework.m5_wall_clock import wall_bulk_pressure  # noqa: E402

A3 = HERE.parent / "0003" / "bag_bvp.py"
SPEC = importlib.util.spec_from_file_location("p250_bag_attempt3", A3)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {A3}")
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)

SIGMA0 = json.loads(
    (HERE.parent / "0005" / "sigma0_results_repaired.json").read_text()
)["headline"]["sigma0_best"]


def _load_reference() -> tuple[list[CubicSpline], np.ndarray]:
    path = HERE.parent / "0005" / "profile_L12_anchor.csv"
    data = np.genfromtxt(path, delimiter=",", names=True)
    splines = [CubicSpline(data["x"], data[name]) for name in ("m", "c", "b", "f")]
    b_star = np.array(
        [0.0, float(BASE.BULK_B0[0]), float(BASE.BULK_B0[1]), float(BASE.BULK_B0[2])]
    )
    return splines, b_star


REFERENCE, B_STAR = _load_reference()
A_STATE = np.array([1.0, 0.0, 0.0, 0.0])


def _reference(z: np.ndarray, bulk: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return a centered B-to-A reference and its z derivative."""
    x = np.clip(-np.asarray(z, dtype=float), -12.0, 12.0)
    values = np.stack([s(x) for s in REFERENCE])
    derivatives_x = np.stack([s(x, 1) for s in REFERENCE])
    scale = (bulk - A_STATE) / (B_STAR - A_STATE)
    u = A_STATE[:, None] + (values - A_STATE[:, None]) * scale[:, None]
    return u, -derivatives_x * scale[:, None]


def _dynamical_matrix(state: np.ndarray, omega_sq: float) -> np.ndarray:
    hessian = np.array(
        [[term(*state, omega_sq) for term in row] for row in BASE._H_terms],
        dtype=float,
    )
    return BASE.K2inv @ hessian


def _positive_matrix_sqrt(matrix: np.ndarray) -> np.ndarray:
    root = sqrtm(matrix)
    if np.max(np.abs(np.imag(root))) > 1e-10:
        raise RuntimeError("asymptotic matrix square root is not real")
    root = np.asarray(np.real(root), dtype=float)
    if np.linalg.norm(root @ root - matrix, ord=np.inf) > 2e-10:
        raise RuntimeError("asymptotic matrix square root failed its residual")
    return root


def solve_phase_fixed(
    omega_sq: float,
    bulk: np.ndarray,
    radius_guess: float,
    *,
    match_exponent: float = 14.0,
    tol: float = 1e-8,
    n_initial: int = 1201,
) -> tuple[object, float, float]:
    """Solve the radial wall in z=r-R with R as a phase-fixed parameter."""
    half, exterior, _, _ = BASE.window_extents(
        omega_sq, bulk, target=match_exponent
    )
    z = np.linspace(-half, exterior, n_initial)
    reference, reference_prime = _reference(z, bulk)
    y0 = np.vstack((reference, reference_prime, np.zeros_like(z)))
    sqrt_b = _positive_matrix_sqrt(_dynamical_matrix(bulk, omega_sq))
    sqrt_a = _positive_matrix_sqrt(_dynamical_matrix(A_STATE, omega_sq))

    def fun(z_value, y, parameters):
        radius = parameters[0]
        r = radius + z_value
        kernel = BASE._kernel(y[0], y[1], y[2], y[3], omega_sq)
        acceleration = BASE.K2inv @ np.stack(kernel[1:5]) - 2.0 * y[4:8] / r
        ref, ref_prime = _reference(z_value, bulk)
        phase_derivative = np.sum((y[:4] - ref) * ref_prime, axis=0)
        return np.vstack((y[4:8], acceleration, phase_derivative))

    def boundary(ya, yb, parameters):
        radius = parameters[0]
        inner_log_derivative = sqrt_b - np.eye(4) / (radius - half)
        outer_log_derivative = sqrt_a + np.eye(4) / (radius + exterior)
        return np.concatenate(
            (
                ya[4:8] - inner_log_derivative @ (ya[:4] - bulk),
                yb[4:8] + outer_log_derivative @ (yb[:4] - A_STATE),
                [ya[8], yb[8]],
            )
        )

    solution = solve_bvp(
        fun,
        boundary,
        z,
        y0,
        p=[radius_guess],
        tol=tol,
        max_nodes=200_000,
        verbose=0,
    )
    boundary_residual = float(
        np.max(np.abs(boundary(solution.y[:, 0], solution.y[:, -1], solution.p)))
    )
    return solution, half, exterior, boundary_residual


def _compensated_trapezoid(x: np.ndarray, y: np.ndarray) -> float:
    return math.fsum(
        (x[i + 1] - x[i]) * (y[i] + y[i + 1]) / 2
        for i in range(len(x) - 1)
    )


def analyze(solution, half: float, exterior: float, boundary_residual: float,
            omega_sq: float, bulk: np.ndarray, pressure: float) -> dict:
    radius_phase = float(solution.p[0])
    z = np.linspace(-half, exterior, 30_001)
    y = solution.sol(z)
    r = radius_phase + z
    kernel = BASE._kernel(y[0], y[1], y[2], y[3], omega_sq)
    potential = np.asarray(kernel[0], dtype=float)
    gradient = (y[4] ** 2 + 2 * y[5] ** 2 + 2 * y[6] ** 2) / 4 + y[7] ** 2 / 2
    inertia_density = y[3] ** 2 + 4 * y[2] ** 2

    roots = [
        float(value)
        for value in CubicSpline(z, y[0]).solve(0.5)
        if -half < float(value) < exterior
    ]
    if len(roots) != 1:
        raise RuntimeError(f"expected one m=1/2 crossing, found {roots}")
    radius_m = radius_phase + roots[0]
    radius_t = float(r[np.argmax(gradient)])

    inner_radius = radius_phase - half
    iota_bulk = float(bulk[3] ** 2 + 4 * bulk[2] ** 2)
    eomega_inner = -4 * math.pi * pressure * inner_radius ** 3 / 3
    inertia_inner = 4 * math.pi * iota_bulk * inner_radius ** 3 / 3
    virial_inner = -pressure * inner_radius ** 3
    virial_norm_inner = pressure * inner_radius ** 3

    eomega_window = 4 * math.pi * _compensated_trapezoid(
        r, r ** 2 * (gradient + potential)
    )
    inertia_window = 4 * math.pi * _compensated_trapezoid(
        r, r ** 2 * inertia_density
    )
    virial_window = _compensated_trapezoid(
        r, r ** 2 * (gradient + 3 * potential)
    )
    virial_norm_window = _compensated_trapezoid(
        r, r ** 2 * (gradient + 3 * np.abs(potential))
    )

    eomega = eomega_inner + eomega_window
    inertia = inertia_inner + inertia_window
    omega = math.sqrt(omega_sq)
    charge = omega * inertia
    energy = eomega + omega * charge
    virial = virial_inner + virial_window
    virial_norm = virial_norm_inner + virial_norm_window

    return {
        "R_phase": radius_phase,
        "R_m": radius_m,
        "R_T": radius_t,
        "radius_estimator_spread": abs(radius_m - radius_t),
        "p_exact": pressure,
        "iota_B": iota_bulk,
        "chi": radius_m * pressure / (2 * SIGMA0),
        "E_omega": eomega,
        "I": inertia,
        "Q": charge,
        "E": energy,
        "virial": virial,
        "virial_relative": abs(virial) / max(virial_norm, 1e-300),
        "solver_status": int(solution.status),
        "solver_message": solution.message,
        "nodes": int(solution.x.size),
        "residual_rms": float(np.sqrt(np.mean(solution.rms_residuals ** 2))),
        "residual_max": float(np.max(solution.rms_residuals)),
        "boundary_residual": boundary_residual,
        "inner_endpoint_deviation": float(np.max(np.abs(y[:4, 0] - bulk))),
        "outer_endpoint_deviation": float(np.max(np.abs(y[:4, -1] - A_STATE))),
    }


def solve_row(delta: float, *, match_exponent=14.0, tol=1e-8,
              n_initial=1201) -> tuple[dict, object, float, float]:
    omega_sq_mp = BASE.W2_STAR + BASE.mp.mpf(str(delta))
    omega_sq = float(omega_sq_mp)
    branch, branch_residual, branch_det = BASE.branch_solve(
        omega_sq_mp, BASE.BULK_B0
    )
    bulk = np.array([0.0, *branch], dtype=float)
    pressure_expr = wall_bulk_pressure(*bulk, omega_sq)
    pressure = float(pressure_expr)
    if not pressure > 0:
        raise RuntimeError(f"nonpositive pressure at delta={delta}: {pressure}")
    radius_guess = 2 * SIGMA0 / pressure
    solution, half, exterior, boundary_residual = solve_phase_fixed(
        omega_sq,
        bulk,
        radius_guess,
        match_exponent=match_exponent,
        tol=tol,
        n_initial=n_initial,
    )
    row = analyze(
        solution, half, exterior, boundary_residual, omega_sq, bulk, pressure
    )
    row.update(
        {
            "delta": delta,
            "omega_sq": omega_sq,
            "omega": math.sqrt(omega_sq),
            "bulk": bulk.tolist(),
            "branch_residual": float(branch_residual),
            "branch_jacobian_determinant": float(branch_det),
            "match_exponent": match_exponent,
            "tol": tol,
        }
    )
    return row, solution, half, exterior


def main() -> int:
    rows = []
    solutions = []
    print("== P250 attempt 0006: phase-fixed Robin bag replay ==")
    print(f"sigma_0 = {SIGMA0:.15f}")
    for index in range(1, 8):
        row, solution, half, exterior = solve_row(index / 1000)
        rows.append(row)
        solutions.append((solution, half, exterior))
        print(
            f"delta={row['delta']:.3f} R={row['R_m']:.9f} "
            f"p={row['p_exact']:.12g} chi={row['chi']:.9f} "
            f"Q={row['Q']:.9e} vir={row['virial_relative']:.3e} "
            f"res={row['residual_max']:.3e} status={row['solver_status']} "
            f"nodes={row['nodes']}"
        )

    sensitivities = {"tolerance": [], "matching": []}
    delta_mid = 0.004
    baseline = rows[3]
    for tolerance in (1e-6, 1e-8, 1e-10):
        row, _, _, _ = solve_row(delta_mid, tol=tolerance)
        sensitivities["tolerance"].append(
            {
                "tol": tolerance,
                "R_m": row["R_m"],
                "E": row["E"],
                "Q": row["Q"],
                "status": row["solver_status"],
                "residual_max": row["residual_max"],
            }
        )
    for exponent in (12.0, 14.0, 16.0):
        row, _, _, _ = solve_row(delta_mid, match_exponent=exponent)
        sensitivities["matching"].append(
            {
                "match_exponent": exponent,
                "R_m": row["R_m"],
                "E": row["E"],
                "Q": row["Q"],
                "status": row["solver_status"],
                "residual_max": row["residual_max"],
            }
        )

    envelope = []
    for i in range(1, len(rows) - 1):
        derivative = (rows[i + 1]["E"] - rows[i - 1]["E"]) / (
            rows[i + 1]["Q"] - rows[i - 1]["Q"]
        )
        relative = abs(derivative - rows[i]["omega"]) / rows[i]["omega"]
        envelope.append(
            {
                "delta": rows[i]["delta"],
                "dE_dQ": derivative,
                "omega": rows[i]["omega"],
                "relative_error": relative,
            }
        )

    pressure_regression = []
    for row in rows:
        endpoint_product = row["delta"] * row["iota_B"] / 2
        pressure_regression.append(
            {
                "delta": row["delta"],
                "p_exact": row["p_exact"],
                "p_endpoint_product": endpoint_product,
                "relative_endpoint_error": (endpoint_product - row["p_exact"])
                / row["p_exact"],
            }
        )

    checks = {
        "all_solver_status_zero": all(row["solver_status"] == 0 for row in rows),
        "all_residual_max_below_1p1e_8": all(
            row["residual_max"] < 1.1e-8 for row in rows
        ),
        "all_boundary_residual_below_1e_10": all(
            row["boundary_residual"] < 1e-10 for row in rows
        ),
        "radius_strictly_decreases": all(
            rows[i]["R_m"] > rows[i + 1]["R_m"] for i in range(len(rows) - 1)
        ),
        "charge_strictly_decreases": all(
            rows[i]["Q"] > rows[i + 1]["Q"] for i in range(len(rows) - 1)
        ),
        "chi_tends_to_one": abs(rows[0]["chi"] - 1)
        < abs(rows[-1]["chi"] - 1)
        and abs(rows[0]["chi"] - 1) < 0.002,
        "virial_relative_below_2e_5": all(
            row["virial_relative"] < 2e-5 for row in rows
        ),
        "envelope_relative_below_2e_4": max(
            item["relative_error"] for item in envelope
        )
        < 2e-4,
        "sensitivity_status_zero": all(
            item["status"] == 0
            for family in sensitivities.values()
            for item in family
        ),
        "exact_pressure_differs_from_endpoint_approximation": all(
            item["relative_endpoint_error"] > 0
            for item in pressure_regression
        ),
        "rung1_critical_energy_relative_below_1p6e_3": abs(
            16 * math.pi * SIGMA0 ** 3 / (3 * rows[0]["p_exact"] ** 2)
            - rows[0]["E_omega"]
        )
        / rows[0]["E_omega"]
        < 1.6e-3,
    }

    output = {
        "method": "moving-coordinate phase-fixed radial BVP with asymptotic Robin data",
        "sigma0": SIGMA0,
        "omega_sq_star": float(BASE.W2_STAR),
        "rungs": rows,
        "envelope": envelope,
        "pressure_regression": pressure_regression,
        "sensitivities": sensitivities,
        "checks": checks,
    }
    (HERE / "bag_results_repaired.json").write_text(json.dumps(output, indent=2))

    for tag, index in (("thin", 0), ("thick", len(rows) - 1)):
        solution, half, exterior = solutions[index]
        z = np.linspace(-half, exterior, 2001)
        y = solution.sol(z)
        r = solution.p[0] + z
        data = np.column_stack((r, y[:8].T))
        np.savetxt(
            HERE / f"profile_{tag}.csv",
            data,
            delimiter=",",
            header="r,m,c,b,f,dm,dc,db,df",
            comments="",
        )

    for name, passed in checks.items():
        print(f"[{'PASS' if passed else 'FAIL'}] {name}")
    print(f"max envelope relative error = {max(x['relative_error'] for x in envelope):.3e}")
    print("wrote bag_results_repaired.json and endpoint profiles")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    transcript = io.StringIO()
    with contextlib.redirect_stdout(transcript):
        status = main()
    rendered = transcript.getvalue()
    print(rendered, end="")
    (HERE / "run_stdout.txt").write_text(rendered)
    raise SystemExit(status)

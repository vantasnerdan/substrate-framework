"""Run only when authorized: scalar FULL-functional differentiation, attempt0013.

PYTHONPATH=src .venv/bin/python proposals/P254-ni-collective-response/attempts/0013/static_pilot.py
Inputs and steps are frozen in fixture.py and registration.json; no fitted terms.
"""
from __future__ import annotations

import json
import sys

import numpy as np

from fixture import make_fixture
from substrate_framework.variational_fermion_static import analytic_hessian, solve_minimum

STEPS = (.02, .01, .005)
MUTATIONS = ("wrongGC", "missingW", "droppedUcross")


def _json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def run_ensemble(fixture, fixed_m):
    """Differentiate actual optimized energies; analytic prediction enters afterward."""
    states = {}

    def energy(u):
        key = tuple(float(x) for x in u)
        if key not in states:
            states[key] = solve_minimum(fixture, key, fixed_m)
        return states[key]["energy"]

    center = energy((0., 0.))
    differences = []
    for h in STEPS:
        hessian = np.zeros((2, 2))
        for a in range(2):
            offset = np.zeros(2)
            offset[a] = h
            hessian[a, a] = (energy(offset) - 2 * center + energy(-offset)) / h**2
        hessian[0, 1] = hessian[1, 0] = (
            energy((h, h)) - energy((h, -h))
            - energy((-h, h)) + energy((-h, -h))) / (4 * h**2)
        differences.append({"h": h, "hessian": hessian})

    expected = analytic_hessian(fixture, fixed_m)
    previous_error = None
    for index, row in enumerate(differences):
        error = row["hessian"] - expected
        norm = float(np.max(np.abs(error)))
        allowance = max(2e-5, 2 * row["h"]**2)
        row.update({"error_matrix": error, "error_max": norm,
                    "error_frobenius": float(np.linalg.norm(error)),
                    "registered_allowance": allowance, "agrees": bool(norm <= allowance)})
        if previous_error is not None:
            ratio = previous_error / norm if norm > 0 else None
            row["previous_error_over_error"] = ratio
            row["observed_order"] = float(np.log2(ratio)) if ratio and ratio > 0 else None
            row["hessian_change_from_previous"] = (
                row["hessian"] - differences[index - 1]["hessian"])
        previous_error = norm

    finest = differences[-1]
    sensitivities = {}
    for mutation in MUTATIONS:
        altered = analytic_hessian(fixture, fixed_m, mutation)
        error = finest["hessian"] - altered
        norm = float(np.max(np.abs(error)))
        sensitivities[mutation] = {
            "hessian": altered, "difference_from_unmutated": altered - expected,
            "error_against_independent_scalar_hessian": error,
            "error_max": norm,
            "detected_at_registered_allowance": bool(norm > finest["registered_allowance"]),
        }

    records = []
    for state in states.values():
        records.append({key: value for key, value in state.items() if key != "gamma"})
    minimum_pass = all(s["accepted_minimum"] and s["solver_success"] for s in records)
    curvature_pass = all(row["agrees"] for row in differences)
    mutation_pass = all(s["detected_at_registered_allowance"] for s in sensitivities.values())
    return {
        "ensemble": "fixedN-fixedM" if fixed_m else "fixedN-relaxedM",
        "constraints": {"N": fixture["N0"], "M": fixture["M0"] if fixed_m else None},
        "reference_minimum": records[0], "scalar_minima": records,
        "analytic_hessian": expected, "finite_difference": differences,
        "sensitivities": sensitivities,
        "all_scalar_minima_converged": minimum_pass,
        "all_curvatures_within_registered_allowance": curvature_pass,
        "all_registered_mutations_detected": mutation_pass,
        "passed": bool(minimum_pass and curvature_pass and mutation_pass),
    }


def main():
    fixture = make_fixture()
    ensembles = [run_ensemble(fixture, fixed_m) for fixed_m in (False, True)]
    result = {
        "attempt": "0013", "oracle": "independent_minimized_full_scalar_functional",
        "scope": {
            "model": "Four total one-body orbitals, two per collinear spin; arbitrary dimensionless units",
            "admissible_space": "All Hermitian within-spin 2x2 density blocks with 0<gamma<I",
            "thermodynamics": "Isothermal FD mean constraints; NOT exact finite canonical Gibbs",
            "interaction": "Frozen positive density Hessian U including cross-spin coupling",
            "minimum_justification": "Strictly convex spectral entropy plus PSD density quadratic; feasible stationarity",
            "independence": "FD selfconsistency and full scalar energy do not call the analytic response/Schur matrix",
            "analytic_comparison": "Gibbs divided differences, constraint Schur projection, induced-density U loop, plus <W>",
            "does_not_license": ["Actual Ni or H calibration", "SG15/PBE functional agreement",
                                 "Nonzero-q native transfer", "Physical loss or nuclear/product-energy rates",
                                 "Isolated dynamics from isothermal population response"],
            "mutation_scope": "Only analytic comparison changed; independent scalar functional remains frozen",
            "wrongGC_scope": "Drop both number and (where present) magnetization constraints",
        },
        "solver": {"routine": "scipy.optimize.root", "method": "hybr", "xtol": 1e-11,
                   "precision": "float64", "density_constraint_residual_limit": 1e-10,
                   "stationarity_frobenius_limit": 1e-9},
        "fixture": fixture, "ensembles": ensembles,
        "ensemble_hessian_difference_fixedM_minus_relaxedM": (
            ensembles[1]["analytic_hessian"] - ensembles[0]["analytic_hessian"]),
        "passed": all(row["passed"] for row in ensembles),
    }
    print(json.dumps(result, default=_json_default, indent=2, allow_nan=False))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

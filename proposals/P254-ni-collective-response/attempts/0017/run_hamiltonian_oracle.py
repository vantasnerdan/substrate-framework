#!/usr/bin/env python3
"""Corroborate full target-PW H against every retained Native target eigenpair.

Run with PYTHONPATH=src; --output must name a new JSON receipt. Its companion
NPZ stores complete per-state residual/error vectors, not dense H copies.
Agreement corroborates retained Native eigenpairs, not independent full-space
Native H-action equivalence, material convergence or any coherence/loss claim.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback
from typing import Any

import numpy as np

import substrate_framework.qe_pw_hamiltonian as pw
from substrate_framework.qe_response_bundle import load_bundle
from substrate_framework.qe_response_fft import NativeRhoFFT


ATTEMPT = Path(__file__).resolve().parent
RESIDUAL_LIMIT_RY = 1e-4
PROJECTION_LIMIT_RY = 1e-4
NORMALIZATION_LIMIT = 1e-8
HERMITICITY_FACTOR = 1e-10
MUTATIONS = ("omitted_local", "omitted_nonlocal", "wrong_fft_normalization")


def digest(path: Path, data: bytes | None = None) -> dict[str, str]:
    raw = path.read_bytes() if data is None else data
    return {"path": str(path), "bytes": str(len(raw)),
            "sha256": hashlib.sha256(raw).hexdigest()}


def source_receipts() -> tuple[dict[str, Any], dict[str, Any]]:
    """Actually read every preregistered source and bind the running consumers."""
    inventory_path = ATTEMPT / "source-inventory.json"
    registration_path = ATTEMPT / "registration.json"
    inventory = json.loads(inventory_path.read_bytes())
    registration = json.loads(registration_path.read_bytes())
    receipts = []
    for record in inventory["raw_byte_inventory"] + inventory["consumer_sources"]:
        actual = digest(Path(record["path"]))
        if any(actual[key] != record[key] for key in ("bytes", "sha256")):
            raise ValueError(f"preregistered source changed: {record['path']}")
        receipts.append(actual)
    return inventory, {
        "preregistration_commit": "c58135d1",
        "registration": digest(registration_path),
        "source_inventory": digest(inventory_path),
        "registered_sources": receipts,
        "new_implementation": [digest(Path(pw.__file__)), digest(Path(__file__).resolve())],
        "frozen_acceptance": registration["acceptance"],
    }


def error_diagnostics(
    action: np.ndarray, psi: np.ndarray, eigenvalues: np.ndarray,
    prefix: str, payload: dict[str, np.ndarray],
) -> dict[str, Any]:
    """Keep raw errors and report exactly the registered eigenpair comparisons."""
    residual = action - psi * eigenvalues[None, :]
    projection_error = psi.conj().T @ action - np.diag(eigenvalues)
    residual_norms = np.linalg.norm(residual, axis=0)
    projection_maxima = np.max(np.abs(projection_error), axis=0)
    if not (np.isfinite(residual).all() and np.isfinite(projection_error).all()):
        raise ValueError(f"{prefix}: nonfinite action/error; no repair is permitted")
    payload.update({prefix + "residual_vectors_ry": residual,
                    prefix + "residual_2norm_ry": residual_norms,
                    prefix + "projection_energy_error_ry": projection_error,
                    prefix + "projection_maxabs_per_band_ry": projection_maxima})
    return {
        "residual_2norm_per_band_ry": residual_norms.tolist(),
        "max_band_residual_2norm_ry": float(np.max(residual_norms)),
        "projection_energy_maxabs_per_band_ry": projection_maxima.tolist(),
        "max_projection_energy_error_ry": float(np.max(projection_maxima)),
        "rayleigh_energy_error_per_band_ry": projection_error.diagonal().real.tolist(),
        "residual_pass": bool(np.all(residual_norms <= RESIDUAL_LIMIT_RY)),
        "projection_pass": bool(np.all(projection_maxima <= PROJECTION_LIMIT_RY)),
    }


def evaluate_state(bundle: Any, state: Any, label: str,
                   payload: dict[str, np.ndarray]) -> dict[str, Any]:
    meta = state.metadata
    n, nbnd = (meta[key] for key in ("ngw_target", "nbnd"))
    prefix = f"{label}_spin{meta['spin']}_ik{meta['ik']}_"
    psi = state.array("evq", (meta["npwx"], nbnd), "complex128")[:n]
    eigenvalues = state.array("eigenvalues_target", (nbnd,), "float64")
    h = pw.assemble_target_hamiltonian(bundle, state)
    if h.shape != (n, n) or h.dtype != np.dtype("complex128"):
        raise ValueError("assembler did not return the full complex target-PW square")
    hermiticity_rows = np.max(np.abs(h - h.conj().T), axis=1)
    hscale = float(np.max(np.abs(h)))
    hermiticity_limit = HERMITICITY_FACTOR * max(1., hscale)
    gram_error = psi.conj().T @ psi - np.eye(nbnd)
    normalization_error = float(np.max(np.abs(gram_error)))
    baseline = error_diagnostics(h @ psi, psi, eigenvalues, prefix, payload)
    residual = payload[prefix + "residual_vectors_ry"]
    retained_residual = psi.conj().T @ residual
    payload.update({prefix + "gram_error": gram_error,
                    prefix + "native_eigenvalues_target_ry": eigenvalues,
                    prefix + "igk_target_one_based": state.array("igk_target"),
                    prefix + "hermiticity_row_maxabs_ry": hermiticity_rows,
                    prefix + "retained_projected_residual_ry": retained_residual,
                    prefix + "complement_residual_vectors_ry": residual - psi @ retained_residual})
    mutations = {}
    for name, action in pw.hamiltonian_sensitivity_actions(bundle, state, psi).items():
        errors = error_diagnostics(action, psi, eigenvalues, prefix + name + "_", payload)
        errors["exposed"] = not (errors["residual_pass"] and errors["projection_pass"])
        mutations[name] = errors
    hermitian = bool(np.all(hermiticity_rows <= hermiticity_limit))
    normalized = normalization_error <= NORMALIZATION_LIMIT
    return {
        "document": state.path, "npz_prefix": prefix,
        # Native identities/dimensions remain exact decimal strings in the receipt.
        "native_dimensions_and_identities": {key: str(meta[key]) for key in
            ("spin", "ik", "ikk", "ikq", "ngw_source", "ngw_target", "npwx", "nbnd", "nkb")},
        "native_projector_dimensions": {
            "nh": [str(int(value)) for value in state.array("nh")],
            "ityp_one_based": [str(int(value)) for value in state.array("ityp")],
            "deeq_shape": [str(value) for value in state.array("deeq").shape],
            "vkb_target_shape": [str(value) for value in state.array("vkb_target").shape]},
        "hamiltonian_shape": [str(n), str(n)],
        "unretained_complement_dimension": str(n - nbnd),
        "basis": "exact igk_target order; leading physical evq/vkb_target rows; no PW truncation",
        "dtype": str(h.dtype), "maxabs_h_ry": hscale,
        "hermiticity_maxabs_ry": float(np.max(hermiticity_rows)),
        "hermiticity_allowance_ry": hermiticity_limit, "hermiticity_pass": hermitian,
        "normalization_maxabs": normalization_error, "normalization_pass": normalized,
        "baseline": baseline, "mutations": mutations,
        "agreement": bool(hermitian and normalized and baseline["residual_pass"]
                          and baseline["projection_pass"]),
    }


def evaluate_bundle(label: str, path: Path, inventory: dict[str, Any],
                    payload: dict[str, np.ndarray]) -> dict[str, Any]:
    bundle = load_bundle(path)
    manifest = digest(bundle.root / "manifest.json", bundle.buffers["manifest.json"])
    expected = {record["sha256"]: record for record in inventory["settled_native_inputs"]}
    if (manifest["sha256"] not in expected or
            manifest["bytes"] != expected[manifest["sha256"]]["bytes"]):
        raise ValueError("input is not an immutable preregistered generation-003 bundle")
    expected_branch = "gamma" if label == "gamma" else "half-b1"
    if expected_branch not in Path(expected[manifest["sha256"]]["path"]).parts:
        raise ValueError("Gamma/half inputs exchanged relative to preregistration")
    bg = bundle.background
    fft = NativeRhoFFT.from_bundle(bundle)
    q = bg.array("q_cart_2pi_alat", (3,), "float64")
    reciprocal = bg.array("bg", (3, 3), "float64")
    q_error = (float(np.max(np.abs(q))) if label == "gamma" else
               float(min(np.max(np.abs(q - .5 * reciprocal[:, 0])),
                         np.max(np.abs(q + .5 * reciprocal[:, 0])))))
    if q_error > 1e-13:
        raise ValueError(f"{label}: Native q is not the requested Gamma/half-b1")
    nksq = bg.metadata["nksq"]
    if len(bundle.states) != nksq or {s.metadata["ik"] for s in bundle.states} != set(range(1, nksq + 1)):
        raise ValueError("closed bundle does not contain every Native electronic state")
    results = []
    for state in bundle.states:
        try:
            results.append(evaluate_state(bundle, state, label, payload))
        except Exception as error:
            # Preserve the first failure and continue collecting all other states.
            results.append({"document": state.path, "agreement": False,
                            "native_dimensions_and_identities": {
                                key: str(state.metadata[key]) for key in
                                ("spin", "ik", "ikk", "ikq", "ngw_source", "ngw_target",
                                 "npwx", "nbnd", "nkb")},
                            "error_type": type(error).__name__, "error": str(error),
                            "traceback": traceback.format_exc()})
    return {
        "bundle": str(bundle.root), "manifest": manifest,
        "native_input_closure": json.loads(bundle.buffers["manifest.json"])["files"],
        "native_nksq": str(nksq), "reported_states": str(len(results)),
        "physical_fft": [str(n) for n in fft.physical],
        "storage_fft": [str(n) for n in fft.storage],
        "native_nnr": str(fft.nnr), "q_cart_2pi_alat": q.tolist(),
        "states": results, "agreement": all(state["agreement"] for state in results),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gamma", type=Path, required=True, help="closed Native Gamma q directory")
    parser.add_argument("--half", type=Path, required=True, help="closed Native half-b1 q directory")
    parser.add_argument("--output", type=Path, required=True, help="new JSON; companion .npz keeps raw errors")
    args = parser.parse_args()
    archive = args.output.with_suffix(".npz")
    if args.output.suffix.lower() != ".json" or args.output.exists() or archive.exists():
        raise ValueError("preserve receipts: choose a new .json path and absent companion .npz")
    started = time.monotonic()
    payload: dict[str, np.ndarray] = {}
    report: dict[str, Any] = {
        "schema_version": 1, "started_at": datetime.now(timezone.utc).isoformat(),
        "object": "full ngw_target-square NC target Hamiltonian in Ry",
        "source_equation": "diag(g2kin_target) + FFT_forward(vrs_spin)[G_i-G_j] + B D Bdagger",
        "allowances": {"max_band_residual_2norm_ry": RESIDUAL_LIMIT_RY,
                       "max_projection_energy_error_ry": PROJECTION_LIMIT_RY,
                       "maxabs_gram_error": NORMALIZATION_LIMIT,
                       "hermiticity": "1e-10*max(1,maxabs(H)) Ry"},
        "environment": {"command": [sys.executable, *sys.argv], "python": platform.python_version(),
                        "numpy": np.__version__, "cpu_affinity": sorted(os.sched_getaffinity(0)),
                        "threads": {key: os.environ.get(key) for key in
                            ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")}},
        "results": {},
        "independent_full_space_native_h_action_proof": False,
        "scope": "Native retained-eigenpair corroboration only. The complete PW complement is constructed, but independent Native H actions on its unretained vectors are not exported. No material/coherence/loss/nuclear claim.",
        "missing_native_proof_field": "independent full-PW Native H-action probes on the unretained complement",
    }
    try:
        inventory, report["sources"] = source_receipts()
        for label, path in (("gamma", args.gamma), ("half", args.half)):
            try:
                report["results"][label] = evaluate_bundle(label, path, inventory, payload)
            except Exception as error:
                report["results"][label] = {"agreement": False, "error_type": type(error).__name__,
                                            "error": str(error), "traceback": traceback.format_exc()}
    except Exception as error:
        report["source_error"] = {"error_type": type(error).__name__, "error": str(error),
                                  "traceback": traceback.format_exc()}
    all_states = [state for result in report["results"].values() for state in result.get("states", [])]
    exposure = {name: [state["npz_prefix"] for state in all_states
                       if state.get("mutations", {}).get(name, {}).get("exposed", False)]
                for name in MUTATIONS}
    report["mutation_exposure"] = {
        name: {"exposed_state_prefixes": prefixes, "exposed": bool(prefixes),
               "interpretation": "exposed" if prefixes else "not exposed; no sensitivity credit"}
        for name, prefixes in exposure.items()}
    report["sensitivity_pass"] = any(exposure.values())
    report["retained_eigenpair_agreement"] = bool(len(report["results"]) == 2 and
        all(result["agreement"] for result in report["results"].values()))
    passed = report["retained_eigenpair_agreement"] and report["sensitivity_pass"]
    report["acceptance_pass"] = passed
    report["elapsed_seconds"] = time.monotonic() - started
    with archive.open("xb") as stream:
        np.savez_compressed(stream, **payload)
    report["diagnostics_archive"] = digest(archive)
    text = json.dumps(report, indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as stream:
        stream.write(text)
    print(text, end="")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Compare factorized arbitrary-PW bare ionic actions with new closed Native data.

PYTHONPATH=src; --output must name a new JSON receipt (with a new companion NPZ).
All three inputs must be distinct actual Gamma/+half-b1/-half-b1 bundles carrying
the 0018 extension. Old generation-003 data are not a substitute. No Native job,
fit, hermitization, precision change or independent complement oracle is run.
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

import substrate_framework.qe_ionic_operator as ionic
import substrate_framework.qe_response_bundle as decoder
import substrate_framework.qe_response_fft as fft_module
import substrate_framework.qe_response_schema as schema


ATTEMPT = Path(__file__).resolve().parent
FACTOR = 1e-9
MUTATIONS = ("dropped_local", "dropped_nonlocal", "wrong_target_derivative_sign",
             "conjugated_projector_phases", "conjugated_displacement")


def binding(path: Path, data: bytes | None = None) -> dict[str, str]:
    raw = path.read_bytes() if data is None else data
    return {"path": str(path), "bytes": str(len(raw)),
            "sha256": hashlib.sha256(raw).hexdigest()}


def failure(error: Exception) -> dict[str, Any]:
    return {"agreement": False, "error_type": type(error).__name__, "error": str(error),
            "traceback": traceback.format_exc()}


def source_bindings() -> dict[str, Any]:
    inventory_path = ATTEMPT / "source-inventory.json"
    inventory_bytes = inventory_path.read_bytes()
    inventory = json.loads(inventory_bytes)
    registered = []
    for record in inventory["sources"]:
        actual = binding(Path(record["path"]))
        unchanged = all(actual[key] == record[key] for key in ("bytes", "sha256"))
        # Native exporters are deliberately extended after registration; QE is immutable.
        if str(record["path"]).startswith("/opt/qebuild/") and not unchanged:
            raise ValueError(f"immutable QE authority changed: {record['path']}")
        registered.append({"actual": actual, "registered": record,
                           "unchanged_since_registration": unchanged})
    native_dir = Path("/home/administrator/projects/catalyst-nc-operator/native/nc-response")
    return {
        "preregistration_commit": "b3f358db",
        "registration": binding(ATTEMPT / "registration.json"),
        "operator_contract": binding(ATTEMPT / "operator-contract.json"),
        "source_inventory": binding(inventory_path, inventory_bytes),
        "registered_source_bindings": registered,
        "actual_native_exporter_sources": [binding(path) for path in sorted(native_dir.glob("*.f90"))],
        "actual_consumer_sources": [binding(Path(module.__file__).resolve()) for module in
                                    (ionic, decoder, fft_module, schema)]
                                   + [binding(Path(__file__).resolve())],
    }


def comparison(action: np.ndarray, native: np.ndarray,
               allowance: float) -> tuple[np.ndarray, dict[str, Any]]:
    if action.shape != native.shape or not np.isfinite(action).all():
        raise ValueError("action must retain Native physical target rows and be finite")
    error = action - native
    if not np.isfinite(error).all():
        raise ValueError("nonfinite error; no repair is permitted")
    per_band = np.max(np.abs(error), axis=0)
    maximum = float(np.max(per_band))
    return error, {"maxabs_error_ry_per_bohr": maximum,
                   "maxabs_error_per_band_ry_per_bohr": per_band.tolist(),
                   "allowance_ry_per_bohr": allowance, "pass": maximum <= allowance}


def pattern_result(operator: ionic.BareIonicOperator, state: Any, bare: Any,
                   prefix: str, payload: dict[str, np.ndarray]) -> dict[str, Any]:
    meta = state.metadata
    ns, nt, npwx, nbnd = (meta[key] for key in ("ngw_source", "ngw_target", "npwx", "nbnd"))
    u = bare.array("Cartesian_complex_pattern", (6,), "complex128")
    psi = state.array("evc", (npwx, nbnd), "complex128")[:ns]
    native = bare.array("fullPWbareRHS_dvpsi_before_screening", (npwx, nbnd), "complex128")[:nt]
    scale = float(np.max(np.abs(native)))
    allowance = FACTOR * max(1., scale)
    error, baseline = comparison(operator.apply(u, psi), native, allowance)
    payload[prefix + "error_vectors_ry_per_bohr"] = error
    payload[prefix + "actual_cartesian_complex_pattern"] = u
    mutations = {}
    for name, action in operator.sensitivity_actions(u, psi).items():
        _, diagnostic = comparison(action, native, allowance)
        diagnostic["exposed"] = not diagnostic["pass"]
        diagnostic["sensitivity_credit"] = bool(baseline["pass"] and diagnostic["exposed"])
        mutations[name] = diagnostic
    return {
        "bare_document": bare.path, "npz_prefix": prefix,
        "native_identities": {key: str(bare.metadata[key]) for key in
                              ("spin", "ik", "mode", "irrep", "perturbation_index")},
        "physical_action_shape": [str(nt), str(nbnd)],
        "maxabs_native_ry_per_bohr": scale, "baseline": baseline,
        "mutations": mutations, "agreement": baseline["pass"],
    }


def evaluate_state(bundle: Any, state: Any, label: str,
                   payload: dict[str, np.ndarray]) -> dict[str, Any]:
    meta = state.metadata
    prefix = f"{label}_spin{meta['spin']}_ik{meta['ik']}_"
    result: dict[str, Any] = {
        "document": state.path, "npz_prefix": prefix,
        "native_dimensions_and_identities": {key: str(meta[key]) for key in
            ("spin", "ik", "ikk", "ikq", "ngw_source", "ngw_target", "npwx", "nbnd", "nkb")},
        "projector_shapes": {name: [str(n) for n in state.array(name).shape]
                             for name in ("vkb_source", "vkb_target", "deeq", "nh", "ityp")},
        "patterns": [], "agreement": False,
    }
    operator = ionic.BareIonicOperator.from_bundle(bundle, state)
    records = [bare for bare in bundle.bare if
               (bare.metadata["spin"], bare.metadata["ik"]) == (meta["spin"], meta["ik"])]
    if len(records) != 6 or {bare.metadata["mode"] for bare in records} != set(range(1, 7)):
        raise ValueError(f"{state.path}: all six actual Native patterns are required")
    for bare in records:
        pattern_prefix = prefix + f"mode{bare.metadata['mode']}_"
        try:
            result["patterns"].append(pattern_result(operator, state, bare, pattern_prefix, payload))
        except Exception as error:
            result["patterns"].append({"bare_document": bare.path, "npz_prefix": pattern_prefix,
                                       **failure(error)})
    result["agreement"] = all(item["agreement"] for item in result["patterns"])
    return result


def gamma_hermiticity(bundle: Any, payload: dict[str, np.ndarray]) -> dict[str, Any]:
    """Use a complete source basis, not a 20-band projection, for each real mode."""
    state = bundle.states[0]
    meta = state.metadata
    n = meta["ngw_source"]
    if (n != meta["ngw_target"] or not np.array_equal(state.array("igk_source"), state.array("igk_target"))
            or not np.array_equal(state.array("k_source"), state.array("k_target"))):
        raise ValueError("Gamma Hermiticity needs identical actual source/target PW bases")
    operator = ionic.BareIonicOperator.from_bundle(bundle, state)
    modes = []
    # One dense Gamma invariant matrix at a time; no six-matrix archive/copies.
    for mode in range(6):
        u = np.zeros(6, dtype=np.complex128)
        u[mode] = 1.
        matrix = np.empty((n, n), dtype=np.complex128)
        for start in range(0, n, 32):
            stop = min(n, start + 32)
            basis = np.zeros((n, stop - start), dtype=np.complex128)
            basis[np.arange(start, stop), np.arange(stop - start)] = 1.
            matrix[:, start:stop] = operator.apply(u, basis)
        scale = float(np.max(np.abs(matrix)))
        allowance = FACTOR * max(1., scale)
        row_maxima = np.array([np.max(np.abs(matrix[row] - matrix[:, row].conj()))
                               for row in range(n)])
        payload[f"gamma_full_operator_cartesian{mode + 1}_hermiticity_row_maxabs"] = row_maxima
        modes.append({"cartesian_mode": str(mode + 1), "full_operator_shape": [str(n), str(n)],
                      "maxabs_operator_ry_per_bohr": scale,
                      "maxabs_antihermitian_ry_per_bohr": float(np.max(row_maxima)),
                      "allowance_ry_per_bohr": allowance,
                      "pass": bool(np.all(row_maxima <= allowance))})
    return {"state_document": state.path, "complete_basis_dimension": str(n),
            "native_retained_band_count": str(meta["nbnd"]),
            "invariant_only_not_independent_native_complement_action": True,
            "real_cartesian_modes": modes, "pass": all(item["pass"] for item in modes)}


def evaluate_bundle(label: str, bundle: Any,
                    payload: dict[str, np.ndarray]) -> dict[str, Any]:
    bg = bundle.background
    q = bg.array("q_cart_2pi_alat", (3,), "float64")
    reciprocal = bg.array("bg", (3, 3), "float64")
    expected = np.zeros(3) if label == "gamma" else (.5 if label == "plus" else -.5) * reciprocal[:, 0]
    if np.max(np.abs(q - expected)) > 1e-13 or bg.metadata["lgamma"] != (label == "gamma"):
        raise ValueError(f"{label}: actual Native q/flag is not the requested separate branch")
    nksq = bg.metadata["nksq"]
    if (len(bundle.states) != nksq
            or {state.metadata["ik"] for state in bundle.states} != set(range(1, nksq + 1))):
        raise ValueError("every actual Native electronic state is required")
    results = []
    for state in bundle.states:
        try:
            results.append(evaluate_state(bundle, state, label, payload))
        except Exception as error:
            results.append({"document": state.path, **failure(error)})
    invariant = None
    if label == "gamma":
        try:
            invariant = gamma_hermiticity(bundle, payload)
        except Exception as error:
            invariant = {"pass": False, **failure(error)}
    fft = fft_module.NativeRhoFFT.from_bundle(bundle)
    return {
        "bundle": str(bundle.root),
        "actual_input_byte_bindings": [binding(bundle.root / name, raw)
                                       for name, raw in bundle.buffers.items()],
        "q_cart_2pi_alat": q.tolist(), "native_nksq": str(nksq),
        "reported_states": str(len(results)), "native_nnr": str(fft.nnr),
        "physical_fft": [str(n) for n in fft.physical],
        "storage_fft": [str(n) for n in fft.storage],
        "states": results, "gamma_full_operator_hermiticity": invariant,
        "agreement": bool(all(state["agreement"] for state in results)
                          and (invariant is None or invariant["pass"])),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for label in ("gamma", "plus", "minus"):
        parser.add_argument("--" + label, type=Path, required=True, help="new actual closed Native branch")
    parser.add_argument("--output", type=Path, required=True, help="new JSON and companion NPZ")
    args = parser.parse_args()
    archive = args.output.with_suffix(".npz")
    if args.output.suffix.lower() != ".json" or args.output.exists() or archive.exists():
        raise ValueError("choose a new .json path and absent companion .npz; preserve first failures")
    started = time.monotonic()
    payload: dict[str, np.ndarray] = {}
    report: dict[str, Any] = {
        "schema_version": 1, "started_at": datetime.now(timezone.utc).isoformat(),
        "object": "factorized complete source-to-target PW bare Cartesian ionic operator",
        "scope": "Selected ordinaryNC/collinearPBE/noNLCC/noHubbard/twoatom/serial24cubed backgrounds. Engineering-only retained-Native-bare-RHS corroboration and full-operator Gamma invariant; no independent Native arbitrary-complement action, material/coherence/loss/nuclear result.",
        "comparator": "fullPWbareRHS_dvpsi_before_screening on every physical target row, every state and all6 actual complexNative patterns",
        "bare_action_allowance": "1e-9*max(1,maxabsNative) Ry/bohr",
        "gamma_invariant_allowance": "1e-9*max(1,maxabsFullOperator) Ry/bohr; real Cartesian unit displacements",
        "source_equation": "vloc_bare_cartesian*u by actual normalized FFT + i beta_t D beta_s_dagger (p_s.u psi) - i (p_t.u) beta_t D beta_s_dagger psi; p=tpiba*(k+G)",
        "independent_full_complement_native_action_proof": False,
        "environment": {"command": [sys.executable, *sys.argv], "python": platform.python_version(),
                        "numpy": np.__version__,
                        "cpu_affinity": [str(cpu) for cpu in sorted(os.sched_getaffinity(0))],
                        "threads": {key: os.environ.get(key) for key in
                                    ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                     "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")}},
        "results": {},
    }
    try:
        report["sources"] = source_bindings()
        paths = {label: getattr(args, label).resolve() for label in ("gamma", "plus", "minus")}
        if len(set(paths.values())) != 3:
            raise ValueError("Gamma/+q/-q must be three distinct actual closed bundles")
        bundles = {label: decoder.load_bundle(path) for label, path in paths.items()}
        report["actual_input_byte_bindings"] = {
            label: [binding(bundle.root / name, raw) for name, raw in bundle.buffers.items()]
            for label, bundle in bundles.items()
        }
        manifests = {hashlib.sha256(bundle.buffers["manifest.json"]).hexdigest()
                     for bundle in bundles.values()}
        if len(manifests) != 3:
            raise ValueError("copied identical bundle manifests cannot supply separate actual q branches")
        for label, bundle in bundles.items():
            try:
                report["results"][label] = evaluate_bundle(label, bundle, payload)
            except Exception as error:
                report["results"][label] = failure(error)
    except Exception as error:
        report["input_or_source_error"] = failure(error)
    patterns = [pattern for result in report["results"].values()
                for state in result.get("states", []) for pattern in state.get("patterns", [])]
    report["mutation_exposure"] = {}
    for name in MUTATIONS:
        exposed = [pattern["npz_prefix"] for pattern in patterns
                   if pattern.get("mutations", {}).get(name, {}).get("sensitivity_credit", False)]
        report["mutation_exposure"][name] = {
            "credited_pattern_prefixes": exposed, "exposed": bool(exposed),
            "interpretation": "exposed against same Native comparator" if exposed
                              else "not exposed; no sensitivity credit",
        }
    report["all_actual_pattern_actions_and_gamma_invariant_pass"] = bool(
        len(report["results"]) == 3 and all(item["agreement"] for item in report["results"].values()))
    # Each required error family earns its own credit; unexposed families do not disappear.
    report["required_mutation_families_pass"] = all(
        report["mutation_exposure"][name]["exposed"] for name in MUTATIONS[:4])
    passed = (report["all_actual_pattern_actions_and_gamma_invariant_pass"]
              and report["required_mutation_families_pass"])
    report["acceptance_pass"] = passed
    report["elapsed_seconds"] = time.monotonic() - started
    with archive.open("xb") as stream:
        np.savez_compressed(stream, **payload)
    report["diagnostics_archive"] = binding(archive)
    text = json.dumps(report, indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as stream:
        stream.write(text)
    print(text, end="")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())

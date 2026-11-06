#!/usr/bin/env python3
"""Compare fresh closed Gamma/half-b1 native bundles with the spatial kernel.

Run from the repository with PYTHONPATH=src. JSON records full 3x3 density-
potential matrices; the companion NPZ preserves every input, native/Python
output, error and scientific mutation field, rather than only error summaries.
No finite-z, coherence, loss, full-material accuracy or nuclear claim is made.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from substrate_framework.qe_response_bundle import load_bundle
from substrate_framework.qe_response_fft import NativeRhoFFT
from substrate_framework.spin_pbe_kernel import (
    apply_spin_pbe_kernel, density_potential_pairing, kernel_sensitivity_actions,
)

_MODES = np.array([[1, 0, 1], [0, 1, -1], [0, 1, 0]], dtype=np.int32)
_COEFFICIENTS = np.array([[1 + .3j, 1 + .3j, 0], [-.4 + .2j, 0, -.4 + .2j]], dtype=np.complex128)


def complex_json(value: Any) -> dict[str, Any]:
    array = np.asarray(value)
    return {"real": array.real.tolist(), "imag": array.imag.tolist()}


def evaluate(label: str, path: Path, payload: dict[str, np.ndarray]) -> dict[str, Any]:
    bundle = load_bundle(path)
    bg = bundle.background
    fft = NativeRhoFFT.from_bundle(bundle)
    if not bundle.states or not bundle.bare:
        raise ValueError(f"{label}: complete fresh export needs state and bare-RHS documents")
    # Full requested HXC oracle, not an lrpa/lnolr restricted alternative.
    if bg.metadata["lrpa"] or bg.metadata["lnolr"]:
        raise ValueError("oracle requires full Hartree+PBE with lrpa=lnolr=false")
    q = bg.array("q_cart_2pi_alat", (3,), "float64")
    reciprocal = bg.array("bg", (3, 3), "float64")
    if label == "gamma":
        q_valid = bool(np.max(np.abs(q)) < 1e-13)
    else:
        q_valid = bool(min(np.max(np.abs(q - .5 * reciprocal[:, 0])),
                           np.max(np.abs(q + .5 * reciprocal[:, 0]))) < 1e-13)
    if not q_valid:
        raise ValueError(f"{label}: bundle q does not match the frozen requested wavevector")
    modes = bg.array("kernel_probe_modes", (3, 3), "int32")
    coefficients = bg.array("kernel_probe_coefficients", (2, 3), "complex128")
    if not np.array_equal(modes, _MODES) or not np.array_equal(coefficients, _COEFFICIENTS):
        raise ValueError("native probe definitions differ from the preregistered three fields")
    inputs = bg.array("kernel_probe_in", (fft.nnr, 2, 3), "complex128")
    native = bg.array("kernel_probe_out", inputs.shape, "complex128")
    # Direct trigonometric construction is independent of the native FFT and
    # of the FFT implementation used in apply_spin_pbe_kernel.
    grid = np.indices(fft.physical)
    coordinates = grid / np.array(fft.physical)[:, None, None, None]
    expected = np.zeros(inputs.shape, dtype=np.complex128)
    for p in range(3):
        phase = np.exp(2j * np.pi * np.einsum("c,cijk->ijk", modes[:, p], coordinates))
        fft.physical_view(expected)[:, :, :, :, p] = 1e-4 * phase[:, :, :, None] * coefficients[:, p]
    input_error = float(np.max(np.abs(inputs - expected)))
    if input_error > 1e-15:
        raise ValueError(f"native probe input normalization/phase mismatch: {input_error}")
    python = np.stack([apply_spin_pbe_kernel(bundle, inputs[:, :, p]) for p in range(3)], axis=2)
    error = python - native
    limits = 1e-9 * np.maximum(1., np.max(np.abs(native), axis=(0, 1)))
    maxima = np.max(np.abs(error), axis=(0, 1))
    native_matrix = np.array([[density_potential_pairing(bundle, inputs[:, :, i], native[:, :, j])
                               for j in range(3)] for i in range(3)])
    python_matrix = np.array([[density_potential_pairing(bundle, inputs[:, :, i], python[:, :, j])
                               for j in range(3)] for i in range(3)])
    prefix = f"{label}_"
    payload.update({prefix + "input": inputs, prefix + "native_output": native,
                    prefix + "python_output": python, prefix + "error": error,
                    prefix + "native_pairing_matrix": native_matrix,
                    prefix + "python_pairing_matrix": python_matrix,
                    prefix + "pairing_matrix_error": python_matrix - native_matrix})
    mutation_fields: dict[str, list[np.ndarray]] = {}
    for p in range(3):
        for name, field in kernel_sensitivity_actions(bundle, inputs[:, :, p]).items():
            mutation_fields.setdefault(name, []).append(field)
    sensitivity = {}
    for name, fields in mutation_fields.items():
        action = np.stack(fields, axis=2)
        difference = action - native
        separation = np.max(np.abs(difference), axis=(0, 1))
        payload[prefix + name] = action
        payload[prefix + name + "_error"] = difference
        sensitivity[name] = {
            "max_abs_error_ry_per_probe": separation.tolist(),
            "exposed_by_probe": (separation > limits).tolist(),
            "interpretation": "Only probes exceeding the frozen native comparison allowance expose this mutation.",
        }
    return {
        "bundle": str(bundle.root),
        "closed_manifest_sha256": hashlib.sha256(bundle.buffers["manifest.json"]).hexdigest(),
        "q_cart_2pi_alat": q.tolist(), "physical_fft": fft.physical,
        "storage_fft": fft.storage, "nnr": fft.nnr,
        "native_input_max_abs_error_bohr_minus3": input_error,
        "native_output_max_abs_ry_per_probe": np.max(np.abs(native), axis=(0, 1)).tolist(),
        "allowance_ry_per_probe": limits.tolist(),
        "max_abs_error_ry_per_probe": maxima.tolist(),
        "max_abs_error_ry_per_spin_probe": np.max(np.abs(error), axis=0).tolist(),
        "agreement_per_probe": (maxima <= limits).tolist(),
        "native_pairing_matrix": complex_json(native_matrix),
        "python_pairing_matrix": complex_json(python_matrix),
        "pairing_matrix_error": complex_json(python_matrix - native_matrix),
        "pairing_units": "Ry",
        "pairing_convention": "omega/Ngrid sum_r,spin conj(delta_n)*delta_v; no extra spin factor",
        "sensitivity": sensitivity,
        "electronic_state_documents": [document.path for document in bundle.states],
        "bare_operator_documents": [document.path for document in bundle.bare],
        "full_verified_files": sorted(bundle.buffers),
        "operator_agreement": bool(np.all(maxima <= limits)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gamma", type=Path, required=True, help="closed native Gamma q directory")
    parser.add_argument("--half", type=Path, required=True, help="closed native half-b1 q directory")
    parser.add_argument("--output", type=Path, required=True, help="new JSON receipt path; companion .npz stores full fields")
    args = parser.parse_args()
    archive = args.output.with_suffix(".npz")
    if args.output.exists() or archive.exists():
        raise ValueError("preserve existing receipts: choose new JSON/NPZ output paths")
    payload: dict[str, np.ndarray] = {}
    results = {label: evaluate(label, path, payload)
               for label, path in (("gamma", args.gamma), ("half", args.half))}
    passed = all(result["operator_agreement"] for result in results.values())
    np.savez_compressed(archive, **payload)
    report = {
        "schema_version": 1,
        "object": "native spatial collinear Hartree+PBE-GGA density-to-potential action",
        "frozen_allowance": "1e-9*max(1,maxabs(native_output)) Ry, separately for each fixed probe",
        "results": results, "operator_agreement": passed,
        "full_fields_archive": {"path": str(archive), "bytes": str(archive.stat().st_size),
                                "sha256": hashlib.sha256(archive.read_bytes()).hexdigest()},
        "scope": "Registered engineering background; kernel agreement does not establish finite-z response, material convergence, coherence, loss or nuclear bridge.",
    }
    text = json.dumps(report, indent=2, allow_nan=False) + "\n"
    args.output.write_text(text)
    print(text, end="")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())

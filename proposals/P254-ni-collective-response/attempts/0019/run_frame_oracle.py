#!/usr/bin/env python3
"""Construct actual full-PW/full-FD frames and save bounded engineering evidence.

Run with shared SF venv, CPU6 affinity and BLAS1; no Native process is launched.
--output must be a new JSON path with an absent companion density-only NPZ.
The saved-SCF density difference is measured, never used as a rejection gate.
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
import scipy
import sympy as sp

import substrate_framework.pw_response_data as interface
import substrate_framework.qe_ionic_operator as ionic
import substrate_framework.qe_pw_hamiltonian as hamiltonian
import substrate_framework.qe_pw_response_frames as frames
import substrate_framework.qe_response_bundle as decoder
import substrate_framework.qe_response_fft as fft_module
import substrate_framework.qe_response_schema as schema


ATTEMPT = Path(__file__).resolve().parent


def binding(path: Path, data: bytes | None = None) -> dict[str, str]:
    raw = path.read_bytes() if data is None else data
    return {"path": str(path), "bytes": str(len(raw)), "sha256": hashlib.sha256(raw).hexdigest()}


def source_bindings() -> dict[str, Any]:
    inventory = json.loads((ATTEMPT / "source-inventory.json").read_bytes())
    registered = []
    for record in inventory["sources"]:
        actual = binding(Path(record["path"]))
        unchanged = all(actual[key] == record[key] for key in ("bytes", "sha256"))
        if not unchanged:
            raise ValueError(f"settled registered source changed: {record['path']}")
        registered.append({"actual": actual, "registered": record, "unchanged": unchanged})
    return {"preregistration_commit": "f78b9cad", "registration": binding(ATTEMPT / "registration.json"),
            "source_inventory": binding(ATTEMPT / "source-inventory.json"), "registered_sources": registered,
            "actual_consumer_sources": [binding(Path(module.__file__).resolve()) for module in
                                        (interface, ionic, hamiltonian, frames, decoder, fft_module, schema)]
                                       + [binding(Path(__file__).resolve())]}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def evaluate_models(models: dict[str, interface.PWResponseModel],
                    bundles: dict[str, decoder.QEResponseBundle],
                    payload: dict[str, np.ndarray]) -> dict[str, Any]:
    """Independently check frame shapes, shared sources, FD, Parseval and N."""
    _require(set(models) == {"gamma", "half"}, "frame builder must return Gamma and half models")
    gamma, half = models["gamma"], models["half"]
    fft = gamma.fft
    _require(len(gamma.channels) == len(half.channels) == len(bundles["gamma"].states),
             "all actual source k/spin channels must participate")
    density = np.zeros((fft.nnr, 2))
    spin_counts = np.zeros(2)
    channel_reports = []
    mu = float(gamma.metadata["common_mu_Ry"])
    width = float(gamma.metadata["actual_degauss_Ry"])
    for index, (channel, other) in enumerate(zip(gamma.channels, half.channels)):
        _require(channel.identity == other.identity and channel.spin in (0, 1)
                 and channel.spin == other.spin, "explicit source k/spin identity mismatch")
        for name in ("source_energies", "weighted_occupations", "source_grid"):
            _require(getattr(channel, name) is getattr(other, name), f"source arrays not shared: {name}")
        _require(channel.g_plus is channel.g_minus and channel.plus_vectors is channel.minus_vectors,
                 "Gamma must legitimately reuse its actual same-q0 operator and eigenspace")
        _require(other.g_plus is not other.g_minus, "half requires distinct actual +/-q ionic actions")
        active = channel.source_energies.size
        _require(channel.source_grid.shape == (fft.nnr, active)
                 and channel.weighted_occupations.shape == (active,), "source dimensions disagree")
        weight = gamma.metadata["source_channels"][index]["quadrature_weight"]
        # logaddexp avoids positive-exponent overflow, retaining subnormal tails.
        expected = weight * np.exp(-np.logaddexp(0., (channel.source_energies - mu) / width))
        full_fd = np.exp(-np.logaddexp(0., (channel.plus_energies - mu) / width))
        _require(np.array_equal(channel.source_energies, channel.plus_energies[full_fd != 0]),
                 "every representably nonzero full-FD source must participate")
        fd_error = float(np.max(np.abs(channel.weighted_occupations - expected), initial=0))
        _require(fd_error <= 16 * np.finfo(float).eps * max(1., weight), "common-mu FD/no-extra2 mismatch")
        _require(np.isfinite(channel.source_grid).all() and np.all(channel.weighted_occupations >= 0),
                 "nonfinite source or negative physical occupation")
        grid_norm = np.mean(np.abs(fft.physical_view(channel.source_grid)) ** 2, axis=(0, 1, 2))
        norm_error = float(np.max(np.abs(grid_norm - 1), initial=0))
        _require(norm_error <= 1e-10, "independent source FFT Parseval check failed")
        density[:, channel.spin] += np.sum(np.abs(channel.source_grid) ** 2
                                           * channel.weighted_occupations[None, :], axis=1) / fft.omega
        spin_counts[channel.spin] += channel.weighted_occupations.sum()
        branches = {}
        for model_name, current in (("gamma", channel), ("half", other)):
            for side in ("plus", "minus"):
                energies = getattr(current, side + "_energies")
                vectors = getattr(current, side + "_vectors")
                slots = getattr(current, side + "_slots")
                action = getattr(current, "g_" + side)
                n = energies.size
                label = "gamma" if model_name == "gamma" else side
                branch_record = gamma.metadata["source_channels"][index]["branches"][label]
                actual_state = bundles[label].documents[branch_record["state_document"]]
                _require(n == actual_state.metadata["ngw_target"], "target PW spectrum was truncated")
                actual_fft = fft_module.NativeRhoFFT.from_bundle(bundles[label])
                native_slots = actual_fft.nl[actual_state.array("igk_target").astype(np.int64) - 1]
                shift = np.asarray(branch_record["source_reciprocal_shift"], dtype=np.int64)
                native_coordinates = np.asarray(np.unravel_index(native_slots, fft.storage, order="F"))
                shifted_coordinates = tuple((native_coordinates[i] + shift[i]) % fft.physical[i]
                                            for i in range(3))
                expected_slots = np.ravel_multi_index(shifted_coordinates, fft.storage, order="F")
                _require(np.array_equal(slots, expected_slots),
                         "target slots must transform actual Native map into the common source envelope")
                _require(vectors.shape == (n, n) and slots.shape == (n,)
                         and action.shape == (6, n, active), "complete target/g dimensions disagree")
                _require(np.unique(slots).size == n and np.all((slots >= 0) & (slots < fft.nnr)),
                         "target FFT slots invalid or duplicated")
                if np.any(shift):
                    # Independent direct Fourier phase comparison, not the production Miller map.
                    points = np.array([[.13, .21, .37], [.07, .42, .31], [.29, .17, .11]])
                    miller = bundles[label].background.array("mill")[:, actual_state.array("igk_target") - 1]
                    coefficients = np.exp(.37j * np.arange(n)) / n
                    actual_wave = coefficients @ np.exp(2j * np.pi * miller.T @ points.T)
                    common_wave = coefficients @ np.exp(2j * np.pi * (miller + shift[:, None]).T @ points.T)
                    expected_wave = actual_wave * np.exp(2j * np.pi * points @ shift)
                    gauge_error = float(np.max(np.abs(common_wave - expected_wave)))
                    _require(gauge_error <= 1e-10, "independent Bloch-envelope phase check failed")
                    for wrong_shift in (np.zeros(3), -shift):
                        wrong_wave = coefficients @ np.exp(2j * np.pi *
                                                           (miller + wrong_shift[:, None]).T @ points.T)
                        _require(np.max(np.abs(wrong_wave - expected_wave)) > 1e-10,
                                 "omitted/opposite reciprocal-envelope mutation was not exposed")
                _require(np.isfinite(energies).all() and np.isfinite(vectors).all()
                         and np.isfinite(action).all() and np.all(np.diff(energies) >= 0),
                         "target spectra/actions must be finite and complete spectra sorted")
                eigen_norm = float(np.max(np.abs(np.sum(np.abs(vectors) ** 2, axis=0) - 1)))
                _require(eigen_norm <= 64 * np.finfo(float).eps * n, "target eigenvector norm failed")
                branches[model_name + "_" + side] = {"complete_dimension": str(n),
                    "ionic_shape": [str(value) for value in action.shape], "column_norm_maxabs_error": eigen_norm}
        channel_reports.append({"identity": channel.identity, "spin_zero_based": str(channel.spin),
                                "active_dimension": str(active), "FD_maxabs_error": fd_error,
                                "source_grid_Parseval_maxabs_error": norm_error, "branches": branches})
    integral = fft.omega * np.mean(fft.physical_view(density), axis=(0, 1, 2))
    nelec = float(bundles["gamma"].background.metadata["nelec"])
    _require(abs(float(spin_counts.sum()) - nelec) <= 1e-10, "independent total-N FD check failed")
    _require(np.max(np.abs(integral - spin_counts)) <= 1e-10, "independent density/spin integral check failed")
    _require(np.allclose(half.q_cart_2pi_alat, bundles["plus"].background.array("q_cart_2pi_alat"),
                        rtol=0, atol=2e-12), "half envelope q does not match actual positive branch")
    actual = bundles["gamma"].background.array("rho_total_magnetization")
    saved = np.column_stack(((actual[:, 0] + actual[:, 1]) / 2, (actual[:, 0] - actual[:, 1]) / 2))
    payload.update(reference_density_up_down_bohr_minus3=density,
                   saved_SCF_density_up_down_bohr_minus3=saved,
                   reference_minus_saved_up_down_bohr_minus3=density - saved,
                   q_gamma_2pi_alat=gamma.q_cart_2pi_alat, q_half_2pi_alat=half.q_cart_2pi_alat)
    return {"models": {name: model.metadata for name, model in models.items()}, "channels": channel_reports,
            "independent_spin_occupancies_electrons": spin_counts.tolist(),
            "independent_density_integral_up_down_electrons": integral.tolist(),
            "independent_total_N_error_electrons": abs(float(spin_counts.sum()) - nelec),
            "normalization_total_N_and_explicit_branch_assertions_pass": True,
            "shared_sources_and_Gamma_same_q_operator_pass": True,
            "saved_density_difference_is_rejection_gate": False}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for label in ("gamma", "plus", "minus"):
        parser.add_argument("--" + label, type=Path, required=True, help="actual closed Native per-q bundle")
    parser.add_argument("--output", type=Path, required=True, help="new JSON and compact density NPZ")
    args = parser.parse_args()
    archive = args.output.with_suffix(".npz")
    if args.output.suffix.lower() != ".json" or args.output.exists() or archive.exists():
        raise ValueError("choose a new .json path and absent companion .npz; preserve first failures")
    started = time.monotonic()
    payload: dict[str, np.ndarray] = {}
    report: dict[str, Any] = {"schema_version": 1, "started_at": datetime.now(timezone.utc).isoformat(),
        "object": "complete H0 spectra/common-mu full-FD explicit-two-branch ionic response frames",
        "scope": "Engineering applicability only; full-FD spectral H0 density is not saved Native SCF rho. No response evaluation, Native execution, convergence, coherence-domain, physical-loss or nuclear-rate claim.",
        "source_convention": "single-spin f_FD times actual quadrature weight; no extra2; unnormalised periodic PW grid; physical density 1/omega",
        "energy_units": "Ry", "FD_total_N_absolute_allowance_electrons": 1e-10,
        "environment": {"command": [sys.executable, *sys.argv], "python": platform.python_version(),
                        "numpy": np.__version__, "scipy": scipy.__version__,
                        "cpu_affinity": [str(cpu) for cpu in sorted(os.sched_getaffinity(0))],
                        "threads": {name: os.environ.get(name) for name in
                                    ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                                     "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")}}, "acceptance_pass": False}
    try:
        report["sources"] = source_bindings()
        k, g, shift, x = sp.symbols("k g shift x", real=True)
        residual = sp.expand(((k + shift) + (g - shift)) * x - (k + g) * x)
        wrong = sp.expand(((k + shift) + (g + shift)) * x - (k + g) * x)
        _require(residual == 0 and wrong != 0, "exact Bloch reciprocal-gauge identity/mutation failed")
        report["exact_reciprocal_gauge_oracle"] = {
            "correct_phase_residual": str(residual), "wrong_sign_residual": str(wrong),
            "repair_registration": binding(ATTEMPT / "repair-registration-001.json")}
        paths = {label: getattr(args, label).resolve() for label in ("gamma", "plus", "minus")}
        bundles = {label: decoder.load_bundle(path) for label, path in paths.items()}
        report["actual_input_byte_bindings"] = {
            label: [binding(bundle.root / name, raw) for name, raw in bundle.buffers.items()]
            for label, bundle in bundles.items()}
        build_started = time.monotonic()
        models = frames.build_pw_response_models(bundles["gamma"], bundles["plus"], bundles["minus"])
        report["frame_build_seconds"] = time.monotonic() - build_started
        report["results"] = evaluate_models(models, bundles, payload)
        report["acceptance_pass"] = True
    except Exception as error:
        report["failure"] = {"error_type": type(error).__name__, "error": str(error),
                             "traceback": traceback.format_exc()}
    report["elapsed_seconds"] = time.monotonic() - started
    report["registered_180s_elapsed_expectation_met"] = report["elapsed_seconds"] <= 180
    with archive.open("xb") as stream:
        np.savez_compressed(stream, **payload)
    report["density_diagnostics_archive"] = binding(archive)
    text = json.dumps(report, indent=2, allow_nan=False) + "\n"
    with args.output.open("x") as stream:
        stream.write(text)
    print(json.dumps({"output": str(args.output), "acceptance_pass": report["acceptance_pass"],
                      "elapsed_seconds": report["elapsed_seconds"], "failure": report.get("failure")}, indent=2))
    return 0 if report["acceptance_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

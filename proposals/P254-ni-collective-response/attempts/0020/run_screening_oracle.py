#!/usr/bin/env python3
"""Finite cross-spin screening oracle, not Native data or a material claim.

Run from the repository with PYTHONPATH=src and --output NEW.json. Only the
immutable 0019 fixture constructor is reused. The comparator assembles and
solves the coupled, original-PW-basis density-matrix EOM, including both spin
sectors in ONE linear system. It never calls the response or screening solver.
Every case retains failures independently; existing receipts are never replaced.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from time import monotonic, sleep
import traceback

import numpy as np
import scipy

import substrate_framework.pw_independent_response as independent_module
import substrate_framework.pw_response_data as data_module
import substrate_framework.pw_screened_response as screened_module
import substrate_framework.qe_response_fft as fft_module
from substrate_framework.pw_independent_response import PWIndependentResponse
from substrate_framework.pw_response_data import PWResponseValue
from substrate_framework.pw_screened_response import ScreenedResponseError, solve_screened_response


ATTEMPT = Path(__file__).resolve().parent
FIXTURE_SOURCE = ATTEMPT.parent / "0019" / "run_two_branch_oracle.py"
ALLOWANCE = 1e-10
FREQUENCIES = (0.02 + 0.02j, 0.05 + 0.02j)


def load_fixture_source():
    spec = importlib.util.spec_from_file_location("p254_0019_immutable_fixture", FIXTURE_SOURCE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def binding(path):
    data = path.read_bytes()
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def maxabs(value):
    return float(np.max(np.abs(value), initial=0.0))


def encoded(value):
    value = np.asarray(value)
    return {"real": value.real.tolist(), "imag": value.imag.tolist()}


def compare(actual, expected):
    result = {}
    for name in ("density", "ionic_gradient"):
        reference = getattr(expected, name)
        error = maxabs(getattr(actual, name) - reference)
        scale = max(1.0, maxabs(reference))
        result[name] = {"maxabs_error": error, "scale": scale,
                        "allowance": ALLOWANCE * scale, "agreement": error <= ALLOWANCE * scale}
    result["agreement"] = all(result[name]["agreement"] for name in ("density", "ionic_gradient"))
    return result


def kernel_coefficients(f):
    count = f.envelopes.shape[1]
    position = 2 * np.pi * np.arange(count) / count
    coefficients = np.empty((count, 2, 2))
    coefficients[:, 0, 0] = 0.62 + 0.19 * np.cos(position)
    coefficients[:, 1, 1] = 0.39 + 0.11 * np.sin(position)
    coefficients[:, 0, 1] = 0.27 + 0.08 * np.cos(position + 0.4)
    coefficients[:, 1, 0] = coefficients[:, 0, 1]
    return coefficients


def consumer_kernel(coefficients, nnr):
    def action(density):
        result = np.zeros((nnr, 2), dtype=complex)
        count = coefficients.shape[0]
        result[:count] = np.einsum("rst,rt->rs", coefficients, density[:count])
        return result
    return action


def direct_eom(f, z, coefficients, u, check_deadline):
    """Solve z*dR-[H,dR]-[V(K n[dR]),rho0]=[g*u,rho0].

    Full original-basis matrices include every k block, not only the driven
    sector. Local-potential and density maps are literal Fourier-sum loops;
    neither PWIndependentResponse nor its contractions enter this comparator.
    """
    spins, nk, npw, _ = f.h.shape
    ngrid = f.envelopes.shape[1]
    size = nk * npw
    dimension = spins * size * size
    h = np.zeros((spins, size, size), dtype=complex)
    rho = np.zeros_like(h)
    ionic = np.zeros((spins, 6, size, size), dtype=complex)
    for s in range(spins):
        for k in range(nk):
            source = slice(k * npw, (k + 1) * npw)
            target_k = (k + f.step) % nk
            target = slice(target_k * npw, (target_k + 1) * npw)
            h[s, source, source] = f.h[s, k]
            rho[s, source, source] = (f.vectors[s, k] * f.weights[s, k]) @ f.vectors[s, k].conj().T
            ionic[s, :, target, source] = f.gq[s, :, k]

    def density_of(delta):
        density = np.zeros((f.model.fft.nnr, spins), dtype=complex)
        for s in range(spins):
            for k in range(nk):
                target_k = (k + f.step) % nk
                block = delta[s, target_k * npw:(target_k + 1) * npw,
                              k * npw:(k + 1) * npw]
                for r in range(ngrid):
                    density[r, s] += (f.envelopes[target_k, r] @ block
                                      @ f.envelopes[k, r].conj()) / f.model.fft.omega
        return density

    def induced_drive(delta):
        density = density_of(delta)
        potential = np.zeros((ngrid, spins), dtype=complex)
        # Independent explicit spin contraction, NOT the consumer kernel call.
        for r in range(ngrid):
            for s in range(spins):
                for t in range(spins):
                    potential[r, s] += coefficients[r, s, t] * density[r, t]
        drive = np.zeros_like(h)
        for s in range(spins):
            for k in range(nk):
                target_k = (k + f.step) % nk
                block = np.zeros((npw, npw), dtype=complex)
                for r in range(ngrid):
                    block += potential[r, s] * np.outer(
                        f.envelopes[target_k, r].conj(), f.envelopes[k, r]) / ngrid
                drive[s, target_k * npw:(target_k + 1) * npw,
                      k * npw:(k + 1) * npw] = block
        return drive

    def action(delta):
        potential = induced_drive(delta)
        return z * delta - h @ delta + delta @ h - potential @ rho + rho @ potential

    ionic_drive = np.zeros_like(h)
    for a in range(6):
        ionic_drive += u[a] * ionic[:, a]
    rhs = ionic_drive @ rho - rho @ ionic_drive
    operator = np.empty((dimension, dimension), dtype=complex)
    unit = np.zeros_like(h)
    for column in range(dimension):
        check_deadline()
        unit.flat[column] = 1
        operator[:, column] = action(unit).reshape(-1)
        unit.flat[column] = 0
    check_deadline()
    delta = np.linalg.solve(operator, rhs.reshape(-1)).reshape(h.shape)
    check_deadline()
    density = density_of(delta)
    gradient = np.zeros(6, dtype=complex)
    for s in range(spins):
        for a in range(6):
            gradient[a] += np.trace(ionic[s, a].conj().T @ delta[s])
    residual = maxabs(action(delta) - rhs) / max(1.0, maxabs(rhs))
    stationary = maxabs(h @ rho - rho @ h)
    return PWResponseValue(density, gradient), {
        "eom_scaled_residual": residual, "stationary_reference_maxabs": stationary,
        "dimension": dimension, "full_original_basis": True,
        "agreement": residual <= ALLOWANCE and stationary <= ALLOWANCE,
    }


def json_safe(value):
    """Keep even nonfinite failure diagnostics bankable as strict JSON."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return {"nonfinite_float": str(value)}
    if isinstance(value, np.generic):
        return value.item()
    return value


def failure(error):
    record = {"type": type(error).__name__, "message": str(error),
              "traceback": traceback.format_exc()}
    if isinstance(error, ScreenedResponseError):
        record.update(reason=error.reason, diagnostics=json_safe(error.diagnostics))
    return record


def run_suite(report, fixture_source, deadline):
    def check_deadline():
        if monotonic() >= deadline:
            raise TimeoutError("finite-oracle monotonic wall cap reached")

    def solve(response, kernel, u, **options):
        check_deadline()
        settings = {"rtol": ALLOWANCE, "max_matvec": 60,
                    "wall_seconds": deadline - monotonic()}
        settings.update(options)
        return solve_screened_response(response, kernel, u, **settings)

    cases = report["cases"] = []

    def capture(name, call):
        record = {"name": name}
        cases.append(record)
        try:
            check_deadline()
            record.update(call())
        except Exception as error:
            record.update(agreement=False, error=failure(error))

    def completed(value):
        return {"density": encoded(value.density), "ionic_gradient": encoded(value.ionic_gradient),
                "induced_potential": encoded(value.induced_potential), "diagnostics": value.diagnostics}

    def physical_case(f, coefficients, z, u):
        response = PWIndependentResponse(f.model, z)
        actual = solve(response, consumer_kernel(coefficients, f.model.fft.nnr), u)
        expected, eom = direct_eom(f, z, coefficients, u, check_deadline)
        comparison = compare(actual, expected)
        integral = f.model.fft.omega / f.envelopes.shape[1] * np.sum(
            actual.density[:f.envelopes.shape[1]], axis=0)
        return {"z_ry": [z.real, z.imag], "drive": encoded(u), "actual": completed(actual),
                "expected": {"density": encoded(expected.density), "ionic_gradient": encoded(expected.ionic_gradient)},
                "comparison": comparison, "direct_EOM": eom,
                "per_spin_density_integral": encoded(integral),
                "agreement": comparison["agreement"] and eom["agreement"]
                and (f.step != 0 or maxabs(integral) <= ALLOWANCE)}

    def zero_kernel_case(f, z, u):
        response = PWIndependentResponse(f.model, z)
        coefficients = np.zeros_like(kernel_coefficients(f))
        actual = solve(response, consumer_kernel(coefficients, f.model.fft.nnr), u)
        isolated = response.apply(np.zeros((f.model.fft.nnr, 2), dtype=complex), u)
        expected, eom = direct_eom(f, z, coefficients, u, check_deadline)
        comparison, isolated_comparison = compare(actual, expected), compare(actual, isolated)
        return {"actual": completed(actual), "comparison_to_EOM": comparison,
                "comparison_to_isolated": isolated_comparison, "direct_EOM": eom,
                "agreement": comparison["agreement"] and isolated_comparison["agreement"]
                and eom["agreement"] and maxabs(actual.induced_potential) == 0}

    def mutation_case(f, z, u, mutation):
        coefficients = kernel_coefficients(f)
        expected, eom = direct_eom(f, z, coefficients, u, check_deadline)
        altered = coefficients.copy()
        if mutation == "wrong_screening_sign":
            altered *= -1
        elif mutation == "dropped_cross_spin":
            altered[:, 0, 1] = 0
            altered[:, 1, 0] = 0
        elif mutation == "mixed_spin_basis":
            # Wrong total/magnetization conversion on already up/down input.
            altered = altered @ np.array([[0.5, 0.5], [0.5, -0.5]])
        else:
            raise ValueError(mutation)
        mutant = solve(PWIndependentResponse(f.model, z), consumer_kernel(altered, f.model.fft.nnr), u)
        comparison = compare(mutant, expected)
        exposed = not comparison["density"]["agreement"] and not comparison["ionic_gradient"]["agreement"]
        return {"mutation": mutation, "mutant": completed(mutant), "comparison": comparison,
                "direct_EOM": eom, "exposed": exposed,
                "agreement": exposed and eom["agreement"]}

    for step in (0, 1):
        f = fixture_source.fixture(step)
        coefficients = kernel_coefficients(f)
        rng = np.random.default_rng(202020 + step)
        u = 0.23 * (rng.normal(size=6) + 1j * rng.normal(size=6))
        for index, z in enumerate(FREQUENCIES):
            capture(f"qstep{step}_frequency{index}_complex_drive",
                    lambda: physical_case(f, coefficients, z, u))
        capture(f"qstep{step}_zero_kernel", lambda: zero_kernel_case(f, FREQUENCIES[0], u))
        for mutation in ("wrong_screening_sign", "dropped_cross_spin", "mixed_spin_basis"):
            capture(f"qstep{step}_{mutation}", lambda: mutation_case(f, FREQUENCIES[0], u, mutation))

    # Operational checks use the same finite physical fixture, not mock Native data.
    f = fixture_source.fixture(1)
    response = PWIndependentResponse(f.model, FREQUENCIES[0])
    kernel = consumer_kernel(kernel_coefficients(f), f.model.fft.nnr)
    u = np.array([0.3 + 0.2j, -0.1j, 0.4, -0.2 + 0.3j, 0.1, 0.17j])

    def expected_failure(name, call, reasons):
        try:
            value = call()
        except ScreenedResponseError as error:
            details = failure(error)
            return {"expected_failure": name, "error": details,
                    "agreement": error.reason in reasons
                    and error.diagnostics["completed_outputs"] == 0
                    and error.diagnostics["response_actions"] <= error.diagnostics["max_matvec"]}
        return {"expected_failure": name, "unexpected_result": completed(value), "agreement": False}

    capture("physical_nonconvergence_work_cap", lambda: expected_failure(
        "hard work cap before convergence", lambda: solve(response, kernel, u, max_matvec=4), {"work_limit"}))
    capture("unattainable_accuracy_not_partial_success", lambda: expected_failure(
        "unattainable accuracy", lambda: solve(response, kernel, u, rtol=1e-30),
        {"nonconvergence", "work_limit"}))

    def slow_kernel(density):
        sleep(0.03)
        return kernel(density)

    capture("wall_cap_during_physical_action", lambda: expected_failure(
        "wall cap", lambda: solve(response, slow_kernel, u, wall_seconds=0.02), {"wall_limit"}))

    def zero_case(empty_occupations):
        model = f.model
        if empty_occupations:
            model = replace(model, channels=tuple(replace(channel,
                weighted_occupations=np.zeros_like(channel.weighted_occupations)) for channel in model.channels))
        value = solve(PWIndependentResponse(model, FREQUENCIES[0]), kernel,
                      u if empty_occupations else np.zeros(6), max_matvec=2)
        return {"actual": completed(value), "agreement": maxabs(value.density) == 0
                and maxabs(value.ionic_gradient) == 0 and maxabs(value.induced_potential) == 0
                and value.diagnostics["response_actions"] == 2
                and value.diagnostics["gmres_matvecs"] == 0
                and value.diagnostics["relative_residual"] == 0}

    class FinalDriveFault:
        """Deliberate action inconsistency, solely testing the final audit."""
        model = response.model
        z = response.z

        def __init__(self):
            self.ionic_calls = 0

        def apply(self, potential, displacement):
            if np.any(displacement):
                self.ionic_calls += 1
                if self.ionic_calls > 1:
                    displacement = 1.01 * displacement
            return response.apply(potential, displacement)

    capture("independent_final_residual_rejects_changed_action", lambda: expected_failure(
        "final ionic drive changed by one percent after the bare action",
        lambda: solve(FinalDriveFault(), consumer_kernel(
            np.zeros_like(kernel_coefficients(f)), f.model.fft.nnr), u),
        {"nonconvergence"}))

    capture("zero_drive", lambda: zero_case(False))
    capture("zero_response_nonzero_drive", lambda: zero_case(True))

    def invalid_tolerance(value):
        try:
            solve(response, kernel, u, rtol=value)
        except ValueError as error:
            return {"rejected": True, "exception": str(error), "agreement": True}
        return {"rejected": False, "agreement": False}

    for label, tolerance in (("zero", 0.0), ("negative", -1e-4), ("nan", float("nan"))):
        capture(f"reject_{label}_rtol", lambda: invalid_tolerance(tolerance))
    report["agreement"] = bool(cases) and all(case["agreement"] for case in cases)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new JSON receipt; never overwrite")
    parser.add_argument("--wall-seconds", type=float, default=120.0)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("--output must name a new file")
    if not np.isfinite(args.wall_seconds) or args.wall_seconds <= 0:
        parser.error("--wall-seconds must be finite and positive")
    started = monotonic()
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(), "preregistration_commit": "eacf18e0",
        "scaled_allowance": ALLOWANCE, "solver_rtol": ALLOWANCE, "wall_cap_seconds": args.wall_seconds,
        "input_scope": "synthetic finite original-PW-basis four-k systems; no Native execution",
        "equation": "z*dR-[H,dR]-[V(K n[dR]),rho0]=[g*u,rho0]",
        "screening_equation": "(I-R0_nn(z) K_q)n=R0_ng(z)u",
        "versions": {"numpy": np.__version__, "scipy": scipy.__version__},
        "not_claimed": ["material accuracy", "physical damping or loss", "coherence domain",
                        "nuclear rates", "selfconsistent Native material reference"],
        "source_bindings": [], "cases": [], "agreement": False,
    }
    try:
        paths = [Path(module.__file__).resolve() for module in
                 (independent_module, screened_module, data_module, fft_module)]
        paths += [Path(__file__).resolve(), FIXTURE_SOURCE, ATTEMPT / "registration.json",
                  ATTEMPT / "source-inventory.json"]
        for path in paths:
            report["source_bindings"].append(binding(path))
        run_suite(report, load_fixture_source(), started + args.wall_seconds)
    except Exception as error:
        report["agreement"] = False
        report["error"] = failure(error)
    report["elapsed_seconds"] = monotonic() - started
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"output": str(args.output), "agreement": report["agreement"]}))
    return 0 if report["agreement"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

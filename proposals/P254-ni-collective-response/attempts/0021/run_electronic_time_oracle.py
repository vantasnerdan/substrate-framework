#!/usr/bin/env python3
"""Finite original-basis commutator oracle for the full-PW time generator.

PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python proposals/P254-ni-collective-response/attempts/0021/run_electronic_time_oracle.py \
  --output NEW.json

Root supplies the registered external 120-second wall cap / six-CPU affinity.
Only immutable 0019 mathematical finite EigenKets are inputs, NOT Native time
trajectories. No new random Hamiltonian fixture or electronic approximation is
introduced. All random drives and initial amplitudes have fixed seeds.
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
import traceback

import numpy as np
from scipy.linalg import expm

import substrate_framework.pw_electronic_dynamics as dynamics_module
import substrate_framework.pw_independent_response as response_module
import substrate_framework.pw_response_data as data_module
import substrate_framework.pw_transition_space as space_module
import substrate_framework.qe_response_fft as fft_module
from substrate_framework.pw_electronic_dynamics import PWElectronicDynamics
from substrate_framework.pw_independent_response import PWIndependentResponse
from substrate_framework.pw_response_data import PWResponseValue


ATTEMPT = Path(__file__).resolve().parent
FIXTURE_PATH = ATTEMPT.parent / "0019" / "run_two_branch_oracle.py"
FIXTURE_SHA256 = "58bc1b7812d81d81ba73aeeffeac117b3711587f6d5c1094bce48f029e5b28db"
ALLOWANCE = 1e-10
FREQUENCIES = (0.02 + 0.02j, 0.05 + 0.02j)
TIMES = (0.07, 0.31, 1.0)


def binding(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def load_fixture():
    if binding(FIXTURE_PATH)["sha256"] != FIXTURE_SHA256:
        raise ValueError("immutable 0019 finite fixture source hash mismatch")
    spec = importlib.util.spec_from_file_location("p254_0019_immutable_fixture", FIXTURE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def maxabs(value):
    return float(np.max(np.abs(value), initial=0.))


def array_comparison(actual, expected):
    error = maxabs(actual - expected)
    scale = max(1., maxabs(expected))
    return {"maxabs_error": error, "scale": scale, "allowance": ALLOWANCE * scale,
            "agreement": error <= ALLOWANCE * scale}


def value_comparison(actual, expected):
    result = {name: array_comparison(getattr(actual, name), getattr(expected, name))
              for name in ("density", "ionic_gradient")}
    result["agreement"] = all(item["agreement"] for item in result.values())
    return result


def record(checks, name, comparison, **context):
    checks.append({"name": name, **context, **comparison})


def original_state(f, dynamic, state, fixture):
    """Map spectral ket/bra amplitudes to ORIGINAL PW-basis directed matrices.

    No production FFT, projection, contraction, or derivative is used here.
    The source's occupation weight multiplies its outer-product contribution.
    """
    size = fixture.NK * fixture.NPW
    delta = np.zeros((2, size, size), dtype=complex)
    for block in dynamic.state_layout:
        spin, k = divmod(block.channel_index, fixture.NK)
        amplitude = state[block.start:block.stop].reshape(block.shape)
        source = f.vectors[spin, k]
        weight = f.weights[spin, k]
        if block.branch == "plus":
            target = (k + f.step) % fixture.NK
            delta[spin, fixture.sector(target), fixture.sector(k)] += (
                f.vectors[spin, target] @ (amplitude * weight) @ source.conj().T)
        else:
            target = (k - f.step) % fixture.NK
            delta[spin, fixture.sector(k), fixture.sector(target)] += (
                (source * weight) @ amplitude.T @ f.vectors[spin, target].conj().T)
    return delta


def original_observables(f, delta, fixture):
    """Literal Fourier sums and full original-basis trace, independent of FFT maps."""
    density = np.zeros((f.model.fft.nnr, 2), dtype=complex)
    gradient = np.zeros(6, dtype=complex)
    for spin in range(2):
        for k in range(fixture.NK):
            target = (k + f.step) % fixture.NK
            block = delta[spin, fixture.sector(target), fixture.sector(k)]
            density[:fixture.NGRID, spin] += np.einsum(
                "ri,ij,rj->r", f.envelopes[target], block,
                f.envelopes[k].conj()) / f.model.fft.omega
            for a in range(6):
                gradient[a] += np.trace(f.gq[spin, a, k].conj().T @ block)
    return PWResponseValue(density, gradient)


def original_operators(f, potential, u, fixture):
    hs, rhos, forces = [], [], []
    for spin in range(2):
        h = fixture.block_matrix(f.h[spin])
        rho = fixture.block_matrix(np.array([
            (f.vectors[spin, k] * f.weights[spin, k]) @ f.vectors[spin, k].conj().T
            for k in range(fixture.NK)]))
        drive, _ = fixture.full_drive(f, spin, potential, u)
        hs.append(h)
        rhos.append(rho)
        forces.append(-1j * (drive @ rho - rho @ drive))
    return np.array(hs), np.array(rhos), np.array(forces)


def original_derivative(hs, delta, forcing):
    return -1j * (hs @ delta - delta @ hs) + forcing


def state_generator(dynamic, zero_v, zero_u):
    unit = np.eye(dynamic.state_size, dtype=complex)
    return np.column_stack([dynamic.derivative(column, zero_v, zero_u) for column in unit])


def original_generator(h):
    """Column-by-column commutator, not a spectral or Kronecker construction."""
    n = h.shape[0]
    columns = []
    for j in range(n):
        for i in range(n):
            unit = np.zeros_like(h)
            unit[i, j] = 1
            columns.append((-1j * (h @ unit - unit @ h)).reshape(-1, order="F"))
    return np.column_stack(columns)


def driven_evolution(generator, forcing, initial, z, time):
    """Exact finite exponential with a scalar exp(-i*z*t) drive coordinate."""
    size = initial.size
    augmented = np.zeros((size + 1, size + 1), dtype=complex)
    augmented[:size, :size] = generator
    augmented[:size, size] = forcing
    augmented[size, size] = -1j * z
    return (expm(time * augmented) @ np.r_[initial, 1.])[:size]


def driven_resolvent(generator, forcing, z):
    return np.linalg.solve(-1j * z * np.eye(generator.shape[0]) - generator, forcing)


def run_fixture(step, fixture, checks):
    f = fixture.fixture(step)
    dynamic = PWElectronicDynamics(f.model)
    rng = np.random.default_rng(210021 + step)
    potential = 0.31 * fixture.complex_random(rng, (f.model.fft.nnr, 2))
    u = 0.23 * fixture.complex_random(rng, (6,))
    initial = 0.17 * fixture.complex_random(rng, (dynamic.state_size,))
    zero = dynamic.zero_state()
    zero_v, zero_u = np.zeros_like(potential), np.zeros_like(u)
    generator = state_generator(dynamic, zero_v, zero_u)
    original = original_state(f, dynamic, initial, fixture)
    hs, rhos, forcing = original_operators(f, potential, u, fixture)
    context = {"q_step": step}
    record(checks, "packing_and_state_observables", value_comparison(
        dynamic.observables(initial), original_observables(f, original, fixture)), **context)
    record(checks, "original_basis_commutator_derivative", array_comparison(
        original_state(f, dynamic, dynamic.derivative(initial, potential, u), fixture),
        original_derivative(hs, original, forcing)), **context)
    record(checks, "stationary_reference", array_comparison(hs @ rhos - rhos @ hs,
                                                           np.zeros_like(hs)), **context)
    alpha, beta = 0.7 + 0.4j, -0.2 + 0.6j
    second = fixture.complex_random(rng, initial.shape)
    record(checks, "complex_linear_derivative", array_comparison(
        dynamic.derivative(alpha * initial + beta * second, alpha * potential, beta * u),
        alpha * dynamic.derivative(initial, potential, zero_u)
        + beta * dynamic.derivative(second, zero_v, u)), **context)
    left, right = dynamic.observables(initial), dynamic.observables(second)
    record(checks, "complex_linear_observables", value_comparison(
        dynamic.observables(alpha * initial + beta * second), PWResponseValue(
            alpha * left.density + beta * right.density,
            alpha * left.ionic_gradient + beta * right.ionic_gradient)), **context)
    up, down = potential.copy(), potential.copy()
    up[:, 1] = 0
    down[:, 0] = 0
    drives = {"mixed_complex": (potential, u), "ionic_complex": (zero_v, u),
              "local_up": (up, zero_u), "local_down": (down, zero_u)}
    if step == 0:
        constant = np.empty_like(potential)
        constant[:] = [0.37 + 0.19j, -0.21 + 0.43j]
        drives["constant_per_spin_gauge"] = constant, zero_u
    for z in FREQUENCIES:
        frequency = {**context, "z_ry": [z.real, z.imag]}
        for name, (v, displacement) in drives.items():
            force = dynamic.derivative(zero, v, displacement)
            solved = driven_resolvent(generator, force, z)
            actual = dynamic.observables(solved)
            expected, residuals = fixture.direct_eom(f, z, v, displacement)
            record(checks, "drive_resolvent_original_EOM", value_comparison(actual, expected),
                   drive=name, **frequency)
            record(checks, "direct_EOM_residuals", {"values": residuals,
                   "agreement": all(value <= ALLOWANCE for value in residuals.values())},
                   drive=name, **frequency)
            record(checks, "drive_resolvent_settled_R0", value_comparison(
                actual, PWIndependentResponse(f.model, z).apply(v, displacement)),
                drive=name, **frequency)
            if name == "constant_per_spin_gauge":
                record(checks, "constant_per_spin_zero_observable", value_comparison(
                    actual, PWResponseValue(np.zeros_like(actual.density),
                                            np.zeros_like(actual.ionic_gradient))), **frequency)
            if step == 0:
                record(checks, "gamma_per_spin_particle_conservation", array_comparison(
                    f.model.fft.omega / fixture.NGRID * actual.density[:fixture.NGRID].sum(axis=0),
                    np.zeros(2)), drive=name, **frequency)
        force = dynamic.derivative(zero, potential, u)
        original_generators = [original_generator(h) for h in hs]
        for time in TIMES:
            evolved = driven_evolution(generator, force, initial, z, time)
            direct = np.array([
                driven_evolution(op, forcing[spin].reshape(-1, order="F"),
                                 original[spin].reshape(-1, order="F"), z, time)
                .reshape(hs[spin].shape, order="F")
                for spin, op in enumerate(original_generators)])
            record(checks, "time_domain_original_state", array_comparison(
                original_state(f, dynamic, evolved, fixture), direct), tau_ry=time, **frequency)
            record(checks, "time_domain_density_and_gradient", value_comparison(
                dynamic.observables(evolved), original_observables(f, direct, fixture)),
                tau_ry=time, **frequency)
            if step == 0:
                record(checks, "time_domain_spin_trace_conservation", array_comparison(
                    np.trace(direct, axis1=1, axis2=2), np.trace(original, axis1=1, axis2=2)),
                    tau_ry=time, **frequency)
    if step == 1:
        mutation_checks(f, dynamic, generator, potential, u, fixture, checks)
    invalid_calls = {
        "wrong_state_shape": lambda: dynamic.observables(zero[:-1]),
        "nonfinite_state": lambda: dynamic.derivative(zero * np.nan, zero_v, zero_u),
        "wrong_potential_shape": lambda: dynamic.derivative(zero, zero_v[:, 0], zero_u),
        "wrong_displacement_shape": lambda: dynamic.derivative(zero, zero_v, zero_u[:-1]),
    }
    for name, call in invalid_calls.items():
        try:
            call()
        except ValueError:
            rejected = True
        else:
            rejected = False
        record(checks, name, {"agreement": rejected}, **context)
    return {"q_step": step, "state_size": dynamic.state_size,
            "state_layout": [vars(block) for block in dynamic.state_layout],
            "input_scope": "immutable 0019 synthetic finite full EigenKets; not Native"}


def mutation_checks(f, dynamic, generator, potential, u, fixture, checks):
    z = FREQUENCIES[0]
    forcing = dynamic.derivative(dynamic.zero_state(), potential, u)
    reference, _ = fixture.direct_eom(f, z, potential, u)
    minus = np.concatenate([np.arange(b.start, b.stop) for b in dynamic.state_layout
                            if b.branch == "minus"])
    wrong_generator, wrong_force = generator.copy(), forcing.copy()
    wrong_generator[minus] *= -1
    wrong_force[minus] *= -1
    mutants = {"wrong_minus_EOM_sign": dynamic.observables(
        driven_resolvent(wrong_generator, wrong_force, z))}
    wrong_force = forcing.copy()
    # Wrongly conjugate the stored bra's complex drive (undoing the required
    # drive conjugation pair); include its i factor explicitly.
    wrong_force[minus] = 1j * (forcing[minus] / 1j).conj()
    mutants["wrong_minus_drive_conjugation"] = dynamic.observables(
        driven_resolvent(generator, wrong_force, z))
    solved = driven_resolvent(generator, forcing, z)
    wrong_state = solved.copy()
    wrong_state[minus] = wrong_state[minus].conj()
    mutants["wrong_minus_observable_conjugation"] = dynamic.observables(wrong_state)
    wrong_map = replace(f.model, channels=tuple(replace(
        c, plus_slots=c.minus_slots, minus_slots=c.plus_slots) for c in f.model.channels))
    wrong_weights = replace(f.model, channels=tuple(replace(
        c, weighted_occupations=np.ones_like(c.weighted_occupations)) for c in f.model.channels))
    for name, model in (("swapped_plus_minus_FFT_maps", wrong_map),
                        ("missing_source_weights", wrong_weights)):
        mutant = PWElectronicDynamics(model)
        mutants[name] = mutant.observables(driven_resolvent(
            generator, mutant.derivative(mutant.zero_state(), potential, u), z))
    correct = dynamic.observables(solved)
    mutants["missing_inverse_volume"] = PWResponseValue(
        correct.density * f.model.fft.omega, correct.ionic_gradient)
    mutants["extra_spin_factor_two"] = PWResponseValue(
        correct.density * 2, correct.ionic_gradient * 2)
    for name, value in mutants.items():
        comparison = value_comparison(value, reference)
        record(checks, name, {"agreement": not comparison["agreement"],
                             "exposed": not comparison["agreement"], "comparison": comparison},
               kind="intentional_mutation", q_step=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new JSON receipt; never overwrite")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("--output must name a new file")
    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "preregistration_commit": "2cb8638a", "scaled_allowance": ALLOWANCE,
        "equation": "i*d(delta_rho)/d(tau)=[H0,delta_rho]+[v+g*u,rho0]",
        "time_convention": "Ry-energy time, hbar=1, exp(-i*z*tau)",
        "state_normalization": "unweighted spectral plus ket and conjugated minus ket amplitudes; "
                               "weighted occupations enter observables only; no state cutoff",
        "state_packing": "channel order, plus then minus, (target,source) flattened C order",
        "source_bindings": [], "checks": [], "fixtures": [], "errors": [],
        "expected_immutable_fixture_sha256": FIXTURE_SHA256,
        "not_claimed": ["Native full-time comparison", "physical damping or loss", "screening",
                        "ionic inertia", "coupled equilibrium noise", "material coherence domain"],
    }
    try:
        paths = [Path(module.__file__) for module in
                 (dynamics_module, response_module, data_module, space_module, fft_module)]
        paths += [Path(__file__), ATTEMPT / "registration.json", ATTEMPT / "source-inventory.json",
                  FIXTURE_PATH]
        report["source_bindings"] = [binding(path) for path in paths]
        fixture = load_fixture()
        for step in (0, 1):
            try:
                report["fixtures"].append(run_fixture(step, fixture, report["checks"]))
            except Exception as error:
                report["errors"].append({"q_step": step, "type": type(error).__name__,
                                         "message": str(error), "traceback": traceback.format_exc()})
    except Exception as error:
        report["errors"].append({"type": type(error).__name__, "message": str(error),
                                 "traceback": traceback.format_exc()})
    report["failures"] = [check for check in report["checks"] if not check["agreement"]]
    report["agreement"] = (len(report["fixtures"]) == 2 and not report["errors"]
                            and not report["failures"] and bool(report["checks"]))
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"output": str(args.output), "agreement": report["agreement"],
                      "checks": len(report["checks"]), "failures": len(report["failures"]),
                      "errors": len(report["errors"])}))
    return 0 if report["agreement"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

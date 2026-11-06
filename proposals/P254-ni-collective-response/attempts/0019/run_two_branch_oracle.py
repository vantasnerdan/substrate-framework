#!/usr/bin/env python3
"""New tiny mathematical four-k oracle for the isolated full-PW two-branch path.

Run with PYTHONPATH=src and --output NEW.json. This is not Native input and
runs no Native calculation. A full block-H density-matrix Liouvillian solve,
in the original non-diagonal PW basis, is the independent comparator. Only
fixture preparation diagonalizes H; the comparator never uses its spectrum
or the response's resolvents. The finite-q shift is one step on a four-k ring,
so plus and minus targets differ. q=0 is a separate conservation fixture.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import traceback
from typing import Any

import numpy as np

import substrate_framework.pw_independent_response as response_module
import substrate_framework.pw_response_data as data_module
import substrate_framework.qe_response_fft as fft_module
from substrate_framework.pw_independent_response import PWIndependentResponse
from substrate_framework.pw_response_data import PWResponseChannel, PWResponseModel, PWResponseValue
from substrate_framework.qe_response_fft import NativeRhoFFT


ATTEMPT = Path(__file__).resolve().parent
ALLOWANCE = 1e-10
NK, NPW, NGRID = 4, 3, 7
FREQUENCIES = (0.02 + 0.02j, 0.05 + 0.02j)


@dataclass(frozen=True)
class Fixture:
    model: PWResponseModel
    h: np.ndarray                    # spin,k,PW,PW, original non-diagonal basis
    vectors: np.ndarray
    weights: np.ndarray
    envelopes: np.ndarray            # k,physical_grid,PW, literal Fourier sums
    gq: np.ndarray                   # spin,ionic,k,target_PW,source_PW
    step: int


def complex_random(rng, shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


def fixture(step: int) -> Fixture:
    """Complete three-PW sectors; no truncation of any fixture eigenspace."""
    rng = np.random.default_rng(191021 + step)
    slots = np.array([[0, 1, 6], [2, 0, 6], [5, 1, 0], [0, 5, 2]])
    envelopes = np.exp(2j * np.pi * np.arange(NGRID)[None, :, None]
                       * slots[:, None, :] / NGRID)
    geometry = NativeRhoFFT((NGRID, 1, 1), (NGRID + 2, 1, 1), NGRID + 2,
                            np.arange(NGRID), np.zeros((3, NGRID)), 3.7)
    h = np.empty((2, NK, NPW, NPW), dtype=complex)
    vectors = np.empty_like(h)
    energies = np.empty((2, NK, NPW))
    weights = np.empty_like(energies)
    gq = 0.27 * complex_random(rng, (2, 6, NK, NPW, NPW))
    quadrature = np.array([0.13, 0.21, 0.29, 0.37])
    channels = []
    for spin in range(2):
        for k in range(NK):
            raw = 0.22 * complex_random(rng, (NPW, NPW))
            h[spin, k] = raw + raw.conj().T + np.diag(
                np.array([-0.6, 0.7, 1.9]) + 0.13 * k + 0.31 * spin)
            energies[spin, k], vectors[spin, k] = np.linalg.eigh(h[spin, k])
            weights[spin, k] = quadrature[k] / (1 + np.exp((energies[spin, k] - 0.4) / 0.63))
            if step == 0:
                gq[spin, :, k] += gq[spin, :, k].conj().transpose(0, 2, 1)
    for spin in range(2):
        for k in range(NK):
            plus, minus = (k + step) % NK, (k - step) % NK
            source = np.zeros((geometry.nnr, NPW), dtype=complex)
            source[:NGRID] = envelopes[k] @ vectors[spin, k]
            channels.append(PWResponseChannel(
                identity=f"mathematical-spin{spin}-k{k}", spin=spin,
                source_energies=energies[spin, k], weighted_occupations=weights[spin, k],
                source_grid=source,
                plus_energies=energies[spin, plus], plus_vectors=vectors[spin, plus],
                plus_slots=slots[plus],
                minus_energies=energies[spin, minus], minus_vectors=vectors[spin, minus],
                minus_slots=slots[minus],
                g_plus=gq[spin, :, k] @ vectors[spin, k],
                g_minus=gq[spin, :, minus].conj().transpose(0, 2, 1) @ vectors[spin, k]))
    model = PWResponseModel(geometry, np.array([step / NK, 0., 0.]), tuple(channels),
                            {"fixture": "synthetic complete finite four-k ring; NOT Native"})
    return Fixture(model, h, vectors, weights, envelopes, gq, step)


def block_matrix(blocks: np.ndarray) -> np.ndarray:
    result = np.zeros((NK * NPW, NK * NPW), dtype=complex)
    for k, block in enumerate(blocks):
        sector = slice(k * NPW, (k + 1) * NPW)
        result[sector, sector] = block
    return result


def sector(k: int) -> slice:
    return slice(k * NPW, (k + 1) * NPW)


def full_drive(f: Fixture, spin: int, potential: np.ndarray, u: np.ndarray):
    """Assemble an actual directed q operator and its paired ionic observables."""
    size = NK * NPW
    drive = np.zeros((size, size), dtype=complex)
    ionic = np.zeros((6, size, size), dtype=complex)
    for k in range(NK):
        target = (k + f.step) % NK
        local = f.envelopes[target].conj().T @ (
            potential[:NGRID, spin, None] * f.envelopes[k]) / NGRID
        ionic[:, sector(target), sector(k)] = f.gq[spin, :, k]
        drive[sector(target), sector(k)] = local + np.einsum("a,aij->ij", u, f.gq[spin, :, k])
    return drive, ionic


def liouvillian(h: np.ndarray, z: complex) -> np.ndarray:
    """Column-by-column EOM action; no Kronecker/vec convention is assumed."""
    n = h.shape[0]
    operator = np.empty((n * n, n * n), dtype=complex)
    for j in range(n):
        for i in range(n):
            unit = np.zeros((n, n), dtype=complex)
            unit[i, j] = 1
            operator[:, i + n * j] = (z * unit - h @ unit + unit @ h).reshape(-1, order="F")
    return operator


def direct_eom(f: Fixture, z: complex, potential: np.ndarray, u: np.ndarray):
    """Solve z*dR-[H,dR]=[V,R0] directly in the full original block basis."""
    density = np.zeros((f.model.fft.nnr, 2), dtype=complex)
    gradient = np.zeros(6, dtype=complex)
    residuals = []
    stationary = []
    traces = []
    for spin in range(2):
        h = block_matrix(f.h[spin])
        rho = block_matrix(np.array([
            (f.vectors[spin, k] * f.weights[spin, k]) @ f.vectors[spin, k].conj().T
            for k in range(NK)]))
        drive, ionic = full_drive(f, spin, potential, u)
        rhs = drive @ rho - rho @ drive
        solution = np.linalg.solve(liouvillian(h, z), rhs.reshape(-1, order="F"))
        delta = solution.reshape(h.shape, order="F")
        residuals.append(maxabs(z * delta - h @ delta + delta @ h - rhs) / max(1., maxabs(rhs)))
        stationary.append(maxabs(h @ rho - rho @ h))
        traces.append(abs(np.trace(delta)))
        for k in range(NK):
            target = (k + f.step) % NK
            density[:NGRID, spin] += np.einsum(
                "ri,ij,rj->r", f.envelopes[target], delta[sector(target), sector(k)],
                f.envelopes[k].conj()) / f.model.fft.omega
        for a in range(6):
            gradient[a] += np.trace(ionic[a].conj().T @ delta)
    return PWResponseValue(density, gradient), {
        "eom_scaled_residual": max(residuals), "stationary_reference_maxabs": max(stationary),
        "trace_response_maxabs": max(traces)}


def plus_only_mutant(f: Fixture, z: complex, potential: np.ndarray, u: np.ndarray):
    """Intentional one-branch defect, using small original-basis ket solves."""
    density = np.zeros((f.model.fft.nnr, 2), dtype=complex)
    gradient = np.zeros(6, dtype=complex)
    for spin in range(2):
        for k in range(NK):
            target = (k + f.step) % NK
            source = f.vectors[spin, k]
            g = f.gq[spin, :, k] @ source
            local = f.envelopes[target].conj().T @ (
                potential[:NGRID, spin, None] * (f.envelopes[k] @ source)) / NGRID
            drive = local + np.einsum("a,aij->ij", u, g)
            ket = np.column_stack([
                np.linalg.solve((f.model.channels[spin * NK + k].source_energies[n] + z)
                                * np.eye(NPW) - f.h[spin, target], drive[:, n])
                for n in range(NPW)])
            density[:NGRID, spin] += np.einsum(
                "rn,rn,n->r", (f.envelopes[k] @ source).conj(),
                f.envelopes[target] @ ket, f.weights[spin, k]) / f.model.fft.omega
            gradient += np.einsum("atn,tn,n->a", g.conj(), ket, f.weights[spin, k])
    return PWResponseValue(density, gradient)


def maxabs(array) -> float:
    return float(np.max(np.abs(array), initial=0.))


def compare(actual: PWResponseValue, expected: PWResponseValue) -> dict[str, Any]:
    details = {}
    for name in ("density", "ionic_gradient"):
        reference = getattr(expected, name)
        scale = max(1., maxabs(reference))
        error = maxabs(getattr(actual, name) - reference)
        details[name] = {"maxabs_error": error, "scale": scale,
                         "allowance": ALLOWANCE * scale, "agreement": error <= ALLOWANCE * scale}
    details["agreement"] = all(details[name]["agreement"] for name in ("density", "ionic_gradient"))
    return details


def mutation_records(f, z, potential, u, correct, direct):
    wrong_pair = replace(f.model, channels=tuple(
        replace(c, g_minus=c.g_plus) for c in f.model.channels))
    wrong_map = replace(f.model, channels=tuple(
        replace(c, plus_slots=c.minus_slots, minus_slots=c.plus_slots) for c in f.model.channels))
    mutants = {
        "plus_branch_only": plus_only_mutant(f, z, potential, u),
        "wrong_minus_operator_pairing": PWIndependentResponse(wrong_pair, z).apply(potential, u),
        "swapped_plus_minus_fft_maps": PWIndependentResponse(wrong_map, z).apply(potential, u),
        "opposite_physical_sign": PWResponseValue(-correct.density, -correct.ionic_gradient),
        "extra_spin_factor_two": PWResponseValue(2 * correct.density, 2 * correct.ionic_gradient),
        "missing_inverse_volume": PWResponseValue(f.model.fft.omega * correct.density,
                                                   correct.ionic_gradient),
        "swapped_spin_density": PWResponseValue(correct.density[:, ::-1], correct.ionic_gradient),
    }
    records = {}
    for name, value in mutants.items():
        comparison = compare(value, direct)
        records[name] = {"exposed": not comparison["agreement"], "comparison": comparison}
    return records


def invalid_input_records(model, z, potential, u):
    cases = {
        "real_frequency": lambda: PWIndependentResponse(model, z.real),
        "lower_half_plane": lambda: PWIndependentResponse(model, z.conjugate()),
        "nonfinite_frequency": lambda: PWIndependentResponse(model, complex(0., float("inf"))),
        "wrong_potential_shape": lambda: PWIndependentResponse(model, z).apply(potential[:, 0], u),
        "wrong_displacement_shape": lambda: PWIndependentResponse(model, z).apply(potential, u[:5]),
        "nonfinite_potential": lambda: PWIndependentResponse(model, z).apply(potential * np.nan, u),
    }
    records = {}
    for name, call in cases.items():
        try:
            call()
        except ValueError:
            records[name] = True
        else:
            records[name] = False
    return records


def run_fixture(step: int) -> dict[str, Any]:
    f = fixture(step)
    rng = np.random.default_rng(5519 + step)
    potential = 0.31 * complex_random(rng, (f.model.fft.nnr, 2))
    u = 0.23 * complex_random(rng, (6,))
    up = potential.copy()
    up[:, 1] = 0
    down = potential.copy()
    down[:, 0] = 0
    zero_v, zero_u = np.zeros_like(potential), np.zeros_like(u)
    drives = {"mixed_complex": (potential, u), "local_up_only": (up, zero_u),
              "local_down_only": (down, zero_u), "ionic_complex_only": (zero_v, u)}
    if step == 0:
        constant = np.empty_like(potential)
        constant[:] = [0.37 + 0.19j, -0.21 + 0.43j]
        drives["constant_per_spin"] = constant, zero_u
    cases = []
    mutations = {}
    for frequency_index, z in enumerate(FREQUENCIES):
        response = PWIndependentResponse(f.model, z)
        for name, (v, displacement) in drives.items():
            actual = response.apply(v, displacement)
            expected, eom = direct_eom(f, z, v, displacement)
            comparison = compare(actual, expected)
            record = {"drive": name, "z_ry": [z.real, z.imag],
                      "comparison": comparison, "direct_EOM": eom}
            record["agreement"] = comparison["agreement"] and all(
                value <= ALLOWANCE for value in eom.values())
            if step == 0:
                integral = f.model.fft.omega / NGRID * np.sum(actual.density[:NGRID], axis=0)
                record["per_spin_particle_response_maxabs"] = maxabs(integral)
                record["agreement"] &= maxabs(integral) <= ALLOWANCE
            if name == "constant_per_spin":
                record["zero_response_maxabs"] = max(maxabs(actual.density), maxabs(actual.ionic_gradient))
                record["agreement"] &= record["zero_response_maxabs"] <= ALLOWANCE
            if name.startswith("local_"):
                inactive = 1 if name == "local_up_only" else 0
                record["inactive_spin_density_maxabs"] = maxabs(actual.density[:, inactive])
                record["agreement"] &= record["inactive_spin_density_maxabs"] <= ALLOWANCE
            cases.append(record)
            if step == 1 and frequency_index == 0 and name == "mixed_complex":
                mutations = mutation_records(f, z, v, displacement, actual, expected)
        alpha, beta = 0.7 + 0.4j, -0.2 + 0.6j
        left = response.apply(alpha * potential, alpha * u)
        right = response.apply(beta * up, zero_u)
        combined = response.apply(alpha * potential + beta * up, alpha * u)
        linearity = compare(combined, PWResponseValue(left.density + right.density,
                                                     left.ionic_gradient + right.ionic_gradient))
        cases.append({"drive": "complex_linearity", "z_ry": [z.real, z.imag],
                      "comparison": linearity, "agreement": linearity["agreement"]})
    targets_distinct = float(min(
        np.linalg.norm(f.h[spin, (k + step) % NK] - f.h[spin, (k - step) % NK])
        for spin in range(2) for k in range(NK)))
    invalid = invalid_input_records(f.model, FREQUENCIES[0], potential, u)
    return {"fixture": "synthetic finite four-k ring; NOT Native", "q_step": step,
            "k_sectors": NK, "full_PW_dimension_per_sector": NPW,
            "grid_physical": list(f.model.fft.physical), "grid_storage": list(f.model.fft.storage),
            "omega_bohr3": f.model.fft.omega, "explicit_spins": 2,
            "plus_minus_target_H_min_difference": float(targets_distinct),
            "finite_q_distinct_targets": step == 1 and targets_distinct > ALLOWANCE,
            "cases": cases, "mutations": mutations, "invalid_inputs_rejected": invalid,
            "agreement": all(case["agreement"] for case in cases) and all(invalid.values())
                         and all(record["exposed"] for record in mutations.values())
                         and (step == 0 or targets_distinct > ALLOWANCE)}


def binding(path: Path) -> dict[str, str]:
    data = path.read_bytes()
    return {"path": str(path), "bytes": str(len(data)), "sha256": hashlib.sha256(data).hexdigest()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new JSON receipt; never overwrite")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("--output must name a new file")
    report: dict[str, Any] = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "preregistration_commit": "f78b9cad", "scaled_allowance": ALLOWANCE,
        "input_scope": "new synthetic finite systems only; no Native execution or material inference",
        "independent_equation": "z*dR-[H,dR]=[V,R0], full original-basis block-H linear solve",
        "observable_convention": "physical response to +v+g*u; delta<g_-q>, not negative force",
        "source_bindings": [binding(Path(module.__file__).resolve()) for module in
                            (response_module, data_module, fft_module)]
                           + [binding(Path(__file__).resolve()), binding(ATTEMPT / "registration.json"),
                              binding(ATTEMPT / "source-inventory.json")],
        "not_claimed": ["actual Native frame consumption", "physical eta or loss", "screening",
                        "coherence-domain solution", "static thermal population response"],
    }
    try:
        report["finite_q"] = run_fixture(1)
        report["q0_conservation"] = run_fixture(0)
        report["agreement"] = report["finite_q"]["agreement"] and report["q0_conservation"]["agreement"]
    except Exception as error:
        report["agreement"] = False
        report["error"] = {"type": type(error).__name__, "message": str(error),
                           "traceback": traceback.format_exc()}
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"output": str(args.output), "agreement": bool(report["agreement"])}))
    return 0 if report["agreement"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

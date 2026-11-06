#!/usr/bin/env python3
"""Exact finite real-space/BvK/Fock witnesses for prepared bilinear noise.

Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python proposals/P254-ni-collective-response/attempts/0023/run_quantum_noise_oracle.py
--output PATH. No physical material, bath, interacting-equilibrium or Gaussian
higher-moment claim follows from these finite fixtures. Leaves execute before
composed production checks. No output is a pass until this program executes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import sympy as sp

from substrate_framework.finite_fermion_state import (
    annihilators, fermi_dirac, fock_gibbs, second_quantize,
)
from substrate_framework.quasifree_spectral_noise import (
    ordered_transition_correlation, paired_correlations,
)


TOLERANCE = 2e-12


def exact(value):
    return sp.simplify(sp.expand_complex(value))


def clean(matrix):
    return sp.ImmutableMatrix(matrix.applyfunc(exact))


def require(condition, label, receipt):
    if not condition:
        raise AssertionError(label)
    receipt["checks"].append(label)


def zero_matrix(matrix):
    return all(exact(value) == 0 for value in matrix)


def numeric(matrices):
    return np.asarray([np.asarray(matrix.evalf(17), dtype=complex)
                       for matrix in matrices])


def close_to(actual, expected, label, receipt):
    error = float(np.max(np.abs(actual - expected), initial=0.0))
    require(error <= TOLERANCE, label, receipt)
    receipt["errors"][label] = error


def bind_sources(receipt):
    attempt = Path(__file__).resolve().parent
    inventory = json.loads((attempt / "source-inventory.json").read_text())
    recorded = {Path(item["path"]).name: item for item in inventory["files"]}
    paths = [ROOT / "src/substrate_framework" / name for name in (
        "finite_fermion_state.py", "finite_fermion_response.py",
        "quasifree_spectral_noise.py")]
    paths += [Path(__file__).resolve(), attempt / "registration.json",
              attempt / "source-inventory.json"]
    for path in paths:
        data = path.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        receipt["source_bindings"].append({"path": str(path.relative_to(ROOT)),
                                           "bytes": len(data), "sha256": digest})
        if path.name in ("finite_fermion_state.py", "finite_fermion_response.py"):
            require(digest == recorded[path.name]["sha256"],
                    "settled source bytes: " + path.name, receipt)


def fixture_leaves(receipt):
    """Construct observables in real space before any momentum contraction."""
    cells = 3
    omega = (-1 + sp.I * sp.sqrt(3)) / 2
    fourier = clean(sp.Matrix(cells, cells,
                             lambda r, k: omega ** (r * k) / sp.sqrt(cells)))
    require(zero_matrix(fourier.H * fourier - sp.eye(cells)),
            "unitary three-cell Fourier transform", receipt)
    c = annihilators(cells)
    real_c = tuple(clean(sum((fourier[r, k] * c[k] for k in range(cells)),
                            sp.zeros(8))) for r in range(cells))
    for r in range(cells):
        for s in range(cells):
            require(zero_matrix(real_c[r] * real_c[s].H
                                + real_c[s].H * real_c[r]
                                - int(r == s) * sp.eye(8)),
                    f"real-space CAR {r},{s}", receipt)
            require(zero_matrix(real_c[r] * real_c[s]
                                + real_c[s] * real_c[r]),
                    f"real-space annihilator CAR {r},{s}", receipt)
    energies = (sp.S.Zero, sp.S.Zero, sp.S.One)
    beta = sp.log(2)  # Auxiliary fixture parameter, not a material temperature.
    rho, partition, many_energies = fock_gibbs(energies, beta)
    occupations = fermi_dirac(energies, beta)
    require(occupations == (sp.Rational(1, 2), sp.Rational(1, 2),
                            sp.Rational(1, 3)), "raw FD probabilities", receipt)
    require(sp.trace(rho) == 1 and partition == 6,
            "normalized full eight-state grandcanonical trace", receipt)
    for k in range(cells):
        require(exact(sp.trace(rho * c[k].H * c[k])) == occupations[k],
                f"Fock marginal equals raw FD {k}", receipt)

    def real_observables(q):
        density = sp.zeros(cells)
        hopping = sp.zeros(cells)
        for r in range(cells):
            density[r, r] = omega ** (q * r) / sp.sqrt(cells)
            hopping[(r + 1) % cells, r] = omega ** (q * r) / sp.sqrt(cells)
        matrices = (clean(density), clean((1 + sp.I) * density
                                         + (2 - sp.I) * hopping / 3))
        fock = []
        for matrix in matrices:
            # Direct real-space bilinear construction, NOT momentum Wick weights.
            observable = sp.zeros(8)
            for r in range(cells):
                for s in range(cells):
                    observable += matrix[r, s] * real_c[r].H * real_c[s]
            fock.append(clean(observable))
        momentum = tuple(clean(fourier.H * matrix * fourier) for matrix in matrices)
        for index, (one_body, full_fock) in enumerate(zip(momentum, fock)):
            require(zero_matrix(second_quantize(one_body) - full_fock),
                    f"real-space to momentum Fock map q={q},a={index}", receipt)
            for target in range(cells):
                for source in range(cells):
                    if target != (source + q) % cells:
                        require(exact(one_body[target, source]) == 0,
                                f"directed momentum support q={q},a={index},{target},{source}",
                                receipt)
        return momentum, tuple(fock)

    plus, fock_plus = real_observables(1)
    gamma, fock_gamma = real_observables(0)
    minus = tuple(clean(matrix.H) for matrix in plus)
    fock_minus = tuple(clean(matrix.H) for matrix in fock_plus)
    require(not zero_matrix(plus[0] - minus[0]),
            "q and minus-q genuinely distinct", receipt)
    require(exact(plus[0][1, 0]) != 0 and energies[0] == energies[1],
            "non-diagonal zero-gap coherence present", receipt)
    for source in range(cells):
        require(exact(plus[0][(source + 1) % cells, source]) == 1 / sp.sqrt(cells),
                f"explicit inverse-square-root cell normalization {source}", receipt)
    total_number = sum((op.H * op for op in c), sp.zeros(8))
    require(zero_matrix(fock_gamma[0] - total_number / sp.sqrt(cells)),
            "Gamma density is total-number over square-root cells", receipt)
    receipt["fixture"] = {
        "cells": cells, "explicit_spin_channels": 1, "fock_states": 8,
        "energies": list(map(str, energies)), "auxiliary_beta": str(beta),
        "raw_probabilities": list(map(str, occupations)),
        "normalization": "Real-space operators contain 1/sqrt(3); production g omits it and weight=1/3",
        "real_space_observables": "D_q=sum_R exp(i q R)n_R/sqrt(3); second=(1+i)D_q+(2-i)sum_R exp(i q R)c_(R+1)^dag c_R/(3sqrt(3)); minus=plus^dag",
    }
    return energies, occupations, rho, many_energies, plus, minus, gamma, fock_plus, fock_minus, fock_gamma


def connected_fock(left, right, time_value, rho, many_energies, reverse=False):
    """Exact connected trace over ALL Fock states, with explicit H0 evolution."""
    total = sp.S.Zero
    for i in range(rho.rows):
        for j in range(rho.rows):
            if reverse:
                term = right[i, j] * left[j, i]
                gap = many_energies[j] - many_energies[i]
            else:
                term = left[i, j] * right[j, i]
                gap = many_energies[i] - many_energies[j]
            if term != 0:
                total += rho[i, i] * term * sp.exp(sp.I * gap * time_value)
    mean_left = sp.trace(rho * left)
    mean_right = sp.trace(rho * right)
    return exact(total - mean_left * mean_right)


def mutation_witness(receipt, name, mutated, reference):
    difference = np.abs(mutated - reference)
    index = np.unravel_index(int(np.argmax(difference)), difference.shape)
    distance = float(difference[index])
    require(distance > 1e-8, "reject mutation: " + name, receipt)
    receipt["mutation_witnesses"][name] = {
        "index": list(map(int, index)), "distance": distance,
        "mutated": [float(mutated[index].real), float(mutated[index].imag)],
        "full_fock_reference": [float(reference[index].real), float(reference[index].imag)],
    }


def composed_checks(fixture, receipt):
    energies, occupations, rho, many, plus, minus, gamma, fp, fm, fg = fixture
    exact_times = (sp.S.Zero, sp.pi / 2, -sp.pi / 2)
    times = np.asarray([float(t) for t in exact_times])
    e = np.asarray(energies, dtype=float)
    f = np.asarray(occupations, dtype=float)
    # Undo ONLY the operator normalization established by the real-space leaves.
    ops = {name: numeric([clean(sp.sqrt(3) * m) for m in matrices])
           for name, matrices in (("plus", plus), ("minus", minus), ("gamma", gamma))}

    def traces(observables, reverse=False):
        return tuple(sp.ImmutableMatrix(2, 2, lambda a, b:
                     connected_fock(observables[a].H, observables[b], t, rho,
                                    many, reverse)) for t in exact_times)

    expected_exact = {name: traces(observables) for name, observables in
                      (("plus", fp), ("minus", fm), ("gamma", fg))}
    expected = {name: numeric(matrices) for name, matrices in expected_exact.items()}
    receipt["exact_ordered_traces"] = {
        name: [[[str(value) for value in row] for row in matrix.tolist()]
               for matrix in matrices] for name, matrices in expected_exact.items()}
    receipt["times"] = list(map(str, exact_times))

    def production(name="plus", **overrides):
        kwargs = dict(source_energies=e, target_energies=e,
                      source_occupations=f, target_occupations=f,
                      operators=ops[name], times=times, weight=1 / 3)
        kwargs.update(overrides)
        return ordered_transition_correlation(**kwargs)

    actual = {name: production(name) for name in ops}
    for name in ops:
        close_to(actual[name], expected[name], "full-Fock ordered " + name, receipt)
    # Direct real-space Fock number variance: f(1-f), never a fixed-exact-N state.
    require(expected_exact["gamma"][0][0, 0] == sp.Rational(13, 54),
            "Gamma population variance survives exactly", receipt)
    require(any(exact(sp.im(value)) != 0 for value in expected_exact["plus"][0]),
            "complex cross-channel witness nonzero", receipt)

    ordered_checks(production, actual, expected, e, f, ops, times, receipt)
    paired_checks(production, actual, expected, numeric(traces(fp, reverse=True)),
                  times, receipt)


def ordered_checks(production, actual, expected, e, f, ops, times, receipt):
    rectangular = sum((production(source_energies=e[selection],
                                  source_occupations=f[selection],
                                  operators=ops["plus"][:, :, selection])
                       for selection in (slice(0, 2), slice(2, 3))),
                      np.zeros_like(actual["plus"]))
    close_to(rectangular, expected["plus"], "rectangular directed source-sector sum", receipt)
    # Zero-operator spectator padding makes the kernel cross its 4096 boundary.
    size = 65
    locations = (0, 32, 64)
    padded = np.zeros((2, size, size), dtype=complex)
    pe = np.zeros(size)
    pf = np.full(size, 0.5)
    for i, target in enumerate(locations):
        pe[target], pf[target] = e[i], f[i]
        for j, source in enumerate(locations):
            padded[:, target, source] = ops["plus"][:, i, j]
    close_to(production(source_energies=pe, target_energies=pe,
                        source_occupations=pf, target_occupations=pf,
                        operators=padded), expected["plus"],
             "block-boundary contraction with invisible spectators", receipt)


    mutation_witness(receipt, "missing BvK weight", production(weight=1), expected["plus"])
    mutation_witness(receipt, "spurious explicit-spin factor two", production(weight=2 / 3),
                     expected["plus"])
    mutation_witness(receipt, "weighted occupations as Bernoulli",
                     production(source_occupations=f / 3, target_occupations=f / 3),
                     expected["plus"])
    mutation_witness(receipt, "missing Pauli holes", production(target_occupations=np.zeros(3)),
                     expected["plus"])
    equal_gap = e[:, None] == e[None, :]
    for name, mask in (("all zero gaps discarded", equal_gap),
                       ("non-diagonal degenerate coherence discarded",
                        equal_gap & ~np.eye(3, dtype=bool))):
        changed = ops["plus"].copy()
        changed[:, mask] = 0
        mutation_witness(receipt, name, production(operators=changed), expected["plus"])
    changed = ops["gamma"].copy()
    changed[:, np.eye(3, dtype=bool)] = 0
    mutation_witness(receipt, "Gamma population fluctuations discarded",
                     production("gamma", operators=changed), expected["gamma"])
    mutation_witness(receipt, "ordered evolution time sign", production(times=-times),
                     expected["plus"])
    mutation_witness(receipt, "ordered operator conjugation reversed",
                     production(operators=ops["plus"].conj()), expected["plus"])

    invalid = [dict(source_occupations=[-0.1, 0.5, 0.5]),
               dict(target_occupations=[0.5, 1.1, 0.5]),
               dict(source_energies=[0, 1]), dict(times=[np.nan]),
               dict(weight=-1), dict(weight=np.inf),
               dict(operators=np.full((2, 3, 3), complex(np.nan, 0))),
               dict(operators=np.zeros((2, 3, 2)))]
    for index, change in enumerate(invalid):
        try:
            production(**change)
        except ValueError:
            receipt["checks"].append(f"invalid input rejected {index}")
        else:
            raise AssertionError(f"invalid input accepted {index}")


def paired_checks(production, actual, expected, reverse, times, receipt):
    opposite_negative = production("minus", times=-times)
    close_to(opposite_negative.swapaxes(1, 2), reverse,
             "opposite-sector negative-time TRANSPOSE equals reverse Fock trace", receipt)
    paired = paired_correlations(actual["plus"], opposite_negative)
    close_to(paired["symmetrized"], (expected["plus"] + reverse) / 2,
             "full-Fock symmetrized", receipt)
    close_to(paired["commutator"], expected["plus"] - reverse,
             "full-Fock commutator", receipt)
    mutation_witness(receipt, "paired adjoint instead of transpose",
                     (actual["plus"] + opposite_negative.conj().swapaxes(1, 2)) / 2,
                     (expected["plus"] + reverse) / 2)
    mutation_witness(receipt, "paired opposite time not reversed",
                     paired_correlations(actual["plus"], actual["minus"])["symmetrized"],
                     (expected["plus"] + reverse) / 2)
    mutation_witness(receipt, "paired transpose omitted",
                     (actual["plus"] + opposite_negative) / 2,
                     (expected["plus"] + reverse) / 2)
    for bad in (np.zeros((3, 2)), np.full((3, 2, 2), np.nan)):
        try:
            paired_correlations(actual["plus"], bad)
        except ValueError:
            receipt["checks"].append("invalid paired input rejected")
        else:
            raise AssertionError("invalid paired input accepted")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    receipt = {"status": "running", "checks": [], "errors": {},
               "mutation_witnesses": {}, "source_bindings": [],
               "scope": "Exact finite independent-H0 grandcanonical fixtures corroborate this construction, not a general interacting-equilibrium theorem, physical temperature, irreversible bath, Gaussian higher moments, or nuclear/coherence verdict.",
               "numeric_tolerance": TOLERANCE}
    try:
        bind_sources(receipt)
        fixture = fixture_leaves(receipt)
        composed_checks(fixture, receipt)
        receipt["status"] = "passed"
    except Exception as error:
        receipt["status"] = "failed"
        receipt["error"] = f"{type(error).__name__}: {error}"
        raise
    finally:
        receipt["elapsed_seconds"] = time.monotonic() - started
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": receipt["status"], "checks": len(receipt["checks"]),
                      "output": str(args.output)}))


if __name__ == "__main__":
    main()

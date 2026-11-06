#!/usr/bin/env python3
"""Independent finite real-space proof for the new PW preparation adapter only.

No Native/material calculation, frame reconstruction, or peer correlation module
is used. Plane waves and finite BvK matrix elements are summed explicitly in real
space; production uses FFT source projection. Run with --output receipt.json.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import replace
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "src"))
from substrate_framework.pw_prepared_noise import prepare_pw_noise
from substrate_framework.pw_response_data import PWResponseChannel, PWResponseModel
from substrate_framework.qe_response_fft import NativeRhoFFT


def _waves(points, miller, coefficients):
    """Direct periodic plane-wave sums, with no FFT/projector implementation."""
    values = np.zeros((len(points), coefficients.shape[1]), complex)
    for r, point in enumerate(points):
        for g, reciprocal in enumerate(miller):
            values[r] += np.exp(2j * np.pi * np.dot(reciprocal, point)) * coefficients[g]
    return values


def _unitary(rng, size):
    values = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    return np.linalg.qr(values)[0]


def _fixture():
    physical, storage = (5, 4, 3), (6, 5, 3)
    sites = list(product(*(range(n) for n in physical)))
    points = np.asarray(sites) / np.asarray(physical)
    nnr = math.prod(storage)
    fft = NativeRhoFFT(physical, storage, nnr, np.array([], dtype=int),
                       np.empty((3, 0)), 13.7)
    rng = np.random.default_rng(2540023)
    source_miller = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    miller = {"plus": np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0], [1, 0, 0]]),
              "minus": np.array([[-1, 0, 0], [1, 0, 0], [0, 0, 0]])}
    source_coefficients = _unitary(rng, 3)
    periodic_source = _waves(points, source_miller, source_coefficients)
    source_grid = np.zeros((nnr, 3), complex)
    for site, values in zip(sites, periodic_source):
        source_grid[np.ravel_multi_index(site, storage, order="F")] = values
    vectors = {name: _unitary(rng, len(gs)) for name, gs in miller.items()}
    slots = {name: np.array([np.ravel_multi_index(tuple(g % physical), storage, order="F")
                             for g in gs]) for name, gs in miller.items()}
    target_grid = {name: _waves(points, miller[name], vectors[name]) for name in miller}
    energies = np.array([0., 1., 744.])
    # Last raw probability survives, but weighting by 1/8 rounds it to zero.
    raw = np.array([0.5, 1 / (1 + math.exp(1)), math.exp(-744.)])
    target_energies = {"plus": np.array([0., 0.4, 2., 744.]),
                       "minus": np.array([0., -0.6, 744.])}
    basis = np.array([[1.2, 0.1, -0.2], [0., 0.9, 0.3], [0.2, 0., 1.1]])
    records, channels = [], []
    for spin in (0, 1):
        for ik, index in enumerate(product(range(2), repeat=3)):
            identity = f"spin{spin + 1}:ik{ik + 1}"
            # Nontrivial reciprocal representatives, deliberately not sorted.
            fractional = np.array(index) / 2 - np.array([1, 0, 1])
            records.append({"identity": identity, "native_spin": str(spin + 1),
                            "native_ik": str(ik + 1), "k_source_2pi_alat": (basis @ fractional).tolist(),
                            "quadrature_weight": 0.125})
            gradients = {name: rng.normal(size=(6, len(gs), 3)) + 1j * rng.normal(size=(6, len(gs), 3))
                         for name, gs in miller.items()}
            channels.append(PWResponseChannel(
                identity, spin, energies, raw / 8, source_grid,
                target_energies["plus"], vectors["plus"], slots["plus"],
                target_energies["minus"], vectors["minus"], slots["minus"],
                gradients["plus"], gradients["minus"]))
    metadata = {"source_channels": records[::-1], "common_mu_Ry": 0., "actual_degauss_Ry": 1.,
                "branch": "finite independent oracle"}
    model = PWResponseModel(fft, basis @ np.array([0.5, 0., 0.]), tuple(channels), metadata)
    return model, basis, points, periodic_source, target_grid, raw


def _cell_overlap(source, target):
    result = np.zeros((target.shape[1], source.shape[1]), complex)
    for t in range(target.shape[1]):
        for s in range(source.shape[1]):
            result[t, s] = sum(target[r, t].conjugate() * source[r, s]
                               for r in range(len(source))) / len(source)
    return result


def _bvk_overlap(points, source, target, q, target_shift=None):
    """Integrate normalized Bloch states and exp(i q.r)/sqrt(8) over 8 cells."""
    q = np.asarray(q)
    target_k = q if target_shift is None else q + np.asarray(target_shift)
    cells = list(product(range(2), repeat=3))
    volume = 13.7
    result = np.zeros((target.shape[1], source.shape[1]), complex)
    for cell in cells:
        for r, point in enumerate(points):
            position = np.asarray(cell) + point
            for t in range(target.shape[1]):
                psi_target = target[r, t] * np.exp(2j * np.pi * np.dot(target_k, position)) / math.sqrt(8 * volume)
                for s in range(source.shape[1]):
                    psi_source = source[r, s] / math.sqrt(8 * volume)
                    probe = np.exp(2j * np.pi * np.dot(q, position)) / math.sqrt(8)
                    result[t, s] += psi_target.conjugate() * probe * psi_source * volume / len(points)
    return result


def _expect_rejected(model, basis, label):
    try:
        prepare_pw_noise(model, basis)
    except ValueError as error:
        return {"case": label, "error": str(error)}
    raise AssertionError(f"accepted invalid preparation: {label}")


def run():
    model, basis, points, source, targets, expected_raw = _fixture()
    branches, metadata = prepare_pw_noise(model, basis)
    errors = {}
    for name, sign in (("plus", 1), ("minus", -1)):
        expected = _cell_overlap(source, targets[name])
        real_space = _bvk_overlap(points, source, targets[name], [sign * 0.5, 0., 0.])
        errors[f"{name}_BvK_matrix"] = float(np.max(np.abs(real_space - expected / math.sqrt(8))))
        errors[f"{name}_forbidden_momentum"] = float(np.max(np.abs(
            _bvk_overlap(points, source, targets[name], [sign * 0.5, 0., 0.], [0., 0.5, 0.]))))
        for channel, block in zip(model.channels, branches[name]):
            assert block["identity"] == channel.identity
            assert block["weight"] == 0.125
            assert np.array_equal(block["operators"][:6], getattr(channel, f"g_{name}"))
            assert np.count_nonzero(block["operators"][7 - channel.spin]) == 0
            np.testing.assert_allclose(block["operators"][6 + channel.spin], expected, rtol=0, atol=2e-14)
            np.testing.assert_allclose(block["source_occupations"][:2], expected_raw[:2], rtol=0, atol=1e-16)
            assert block["source_occupations"][2] == expected_raw[2] > 0
            assert channel.weighted_occupations[2] == 0
            assert block["source_occupations"][0] != channel.weighted_occupations[0]
            for energy, probability in zip(block["target_energies"], block["target_occupations"]):
                expected_probability = math.exp(-energy) / (1 + math.exp(-energy))
                np.testing.assert_allclose(probability, expected_probability, rtol=2e-15, atol=0)
            np.testing.assert_allclose(np.abs(real_space) ** 2,
                                       block["weight"] * np.abs(block["operators"][6 + channel.spin]) ** 2,
                                       rtol=0, atol=2e-14)
        errors[f"{name}_projection"] = float(np.max(np.abs(branches[name][0]["operators"][6] - expected)))
    assert max(errors.values()) < 2e-14
    assert metadata["Ncells"] == 8 and len(metadata["source_channels"]) == 16
    # Gamma must remain admissible without deleting exactly zero gaps.
    gamma, _ = prepare_pw_noise(replace(model, q_cart_2pi_alat=np.zeros(3)), basis)
    assert gamma["plus"][0]["target_energies"][0] == gamma["plus"][0]["source_energies"][0]
    assert gamma["plus"][0]["operators"].shape == (8, 4, 3)
    block = branches["plus"][0]
    good_variance = float(np.sum(block["weight"] * np.abs(block["operators"][6]) ** 2
                                * block["source_occupations"][None, :]
                                * (1 - block["target_occupations"][:, None])))
    wrong_variance = float(np.sum(block["weight"] * np.abs(block["operators"][6]) ** 2
                                 * model.channels[0].weighted_occupations[None, :]
                                 * (1 - block["target_occupations"][:, None])))
    assert good_variance > 0 and abs(good_variance - wrong_variance) > 1e-4
    # Missing mesh cannot be rescued by renormalizing remaining weights to sum one.
    incomplete_metadata = deepcopy(model.metadata)
    missing = model.channels[0].identity
    incomplete_metadata["source_channels"] = [r for r in incomplete_metadata["source_channels"] if r["identity"] != missing]
    incomplete_uniform = replace(model, channels=model.channels[1:], metadata=deepcopy(incomplete_metadata))
    for record in incomplete_metadata["source_channels"]:
        if record["native_spin"] == "1":
            record["quadrature_weight"] = 1 / 7
    rejected = [_expect_rejected(incomplete_uniform, basis, "incomplete mesh with unchanged uniform weights"),
                _expect_rejected(replace(model, channels=model.channels[1:], metadata=incomplete_metadata), basis,
                                 "incomplete mesh with per-spin weight sums still one")]
    duplicate_metadata = deepcopy(model.metadata)
    records = duplicate_metadata["source_channels"]
    records[0]["k_source_2pi_alat"] = records[1]["k_source_2pi_alat"]
    rejected.append(_expect_rejected(replace(model, metadata=duplicate_metadata), basis,
                                     "duplicate source k with unchanged normalized weights"))
    rejected.append(_expect_rejected(replace(model, q_cart_2pi_alat=basis @ np.array([0.25, 0., 0.])), basis,
                                     "incompatible q"))
    off_mesh = deepcopy(model.metadata)
    off_mesh["source_channels"][0]["k_source_2pi_alat"] = (basis @ np.array([0.25, 0., 0.])).tolist()
    rejected.append(_expect_rejected(replace(model, metadata=off_mesh), basis, "off-mesh k"))
    rejected.append(_expect_rejected(model, np.zeros((3, 3)), "singular reciprocal basis"))
    paths = [Path(__file__).resolve(), ROOT / "src/substrate_framework/pw_prepared_noise.py",
             ROOT / "src/substrate_framework/pw_transition_space.py",
             ROOT / "src/substrate_framework/qe_pw_response_frames.py"]
    return {"status": "pass", "scope": "new finite PW preparation adapter only; no actual frame or material calculation",
            "maxabs_errors": errors, "raw_tail": float(expected_raw[-1]), "weighted_tail": 0.,
            "raw_probability_density_variance": good_variance, "wrong_weighted_probability_variance": wrong_variance,
            "rejected": rejected, "metadata": metadata,
            "sources": [{"path": str(path.relative_to(ROOT)), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                        for path in paths]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": receipt["status"], "output": str(args.output)}))

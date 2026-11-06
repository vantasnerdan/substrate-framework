#!/usr/bin/env python3
"""Independent finite oracle for the new checked general-k dyad adapter.

Only the existing synthetic record generator is reused. Literal Fourier
sums and a separately differentiated analytic NC matrix supply expectations.
No material, old oracle, or saved numerical result is executed here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import time
import traceback

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT / "src"), str(HERE.parent / "0027")]
from bloch_oracle_fixture import fixture
from substrate_framework.qe_bloch_grid import build_bloch_grid, UNITARY_POLICY
from substrate_framework.qe_bloch_dyads import BlochDyadCoupling


def binding(path):
    path = Path(path).resolve()
    data = path.read_bytes()
    return dict(path=str(path), bytes=str(len(data)), sha256=hashlib.sha256(data).hexdigest())


def compare(report, name, actual, expected, *, mutant=False):
    actual, expected = np.asarray(actual), np.asarray(expected)
    if actual.shape != expected.shape or not (np.isfinite(actual).all() and np.isfinite(expected).all()):
        raise AssertionError(name + ": different shape or nonfinite comparison")
    error = float(np.max(np.abs(actual - expected), initial=0))
    allowance = 1e-12 + 1e-10 * float(np.max(np.abs(expected), initial=0))
    passed = error > max(allowance, 1e-8) if mutant else error <= allowance
    report["comparisons"].append(dict(name=name, passed=passed, mutation=mutant,
                                      max_abs_error=error, allowance=allowance))
    if not passed:
        raise AssertionError(name)


def symbolic_derivative(report):
    p, q, position = sp.symbols("p q position", real=True)
    amplitude = sp.Symbol("amplitude", complex=True)
    placed = amplitude * sp.exp(-sp.I * (p - q) * position)
    exact = sp.diff(placed, position)
    two_sides = sp.I * q * placed - sp.I * p * placed
    if sp.simplify(exact - two_sides) != 0:
        raise AssertionError("two differentiated projector sides")
    if sp.simplify(exact - sp.I * q * placed) == 0:
        raise AssertionError("omitted target derivative was not exposed")
    report["symbolic"] = dict(two_sided_derivative=True, omitted_side_exposed=True,
                                units="positive dH/dR, displacement in bohr")


def literal_operators(inputs, basis, *, wrong_k=False):
    """Whole small-space matrices from the fixture's explicit Gaussian scalar."""
    ref, miller = inputs.reference, basis.miller
    p = (basis.reciprocal_bohr_inv @ miller).T
    if not wrong_k:
        p = p + basis.k_cart_bohr_inv
    differences = miller[:, :, None] - miller[:, None, :]
    g = np.einsum("ab,bij->aij", basis.reciprocal_bohr_inv, differences)
    v = -.7 * np.exp(-.17 * np.sum(g*g, axis=0))
    # The fixture declares finite density-G support, not a Gaussian on every
    # possible orbital difference. Folded rows can differ outside that support.
    support = {tuple(row) for row in inputs.density.miller.T}
    covered = np.array([tuple(row) in support for row in differences.reshape(3, -1).T])
    v = np.where(covered.reshape(v.shape), v, 0.)
    centered = np.exp(-.2 * np.sum(p*p, axis=1)) / np.sqrt(4*np.pi)
    gradients, one_sided = [], []
    for position in ref.positions_bohr.T:
        local = v * np.exp(-1j * np.einsum("aij,a->ij", g, position))
        beta = centered * np.exp(-1j * (p @ position))
        nonlocal_matrix = 1.3 * beta[:, None] * beta.conj()[None, :]
        for axis in range(3):
            local_gradient = -1j * g[axis] * local
            bra = 1j * nonlocal_matrix * p[None, :, axis]
            ket = -1j * p[:, None, axis] * nonlocal_matrix
            gradients.append(local_gradient + bra + ket)
            one_sided.append(local_gradient + bra)
    return np.asarray(gradients), np.asarray(one_sided)


def one_frame(report, inputs, frame, label, rng):
    c = BlochDyadCoupling(inputs, frame.basis, block_size=2)
    n = c.npw
    left = (rng.normal(size=(n, 3)) + 1j*rng.normal(size=(n, 3))) / np.sqrt(n)
    right = (rng.normal(size=(n, 2)) + 1j*rng.normal(size=(n, 2))) / np.sqrt(n)
    coefficient = rng.normal(size=(3, 2)) + 1j*rng.normal(size=(3, 2))
    dyad = left @ coefficient @ right.conj().T
    density, gradient = c.moments(left, right, coefficient)
    coordinates = np.column_stack([a.ravel(order="F") for a in np.indices(inputs.reference.fft_grid)])
    coordinates = coordinates / np.asarray(inputs.reference.fft_grid)
    waves = np.exp(2j*np.pi*(coordinates @ frame.basis.miller))
    literal_density = np.sum((waves @ left @ coefficient) * (waves @ right).conj(), axis=1) / c.cell_volume_bohr3
    matrices, one_side = literal_operators(inputs, frame.basis)
    expected_gradient = np.einsum("aij,ji->a", matrices, dyad)
    compare(report, label+":density_literal", density, literal_density)
    compare(report, label+":all_ionic_gradients", gradient, expected_gradient)
    compare(report, label+":integrated_number", c.grid_measure_bohr3*np.sum(density), np.trace(dyad))
    direction = rng.normal(size=c.ndof) + 1j*rng.normal(size=c.ndof)
    action = c.ionic_action(direction, left)
    compare(report, label+":ionic_action", action, np.einsum("a,aij->ij", direction, matrices) @ left)
    compare(report, label+":ionic_trace_pair", direction @ gradient, np.trace(coefficient @ right.conj().T @ action))
    field = np.cos(2*np.pi*coordinates[:, 0]) + .3j*np.sin(2*np.pi*coordinates[:, 1])
    local = c.potential_action(field, left)
    expected_local = waves.conj().T @ (field[:, None]*(waves @ left)) / c.nr
    compare(report, label+":potential_action", local, expected_local)
    compare(report, label+":density_transpose", c.grid_measure_bohr3*np.sum(field*density), np.trace(coefficient @ right.conj().T @ local))
    compare(report, label+":complex_bilinear_mutation", np.vdot(direction, gradient), direction @ gradient, mutant=True)
    compare(report, label+":one_projector_side_mutation", np.einsum("aij,ji->a", one_side, dyad), expected_gradient, mutant=True)
    if np.any(frame.basis.k_cart_bohr_inv):
        wrong, _ = literal_operators(inputs, frame.basis, wrong_k=True)
        compare(report, label+":wrong_k_projector_mutation", np.einsum("aij,ji->a", wrong, dyad), expected_gradient, mutant=True)
    # Opposite rectangular coefficient orientation exercises the other bounded leg.
    reverse_density, reverse_gradient = c.moments(right, left, coefficient.conj().T)
    compare(report, label+":adjoint_dyad_density", reverse_density, density.conj())
    compare(report, label+":adjoint_dyad_gradient", reverse_gradient, gradient.conj())
    # Wide source contraction must reproduce the same literal full dyad.
    square = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
    wide_density, wide_gradient = c.prepared_moments(left, square)
    wide_dyad = left @ square @ left.conj().T
    expected_wide_density = np.sum((waves @ left @ square)*(waves @ left).conj(), axis=1)/c.cell_volume_bohr3
    compare(report, label+":wide_source_density", wide_density, expected_wide_density)
    compare(report, label+":wide_source_gradient", wide_gradient, np.einsum("aij,ji->a", matrices, wide_dyad))
    report["frames"].append(dict(label=label, spin=str(frame.ispin), npw=str(n),
                                  k_bohr_inv=frame.basis.k_cart_bohr_inv.tolist(),
                                  fold=frame.basis.slot.reciprocal_shift.tolist()))


def main(output):
    started = time.monotonic()
    report = dict(schema="P254-general-k-dyad-oracle-v1", status="running", comparisons=[], frames=[])
    code = 1
    try:
        symbolic_derivative(report)
        with tempfile.TemporaryDirectory(prefix="p254-bloch-dyad-") as directory:
            inputs, bindings, _, _ = fixture(directory)
            grid = build_bloch_grid(inputs.reference, bindings, symmetry_policy=UNITARY_POLICY)
            rng = np.random.default_rng(20261028)
            for spin in (1, 2):
                folded = next(s for s in grid.slots if s.ispin == spin and np.any(s.reciprocal_shift)
                              and np.any(s.k_cart_bohr_inv))
                one_frame(report, inputs, grid.load_frame(folded), f"spin{spin}-folded", rng)
            gamma = next(s for s in grid.slots if s.ispin == 1 and s.mesh_index == (0, 0, 0))
            one_frame(report, inputs, grid.load_frame(gamma), "gamma-interface", rng)
        report["status"], code = "passed", 0
    except BaseException:
        report["status"] = "failed"
        report["exception"] = traceback.format_exc()
    finally:
        report["elapsed_seconds"] = time.monotonic() - started
        report["sources"] = [binding(p) for p in (
            Path(__file__), ROOT/"src/substrate_framework/qe_bloch_dyads.py",
            HERE.parent/"0027/bloch_oracle_fixture.py")]
        Path(output).write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(dict(status=report["status"], elapsed_seconds=report["elapsed_seconds"], output=str(output))))
    return code


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    raise SystemExit(main(parser.parse_args().output))

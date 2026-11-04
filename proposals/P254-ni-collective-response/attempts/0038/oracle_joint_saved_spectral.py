"""Focused tiny owning oracle for A38 saved moments and joint I/E boundary.

Explicit invocation only: serializes a small spectral fixture through the same
public loader, uses independent expm/time quadrature and direct Fock traces.
Unchanged A31/A36 proofs are carried, not imported or rerun. No material input,
physical context or scientific calculation is executed on import.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[3] / 'src'), str(HERE.parent / '0037'),
               str(HERE.parent / '0031')]

import numpy as np
from numpy.polynomial import Chebyshev, Legendre
from numpy.polynomial.chebyshev import chebval
from numpy.polynomial.legendre import leggauss, legval
from scipy.linalg import expm
from expanded_influence_readout import save_json, save_npz
from integrate_saved_spectral import (adjoint_moments_from_forward,
    integrate_saved_spectral, load_saved_spectral)
from joint_influence_readout import (LABELS, PARTS, field_source_amplitude,
    integrate_field, joint_moments, joint_operators, joint_reconstruction)
from filtered_quadratic_bath import _fourier_legendre_moments
from substrate_framework.streamed_material_state import authenticate, file_binding

TOLERANCE = 2e-11
SENSITIVITY = 1e-5
CONTRACT = dict(path=str(HERE / 'implementation-contract-001.json'), bytes=15306,
    sha256='623ef2c7ebc3fec3b0166783b165458e349abaf6d4d06c4da8bb5c77433f2668')


def agree(name, actual, expected, records):
    error = float(np.linalg.norm(np.asarray(actual) - expected) / max(1.0, np.linalg.norm(expected)))
    records[name] = error
    if not np.isfinite(error) or error > TOLERANCE:
        raise AssertionError((name, error, TOLERANCE))


def expose(name, wrong, expected, records):
    error = float(np.linalg.norm(np.asarray(wrong) - expected) / max(1.0, np.linalg.norm(expected)))
    records[name] = error
    if not np.isfinite(error) or error < SENSITIVITY:
        raise AssertionError((name, error, SENSITIVITY))


def fixture():
    # An explicit nonunitary frame and analytic inverse: no eigensolver needed.
    nilpotent = np.array([[0, .3 + .2j, -.1j], [0, 0, -.25 + .15j], [0, 0, 0]], complex)
    vectors = np.eye(3) + nilpotent
    inverse = np.eye(3) - nilpotent + nilpotent @ nilpotent
    energies = np.array([.1 + .12j, .1 + .12j, 1.3 - .08j], dtype=np.complex128)
    h = (vectors * energies[None, :]) @ inverse
    conversion = np.zeros((8, 8), dtype=float)
    for degree in range(8):
        # Independent NumPy polynomial conversion, not the production recurrence.
        column = Chebyshev.basis(degree).convert(kind=Legendre).coef
        conversion[:len(column), degree] = column
    saved = dict(eigenvalues=energies, right_eigenvectors=vectors,
        inverse_right_eigenvectors=inverse, energy_gaps=energies[:, None] - energies[None, :],
        chebyshev_to_legendre=conversion)
    rng = np.random.default_rng(380071)
    forces = rng.normal(size=(3, 3, 3)) + 1j * rng.normal(size=(3, 3, 3))
    # Intentionally independent, not manufactured F†. The quadrature must follow
    # this actual supplied leg; a force-adjoint substitution must be exposed.
    adjoints = rng.normal(size=(3, 3, 3)) + 1j * rng.normal(size=(3, 3, 3))
    bounds = np.array([[-.2, .17], [.17, .81], [.81, 1.7]])
    coefficients = (rng.normal(size=(3, 8, 3, 2))
                    + 1j * rng.normal(size=(3, 8, 3, 2))) / (1 + np.arange(8)[None, :, None, None])
    return h, saved, forces, adjoints, bounds, coefficients


def direct_quadrature(h, forces, adjoints, bounds, coefficients, order):
    """Independent matrix exponentials at Gauss time nodes; no spectral gaps."""
    nodes, weights = leggauss(order)
    shape = (coefficients.shape[-1], len(h), len(h))
    forward, backward = np.zeros(shape, complex), np.zeros(shape, complex)
    for segment, (left, right) in enumerate(bounds):
        half, midpoint = (right - left) / 2, (right + left) / 2
        for node, weight in zip(nodes, weights):
            time = midpoint + half * node
            plus, minus = expm(1j * h * time), expm(-1j * h * time)
            plus_a, minus_a = expm(1j * h.conj().T * time), expm(-1j * h.conj().T * time)
            value = chebval(node, coefficients[segment])
            evolved, evolved_a = plus @ forces @ minus, plus_a @ adjoints @ minus_a
            forward += half * weight * np.einsum('mp,mij->pij', value, evolved)
            backward += half * weight * np.einsum('mp,mij->pij', value.conj(), evolved_a)
    return forward, backward


def moment_oracle(saved, bounds, records):
    """Every degree, exact zero/complex/signed gaps and all intervals separately."""
    gaps = saved['energy_gaps']
    nodes, weights = leggauss(64)
    for segment, (left, right) in enumerate(bounds):
        half, midpoint = (right - left) / 2, (right + left) / 2
        moments = _fourier_legendre_moments(gaps, half, midpoint, 8)
        expected = np.empty_like(moments)
        for degree in range(8):
            polynomial = legval(nodes, np.eye(8)[degree])
            expected[degree] = half * np.einsum('q,qij->ij', weights * polynomial,
                np.exp(1j * (midpoint + half * nodes[:, None, None]) * gaps.conj()))
        agree('adjoint_moments_all_degrees_interval_' + str(segment),
              adjoint_moments_from_forward(moments), expected, records)
        expose('mutation_missing_transpose_interval_' + str(segment), moments.conj(), expected, records)
        expose('mutation_missing_conjugation_interval_' + str(segment),
               moments.swapaxes(-1, -2), expected, records)


def fock_operators(size):
    annihilators = np.zeros((size, 2 ** size, 2 ** size), complex)
    for mode in range(size):
        for state in range(2 ** size):
            if state & (1 << mode):
                sign = (-1) ** (state & ((1 << mode) - 1)).bit_count()
                annihilators[mode, state ^ (1 << mode), state] = sign
    return annihilators


def lift(matrix, annihilators):
    value = np.zeros(annihilators.shape[1:], complex)
    for row in range(len(matrix)):
        for column in range(len(matrix)):
            value += matrix[row, column] * annihilators[row].conj().T @ annihilators[column]
    return value


def joint_fock_oracle(h, saved, currents, bounds, coefficients, records):
    """Run the same six-channel field routing, operator sum, and 25-pair reader."""
    field = coefficients[:, :, :, :1]
    adjoints = currents.conj().swapaxes(-1, -2)
    projector = np.diag([1., 0., 1.])
    # Nontrivial original polarizations (orthogonal Cartesian rotation).
    polarizations = np.array([[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]])
    matrix = 2 * polarizations
    amplitude = field_source_amplitude(64, 1 / 64, matrix)
    filtered = integrate_field(saved, currents, adjoints, projector, bounds, field, amplitude)
    # Independent scalar contraction is deliberately not production einsum/routing.
    weights = np.empty_like(field)
    for segment in range(len(bounds)):
        for degree in range(8):
            for axis in range(3):
                weights[segment, degree, axis, 0] = sum(
                    field[segment, degree, alpha, 0] * 2 * polarizations[alpha, axis]
                    for alpha in range(3))
    ref, ref_a = direct_quadrature(h, currents, adjoints, bounds, weights, 48)
    retained = np.stack([projector @ value @ projector for value in currents])
    retained_a = np.stack([projector @ value @ projector for value in adjoints])
    ref_r, ref_ra = direct_quadrature(h, retained, retained_a, bounds, weights, 48)
    agree('field_six_channel_forward', filtered['operators'], np.concatenate((ref, ref_r)), records)
    agree('field_six_channel_adjoint', filtered['adjoint_operators'], np.concatenate((ref_a, ref_ra)), records)
    ionic_full = np.array([[.4 + .1j, .3j, -.2], [.6, -.3, .25j], [.1j, .7, .2 - .3j]])
    ionic_retained = projector @ ionic_full @ projector
    ionic = dict(zip(PARTS, (ionic_full, ionic_retained, ionic_full - ionic_retained)))
    ionic_a = {name: value.conj().T for name, value in ionic.items()}
    matrices, matrix_a, _ = joint_operators(ionic, ionic_a, filtered)
    expected_matrices = dict(ionic=ionic_full, field=ref[0], joint=ionic_full + ref[0],
        retained=ionic_retained + ref_r[0], omitted=ionic_full + ref[0] - ionic_retained - ref_r[0])
    annihilators = fock_operators(3)
    index = np.arange(3)
    unitary = np.exp(2j * np.pi * index[:, None] * index[None, :] / 3) / np.sqrt(3)
    populations = np.array([.17, .52, .83])
    rho = (unitary * populations) @ unitary.conj().T
    density = np.eye(8, dtype=complex)
    for mode, population in enumerate(populations):
        c = sum(unitary[row, mode].conjugate() * annihilators[row] for row in range(3))
        density = density @ ((1 - population) * np.eye(8) + (2 * population - 1) * c.conj().T @ c)
    agree('direct_Fock_trace', np.trace(density), 1., records)
    observed_rho = np.array([[np.trace(density @ annihilators[j].conj().T @ annihilators[i])
                             for j in range(3)] for i in range(3)])
    agree('direct_Fock_nonstationary_original_state', observed_rho, rho, records)
    lifted = {name: lift(value, annihilators) for name, value in expected_matrices.items()}
    expected = {(a, b): np.trace(density @ lifted[a] @ lifted[b])
        - np.trace(density @ lifted[a]) * np.trace(density @ lifted[b]) for a in LABELS for b in LABELS}
    bank = dict(projected_frame=np.eye(3), gram=np.eye(3), density_coefficients=rho)
    moments = joint_moments(bank, matrices, matrix_a, 'tiny-coordinate-basis',
        'tiny-original-factors', 'tiny-original-state')
    actual = {(row['left'], row['right']): row['value'] for row in moments['covariance']}
    agree('joint_all25_ordered_complex_Fock_covariances',
        [actual[key] for key in expected], np.array(list(expected.values())), records)
    reconstruction = joint_reconstruction(moments, matrices, matrix_a)
    agree('joint_raw_covariance_reconstruction', reconstruction['raw_joint_minus_ionic_field_both_cross_orders'], 0., records)
    target = expected['joint', 'joint']
    for dropped in (('ionic', 'field'), ('field', 'ionic')):
        wrong = (actual['ionic', 'ionic'] + actual['field', 'field']
                 + sum(actual[key] for key in (('ionic', 'field'), ('field', 'ionic')) if key != dropped))
        expose('mutation_drop_' + '_'.join(dropped), wrong, target, records)
    expose('mutation_independent_marginal_noise_sum', actual['ionic', 'ionic'] + actual['field', 'field'], target, records)
    for name, multiplier in (('extra_N_covariance', np.sqrt(64)), ('missing_N_covariance', 1 / np.sqrt(64)),
                             ('missing_factor_two', .5), ('wrong_field_sign', -1.)):
        wrong_filtered = integrate_field(saved, currents, adjoints, projector, bounds, field, amplitude * multiplier)
        wrong_matrices, wrong_a, _ = joint_operators(ionic, ionic_a, wrong_filtered)
        wrong_moments = joint_moments(bank, wrong_matrices, wrong_a,
            'tiny-coordinate-basis', 'tiny-original-factors', 'tiny-original-state')
        wrong_values = {(row['left'], row['right']): row['value'] for row in wrong_moments['covariance']}
        expose('mutation_' + name, wrong_values['joint', 'joint'], target, records)


def run_oracle(output):
    started = perf_counter()
    directory = Path(output).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    records = {}
    h, coordinates, forces, adjoints, bounds, coefficients = fixture()
    fixture_binding = save_npz(directory / 'small-spectral-fixture.npz',
        {'spectral_' + name: value for name, value in coordinates.items()})
    saved = load_saved_spectral(fixture_binding)
    result = integrate_saved_spectral(saved, forces, adjoints, bounds, coefficients)
    reference = direct_quadrature(h, forces, adjoints, bounds, coefficients, 48)
    coarser = direct_quadrature(h, forces, adjoints, bounds, coefficients, 32)
    for index, name in enumerate(('forward', 'independent_adjoint')):
        agree('time_quadrature_refinement_' + name, coarser[index], reference[index], records)
        agree('serialized_saved_bank_' + name, result[('operators', 'adjoint_operators')[index]], reference[index], records)
    false = direct_quadrature(h, forces, forces.conj().swapaxes(-1, -2), bounds, coefficients, 48)
    expose('mutation_replace_independent_supplied_adjoint', false[1], reference[1], records)
    moment_oracle(saved, bounds, records)
    joint_fock_oracle(h, saved, forces, bounds, coefficients, records)
    contract = json.loads(authenticate(CONTRACT).read_bytes())
    proofs = [row for row in contract['settled_source_inventory'] if Path(row['path']).name
              in ('owning-block-oracles-001.json', 'owning-canonical-field-force-oracle-001.json')]
    for proof in proofs:
        authenticate(proof)
    report = dict(schema='P254-A38-joint-saved-spectral-owning-oracle-v1', status='complete',
        implementation_contract=dict(CONTRACT),
        fixture=fixture_binding, fixture_route='actual serialized NPZ through public load_saved_spectral',
        comparisons_and_exposed_mutations=records, agreement_tolerance=TOLERANCE,
        minimum_mutation_sensitivity=SENSITIVITY, carried_proofs_not_rerun=proofs,
        raw_complex_nonhermitian_H=True, independently_supplied_adjoints=True,
        all_polynomial_degrees=list(range(8)), ordered_joint_pairs=25,
        material_inputs_read=False, material_actions=0, H_basis_spectral_eigensolves=0,
        unchanged_CAR_admissibility='existing small fixture Gram/state eigendiagnostics included in elapsed_seconds',
        elapsed_seconds=perf_counter() - started,
        implementation=[file_binding(path) for path in (Path(__file__), HERE / 'integrate_saved_spectral.py',
            HERE / 'joint_influence_readout.py')],
        scope='New saved-coordinate/moment-adjoint reuse and same-bath joint I/E cross/normalization boundary only; '
            'floating-point synthetic evidence, not material accuracy, enclosure or full physical E/C/N.')
    return save_json(directory / 'result.json', report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run_oracle(args.output), allow_nan=False), flush=True)


if __name__ == '__main__':
    main()

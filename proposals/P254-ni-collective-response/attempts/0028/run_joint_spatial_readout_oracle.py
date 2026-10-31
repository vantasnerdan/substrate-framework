#!/usr/bin/env python3
"""Independent finite full-q oracle; executable-only temporary synthetic inputs.

No production material result is fabricated. Dense real-space SciPy canonical
oscillators plus an exact SymPy character table independently supply expected
values. This script has not been executed by the implementation worker.
"""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]

import numpy as np
from scipy.linalg import expm
import sympy as sp
from substrate_framework.joint_spatial_readout import coherent_local, reconstruct
from substrate_framework.joint_response_observables import spatial_pair
from substrate_framework.spatial_readout_metrics import periodic_distances, selected_weight
from run_joint_spatial_readout import encode, sources
from joint_spatial_receipt_oracle import exercise_receipts
from joint_spatial_source_local_oracle import exercise_source_local


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    checks = {}

    def close(name, actual, expected, tolerance=3e-11):
        error = float(np.max(np.abs(np.asarray(actual)-np.asarray(expected))))
        if not np.isfinite(error) or error > tolerance:
            raise AssertionError(name+': '+str(error))
        checks[name] = dict(max_absolute_error=error, tolerance=tolerance)

    def reject(name, operation):
        try:
            operation()
        except (ValueError, KeyError):
            checks[name] = dict(rejected=True)
        else:
            raise AssertionError(name+' was accepted')

    n, mesh = 4, np.array([4, 1, 1])
    cells = np.array([[r, 0, 0] for r in range(n)])
    q = cells/mesh
    # Exact, independently defined character table: unitary normalization and
    # nonzero seed cell force the sign to matter (cell1 is not its own inverse).
    exact = sp.Matrix(n, n, lambda r, k: sp.exp(2*sp.pi*sp.I*r*k/n)/sp.sqrt(n))
    identity = (exact.conjugate().T*exact).applyfunc(sp.simplify)
    if identity != sp.eye(n):
        raise AssertionError('exact character orthogonality')
    seed = sp.Matrix([sp.exp(-2*sp.pi*sp.I*k/n)/sp.sqrt(n) for k in range(n)])
    delta = (exact*seed).applyfunc(sp.simplify)
    if delta != sp.Matrix([0, 1, 0, 0]):
        raise AssertionError('exact localized kick normalization/sign')
    seed_np = np.array([complex(x.evalf()) for x in seed])
    close('coherent_local_exact_sqrtN_sign', coherent_local(q, seed_np, cells), [0, 1, 0, 0])
    checks['sympy_character_orthogonality'] = True
    if np.linalg.norm(coherent_local(q, seed_np, -cells)-[0, 1, 0, 0]) < 1:
        raise AssertionError('q sign fixture is not discriminating')
    checks['wrong_q_sign_is_detectable'] = True
    ordered_character = seed_np*np.sqrt(n)
    close('ordered_pair_exact_1overN_sign_no_extra_dagger',
          [spatial_pair(q, ordered_character, cell, [0, 0, 0], ncells=n) for cell in cells],
          [0, 1, 0, 0])

    # Two local canonical species (electron/ion origin labels) and ONE global
    # field oscillator. Dense real-space construction, never q evolution.
    d = 2*n+1
    k = np.diag([2.+.7*(i % 2) for i in range(2*n)]+[3.1])
    for r in range(n):
        for target, coefficient in ((2*r+1, .19), (2*((r+1) % n)+1, .31)):
            k[2*r, target] += coefficient
            k[target, 2*r] += coefficient
        k[2*r+1, -1] = k[-1, 2*r+1] = .13/np.sqrt(n)
    j = np.block([[np.zeros((d, d)), np.eye(d)], [-np.eye(d), np.zeros((d, d))]])
    h = np.block([[k, np.zeros((d, d))], [np.zeros((d, d)), np.eye(d)]])
    t, s = .47, .23
    left, right = expm(t*j@h), expm(s*j@h)
    force = np.zeros((n, 2*d))
    force[:, :d] = -k[np.arange(n)*2+1]
    origins = []
    for label in range(3):
        w = np.zeros((2*d, 2*d), complex)
        ix = list(range(label, 2*n, 2)) if label < 2 else [d-1]
        for a in ix:
            omega = 1.1+.3*label
            w[a, a], w[d+a, d+a] = 1/(2*omega), omega/2
            w[a, d+a], w[d+a, a] = .5j, -.5j
        origins.append(w)
    real_ab = [force@left@w@right.T@force.T for w in origins]
    real_ba = [force@left@w.T@right.T@force.T for w in origins]
    fourier = np.array(exact.conjugate().T.evalf()).astype(complex)
    ab_parts = np.array([np.diag(fourier@x@fourier.conj().T) for x in real_ab]).T
    ba_parts = np.array([np.diag(fourier@x@fourier.conj().T) for x in real_ba]).T
    ab = np.column_stack((ab_parts, ab_parts.sum(axis=1)))
    ba = np.column_stack((ba_parts, ba_parts.sum(axis=1)))
    kick = np.zeros(2*d)
    kick[d+3] = 1.7  # ionic momentum in original cell1
    local_force = force@left@kick
    q_force = fourier@local_force
    geometry = np.array([[2., .7, .1], [0., 1.8, .2], [0., 0., 2.3]])
    reconstructed = reconstruct(q, q_force, ab, ba, mesh=mesh, cell_bohr=geometry,
                                kick_cell=[1, 0, 0], reference_cell=[1, 0, 0],
                                probe_position_bohr=np.zeros(3), seed_position_bohr=np.zeros(3))
    close('coherent_dense_real_space', reconstructed['coherent_force_ry_bohr'], local_force)
    expected_ab = np.column_stack([x[:, 1] for x in real_ab]+[sum(real_ab)[:, 1]])
    expected_ba = np.column_stack([x[:, 1] for x in real_ba]+[sum(real_ba)[:, 1]])
    for name, expected in (('AB', expected_ab), ('BA', expected_ba),
                           ('symmetrized', (expected_ab+expected_ba)/2),
                           ('commutator', expected_ab-expected_ba)):
        actual = reconstructed['connected_force_pairs_ry2_bohr2'][name]
        close(name+'_dense_real_space_all_origins', actual, expected)
        close(name+'_origin_sum', actual[:, :3].sum(axis=1), actual[:, 3])
    if np.linalg.norm(expected_ab[:, 2]) == 0:
        raise AssertionError('field origin oracle lacks signal')
    checks['nonzero_global_field_origin'] = True
    if abs(ab[1, 2])+abs(ab[2, 2])+abs(ab[3, 2]) > 1e-11:
        raise AssertionError('one global field leaked out of Gamma')
    close('opposite_ordering_reality', expected_ab, expected_ba.conj())
    perturbed = reconstruct(q, q_force+1j*1e-4, ab, ba, mesh=mesh, cell_bohr=geometry,
                            kick_cell=[1, 0, 0], reference_cell=[1, 0, 0],
                            probe_position_bohr=np.zeros(3), seed_position_bohr=np.zeros(3))
    if perturbed['reality_residuals']['coherent_imaginary_max'] < 1e-4:
        raise AssertionError('imaginary coherent signal was clipped')
    checks['imaginary_signal_retained'] = True

    distances = periodic_distances(cells, [1, 0, 0], mesh, geometry)
    brute = [min(np.linalg.norm(geometry@(cell-[1, 0, 0]-mesh*np.array([a, b, c])))
                 for a in range(-4, 5) for b in range(-4, 5) for c in range(-4, 5)) for cell in cells]
    close('skew_geometry_minimum_images', distances, brute)
    offset = np.array([-1.4, .7, .3])
    shifted = reconstruct(q, q_force, ab, ba, mesh=mesh, cell_bohr=geometry,
        kick_cell=[1, 0, 0], reference_cell=[1, 0, 0],
        probe_position_bohr=offset, seed_position_bohr=np.zeros(3))
    shifted_brute = [min(np.linalg.norm(geometry@(cell-[1, 0, 0]-mesh*np.array([a, b, c]))+offset)
                        for a in range(-4, 5) for b in range(-4, 5) for c in range(-4, 5)) for cell in cells]
    close('distinct_atom_skew_minimum_images',
          shifted['selected_observable_weights']['coherent_total']['distance_bohr_by_cell'], shifted_brute)
    close('distinct_atom_position_does_not_change_force', shifted['coherent_force_ry_bohr'], local_force)
    close('connected_same_probe_offsets_cancel',
          shifted['selected_observable_weights']['connected']['AB']['total']['distance_bohr_by_cell'], brute)
    exact_offset = sp.Matrix([0, 0, 0])-sp.Matrix([1, 0, 0])
    exact_distance = min(sp.sqrt((exact_offset-10*sp.Matrix([a, b, c])).dot(
        exact_offset-10*sp.Matrix([a, b, c])))
        for a in range(-1, 2) for b in range(-1, 2) for c in range(-1, 2))
    atomic = reconstruct([[0, 0, 0]], [1.], np.ones((1, 4)), np.ones((1, 4)),
        mesh=[1, 1, 1], cell_bohr=10*np.eye(3), kick_cell=[0, 0, 0], reference_cell=[0, 0, 0],
        probe_position_bohr=[0, 0, 0], seed_position_bohr=[1, 0, 0])
    close('same_cell_distinct_atom_rms_sympy',
          atomic['selected_observable_weights']['coherent_total']['rms_distance_bohr'], float(exact_distance))
    close('same_cell_same_probe_connected_rms', atomic['selected_observable_weights']['connected']['AB']['total']['rms_distance_bohr'], 0.)
    if exact_distance == 0 or np.allclose(shifted_brute, brute):
        raise AssertionError('omitted intra-cell offset mutation not exposed')
    checks['omitted_atom_offset_mutation_detected'] = True
    localized = selected_weight([0, 2, 0, 0], distances, value_units='Ry/bohr')
    extended = selected_weight([2, 2, 2, 2], distances, value_units='Ry/bohr')
    zero = selected_weight(np.zeros(n), distances, value_units='Ry/bohr')
    close('localized_participation', localized['participation_cells'], 1.)
    close('extended_participation', extended['participation_cells'], 4.)
    close('localized_extent', localized['rms_distance_bohr'], 0.)
    if zero['defined'] or any(zero[key] is not None for key in
            ('participation_cells', 'participation_fraction', 'rms_distance_bohr', 'normalized_weight_by_cell')):
        raise AssertionError('zero signal falsely localized')
    checks['zero_signal_undefined'] = True

    exercise_receipts(q=q, q_force=q_force, ab=ab, ba=ba, local_force=local_force,
                      cells=cells, t=t, s=s, mesh=mesh,
                      expected_ab=expected_ab, expected_ba=expected_ba,
                      encode=encode, close=close, reject=reject)
    exercise_source_local(encode=encode, close=close, reject=reject)
    bindings = sources()
    result = dict(schema='P254-joint-spatial-readout-oracle-v1', status='passed',
        source_bindings=bindings, source_file_count=len(bindings), check_count=len(checks), checks=checks,
        scope='finite4-q independent real-space/SymPy consumer oracle only; no material64-q input available',
        material_spatial_result=False)
    with args.output.open('x') as stream:
        json.dump(encode(result), stream, indent=2, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Focused finite SymPy/dense oracle for new aggregations, not a solver replay.

Fixtures are explicitly synthetic and confined to this oracle. No material
summary is emitted. Settled minimum-image and scalar participation proofs are
carried forward by exact source bindings rather than re-executed here.
"""
import argparse
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[3]/'src'), str(HERE)]
import numpy as np
import sympy as sp
from clamped_scientific_readout import (
    FIELD_ORDER, ION_ORDER, REGISTRATION, authenticate, binding, compare_pair,
    decoded, displacement_readout, encode, maxwell_authority, maxwell_energy,
    readout_time, sources, validate_identity, validate_observation,
)


def fixture(control):
    n, nat = 3, 2
    cells = np.array([[i, 0, 0] for i in range(n)])
    ids = np.array([7, 13])
    mass = dict(variant='original', original_masses_amu=[58., 1.], effective_masses_amu=[58., 1.])
    specification = dict(case='oracle-only', aggregate=dict(path='oracle-only', bytes='0', sha256='oracle-only'))
    identity = dict(case='oracle-only', control=control, mass=mass, phase_space_order=ION_ORDER,
                    field_order=FIELD_ORDER, atom_indices=ids.tolist(), atom_names=['Ni', 'H'],
                    cell_indices=cells.tolist(), cell_bohr=np.eye(3).tolist(),
                    matched_joint_specification=specification,
                    matched_joint_identity=sha256(json.dumps(specification, sort_keys=True).encode()).hexdigest(),
                    aggregate=specification['aggregate'], archive=specification['aggregate'],
                    electrons_held=True, extra_delta_curvature='exact zero',
                    preparation=dict(field_mass=(n if control == 'periodic' else 1)/(8*np.pi)))
    inputs = dict(schema='P254-matched-clamped-inputs-v1', identity=identity, kick_cell=[1, 0, 0],
                  kick_atom_original_one_based=13, kick_axis='y', readout_phase_space_indices=[4, 10],
                  field_index='singleton global field' if control == 'periodic' else 'original cell order, one independent field per region',
                  field_covariance_pair='field[index](t), field at kick region(s); singleton for global field')
    result = dict(case='oracle-only', mass='original', control=control)
    geometry = validate_identity(inputs, result)
    fields = 1 if control == 'periodic' else n
    arrays = dict(baseline_ions=np.zeros((n, 1, 1, 6*nat), complex),
                  baseline_field=np.zeros((fields, 6), complex),
                  local_kick_delta_ions=np.zeros((n, 1, 1, 6*nat), complex),
                  common_kick_delta_ions=np.zeros((n, 1, 1, 6*nat), complex),
                  local_kick_delta_field=np.zeros((fields, 6), complex),
                  common_kick_delta_field=np.zeros((fields, 6), complex),
                  spatial_equal_time_ordered=np.zeros((n, 2, 2), complex),
                  spatial_two_time_ordered=np.zeros((n, 2, 2), complex),
                  spatial_two_time_commutator=np.zeros((n, 2, 2), complex),
                  field_equal_time_ordered=np.zeros((fields, 6, 6), complex),
                  field_two_time_ordered=np.zeros((fields, 6, 6), complex),
                  ion_field_two_time_ordered=np.zeros((n, 2, 6), complex))
    anchor = 0 if fields == 1 else 1
    w = np.diag([.5, .6, .7, 2., 3., 5.]).astype(complex)
    w[:3, 3:] = .5j*np.eye(3)
    w[3:, :3] = -.5j*np.eye(3)
    arrays['field_equal_time_ordered'][anchor] = w
    arrays['field_two_time_ordered'][anchor] = w*(1+.2j)
    arrays['baseline_field'][:, 3:] = np.arange(1, 3*fields+1).reshape(fields, 3)/7
    arrays['local_kick_delta_field'][:, 3:] = .17
    arrays['common_kick_delta_field'][:, 3:] = -.11
    arrays['baseline_ions'][..., :6] = .25
    arrays['baseline_ions'][..., 6:] = 1e9
    arrays['local_kick_delta_ions'][1, 0, 0, 4] = 2+3j
    arrays['common_kick_delta_ions'][:, 0, 0, 4] = (2+3j)/np.sqrt(n)
    arrays['spatial_two_time_ordered'][1] = [[1+2j, 3-4j], [5+6j, 7-8j]]
    return inputs, result, geometry, arrays


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    inventory = sources()  # raw source bytes recorded BEFORE numerical work
    registration = json.loads(REGISTRATION.read_text())
    checks, mutations = {}, {}

    def close(name, actual, expected, tolerance=2e-11):
        error = float(np.max(np.abs(np.asarray(actual)-np.asarray(expected))))
        if not np.isfinite(error) or error > tolerance:
            raise AssertionError(name+': '+str(error))
        checks[name] = dict(max_absolute_error=error, tolerance=tolerance)

    def truth(name, condition):
        if not condition:
            raise AssertionError(name)
        checks[name] = True

    def reject(name, operation):
        try:
            operation()
        except (ValueError, KeyError):
            mutations[name] = dict(rejected=True)
        else:
            raise AssertionError(name+' accepted')

    # Independent exact Maxwell normalization in Ry with E in Hartree units.
    V, E, P, D, W, N = sp.symbols('V E P D W N', positive=True)
    M = V/(8*sp.pi)
    exact = sp.simplify((P**2/(2*M)).subs(P, -2*M*E)-V*E**2/(4*sp.pi))
    truth('sympy_Maxwell_Ry_electric_normalization', exact == 0)
    exact_delta = sp.expand(((P+D)**2+W)/(2*M)-(P**2+W)/(2*M))
    truth('sympy_moving_baseline_cross_term', sp.simplify(exact_delta-(2*P*D+D**2)/(2*M)) == 0)
    global_cov = W/(2*(N*V/(8*sp.pi)))
    independent_cov = N*W/(2*(V/(8*sp.pi)))
    truth('sympy_explicit_volume_and_covariance_multiplicity', sp.simplify(independent_cov/global_cov-N**2) == 0)
    x, y, z = sp.symbols('x y z', real=True)
    u = sp.Matrix([x+sp.I*y, z, 0])
    truth('sympy_displacement_vector_weight', sp.simplify((u.conjugate().T*u)[0]-(x*x+y*y+z*z)) == 0)
    close('string_metadata_real_imag_decoding', decoded({'real': [7, 13], 'imag': [0, 0]}), [7, 13])
    truth('species_string_metadata_decoding', decoded({'real': ['Ni', 'H'], 'imag': ['', '']}).tolist() == ['Ni', 'H'])

    outputs = {}
    for control in ('periodic', 'independent-regions'):
        inputs, result, geometry, arrays = fixture(control)
        authority = maxwell_authority(inputs, geometry)
        validate_observation(arrays, geometry)
        fields, mass = authority['field_count'], authority['mass_per_field']
        mean = arrays['baseline_field']
        actual = maxwell_energy(mean, arrays['field_equal_time_ordered'], geometry, authority)
        # Independent dense quadratic scalar: block-diagonal N-field Hamiltonian,
        # full covariance, and concatenated canonical mean. No reducer is reused.
        h = np.diag(np.tile([0., 0., 0., 1/mass, 1/mass, 1/mass], fields))
        dense_w = np.kron(np.eye(fields), np.diag([.5, .6, .7, 2., 3., 5.]))
        dense_mean = mean.reshape(-1)
        expected_mean = dense_mean@h@dense_mean/2
        expected_covariance = np.trace(h@dense_w)/2
        close(control+'_dense_Maxwell_mean', actual['mean_energy_ry'], expected_mean)
        close(control+'_dense_Maxwell_covariance', actual['covariance_energy_ry'], expected_covariance)
        close(control+'_dense_Maxwell_total', actual['energy_ry'], expected_mean+expected_covariance)
        close(control+'_per_field_sum', np.sum(actual['energy_by_field_ry']), expected_mean+expected_covariance)
        altered = mean.copy()
        altered[:, :3] = 1e7
        close(control+'_excludes_contact_A_squared', maxwell_energy(altered, arrays['field_equal_time_ordered'], geometry, authority)['energy_ry'], actual['energy_ry'])
        complex_mean = mean.copy()
        complex_mean[0, 3] += .23j
        complex_cov = arrays['field_equal_time_ordered'].copy()
        anchor = 0 if fields == 1 else geometry['anchor_index']
        complex_cov[anchor, 3, 3] += .31j
        complex_expected = (complex_mean.reshape(-1)@h@complex_mean.reshape(-1)+np.trace(h@dense_w))/2+fields*.31j/(2*mass)
        close(control+'_retains_complex_no_conjugation', maxwell_energy(complex_mean, complex_cov, geometry, authority)['energy_ry'], complex_expected)
        receipt = dict(time_inverse_ry=2., reference_time_inverse_ry=.5, archive={'fixture': True}, diagnostics={'raw_complex': {'real': 1, 'imag': .3}})
        distances = np.array([1., 0., 1.])
        readout = readout_time(arrays, receipt, inputs, geometry, authority, distances)
        expected_delta = np.sum(2*mean[:, 3:]*.17+.17**2)/(2*mass)
        close(control+'_kick_delta_energy_includes_moving_baseline', readout['states']['local']['maxwell']['change_from_same_time_baseline_ry'], expected_delta)
        close(control+'_preserves_complex_connected_pairs', readout['connected_ordered_complex']['spatial_two_time_ordered'], arrays['spatial_two_time_ordered'])
        truth(control+'_preserves_reported_diagnostics', readout['producer_diagnostics'] == receipt['diagnostics'])
        d = displacement_readout(arrays['local_kick_delta_ions'], geometry, distances)
        close(control+'_all_ion_displacement_weights', d['all_ions']['metric']['squared_weight_by_cell'], [0, 13, 0])
        close(control+'_H_displacement_weights', d['H']['metric']['squared_weight_by_cell'], [0, 13, 0])
        truth(control+'_Ni_zero_is_undefined', not d['Ni']['metric']['defined'])
        common = displacement_readout(arrays['common_kick_delta_ions'], geometry, distances)
        close(control+'_common_local_total_displacement_weight', common['all_ions']['metric']['total_squared_weight'], d['all_ions']['metric']['total_squared_weight'])
        close(control+'_common_species_aggregation', common['H']['metric']['squared_weight_by_cell'], np.full(3, 13/3))
        baseline = readout['states']['baseline']['absolute_displacement']
        close(control+'_excludes_momentum_units', baseline['all_ions']['metric']['squared_weight_by_cell'], np.full(3, 6*.25**2))
        absolute = readout['states']['local']['absolute_displacement']['all_ions']['metric']['total_squared_weight']
        close(control+'_moving_displacement_baseline_separate', absolute, 18*.25**2+13+2*.25*2)
        absent = deepcopy(geometry)
        absent['names'] = np.array(['Ni', 'Ni'])
        truth(control+'_absent_H_is_not_applicable', not displacement_readout(arrays['baseline_ions'], absent, distances)['H']['applicable'])
        bad = deepcopy(inputs)
        bad['identity']['preparation']['field_mass'] *= 3
        reject(control+'_wrong_volume', lambda: maxwell_authority(bad, geometry))
        missing = deepcopy(inputs)
        del missing['identity']['preparation']['field_mass']
        gap = maxwell_authority(missing, geometry)
        truth(control+'_missing_mass_exact_gap', not gap['available'] and 'identity.preparation.field_mass' in gap['gap'])
        still = readout_time(arrays, receipt, missing, geometry, gap, distances)
        truth(control+'_missing_mass_finishes_displacement', still['states']['local']['kick_delta_displacement']['H']['metric']['defined'])
        bad_arrays = {k: v.copy() for k, v in arrays.items()}
        bad_arrays['baseline_field'] = bad_arrays['baseline_field'][:-1]
        reject(control+'_missing_field_mean', lambda: validate_observation(bad_arrays, geometry))
        bad_inputs = deepcopy(inputs)
        bad_inputs['identity']['phase_space_order'] = 'p,u'
        reject(control+'_ionic_phase_order', lambda: validate_identity(bad_inputs, result))
        bad_inputs = deepcopy(inputs)
        bad_inputs['identity']['field_order'] = 'Pi_xyz,a_xyz'
        reject(control+'_field_phase_order', lambda: validate_identity(bad_inputs, result))
        bad_inputs = deepcopy(inputs)
        bad_inputs['readout_phase_space_indices'].reverse()
        reject(control+'_selected_pair_order', lambda: validate_identity(bad_inputs, result))
        bad_inputs = deepcopy(inputs)
        bad_inputs['identity']['cell_indices'].pop(1)
        reject(control+'_missing_region', lambda: validate_identity(bad_inputs, result))
        bad_inputs = deepcopy(inputs)
        bad_inputs['identity']['matched_joint_identity'] = 'mutated'
        reject(control+'_input_identity', lambda: validate_identity(bad_inputs, result))
        bad_arrays = {k: v.copy() for k, v in arrays.items()}
        del bad_arrays['spatial_two_time_ordered']
        reject(control+'_missing_correlations', lambda: validate_observation(bad_arrays, geometry))
        if fields > 1:
            bad_covariance = arrays['field_equal_time_ordered'].copy()
            bad_covariance[0] = bad_covariance[1]
            reject('independent_wrong_covariance_layout', lambda: maxwell_energy(mean, bad_covariance, geometry, authority))
            wrong_sum = np.trace(np.sum(arrays['field_equal_time_ordered'], axis=0)[3:, 3:])/(2*mass)
            truth('independent_anchor_sum_missing_multiplicity_detectable', abs(wrong_sum-expected_covariance) > 1)
        outputs[control] = dict(control=control, case='oracle-only', mass='original', original_identity=inputs['identity'],
            kick=dict(cell=[1, 0, 0], atom_original_one_based=13, axis='y', nominal_added_energy_ry=1.),
            geometry=dict(cells=geometry['cells']), original_result={'fixture': control}, field_authority=authority, observations=[readout])
    pair = compare_pair(outputs['periodic'], outputs['independent-regions'])
    truth('bound_pair_explicit_changed_field_count', pair['periodic_field_authority']['field_count'] == 1 and pair['independent_field_authority']['field_count'] == 3)
    truth('bound_pair_not_full_lattice_ranking', pair['whole_coupled_lattice_ranking'] is False)
    for name, mutate in (
        ('pair_wrong_seed', lambda b: b['kick'].update(atom_original_one_based=7)),
        ('pair_wrong_mass', lambda b: b.update(mass='H2')),
        ('pair_wrong_case', lambda b: b.update(case='other')),
        ('pair_wrong_time', lambda b: b['observations'][0].update(time_inverse_ry=3.)),
    ):
        bad = deepcopy(outputs['independent-regions'])
        mutate(bad)
        reject(name, lambda: compare_pair(outputs['periodic'], bad))
    with TemporaryDirectory(prefix='p254-clamped-readout-oracle-') as folder:
        path = Path(folder)/'raw-source'
        path.write_bytes(b'original')
        record = binding(path)
        path.write_bytes(b'mutated!')
        reject('raw_source_byte_mutation', lambda: authenticate(record))
    if sorted(checks) != sorted(registration['checks']) or sorted(mutations) != sorted(registration['mutation_witnesses']):
        raise AssertionError('oracle coverage differs from preregistration')
    output = dict(schema='P254-clamped-scientific-readout-oracle-v1', status='passed',
                  source_bindings=inventory, checks=checks, mutation_witnesses=mutations,
                  settled_carry_forward=registration['settled_carry_forward'],
                  material_computation=False, independent_expected_values='exact SymPy plus dense canonical scalar; no producer evolution')
    with args.output.open('x') as stream:
        json.dump(encode(output), stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status='passed', checks=len(checks), mutation_witnesses=len(mutations), output=binding(args.output))))


if __name__ == '__main__':
    main()

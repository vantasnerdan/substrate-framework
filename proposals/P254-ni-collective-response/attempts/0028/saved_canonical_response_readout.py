"""Authenticated saved-data API/CLI for original-atom coupled Gamma response.

API: read_saved_canonical_response(mixed_initial=..., mixed_endpoint=...,
mean_initial=..., mean_endpoint=..., metadata=..., energy_endpoint=None), where
EVERY argument is an explicit original {path,bytes,sha256} binding. It returns
NumPy/complex data; clamped_scientific_readout.encode makes lossless JSON
real/imag objects. CLI consumes source_inputs in the parent's registration.

No producer admission, Native runtime, H0, electronic leg access, propagation,
producer writes or work budget is constructed. This consumes the parent's
already admitted immutable boundaries, not a new trajectory qualification.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

import clamped_scientific_readout as clamped
import mixed_excitation_energy
import mixed_response_storage
import saved_canonical_response as algebra
from substrate_framework.material_modes import AMU_RY
from substrate_framework.partitioned_joint_slots import array_identity
from substrate_framework.streamed_material_state import authenticate, file_binding

SCHEMA = 'P254-saved-canonical-Gamma-atom-response-v1'
SCOPE = (
    'Original saved coupled Gamma coefficient, not localized spatial spreading. '
    'Kinetic entries are additive summands of the minimally coupled kinetic mixed '
    'RESPONSE coefficient only, not total energy per atom, finite-dose redistribution, '
    'covariance, heat, dissipation, coherence lifetime or material ranking. '
    'Only the supplied original initial and endpoint times are read; no interpolation '
    'or new time grid. Baseline mean and drift remain separate from excitation response.'
)


def bound_json(binding):
    """Authenticate actual bytes, retaining original JSON integer metadata."""
    return json.loads(authenticate(binding).read_text())


def same_binding(a, b):
    return (a['path'] == b['path'] and a['sha256'] == b['sha256']
            and int(a['bytes']) == int(b['bytes']))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def geometry(inputs):
    """Decode original clamped identity, retaining its source layout and masses."""
    identity = inputs['identity']
    require(inputs.get('schema') == 'P254-matched-clamped-inputs-v1', 'original clamped inputs required')
    require(identity['phase_space_order'] == clamped.ION_ORDER
            and identity['field_order'] == clamped.FIELD_ORDER, 'original atom/xyz and Cartesian field order required')
    require(identity['mass']['variant'] == 'original' and identity['control'] == 'periodic',
            'original masses and one global periodic field required')
    model = sha256(json.dumps(identity['matched_joint_specification'], sort_keys=True).encode()).hexdigest()
    require(model == identity['matched_joint_identity'], 'original model/specification identity differs')
    real = clamped.real_metadata
    ids = real(identity['atom_indices'], 'atom IDs', integer=True)
    names = clamped.decoded(identity['atom_names'])
    positions = real(identity['positions_bohr'], 'Cartesian positions')
    charges = real(identity['valence_charges'], 'valence charges')
    masses = real(identity['preparation']['mass_ry'], 'Cartesian Ry masses')
    original_amu = real(identity['mass']['original_masses_amu'], 'original atomic masses')
    effective_amu = real(identity['mass']['effective_masses_amu'], 'effective atomic masses')
    cells = real(identity['cell_indices'], 'BvK cells', integer=True)
    require(np.array_equal(ids, np.arange(1, 38)) and names.shape == (37,)
            and np.count_nonzero(names == 'Ni') == 31 and np.count_nonzero(names == 'H') == 6,
            'all37 original atom IDs and Ni31/H6 species required')
    require(positions.shape == (3, 37) and charges.shape == (37,)
            and masses.shape == (111,) and original_amu.shape == (37,)
            and np.array_equal(original_amu, effective_amu)
            and np.array_equal(masses, np.repeat(original_amu*AMU_RY, 3)),
            'original geometry/charges/Cartesian mass ordering differs')
    require(np.array_equal(cells, np.indices((4, 4, 4)).reshape(3, -1).T), 'original BvK4^3 cells required')
    # These inputs explicitly declare Cartesian a_xyz,Pi_xyz, not another mode
    # basis. The optional genuine scalar ledger also authenticates this E via
    # its actual model parameter raw-array identity, below.
    polarizations = np.eye(3)
    return dict(atom_ids=ids, species=names, positions=positions, charges=charges,
                masses=masses, polarizations=polarizations, ncells=len(cells), model=model)


def validate_pair(initial, endpoint, mean0, mean1, g):
    require(initial.identity == endpoint.identity == mean0['model_identity']
            == mean1['model_identity'] == g['model'], 'mixed/mean/metadata model identities differ')
    require(initial.time == mean0['time_ry_inverse'] == 0
            and endpoint.time == mean1['time_ry_inverse'] and endpoint.time > initial.time,
            'actual matched initial/endpoint times required')
    for mean in (mean0, mean1):
        require(mean.get('schema') == 'streamed-joint-mean-v1' and mean.get('complete') is True,
                'actual completed original mean boundary required')
    require(mean0['original_electronic_boundary'] == mean1['original_electronic_boundary'],
            'mean original electronic preparation changed')
    slots = lambda mean: [{k: v for k, v in row.items() if k != 'frame'} for row in mean['slots']]
    require(slots(mean0) == slots(mean1), 'original source order, complex D, raw rows/columns or weights changed')
    require(initial.record['provenance']['fixed_raw_D_and_slots'] == slots(mean0),
            'mixed seed differs from the original mean source/D/slot authority')
    for component in algebra.DIRECTIONS:
        require(initial.record['shapes'][component] == [row['shape'] for row in mean0['slots']],
                'Gamma raw PW rows/original source rank differs from actual mean slots')
    for key in ('identity', 'q', 'shapes', 'target_plus', 'target_minus', 'conventions'):
        require(initial.record[key] == endpoint.record[key], 'changed original mixed layout: '+key)
    require(np.array_equal(initial.q, np.zeros(3)), 'only actual Gamma checkpoints are supported')
    require(initial.target_plus == initial.target_minus == tuple(range(len(mean0['slots']))),
            'Gamma must retain original source/target order')
    require(initial.record['provenance']['model'] == g['model'], 'seed original model differs')
    seed = initial.record['provenance']['original_local_pulse']
    require(seed == dict(atom=1, axis='x', amplitude=1., cell=[0, 0, 0]), 'original Ni1/x cell0 J=1 seed required')
    expected = np.zeros(111, complex)
    expected[0] = seed['amplitude']/g['ncells']
    for c in algebra.DIRECTIONS:
        s = initial.small[c]
        require(np.array_equal(s.p, expected if c != 'mixed' else np.zeros(111))
                and not np.any(s.u) and not np.any(s.coordinate) and not np.any(s.momentum),
                'original physical per-cell canonical-momentum-only seed differs')


def energy_comparison(energy, *, mixed_binding, mean_binding, endpoint, g, kinetic):
    """Compare with genuine completed energy, never fitted/reconstructed totals."""
    require(energy.get('schema') == 'P254-native-mixed-material-energy-readout-v1'
            and energy.get('complete') is True and energy.get('status') == 'complete',
            'genuine completed original energy ledger required')
    require(same_binding(energy['input'], mixed_binding)
            and energy['time_ry_inverse'] == endpoint.time
            and energy['model_identity'] == g['model'], 'energy mixed endpoint/model/time differs')
    actual = energy['mean']['actual_stage']
    require(same_binding(actual['boundary'], mean_binding), 'energy actual mean boundary differs')
    params = actual['model_parameters']
    for array, descriptor, label in ((g['masses'], params['reference']['masses'], 'masses'),
            (g['charges'], params['charges'], 'charges'),
            (g['polarizations'], params['polarizations'], 'polarizations')):
        require(array_identity(array) == descriptor, 'energy actual parameter identity differs: '+label)
    require(params['ncells'] == g['ncells'], 'energy actual N differs')
    # Only small arrays are consumed. No density, H0, electronic scalar or leg
    # payload is read. This closes the mean/canonical conversion comparison.
    for key, name in (('u', 'u'), ('p', 'p'), ('coordinate', 'coordinate'),
                      ('momentum', 'momentum'), ('shift', 'shift'), ('velocity', 'velocity')):
        binding = energy['mean']['stage_record']['arrays'][key]
        value = np.load(authenticate(binding), allow_pickle=False)
        expected = endpoint.stages['mean'][name]
        matches = (np.allclose(value.reshape(expected.shape), expected,
                               rtol=64*np.finfo(float).eps, atol=0)
                   if key == 'velocity' else np.array_equal(value.reshape(expected.shape), expected))
        require(matches, 'energy actual saved mean small array differs: '+key)
    pair = energy['components_global_ry']['ionic_kinetic_minimal']
    require(isinstance(pair, list) and len(pair) == 2, 'actual complex kinetic scalar required')
    expected = complex(*pair)
    observed = np.sum(kinetic)
    tolerance = 128*np.finfo(float).eps*max(float(np.sum(np.abs(kinetic))), abs(expected))
    return dict(component='ionic_kinetic_minimal', atom_sum_global_ry=observed,
                genuine_component_global_ry=expected, difference_global_ry=observed-expected,
                roundoff_comparison_absolute_tolerance_ry=tolerance,
                agrees_within_roundoff=bool(abs(observed-expected) <= tolerance),
                comparison_scope='same scalar component; not a total-energy allocation or accuracy certificate')


def read_saved_canonical_response(*, mixed_initial, mixed_endpoint, mean_initial,
                                  mean_endpoint, metadata, energy_endpoint=None):
    """Return authenticated initial/endpoint/difference37-atom and species tables."""
    inputs = bound_json(metadata)
    mean0, mean1 = bound_json(mean_initial), bound_json(mean_endpoint)
    g = geometry(inputs)
    tables, baselines, kinetics, all_stages, layouts = {}, {}, {}, {}, {}
    with mixed_response_storage.MappedMixedState(mixed_initial) as first, \
            mixed_response_storage.MappedMixedState(mixed_endpoint) as last:
        validate_pair(first, last, mean0, mean1, g)
        for label, state, mean in (('initial', first, mean0), ('endpoint', last, mean1)):
            stages = algebra.canonical_stages(mean, state.small, natoms=37, nfield=3)
            stages = algebra.mechanical_stages(stages, masses=g['masses'], charges=g['charges'],
                                               polarizations=g['polarizations'])
            all_stages[label] = stages
            kinetics[label] = algebra.kinetic_atomwise(stages, masses=g['masses'], ncells=g['ncells'])
            tables[label] = algebra.atom_table(stages, kinetics[label], atom_ids=g['atom_ids'],
                                               species=g['species'], positions=g['positions'])
            tables[label]['time_ry_inverse'] = state.time
            baselines[label] = stages['mean']
            layouts[label] = {key: state.record[key] for key in ('identity', 'time', 'kind', 'q',
                                        'shapes', 'target_plus', 'target_minus', 'conventions', 'provenance')}
        delta = {c: {k: all_stages['endpoint'][c][k]-all_stages['initial'][c][k]
                     for k in all_stages['initial'][c]} for c in algebra.STAGES}
        tables['difference'] = algebra.atom_table(delta, kinetics['endpoint']-kinetics['initial'],
            atom_ids=g['atom_ids'], species=g['species'], positions=g['positions'])
        tables['difference']['definition'] = (
            'u,p,v,field and kinetic: endpoint minus initial. Displacement weights: '
            '|u_endpoint-u_initial|^2, not endpoint weight minus initial weight.')
        baselines['drift_endpoint_minus_initial'] = delta['mean']
        comparison = None
        if energy_endpoint is not None:
            # Attach only a small in-memory view for comparison, never a runtime.
            from types import SimpleNamespace
            endpoint = SimpleNamespace(time=last.time, stages=all_stages['endpoint'])
            comparison = energy_comparison(bound_json(energy_endpoint), mixed_binding=mixed_endpoint,
                mean_binding=mean_endpoint, endpoint=endpoint, g=g, kinetic=kinetics['endpoint'])
    sources = [file_binding(m.__file__) for m in (algebra, mixed_excitation_energy,
                                                mixed_response_storage, clamped)]
    sources.append(file_binding(__file__))
    return dict(schema=SCHEMA, complete=True, model_identity=g['model'], q_fractional=[0., 0., 0.],
        source_inputs=dict(mixed_initial=mixed_initial, mixed_endpoint=mixed_endpoint,
            mean_initial=mean_initial, mean_endpoint=mean_endpoint, metadata=metadata,
            energy_endpoint=energy_endpoint), source_bindings=sources,
        original_metadata_identity=inputs['identity'], original_mixed_layout=layouts,
        original_clamped_times_inverse_ry=inputs['times_inverse_ry'],
        actual_read_times_inverse_ry=[first.time, last.time], seed=layouts['initial']['provenance']['original_local_pulse'],
        units=dict(u='bohr', p='bohr^-1', velocity='bohr Ry (d/d time_ry_inverse)',
                   masses='Ry inertial mass', kinetic='global Ry mixed derivative',
                   displacement_weight='bohr^2; magnitude diagnostic only'),
        conventions=dict(stage_order=list(algebra.STAGES),
            derivative='d2/(d lambda d eta), no extra half; independent complex directions, no dagger or reality repair',
            matter='already physical per-cell; no additional sqrt(N) division',
            field='GLOBAL coordinate and momentum; shift=coordinate@polarizations, no sqrt(N) division',
            polarizations=g['polarizations'], polarization_authority='original Cartesian a_xyz,Pi_xyz metadata; genuine energy parameter identity checked when supplied',
            source_positions_shape=[3, 37], response_order='original atom,xyz',
            weights='sum_xyz abs(u_direction)^2; normalized by all37 atoms when nonzero; not covariance',
            kinetic='N*_bilinear(k,k,sum_xyz(a*b/M))/2 by atom; k=p-Z*a'),
        baseline_unexcited_mean=baselines, response=tables, endpoint_energy_comparison=comparison,
        scope=SCOPE, electronic_leg_access_performed=False, trajectory_performed=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration', required=True, type=Path)
    parser.add_argument('--registration-sha256', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    registration_binding = file_binding(args.registration)
    require(registration_binding['sha256'] == args.registration_sha256, 'explicit registration SHA256 differs')
    registration = bound_json(registration_binding)
    require(registration.get('schema') == 'P254-next-physical-readouts-registration-v1', 'parent input registration required')
    bindings = registration['source_inputs']
    result = read_saved_canonical_response(**{key: bindings[key] for key in
        ('mixed_initial', 'mixed_endpoint', 'mean_initial', 'mean_endpoint', 'metadata')},
        energy_endpoint=bindings.get('energy_endpoint'))
    result['registration'] = registration_binding
    require(result['model_identity'] == registration['model_identity'], 'registered original model differs')
    # Exclusive output creation preserves previous original reports. This CLI
    # writes only this requested new report, never any producer artifacts.
    with args.output.open('x') as stream:
        json.dump(clamped.encode(result), stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(file_binding(args.output), sort_keys=True))


if __name__ == '__main__':
    main()

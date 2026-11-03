"""Authenticated actual-result input and compact physical readouts for A29.

No coefficient assembly, Native admission or trajectory replay happens here.
The result's manifest SNAPSHOT, not a subsequently growing producer manifest,
is the authority. Bytes are authenticated before arrays or observables are used.
"""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from clamped_scientific_readout import decoded, verify_nested
from substrate_framework.streamed_material_state import file_binding, authenticate
from reduced_material_readout import original_geometry, spatial_pulse, finite_kinetic_exchange
from reduced_projection_storage import ImmutableItems

CASES = ('loaded-ni31h6', 'pristine-ni32', 'empty-vacancy-ni31', 'dispersed-ni32h6')


def digest(value):
    return sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def physical_context(specification):
    """Implementation revisions do not change the authenticated physical input."""
    keys = ('case', 'aggregate', 'context', 'initial_means', 'positive_time_symmetry',
            'preparation', 'reality_tolerance', 'scalar', 'source_states')
    return {key: specification[key] for key in keys}


def decode(value):
    if isinstance(value, dict):
        if set(value) == {'real', 'imag'}:
            return decoded(value)
        return {key: decode(item) for key, item in value.items()}
    if isinstance(value, list):
        return [decode(item) for item in value]
    return value


def bound_json(binding):
    return json.loads(authenticate(binding).read_text())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def array_archive(binding):
    with np.load(authenticate(binding), allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def compact_spatial(pulse, coordinate):
    displacement = pulse['displacement']
    out = {key: displacement[key] for key in ('definition', 'all_ions', 'Ni', 'H')}
    u = np.asarray(displacement['displacement_complex'])
    out['selected_original_cartesian_coordinate'] = coordinate
    out['selected_displacement_by_cell_complex_bohr'] = u[:, coordinate//3, coordinate % 3]
    return out


def baseline_readout(x, geometry, masses, charges):
    d, n = geometry['d'], geometry['n']
    z = np.repeat(charges, 3)
    u, p = x[:d]/np.sqrt(n), x[d:2*d]/np.sqrt(n)
    a, pi = x[2*d:2*d+3], x[2*d+3:]
    velocity = (p-z*np.tile(a, d//3))/masses
    return dict(displacement_by_atom_bohr=u.reshape(-1, 3),
                canonical_momentum_by_atom_bohr_inverse=p.reshape(-1, 3),
                velocity_by_atom_bohr_ry=velocity.reshape(-1, 3),
                global_field_coordinate=a, global_field_momentum=pi,
                definition='own nonstationary unexcited baseline; never shared across models')


def correlation_readout(record, n):
    if record is None:
        return dict(available=False, reason='no connected correlation saved at this physical time')
    out = decode(record)
    for label in ('equal_time', 'two_time'):
        components = out[label]
        for key in ('ordered_ab', 'ordered_ba', 'symmetrized', 'commutator'):
            array = np.asarray(components[key], dtype=complex)
            require(array.shape == (n, 1, 1) and np.isfinite(array).all(),
                    'saved connected correlation shape/nonfinite: '+label+'.'+key)
            components[key] = array[:, 0, 0]
    out['available'] = True
    out['units'] = 'bohr^2; selected original kick atom/axis displacement'
    return out


def read_finite(path, *, source_reader=None):
    """Read complete boundaries only, reporting an incomplete run as incomplete.

    SHA binding protects rho as well as x; only one selected-time archive is
    resident at a time. No missing q/spin sector or missing time is synthesized.
    """
    result_binding = file_binding(Path(path))
    result = bound_json(result_binding)
    require(result.get('schema') == 'P254-activeband-material-run-v1', 'actual A29 material run required')
    config = result['configuration']
    require(digest(config) == result['identity'], 'finite result/configuration identity differs')
    spec = config['coefficient_specification']
    require(spec['case'] in CASES and not config['baseline_only'], 'original finite-dose case required')
    require(config['bands'] in (8, 12) and config['features'] in (7, 27), 'nested8/12 and7/27 required')
    manifest = bound_json(result['original_coefficient_manifest_snapshot'])
    original_manifest = result['original_input_manifest']
    snapshot_binding = result['original_coefficient_manifest_snapshot']
    require(original_manifest['sha256'] == snapshot_binding['sha256'] and
            int(original_manifest['bytes']) == int(snapshot_binding['bytes']),
            'canonical input manifest differs from the actual saved snapshot')
    authenticate(original_manifest)
    require(manifest['specification'] == spec and manifest['schema'] == 'p254-physical-activeband-cache-v1',
            'result/manifest source specification differs')
    require(manifest['status'] == 'complete', 'finite spatial result requires complete all-q coefficient snapshot')
    source_reader = ImmutableItems({}) if source_reader is None else source_reader
    for item in manifest['items'].values():
        source_reader.authenticate(item['archive'])
    require(spec['case_identity'] == digest(spec['context']), 'original matched preparation identity differs')
    verify_nested(spec['context'])
    authenticate(spec['original_boundary'])
    # A saved source snapshot can replace a changed code path only when its
    # original digest/size are explicitly the same. Historical raw data are
    # never silently rebound to a current file.
    snapshots = {}
    def collect(value):
        if isinstance(value, dict):
            if {'original_binding', 'archived_snapshot'} <= value.keys():
                old, saved = value['original_binding'], value['archived_snapshot']
                require(old['sha256'] == saved['sha256'] and int(old['bytes']) == int(saved['bytes']),
                        'historical source snapshot differs')
                authenticate(saved)
                snapshots[(old['path'], old['sha256'], int(old['bytes']))] = saved
            for item in value.values():
                collect(item)
        elif isinstance(value, list):
            for item in value:
                collect(item)
    collect(manifest)
    collect(result)
    for source in spec['implementation']+spec['physical_implementation']+result['consumer_sources']:
        key = (source['path'], source['sha256'], int(source['bytes']))
        authenticate(snapshots.get(key, source))
    indexing = spec['indexing']
    cells = np.indices((4, 4, 4)).reshape(3, -1).T
    require(spec['mesh'] == [4, 4, 4] and np.array_equal(indexing['q_indices'], cells)
            and np.array_equal(indexing['k_indices'], cells)
            and np.array_equal(indexing['slot_spin'], np.repeat([1, 2], 64))
            and np.array_equal(indexing['slot_k'], np.tile(np.arange(64), 2)),
            'original64 q and128 explicit spin/k slots required')
    for q in range(64):
        expected = np.ravel_multi_index(tuple(((cells+cells[q]) % 4).T), (4, 4, 4))
        require(np.array_equal(indexing['target'][q], np.r_[expected, expected+64]), 'q source/target changed')
        require(indexing['minus'][q] == np.ravel_multi_index(tuple(-cells[q] % 4), (4, 4, 4)), 'minus-q changed')
    for key in ('initial', 'mechanical', 'passive/%02d' % config['bands']):
        authenticate(manifest['items'][key]['archive'])
    for q in range(64):
        require('kernel/%02d' % q in manifest['items'], 'missing physical q kernel')
        for slot in range(128):
            require('edge/%02d/%03d' % (q, slot) in manifest['items'], 'missing q/spin/source edge')
    for slot in range(128):
        require('slot/%03d' % slot in manifest['items'] and 'preparation/%03d' % slot in manifest['items'],
                'missing original source preparation')
    geometry_arrays = array_archive(config['geometry_archive'])
    geometry = original_geometry(geometry_arrays, (4, 4, 4))
    aggregate = bound_json(spec['context']['aggregate'])
    with np.load(authenticate(aggregate['archive']), allow_pickle=False) as original:
        for key in ('atom_indices', 'atom_names', 'cell_bohr', 'valence_charges'):
            require(np.array_equal(geometry_arrays[key], original[key]),
                    'bound geometry differs from original matched aggregate: '+key)
    mechanical = array_archive(manifest['items']['mechanical']['archive'])
    d, m = geometry['d'], config['bands']
    require(mechanical['masses'].shape == (d,), 'original mechanical mass ordering differs')
    selected = np.flatnonzero(geometry['ids'] == config['atom'])
    require(len(selected) == 1 and config['axis'] in 'xyz', 'original seed atom/axis absent')
    coordinate = 3*int(selected[0])+'xyz'.index(config['axis'])
    requested = config['times']
    require(requested == sorted(set(requested)) and np.isfinite(requested).all(), 'unique finite ordered physical times required')
    require(len(result['completed']) <= len(requested), 'unexpected completed times')
    observations = []
    for index, record in enumerate(result['completed']):
        t = record['time_inverse_ry']
        require(t == requested[index], 'completed physical-time prefix differs')
        baseline = array_archive(record['baseline_state'])
        require(set(baseline) == {'rho', 'x'} and baseline['rho'].shape == (128, m, m)
                and baseline['x'].shape == (2*d+6,), 'saved baseline layout differs')
        finite = array_archive(record['finite_state'])
        require(set(finite) == {'rho'} | {'x_%02d' % q for q in range(64)}
                and finite['rho'].shape == (2, 64*m, 64*m), 'saved finite q/spin layout differs')
        for arrays in (baseline, finite):
            require(all(np.isfinite(v).all() for v in arrays.values()), 'nonfinite saved state')
        delta = tuple(finite['x_%02d' % q] for q in range(64))
        require(all(x.shape == (2*d+(6 if q == 0 else 0),) for q, x in enumerate(delta)), 'finite q coordinate layout differs')
        spatial = spatial_pulse(delta, geometry)
        kinetic = finite_kinetic_exchange(baseline['x'], delta, geometry, mechanical['masses'],
                                          geometry_arrays['valence_charges'], float(mechanical['field_mass']))
        energy = decode(record['energy_change_ry'])
        require(set(energy) == {'mechanical', 'band', 'interaction', 'screening', 'total'}, 'saved scalar energy components differ')
        require(all(np.shape(value) == () and np.isfinite(value) for value in energy.values()),
                'finite complex scalar energies required')
        require(('connected_correlation' in record) == config['correlations'],
                'saved correlation/configuration coverage differs')
        if config['correlations']:
            require(record['connected_correlation']['reference_time_inverse_ry'] ==
                    min(t, config['correlation_reference']), 'actual connected reference time differs')
        observations.append(dict(time_inverse_ry=t, baseline_state=record['baseline_state'],
            finite_state=record['finite_state'], spatial=compact_spatial(spatial, coordinate),
            energy_change_ry=energy, finite_kinetic_exchange=kinetic,
            baseline=baseline_readout(baseline['x'], geometry, mechanical['masses'], geometry_arrays['valence_charges']),
            connected_correlation=correlation_readout(record.get('connected_correlation'), 64)))
    require(result['status'] != 'complete' or len(observations) == len(requested), 'false complete physical-time coverage')
    # Only source-item identities, not timings or append-only manifest metadata,
    # participate in the nested-result comparison key.
    source_items = {key: item['archive'] for key, item in manifest['items'].items()
                    if key.startswith(('slot/', 'preparation/', 'edge/', 'kernel/')) or key == 'mechanical'}
    return dict(binding=result_binding, identity=result['identity'], configuration=config,
        manifest_binding=result['original_coefficient_manifest_snapshot'],
        manifest_path=original_manifest['path'],
        specification=spec, manifest=manifest, geometry=geometry, geometry_arrays=geometry_arrays,
        mechanical_arrays=mechanical, coordinate=coordinate, observations=observations,
        model_source_identity=digest(dict(specification=spec, source_items=source_items)),
        status=result['status'], missing_times_inverse_ry=requested[len(observations):],
        seed=dict(atom=config['atom'], axis=config['axis'], cell=[0, 0, 0], amplitude=config['amplitude']))

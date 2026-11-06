"""Read authenticated matched-clamped observations, never evolve a solver.

Spatial weights use displacement ONLY. Maxwell energy is free electric kinetic
energy, not the contact/diamagnetic matter term. Complex defects remain visible.
"""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from substrate_framework.spatial_readout_metrics import periodic_distances, selected_weight
from substrate_framework.physical_field_environment import UniformMaxwellField

HERE = Path(__file__).resolve().parent
REGISTRATION = HERE/'clamped-scientific-readout-registration.json'
ION_ORDER = 'u(original atoms,xyz),p(original atoms,xyz)'
FIELD_ORDER = 'a_xyz,Pi_xyz'
LIMITS = [
    'Clamped electrons: not full coupled lattice dynamics or a material ranking.',
    'Mean coherence/participation is not a connected quantum correlation or entanglement witness.',
    'Connected ordered correlations are selected original atom/axis pairs only; no missing full covariance is synthesized.',
    'Extent is cell-resolved displacement weight about the kick cell in the original finite periodic geometry, not infinite-volume localization or atom-offset extent.',
    'Independent regions remove inter-region ionic coupling AND replace one global field by N independent smaller-volume fields and vacua.',
    'Maxwell-energy change is reversible electromagnetic exchange within the closed model, not loss, heat, dissipation or an energy current.',
    'No numerical-accuracy, time-convergence, lifetime or asymptotic-persistence certification is supplied.'
]


def binding(path):
    path = Path(path).resolve()
    digest = sha256()
    size = 0
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b''):
            digest.update(chunk)
            size += len(chunk)
    return dict(path=str(path), bytes=str(size), sha256=digest.hexdigest())


def authenticate(record):
    actual = binding(record['path'])
    if actual['sha256'] != record['sha256'] or int(actual['bytes']) != int(record['bytes']):
        raise ValueError('raw binding differs: '+record['path'])
    return Path(record['path'])


def bound_json(record):
    return json.loads(authenticate(record).read_text())


def verify_nested(value, seen=None):
    """Authenticate the original nested receipts, not replacement provenance."""
    seen = set() if seen is None else seen
    if isinstance(value, dict):
        if {'path', 'bytes', 'sha256'} <= value.keys():
            key = (value['path'], str(value['bytes']), value['sha256'])
            if key not in seen:
                authenticate(value)
                seen.add(key)
        else:
            for child in value.values():
                verify_nested(child, seen)
    elif isinstance(value, list):
        for child in value:
            verify_nested(child, seen)


def sources():
    registration = json.loads(REGISTRATION.read_text())
    for record in registration['settled_sources']:
        authenticate(record)
    return registration['settled_sources'] + [binding(REGISTRATION)] + [
        binding(HERE/name) for name in registration['new_files']]


def encode(value):
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if isinstance(value, np.ndarray):
        return encode(value.tolist())
    if isinstance(value, (complex, np.complexfloating)):
        return dict(real=float(value.real), imag=float(value.imag))
    if isinstance(value, np.generic):
        return value.item()
    return value


def decoded(value):
    """Producer encodes even string/integer arrays as real/imag objects."""
    if isinstance(value, dict) and set(value) == {'real', 'imag'}:
        real, imag = np.asarray(value['real']), np.asarray(value['imag'])
        if real.shape != imag.shape:
            raise ValueError('real/imag shape differs')
        if real.dtype.kind in 'US':
            if imag.dtype.kind not in 'US' or np.any(imag != ''):
                raise ValueError('invalid string-array imaginary metadata')
            return real
        return real + 1j*imag
    return np.asarray(value)


def real_metadata(value, name, integer=False):
    array = decoded(value)
    if not np.isfinite(array).all() or np.any(np.imag(array) != 0):
        raise ValueError('non-real/nonfinite metadata: '+name)
    array = np.real(array)
    if integer:
        if np.any(array != np.floor(array)):
            raise ValueError('noninteger metadata: '+name)
        return array.astype(int)
    return array


def validate_identity(inputs, result):
    identity = inputs['identity']
    if inputs.get('schema') != 'P254-matched-clamped-inputs-v1':
        raise ValueError('original matched clamped input schema required')
    if identity['phase_space_order'] != ION_ORDER or identity['field_order'] != FIELD_ORDER:
        raise ValueError('original ionic/field phase-space order required')
    if (identity['case'] != result['case'] or identity['control'] != result['control']
            or identity['mass']['variant'] != result['mass']
            or identity.get('electrons_held') is not True
            or identity.get('extra_delta_curvature') != 'exact zero'):
        raise ValueError('result/input physical identity mismatch')
    specification = identity['matched_joint_specification']
    if sha256(json.dumps(specification, sort_keys=True).encode()).hexdigest() != identity['matched_joint_identity']:
        raise ValueError('original matched joint specification identity differs')
    if specification['aggregate'] != identity['aggregate'] or specification['case'] != result['case']:
        raise ValueError('aggregate/case identity mismatch')
    ids = real_metadata(identity['atom_indices'], 'atom indices', integer=True)
    names = decoded(identity['atom_names'])
    cells = real_metadata(identity['cell_indices'], 'cell indices', integer=True)
    cell = real_metadata(identity['cell_bohr'], 'cell geometry')
    if ids.ndim != 1 or len(set(ids)) != len(ids) or names.shape != ids.shape or names.dtype.kind not in 'US':
        raise ValueError('unique original atom identities and species required')
    if cells.ndim != 2 or cells.shape[1] != 3 or not len(cells) or np.any(cells < 0):
        raise ValueError('complete original cell indices required')
    mesh = cells.max(axis=0)+1
    if not np.array_equal(cells, np.indices(tuple(mesh)).reshape(3, -1).T):
        raise ValueError('missing/reordered original regions')
    anchor = real_metadata(inputs['kick_cell'], 'kick cell', integer=True)
    if anchor.shape != (3,) or np.any(anchor < 0) or np.any(anchor >= mesh):
        raise ValueError('kick cell outside original mesh')
    where = np.flatnonzero(ids == inputs['kick_atom_original_one_based'])
    if len(where) != 1 or inputs['kick_axis'] not in ('x', 'y', 'z'):
        raise ValueError('original kick atom/axis absent')
    coordinate = 3*int(where[0])+'xyz'.index(inputs['kick_axis'])
    if inputs['readout_phase_space_indices'] != [coordinate, 3*len(ids)+coordinate]:
        raise ValueError('selected u,p ordering mismatch')
    control = result['control']
    if control not in ('periodic', 'independent-regions'):
        raise ValueError('unsupported clamped control')
    field_index = ('singleton global field' if control == 'periodic' else
                   'original cell order, one independent field per region')
    if (inputs['field_index'] != field_index or inputs['field_covariance_pair'] !=
            'field[index](t), field at kick region(s); singleton for global field'):
        raise ValueError('field indexing/anchor covariance authority differs')
    return dict(ids=ids, names=names, cells=cells, cell=cell, mesh=mesh, anchor=anchor,
                anchor_index=int(np.ravel_multi_index(tuple(anchor), tuple(mesh))),
                n=len(cells), d=3*len(ids), control=control)


def maxwell_authority(inputs, geometry):
    preparation = inputs['identity']['preparation']
    if 'field_mass' not in preparation:
        return dict(available=False, gap='identity.preparation.field_mass absent; cannot authenticate saved field normalization')
    cell = geometry['cell']
    if cell.shape != (3, 3):
        raise ValueError('3x3 original cell required')
    volume = float(abs(np.linalg.det(cell)))
    count = 1 if geometry['control'] == 'periodic' else geometry['n']
    field_volume = volume*(geometry['n'] if count == 1 else 1)
    field = UniformMaxwellField(field_volume, np.eye(3), 'authenticated saved matched-clamped cell and control')
    mass = float(real_metadata(preparation['field_mass'], 'field mass'))
    if not np.isclose(mass, field.mass, rtol=1e-12, atol=0):
        raise ValueError('saved field mass disagrees with physical Maxwell volume/(8*pi)')
    return dict(available=True, field_count=count, volume_per_field_bohr3=field_volume,
                mass_per_field=mass, cell_volume_bohr3=volume,
                normalization='sum_f [Pi_mean_f.T Pi_mean_f + tr(W_PiPi_f)]/(2 M_f); M_f=V_f/(8 pi)',
                independent_covariance_authority='identical uncoupled local evolution/preparation in bound producer; anchor self covariance repeated N times, not sum of cross-region anchor covariances')


def validate_observation(arrays, geometry):
    n, d = geometry['n'], geometry['d']
    fields = 1 if geometry['control'] == 'periodic' else n
    shapes = dict(baseline_ions=(*geometry['mesh'], 2*d), baseline_field=(fields, 6),
                  local_kick_delta_ions=(*geometry['mesh'], 2*d), common_kick_delta_ions=(*geometry['mesh'], 2*d),
                  local_kick_delta_field=(fields, 6), common_kick_delta_field=(fields, 6),
                  spatial_equal_time_ordered=(n, 2, 2), spatial_two_time_ordered=(n, 2, 2),
                  spatial_two_time_commutator=(n, 2, 2), field_equal_time_ordered=(fields, 6, 6),
                  field_two_time_ordered=(fields, 6, 6), ion_field_two_time_ordered=(n, 2, 6))
    if set(arrays) != set(shapes):
        raise ValueError('missing or unexpected saved observables')
    for name, shape in shapes.items():
        if arrays[name].shape != shape or not np.isfinite(arrays[name]).all():
            raise ValueError('saved observable shape/nonfinite: '+name)
    if fields > 1:
        other = np.arange(n) != geometry['anchor_index']
        for name in ('field_equal_time_ordered', 'field_two_time_ordered',
                     'spatial_equal_time_ordered', 'spatial_two_time_ordered',
                     'spatial_two_time_commutator', 'ion_field_two_time_ordered'):
            if np.any(arrays[name][other] != 0):
                raise ValueError('independent-region anchor cross covariance is not zero: '+name)


def displacement_readout(ions, geometry, distances):
    """L2 displacement norm per cell, species and original atom; never momentum."""
    u = ions.reshape(geometry['n'], 2*geometry['d'])[:, :geometry['d']].reshape(geometry['n'], -1, 3)
    def subset(mask):
        if not np.any(mask):
            return dict(applicable=False, reason='species absent in actual saved atom_names')
        amplitudes = np.sqrt(np.sum(np.abs(u[:, mask, :])**2, axis=(1, 2)))
        return dict(applicable=True, original_atom_indices=geometry['ids'][mask],
                    metric=selected_weight(amplitudes, distances, value_units='bohr'))
    return dict(displacement_complex=u, definition='w_R=sum_selected_atoms,xyz |<u_Ratom,xyz>|^2; no inter-atom phase cancellation',
                all_ions=subset(np.ones(len(geometry['ids']), dtype=bool)),
                Ni=subset(geometry['names'] == 'Ni'), H=subset(geometry['names'] == 'H'),
                by_original_atom={str(atom): subset(geometry['ids'] == atom)
                                  for atom in geometry['ids']})


def maxwell_energy(mean, anchor_covariances, geometry, authority):
    if not authority['available']:
        return authority.copy()
    fields = authority['field_count']
    if mean.shape != (fields, 6) or anchor_covariances.shape != (fields, 6, 6):
        raise ValueError('missing physical field/region mean or covariance')
    anchor = 0 if geometry['control'] == 'periodic' else geometry['anchor_index']
    if fields > 1 and np.any(anchor_covariances[np.arange(fields) != anchor] != 0):
        raise ValueError('expected saved anchor cross covariance, not local covariance stack')
    covariance = anchor_covariances[anchor, 3:, 3:]
    mass = authority['mass_per_field']
    # Ordered Pi Pi has no CCR antisymmetric term; keep all numerical complex defects.
    mean_by_field = np.einsum('fi,fi->f', mean[:, 3:], mean[:, 3:])/(2*mass)
    covariance_per_field = np.trace(covariance)/(2*mass)
    by_field = mean_by_field+covariance_per_field
    return dict(available=True, mean_energy_by_field_ry=mean_by_field,
                covariance_energy_per_field_ry=covariance_per_field,
                covariance_multiplicity=fields, energy_by_field_ry=by_field,
                mean_energy_ry=np.sum(mean_by_field), covariance_energy_ry=fields*covariance_per_field,
                energy_ry=np.sum(by_field), energy_imaginary_defect_ry=float(abs(np.imag(np.sum(by_field)))),
                field_mean_complex=mean, anchor_PiPi_ordered_complex=covariance)


def readout_time(arrays, receipt, inputs, geometry, authority, distances):
    validate_observation(arrays, geometry)
    records = {}
    baseline = maxwell_energy(arrays['baseline_field'], arrays['field_equal_time_ordered'], geometry, authority)
    for name in ('baseline', 'local', 'common'):
        ions = arrays['baseline_ions']
        field = arrays['baseline_field']
        if name != 'baseline':
            ions = ions+arrays[name+'_kick_delta_ions']
            field = field+arrays[name+'_kick_delta_field']
        energy = maxwell_energy(field, arrays['field_equal_time_ordered'], geometry, authority)
        if energy['available']:
            energy['change_from_same_time_baseline_ry'] = energy['energy_ry']-baseline['energy_ry']
        records[name] = dict(absolute_displacement=displacement_readout(ions, geometry, distances), maxwell=energy)
        if name != 'baseline':
            records[name]['kick_delta_displacement'] = displacement_readout(arrays[name+'_kick_delta_ions'], geometry, distances)
    correlations = {key: value for key, value in arrays.items() if 'ordered' in key or 'commutator' in key}
    weights = {}
    for name in ('spatial_equal_time_ordered', 'spatial_two_time_ordered', 'spatial_two_time_commutator'):
        weights[name] = {label: selected_weight(arrays[name][:, i, i], distances, value_units=units)
                         for i, label, units in ((0, 'uu', 'bohr^2'), (1, 'pp', 'bohr^-2'))}
    return dict(time_inverse_ry=receipt['time_inverse_ry'], reference_time_inverse_ry=receipt['reference_time_inverse_ry'],
                archive=receipt['archive'], states=records, connected_ordered_complex=correlations,
                selected_correlation_spatial_weights=weights, producer_diagnostics=receipt['diagnostics'])


def load_readout(path, expected_sha):
    record = binding(path)
    if record['sha256'] != expected_sha:
        raise ValueError('explicit original result digest differs')
    result = bound_json(record)
    if (result.get('schema') != 'P254-matched-clamped-controls-v1' or result.get('status') != 'complete'
            or result.get('electrons_held') is not True):
        raise ValueError('completed original matched-clamped result required')
    inputs = bound_json(result['inputs'])
    if result['source_bindings'] != inputs['source_bindings'] or result['oracle'] != inputs['oracle']:
        raise ValueError('producer result/input provenance differs')
    registered = json.loads(REGISTRATION.read_text())['producer_source_bindings']
    if result['source_bindings'] != registered:
        raise ValueError('only registered frozen matched-clamped producer accepted')
    verify_nested(inputs)
    oracle = bound_json(result['oracle'])
    if oracle.get('status') != 'passed' or oracle.get('source_bindings') != result['source_bindings']:
        raise ValueError('original passed producer oracle/source bindings differ')
    geometry = validate_identity(inputs, result)
    aggregate = bound_json(inputs['identity']['aggregate'])
    if aggregate['archive'] != inputs['identity']['archive']:
        raise ValueError('original aggregate archive differs')
    with np.load(authenticate(aggregate['archive']), allow_pickle=False) as original:
        for key in ('atom_indices', 'atom_names', 'cell_bohr', 'positions_bohr', 'valence_charges'):
            if not np.array_equal(decoded(inputs['identity'][key]), original[key]):
                raise ValueError('saved identity does not match original aggregate: '+key)
        if tuple(original['common_force_lattice_ry_bohr2'].shape[:3]) != tuple(geometry['mesh']):
            raise ValueError('saved regions do not cover original force lattice')
        if not np.array_equal(decoded(inputs['identity']['mass']['original_masses_amu']), original['atomic_masses_amu']):
            raise ValueError('original masses differ from bound aggregate')
    authority = maxwell_authority(inputs, geometry)
    distances = periodic_distances(geometry['cells'], geometry['anchor'], geometry['mesh'], geometry['cell'])
    observations = []
    times = inputs['times_inverse_ry']
    if len(times) != len(result['completed_times']) or len(set(times)) != len(times):
        raise ValueError('completed time coverage/uniqueness differs')
    for time, receipt_binding in zip(times, result['completed_times']):
        receipt = bound_json(receipt_binding)
        if (receipt['status'] != 'complete' or receipt['time_inverse_ry'] != time
                or receipt['reference_time_inverse_ry'] != inputs['reference_time_inverse_ry']):
            raise ValueError('saved observation time identity differs')
        with np.load(authenticate(receipt['archive']), allow_pickle=False) as archive:
            arrays = {key: archive[key] for key in archive.files}
        metadata = {key: dict(shape=list(value.shape), dtype=str(value.dtype)) for key, value in arrays.items()}
        if metadata != receipt['arrays']:
            raise ValueError('receipt/NPZ observable metadata differs')
        observation = readout_time(arrays, receipt, inputs, geometry, authority, distances)
        observation['receipt_binding'] = receipt_binding
        observations.append(observation)
    initial = next((item for item in observations if item['time_inverse_ry'] == 0), None)
    for item in observations:
        for name, state in item['states'].items():
            energy = state['maxwell']
            energy['change_from_initial_ry'] = (energy['energy_ry']-initial['states'][name]['maxwell']['energy_ry']
                if initial is not None and energy['available'] else None)
            if initial is None:
                energy['initial_change_gap'] = 'no saved time-zero observation; no synthetic initial observable substituted'
    initial_energies = {key: decoded(value).item() for key, value in inputs['initial_total_energy_ry'].items()}
    return dict(schema='P254-clamped-scientific-readout-v1', status='complete', original_result=record,
                original_input=result['inputs'], producer_source_bindings=result['source_bindings'],
                original_identity=inputs['identity'], case=result['case'], mass=result['mass'], control=result['control'],
                kick=dict(cell=inputs['kick_cell'], atom_original_one_based=inputs['kick_atom_original_one_based'],
                    axis=inputs['kick_axis'], local='delta p=sqrt(2 M omega) at one cell',
                    common='same delta p/sqrt(N) at EVERY cell; same total added energy, not same per-cell amplitude',
                    nominal_added_energy_ry=inputs['kick_added_energy_ry'],
                    actual_initial_added_energy_ry={name: initial_energies[name]-initial_energies['baseline'] for name in ('local', 'common')},
                    initial_total_energy_ry=inputs['initial_total_energy_ry']),
                field_authority=authority, geometry=dict(cell_indices=geometry['cells'], cell_bohr=geometry['cell'],
                    mesh=geometry['mesh'], atom_indices=geometry['ids'], atom_names=geometry['names']),
                observations=observations, limits=LIMITS, whole_coupled_lattice_ranking=False)


def compare_pair(left, right):
    if left['control'] != 'periodic' or right['control'] != 'independent-regions':
        raise ValueError('periodic left and independent-regions right required')
    for key in ('case', 'mass', 'kick', 'geometry'):
        if key == 'kick':
            keys = ('cell', 'atom_original_one_based', 'axis', 'nominal_added_energy_ry')
            a, b = ({k: value[key][k] for k in keys} for value in (left, right))
        else:
            a, b = left[key], right[key]
        if encode(a) != encode(b):
            raise ValueError('comparison case/mass/seed/geometry differs: '+key)
    for key in ('matched_joint_identity', 'aggregate', 'archive', 'mass'):
        if left['original_identity'][key] != right['original_identity'][key]:
            raise ValueError('comparison original physical identity differs: '+key)
    lt = {(x['time_inverse_ry'], x['reference_time_inverse_ry']): x for x in left['observations']}
    rt = {(x['time_inverse_ry'], x['reference_time_inverse_ry']): x for x in right['observations']}
    common = sorted(set(lt) & set(rt))
    if not common:
        raise ValueError('no same-time/reference-time completed observations')
    rows = []
    for time in common:
        states = {}
        for name in ('baseline', 'local', 'common'):
            a, b = lt[time]['states'][name], rt[time]['states'][name]
            metric_key = 'absolute_displacement' if name == 'baseline' else 'kick_delta_displacement'
            groups = {group: dict(periodic=a[metric_key][group], independent_regions=b[metric_key][group])
                      for group in ('all_ions', 'Ni', 'H')}
            states[name] = dict(displacement=groups, periodic_maxwell=a['maxwell'], independent_regions_maxwell=b['maxwell'])
            if a['maxwell']['available'] and b['maxwell']['available']:
                states[name]['independent_minus_periodic_maxwell_ry'] = b['maxwell']['energy_ry']-a['maxwell']['energy_ry']
        rows.append(dict(time_inverse_ry=time[0], reference_time_inverse_ry=time[1], states=states))
    return dict(schema='P254-clamped-bound-pair-comparison-v1', original_results=[left['original_result'], right['original_result']],
                case=left['case'], mass=left['mass'], matched_times=rows,
                unmatched_periodic_times=sorted(set(lt)-set(rt)), unmatched_independent_times=sorted(set(rt)-set(lt)),
                periodic_field_authority=left['field_authority'], independent_field_authority=right['field_authority'],
                changed_control='one field at N Vcell versus N fields each Vcell, distinct vacua plus removed inter-region ionic coupling',
                whole_coupled_lattice_ranking=False, limits=LIMITS)

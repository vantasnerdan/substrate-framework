"""Source-bound saved-control payload accounting; no physical model imports.

Run on CPU7 with BLAS thread environment set to one. Only raw authentication,
JSON metadata and selected NPZ member headers are consumed. Original frame
bindings are recorded, never opened. This is not a filesystem peak or allocation.
"""
import ast
import copy
from hashlib import sha256
import json
import os
from pathlib import Path
import struct
import sys
import zipfile

HERE = Path(__file__).resolve().parent
DIAG = Path('/var/tmp/p254-a28-control-payload-plan-001')
OUTPUT = HERE / 'saved-control-mean-payload-plan-001.json'
INITIALIZATIONS = [
    ('pristine-joint-initialize-002', '25654',
     '1f733a1f56c48cb24fe0438357d0a09cd82a8835cb316668762a23001abc61a1'),
    ('vacancy-joint-initialize-001', '25784',
     '9b032f6b703c700065bf7e5015535ef2367ab01973528b4d337d68ed015db27c'),
    ('dispersed-joint-initialize-001', '33376',
     '12c767ecd67d62afd330d381dc4e3dfbe3920a42d19e5b37de1d73d7e0d9750e'),
]
CONVENTION = dict(path=str(HERE / 'saved_case_mean_execution.py'), bytes='6274',
    sha256='06a1561dbe8fa46658c5667897dbacdd7dee2cd07c32b33fcaa8aa8531ad81f5')
DIMENSIONS = dict(path=str(HERE / 'matched_joint_case.py'), bytes='8457',
    sha256='bdda2e9ba4e75a497eed928443b1c3bb85e22c75874376e35c94b4d24e55c15c')
SLOT_STORAGE = HERE.parents[3] / 'src/substrate_framework/partitioned_joint_slots.py'
LOADED = dict(path=str(HERE / 'matched-working-precision-mean-complete-001.json'),
    bytes='10063', sha256='3048fc8e630facf9403a8f4d58a5cf2d9b960416a736df8e4b7b235a12966b75')
PLAN = dict(path=str(DIAG / 'block-plan-001.json'), bytes='1830',
    sha256='a3f34fab917036f5476bf3ea24fc83cbbe3caf1d83aa269a9b6c1c952cd95c25')
INVENTORY = {}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def raw_binding(path):
    path = Path(path)
    digest, size = sha256(), 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
            size += len(block)
    return dict(path=str(path), bytes=str(size), sha256=digest.hexdigest())


def authenticate(binding):
    actual = raw_binding(binding['path'])
    require(actual['sha256'] == binding['sha256'], 'raw SHA256 mismatch: ' + actual['path'])
    require(actual['bytes'] == str(binding['bytes']), 'raw byte-count mismatch: ' + actual['path'])
    INVENTORY[actual['path']] = actual
    return Path(actual['path'])


def json_source(binding):
    return json.loads(authenticate(binding).read_bytes())


def normalized(binding):
    return dict(path=binding['path'], bytes=str(binding['bytes']), sha256=binding['sha256'])


def same_binding(left, right):
    return normalized(left) == normalized(right)


def header(stream):
    require(stream.read(6) == b'\x93NUMPY', 'NPY magic missing')
    version = tuple(stream.read(2))
    require(version in ((1, 0), (2, 0), (3, 0)), 'unsupported NPY header version')
    width = 2 if version == (1, 0) else 4
    size = struct.unpack('<H' if width == 2 else '<I', stream.read(width))[0]
    require(size < 65536, 'unexpectedly large NPY metadata header')
    data = ast.literal_eval(stream.read(size).decode('utf8' if version == (3, 0) else 'latin1').strip())
    return dict(shape=list(data['shape']), dtype=data['descr'],
                fortran_order=data['fortran_order'])


def selected_headers(archive):
    # No array values are loaded, including the density array.
    names = ('saved_density_bohr3', 'atomic_masses_amu', 'positions_bohr', 'cell_bohr')
    with zipfile.ZipFile(archive['path']) as members:
        result = {}
        for name in names:
            with members.open(name + '.npy') as stream:
                result[name] = header(stream)
    return result


def validate_slots(slots):
    require(len(slots) == 128, 'all 128 original slot rows required')
    require(len({slot['identity'] for slot in slots}) == 128, 'duplicate original slot identity')
    for slot in slots:
        shape = slot['shape']
        require(len(shape) == 2, 'original slot needs row and column dimensions')
        require(all(type(x) is int and x > 0 for x in shape), 'positive integer row/column required')
        require(all(k in slot['frame'] for k in ('path', 'sha256', 'bytes')), 'frame binding missing')


def collect(initialization, geometry=True):
    init = json_source(initialization)
    require(init['status'] == 'complete' and init['mode'] == 'initialize', 'completed original initialization required')
    saved = init['saved']
    boundary = json_source(saved['original_boundary'])
    joint = json_source(saved['joint_boundary'])
    authenticate(saved['preparation'])
    require(boundary['model_identity'] == saved['original_model_identity'], 'original model identity differs')
    require(joint['model_identity'] == saved['model_identity'], 'saved joint model identity differs')
    require(boundary['time_ry_inverse'] == 0, 'original t0 boundary required')
    require(saved['case'] == init['case'] == boundary['evidence']['case'], 'case provenance differs')
    specification = saved['specification']
    old_specification = boundary['evidence'].get('specification')
    if old_specification is not None:
        require(same_binding(specification['aggregate'], old_specification['aggregate']), 'original aggregate changed')
        require(specification['context'] == old_specification['context'], 'original context changed')
    else:
        require(not geometry, 'control original boundary specification missing')
    require(joint['evidence']['specification'] == specification, 'joint saved specification differs')
    json_source(specification['registration'])
    for binding in specification['source_states'].values():
        json_source(binding)
    context = specification['context']
    reference = json_source(context['reference'])
    json_source(context['preparation'])
    json_source(context['background_receipt'])
    authenticate(context['population_ledger'])
    aggregate = json_source(specification['aggregate'])
    require(aggregate['context'] == context, 'aggregate/context binding differs')
    require(aggregate['status'] == 'complete' and aggregate['slot_count'] == 128, 'complete original aggregate required')
    authenticate(aggregate['archive'])
    headers = selected_headers(aggregate['archive'])
    density, mass = headers['saved_density_bohr3'], headers['atomic_masses_amu']
    require(len(density['shape']) == 2 and density['shape'][1] == 2, 'saved spin-density dimensions differ')
    require(len(mass['shape']) == 1, 'atomic mass metadata must be a vector')
    require(headers['positions_bohr']['shape'] == [3, mass['shape'][0]], 'original geometry dimensions differ')
    require(headers['cell_bohr']['shape'] == [3, 3], 'original cell dimensions differ')
    validate_slots(boundary['slots'])
    # Grounded by matched_joint_case.joint_blocks: np.repeat(atomic masses, 3).
    ndof = 3 * mass['shape'][0]
    require(ndof == len(boundary['u_bohr']) == len(boundary['p_canonical']), 'Cartesian ionic dimensions differ')
    require(ndof == len(joint['u_bohr']), 'joint/original ionic dimensions differ')
    xml = next(b for b in reference['background_sources'] if Path(b['path']).name == 'data-file-schema.xml')
    if geometry:
        authenticate(xml)  # Metadata only; no Native context or XML-derived contraction.
    authority = json.loads(aggregate['ionic_reference']['authority'])
    require(authority['coordinate_identity'] == xml['sha256'] + ':original-atom-order', 'ionic geometry identity differs')
    return dict(case=init['case'], initialization=normalized(initialization),
        original_boundary=normalized(saved['original_boundary']),
        joint_boundary=normalized(saved['joint_boundary']),
        original_model_identity=saved['original_model_identity'], model_identity=saved['model_identity'],
        saved_specification=specification,
        aggregate=normalized(specification['aggregate']), archive=normalized(aggregate['archive']),
        aggregate_headers=headers, nr=density['shape'][0], nat=mass['shape'][0], ndof=ndof,
        density_shape=density['shape'],
        density_size=density['shape'][0] * density['shape'][1],
        density_first_dimension=density['shape'][0],
        density_consumer_expression='context.reference.nbar.shape[0]',
        atomic_mass_size=mass['shape'][0],
        slots=boundary['slots'],
        geometry_provenance=dict(original_xml=normalized(xml), raw_xml_authenticated=geometry,
            coordinate_identity=authority['coordinate_identity'],
            coordinate_order=authority['metadata']['coordinate_order'],
            cell_vectors=authority['metadata']['cell_vectors'],
            positions_member='positions_bohr.npy', cell_member='cell_bohr.npy',
            archive=normalized(aggregate['archive']), geometry_policy=aggregate['geometry_policy']))


def direct_payload(case):
    validate_slots(case['slots'])
    frame = sum(16 * slot['shape'][0] * slot['shape'][1] for slot in case['slots'])
    barrier = len(case['slots']) * 16 * (case['nr'] + case['ndof'] + 5)
    barriers = 1 + 1 * (4 + 2)
    return dict(F=frame, M=barrier, barriers=barriers, forecast=frame + barriers * barrier)


def symbolic_payload(case, sp):
    # Independently encoded entrywise symbolic frame sum, and five storage fields
    # grounded by partitioned_joint_slots: density, gradient, ell, one_body, number.
    rows = sp.symbols('r0:128', integer=True, positive=True)
    columns = sp.symbols('c0:128', integer=True, positive=True)
    n, d = sp.symbols('n d', integer=True, positive=True)
    values = {}
    for i, slot in enumerate(case['slots']):
        values[rows[i]], values[columns[i]] = slot['shape']
    frame_expr = sp.Add(*(sp.Integer(16) * rows[i] * columns[i] for i in range(128)))
    fields = [n, d, sp.Integer(3), sp.Integer(1), sp.Integer(1)]
    barrier_expr = sp.Add(*(sp.Integer(16) * field for _ in range(128) for field in fields))
    frame = frame_expr.subs(values, simultaneous=True)
    barrier = barrier_expr.subs({n: case['nr'], d: case['ndof']})
    # Enumerate the initial barrier, four iteration barriers and two final barriers.
    forecast = sp.Add(frame, *(barrier for _ in ['initial', 'iteration1', 'iteration2',
        'iteration3', 'iteration4', 'post_iteration', 'accepted_final']))
    return dict(F=int(frame), M=int(barrier), barriers=7, forecast=int(forecast)), dict(
        frame_symbolic_sum=str(frame_expr), barrier_decomposed_symbolic_sum=str(barrier_expr),
        oracle='SymPy substitution of 128 independent row/column products; per-slot density/gradient/ell/one_body/number sums; seven explicit barriers')


def rejection(name, function):
    try:
        function()
    except (ValueError, KeyError, FileNotFoundError) as error:
        return dict(mutation=name, detected=True, actual_exception=type(error).__name__, reason=str(error))
    raise ValueError('mutation survived: ' + name)


def equal_forecast(candidate, expected):
    require(candidate == expected, 'payload forecast differs from independent symbolic sum')


def mutations(case, direct, oracle):
    records = []
    missing = dict(case['initialization'], path=str(DIAG / 'deliberately-missing-source.json'))
    records.append(rejection('missing original source', lambda: authenticate(missing)))
    slots = copy.deepcopy(case['slots'])
    slots[0]['shape'] = slots[0]['shape'][1:]
    records.append(rejection('missing original row dimension', lambda: validate_slots(slots)))
    records.append(rejection('missing one of 128 original slots', lambda: validate_slots(case['slots'][:-1])))
    altered = copy.deepcopy(case)
    altered['slots'][0]['shape'][0] += 1
    records.append(rejection('one original row changed', lambda:
        equal_forecast(direct_payload(altered)['forecast'], oracle['forecast'])))
    records[-1]['byte_difference'] = str(direct_payload(altered)['forecast'] - oracle['forecast'])
    records.append(rejection('one contraction barrier omitted', lambda:
        equal_forecast(direct['F'] + 6 * direct['M'], oracle['forecast'])))
    records[-1]['missing_bytes'] = str(direct['M'])
    for name, count in [('density', case['nr']), ('gradient', case['ndof']),
                        ('ell', 3), ('one_body', 1), ('number', 1)]:
        omitted = 7 * 128 * 16 * count
        records.append(rejection('omitted ' + name + ' payload term', lambda omitted=omitted:
            equal_forecast(direct['forecast'] - omitted, oracle['forecast'])))
        records[-1]['missing_bytes'] = str(omitted)
    return records


def string_integers(value):
    if type(value) is int:
        return str(value)
    if isinstance(value, dict):
        return {key: string_integers(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [string_integers(item) for item in value]
    return value


def persist_new(path, value):
    with Path(path).open('x') as stream:
        stream.write(json.dumps(string_integers(value), indent=2, sort_keys=True) + '\n')


def main():
    require(os.sched_getaffinity(0) == {7}, 'CPU7-only affinity required')
    for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
        require(os.environ.get(name) == '1', name + '=1 required')
    json_source(PLAN)
    convention_text = authenticate(CONVENTION).read_text()
    convention_function = next(node for node in ast.parse(convention_text).body
                               if isinstance(node, ast.FunctionDef) and node.name == 'resource_plan')
    convention_return = next(node for node in ast.walk(convention_function) if isinstance(node, ast.Return))
    consumer_expressions = {keyword.arg: ast.get_source_segment(convention_text, keyword.value)
                            for keyword in convention_return.value.keywords}
    authenticate(DIMENSIONS)
    INVENTORY[str(SLOT_STORAGE)] = raw_binding(SLOT_STORAGE)
    INVENTORY[str(Path(__file__).resolve())] = raw_binding(Path(__file__).resolve())
    # Only read and register sources before deriving any payload total.
    cases = [collect(dict(path='/var/tmp/p254-a28-matched-' + suffix + '/result.json',
                         bytes=size, sha256=digest)) for suffix, size, digest in INITIALIZATIONS]
    completed = json_source(LOADED)
    require(completed['status'] == 'complete' and completed['case'] == 'loaded-ni31h6', 'completed loaded comparison required')
    loaded = collect(completed['initialization'], geometry=False)
    for name in ('source-access-001.json', 'source-access-002.json', 'source-access-003.json'):
        path = DIAG / name
        INVENTORY[str(path)] = raw_binding(path)
    inventory_path = DIAG / 'source-inventory-001.json'
    persist_new(inventory_path, dict(schema='P254-control-payload-source-inventory-v1',
        status='all_load_bearing_sources_accessed_before_payload_calculation',
        sources=list(INVENTORY.values()),
        scope='Raw authentication plus JSON and selected NPZ headers only. Frames and coefficients are cited from authenticated original boundary metadata, never opened. Historical implementation binding lists are provenance, not a claim those old versions remain current.'))
    import sympy as sp
    reports = []
    for case in cases:
        direct = direct_payload(case)
        oracle, symbolic = symbolic_payload(case, sp)
        require(direct == oracle, 'direct/symbolic payload disagreement: ' + case['case'])
        records = mutations(case, direct, oracle)
        reports.append(dict(**case, payload=direct, independent_symbolic_check=symbolic,
            direct_equals_symbolic=True, sensitive_mutations=records,
            required_output_margin_forecast_bytes=direct['forecast'],
            original_frame_file_bytes_from_boundary=sum(int(slot['frame']['bytes']) for slot in case['slots']),
            margin_scope='additional modeled output payload for one new full frame set and seven contraction barriers only; excludes already-saved initial frames; not free capacity or total usage'))
    loaded_direct = direct_payload(loaded)
    loaded_oracle, _ = symbolic_payload(loaded, sp)
    require(loaded_direct == loaded_oracle, 'loaded direct/symbolic disagreement')
    plan = completed['resource_plan']
    original_file_bytes = sum(int(slot['frame']['bytes']) for slot in loaded['slots'])
    require(original_file_bytes == int(plan['original_frame_bytes']), 'loaded original-file-byte plan differs')
    require(max(int(slot['frame']['bytes']) for slot in loaded['slots']) == int(plan['max_single_frame_bytes']), 'loaded largest-frame plan differs')
    require(loaded_direct['M'] == int(plan['complex_contraction_payload_bytes_per_complete_barrier']), 'loaded barrier-payload plan differs')
    require(loaded_direct['barriers'] == int(plan['maximum_reduction_barriers_before_cache_reuse']), 'loaded barrier count differs')
    require(plan['retained_new_frame_sets'] == 1, 'loaded retained frame count differs')
    comparison = dict(completed_plan=LOADED, initialization=loaded['initialization'],
        original_boundary=loaded['original_boundary'], aggregate=loaded['aggregate'], archive=loaded['archive'],
        dimensions=dict(nr=loaded['nr'], nat=loaded['nat'], ndof=loaded['ndof']),
        payload=loaded_direct, original_frame_file_bytes=original_file_bytes,
        file_bytes_minus_header_excluded_F=original_file_bytes - loaded_direct['F'],
        checks=dict(symbolic_equals_direct=True, original_file_bytes_match=True,
                    largest_frame_file_bytes_match=True, M_matches=True, seven_barriers_match=True,
                    one_retained_frame_set_matches=True),
        convention_distinction='Existing original_frame_bytes sums recorded NPY file bytes. F counts complex frame elements only; their difference is recorded file overhead, not a new payload term. No frame content or Native context was read.')
    report = dict(schema='P254-saved-control-mean-payload-plan-v1', status='complete',
        scope='source-bound payload-only one-accepted-evolve-step/max_iterations4 forecast',
        sources_before_calculation=raw_binding(inventory_path), focused_block_plan=PLAN,
        convention_source=CONVENTION, ionic_dimension_source=DIMENSIONS,
        contraction_fields_source=INVENTORY[str(SLOT_STORAGE)], implementation=raw_binding(__file__),
        cpu_affinity=[7], blas_threads=1, sympy_version=sp.__version__,
        actual_resource_plan_expressions=consumer_expressions,
        formulas=dict(F='16*sum(rows*columns) over all128 original slots',
                      M='128*16*(nr+ndof+5)', barriers='1+1*(4+2)=7', forecast='F+7*M',
                      nr='saved_density_bohr3.shape[0], not spin-flattened size',
                      ndof='3*atomic_masses_amu.size, independently matched to original u/p lengths; atomic array itself is nat, not ndof',
                      five_scalars='ell[3]+one_body[1]+number[1]'),
        cases=reports, loaded_completed_resource_plan_comparison=comparison,
        consumer_scope=dict(owner='Root; informational input for a later bounded C allocation only',
            excludes=['NPY and ZIP headers', 'JSON receipts', 'stage arrays', 'failed and incomplete attempts',
                      'reservations', 'existing initial data', 'other resource consumers', 'filesystem allocation/peak/RSS'],
            no_native_context=True, no_frame_payload_reads=True, no_new_physical_contractions=True,
            no_trajectory=True, no_free_space_measurement=True, no_resource_grant=True,
            no_fit_guarantee=True, no_accuracy_claim=True, no_design_ranking=True,
            capacity='Previously reported rounded shared45GiB at15:44 only; not remeasured, not converted into an available allocation, not subtracted into an invented total-usage forecast.'))
    persist_new(OUTPUT, report)
    print(json.dumps(dict(status='complete', report=raw_binding(OUTPUT), source_inventory=raw_binding(inventory_path),
        cases=[dict(case=c['case'], F=str(c['payload']['F']), M=str(c['payload']['M']),
                    required_output_margin_forecast_bytes=str(c['payload']['forecast'])) for c in reports],
        loaded_comparison=string_integers(comparison)), indent=2))


if __name__ == '__main__':
    main()

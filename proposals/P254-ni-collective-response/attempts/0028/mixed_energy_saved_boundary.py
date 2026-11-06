"""Read-only admission of an original Native mixed physical boundary.

API: ``admit_saved_boundary(producer_directory, owner_sha256, state_binding,
source_archive=None)`` returns original ``producer_directory, owner_binding,
owner, state_binding, state, boundary_kind, authority, source_archive``.
Owner/state JSON and their opaque integers are never reconstructed. An explicit
archive authenticates historical source bytes independently of current code;
otherwise source files must still match live. All requested state payloads and
non-source physical inputs remain authenticated at their original paths.

Admission performs no publication, RHS, history restoration or propagation.
Initial admission verifies the literal raw-zero physical pulse seed. Endpoint
admission verifies all four completed RK4 stages, every source transition and
actual QWork marker, and the contiguous physical predecessor chain. The driver's
post-energy step receipt is deliberately NOT a prerequisite: stage-3 completion
is the producer's physical boundary barrier. Retired historical payloads require
the original storage retirement journals; the requested boundary must be live.

This authenticates boundary authority, not every unrelated borrowed frame in the
owner's inventory. Runtime constructors still authenticate their actual baseline
inputs. The caller must hold the producer's existing exclusive owner lock for a
positive-time readout, and retain it through consumption. Read-only t0 admission
is permitted while the producer runs. Nothing reopens its budget or clock.
"""
from contextlib import ExitStack
import json
from pathlib import Path

import numpy as np

from matched_joint_case import model_identity
from mixed_material_driver_support import initial_small, metadata_geometry, nested_bindings, time_steps
from mixed_payload_output import exact_zero_payload
from mixed_energy_source_archive import read_producer_sources
from mixed_response_edge_transition import read_transition, transition_marker
from mixed_response_rk4 import (
    SCHEMA as STEP_SCHEMA, STAGE_SCHEMA, STAGE_COMPLETE_SCHEMA, _small_outputs)
from mixed_response_storage import (
    COMPONENTS, SMALL_NAMES, MappedMixedState, _binding, _decode_small, _encode_small,
    payload_path)
from run_joint_q_material import admit_oracle, bindings as existing_dependencies
from substrate_framework.partitioned_joint_slots import digest
from substrate_framework.streamed_material_state import authenticate, file_binding

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OWNER_SCHEMA = 'P254-native-mixed-q-material-v1-owner'
NATIVE_SCHEMA = 'P254-native-mixed-material-prepared-v1'
ADAPTER_MODULES = (
    'mixed_coordinate_response', 'mixed_frame_response', 'mixed_material_edge',
    'mixed_material_moments', 'mixed_dyad_grid_moments', 'mixed_dyad_ionic_moments',
    'mixed_dyad_projector_moments', 'mixed_moment_terms', 'mixed_physical_field_dyads',
    'mixed_response_storage', 'mixed_payload_output',
    'substrate_framework.material_joint_response',
    'substrate_framework.material_physical_operators',
    'substrate_framework.qe_bloch_finiteq_ionic', 'substrate_framework.qe_bloch_finiteq',
    'substrate_framework.streamed_joint_reduction', 'matched_joint_case',
    'substrate_framework.qe_material_reference')


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _plain_path(path):
    path = Path(path)
    for part in (path, *path.parents):
        _require(not part.is_symlink(), 'boundary authority cannot traverse a symlink: '+str(path))
    return path.resolve()


class _Authority:
    """Authenticate original bytes once per binding, retaining original integers."""
    def __init__(self):
        self.bindings = {}
        self.records = {}
        self.source_bytes = {}
        self.source_bindings = {}
        self.source_archive = None

    def recover_sources(self, archive, owner_binding, sources):
        archive = self.bind(archive)
        self.source_bytes = read_producer_sources(archive, owner_binding, sources)
        self.source_bindings = {source['path']: source for source in sources}
        self.source_archive = archive

    def bind(self, binding):
        # Historical initialization authorities use decimal-string byte counts.
        # Validate them with the existing raw contract, but retain their original
        # JSON types in the returned evidence (mixed checkpoints remain strict).
        _require(isinstance(binding, dict) and type(binding.get('bytes')) in (int, str),
                 'original raw authority byte count required')
        _binding(dict(binding, bytes=int(binding['bytes'])))
        binding = dict(binding)
        _plain_path(binding['path'])
        previous = self.bindings.get(binding['path'])
        _require(previous is None or previous == binding, 'conflicting authority binding')
        if previous is None:
            if binding['path'] in self.source_bytes:
                _require(binding == self.source_bindings[binding['path']],
                         'historical source binding differs from exact original inventory')
            else:
                authenticate(binding)
            self.bindings[binding['path']] = binding
        return binding

    def read(self, binding):
        binding = self.bind(binding)
        if binding['path'] not in self.records:
            data = (self.source_bytes[binding['path']] if binding['path'] in self.source_bytes
                    else Path(binding['path']).read_bytes())
            self.records[binding['path']] = json.loads(data)
        return self.records[binding['path']]

    def at(self, path):
        return self.bind(file_binding(_plain_path(path)))


def _source_binding(request, path):
    matches = [b for b in request['source_bindings'] if b['path'] == str(path)]
    _require(len(matches) == 1, 'missing or ambiguous original source authority: '+str(path))
    return matches[0]


def _native_source(request, saved, slots):
    adapters = {'adapter': _source_binding(request, HERE/'mixed_material_dynamics.py')}
    for name in ADAPTER_MODULES:
        path = (ROOT/'src'/Path(*name.split('.'))).with_suffix('.py') if '.' in name else HERE/(name+'.py')
        adapters[name] = _source_binding(request, path)
    return dict(model=saved['model_identity'], adapter_sources=adapters,
                completed_mean=request['history'], origin=saved['original_boundary'],
                slots=slots, block_size=request['block_size'])


def _original_authority(request, authority):
    _require(request.get('schema') == OWNER_SCHEMA, 'original Native producer owner required')
    for binding in request['source_bindings']:
        authority.bind(binding)
    initialization = authority.read(request['initialization'])
    _require(initialization.get('schema') == 'P254-matched-joint-case-v1'
             and initialization.get('status') == 'complete'
             and initialization.get('mode') == 'initialize'
             and initialization.get('case') == request['case'], 'original matched initialization required')
    saved = initialization['saved']
    _require(saved['case'] == request['case']
             and saved['model_identity'] == model_identity(saved['specification']),
             'saved model/specification authority differs')
    authority.bind(saved['preparation'])
    electronic = authority.read(saved['original_boundary'])
    initial = authority.read(saved['joint_boundary'])
    history = authority.read(request['history'])
    admit_oracle(authority.read(request['execution_proof']), existing_dependencies())
    _require(initial.get('complete') is True and initial.get('time_ry_inverse') == 0
             and initial.get('model_identity') == saved['model_identity']
             and initial.get('original_electronic_boundary') == saved['original_boundary']
             and initial['slots'] == electronic['slots']
             and electronic.get('time_ry_inverse') == 0,
             'literal original electronic/physical t0 authority differs')
    _require(history.get('schema') == 'streamed-joint-mean-v1'
             and history.get('complete') is True and bool(history.get('history'))
             and history.get('model_identity') == saved['model_identity']
             and history.get('original_electronic_boundary') == saved['original_boundary']
             and history['history'][0]['prior_boundary_binding'] == saved['joint_boundary'],
             'completed original moving mean history required')
    slots = [{k: v for k, v in s.items() if k != 'frame'} for s in initial['slots']]
    _require(slots == [{k: v for k, v in s.items() if k != 'frame'} for s in history['slots']],
             'moving mean changed original fixed D/slot/weight/rank authority')
    for slot in initial['slots']:
        authority.bind(slot['coefficient'])
    return saved, initial, history, slots


def _seed(directory, request, saved, initial, slots, authority):
    binding = authority.at(directory/'initial'/'checkpoint.json')
    record = authority.read(binding)
    shapes, plus, minus = metadata_geometry(initial, request['q'])
    small = _encode_small(initial_small(initial, q=request['q'], atom=request['atom'],
        axis=request['axis'], amplitude=request['amplitude'], cell=request['cell']))
    provenance = dict(model=saved['model_identity'], matched_case=request['case'],
        initialization=request['initialization'], completed_mean=request['history'],
        history_execution_proof=request['execution_proof'],
        physical_inputs=nested_bindings(request['inputs'], request['source_bindings']),
        source_bindings=request['source_bindings'], normalization=request['normalization'],
        original_local_pulse={k: request[k] for k in ('atom', 'axis', 'amplitude', 'cell')},
        electronic_initialization='new literal raw first/mixed zeros; no weighted-response import or D inverse',
        physical_initialization='first per-cell p q/-q kicks; all u/Maxwell/mixed small initially zero',
        fixed_raw_D_and_slots=slots)
    expected = dict(shapes=shapes, small=small, q=request['q'], target_plus=plus,
                    target_minus=minus, identity=saved['model_identity'], provenance=provenance)
    _require(authority.read(authority.at(directory/'initial'/'request.json')) == expected,
             'original physical seed request differs from owner/model/geometry')
    _require(all(record.get(k) == v for k, v in expected.items())
             and record.get('time') == 0
             and record.get('kind') == 'physical-local-pulse-raw-electronic-zero'
             and record.get('legs') == {c: [exact_zero_payload(s) for s in shapes[c]] for c in COMPONENTS},
             'initial state is not the original complete physical raw-zero pulse seed')
    _require(request['normalization'] == dict(
        source=_source_binding(request, HERE/'local_pulse_energy.py'),
        function='fourier_pulse; physical inverse-DFT coefficient = collective/sqrt(N)',
        convention='physical p_plus=J exp(-2 pi i q.cell)/N; p_minus=conj(p_plus)',
        derivative='independent analytic q/-q; eta=lambda*; mixed derivative has no half'),
        'original physical pulse normalization differs')
    return binding, record


def _layout(record, seed):
    for key in ('identity', 'q', 'shapes', 'target_plus', 'target_minus'):
        _require(record[key] == seed[key], 'changed original raw mixed geometry/model: '+key)
    for component in COMPONENTS:
        for name in SMALL_NAMES:
            _require(len(record['small'][component][name]) == len(seed['small'][component][name]),
                     'changed physical matter/global-field dimension')


def _frame_provenance(request, authority):
    provenance = []
    for frame in request['frames']:
        proof = authority.read(frame['proof'])
        descriptor = authority.read(frame['cache_descriptor'])
        complete = authority.read(frame['cache_complete'])
        _require(proof.get('schema') == 'P254-real-frame-consumer-admission-v1'
                 and proof.get('complete') is True and complete.get('complete') is True
                 and complete.get('cache') == frame['cache_descriptor']
                 and descriptor['history'] == request['history']
                 and descriptor['owner'] == proof['owner']
                 and descriptor['original_request'] == frame['original_request'] == proof['original_request']
                 and frame['admitted_request'] == proof['expected_request']
                 and frame['admitted_request']['input'] == frame['original_request']['input'],
                 'original actual frame admission differs')
        provenance.append(dict(proof=frame['proof'], cache=descriptor,
            consumer_request=frame['admitted_request'], original_request=frame['original_request'],
            cache_complete=frame['cache_complete']))
    return provenance


def _baseline(request, history, interval, time, authority):
    stages = history['history']
    stage = stages[interval]
    prior_binding = stage['prior_boundary_binding']
    prior = authority.read(prior_binding)
    endpoint_binding = stages[interval+1]['prior_boundary_binding'] if interval+1 < len(stages) else request['history']
    endpoint = authority.read(endpoint_binding)
    _require(prior.get('model_identity') == history['model_identity']
             and endpoint.get('model_identity') == history['model_identity'], 'foreign mean boundary')
    boundary = endpoint_binding if time == endpoint['time_ry_inverse'] else prior_binding
    return dict(time=float(time), interval=interval, stage=stage, baseline_boundary=boundary,
                model=history['model_identity'], origin=history['original_electronic_boundary'],
                block_size=request['block_size'])


def _stage(directory, *, step_binding, stage, base, current, accumulator,
           dt, identity, seed, authority):
    """Check an existing complete barrier without calling a publication helper."""
    complete = authority.read(authority.at(directory/'complete.json'))
    predecessors = dict(current=current.binding, accumulator=None if accumulator is None else accumulator.binding)
    _require(complete.get('schema') == STAGE_COMPLETE_SCHEMA and complete.get('complete') is True
             and type(complete.get('stage')) is int and complete['stage'] == stage
             and complete.get('step') == step_binding and complete.get('predecessors') == predecessors,
             'incomplete or foreign RK4 physical lineage')
    plan_binding = authority.at(directory/'plan.json')
    planned = authority.read(plan_binding)
    _require(complete.get('plan') == plan_binding and planned.get('schema') == STAGE_SCHEMA
             and type(planned.get('stage')) is int and planned['stage'] == stage
             and planned.get('step') == step_binding and planned.get('predecessors') == predecessors,
             'RK4 prepared stage differs from its complete source barrier')
    authority.bind(planned['prepared'])
    count = len(base._shapes['mixed'])
    _require(len(complete['transitions']) == len(complete['work_receipts']) == count,
             'RK4 source-prefix is not a complete physical boundary')
    transitions = []
    for index in range(count):
        expected = authority.at(directory/f'edge-{index:06d}'/'transition.json')
        _require(complete['transitions'][index] == expected, 'changed or reordered source transition')
        actual, transition = read_transition(expected['path'], plan=plan_binding, stage=stage,
            index=index, base=base, current=current, accumulator=accumulator,
            require_retired=True, live_successors=False)
        _require(actual == expected, 'changed original source transition bytes')
        for binding in transition['barriers'].values():
            authority.bind(binding)
        receipt = authority.read(complete['work_receipts'][index])
        _require(receipt.get('complete') is True and type(receipt.get('index')) is int
                 and receipt['index'] == index and receipt.get('request') == complete['work_request']
                 and set(receipt['arrays']) == {'transition_sha256'}, 'actual complete QWork source receipt required')
        marker_binding = authority.bind(receipt['arrays']['transition_sha256'])
        marker = np.load(marker_binding['path'], allow_pickle=False)
        _require(marker.dtype == np.uint8 and marker.shape == (32,)
                 and np.array_equal(marker, transition_marker(expected)), 'QWork transition marker differs')
        transitions.append(transition)
    small = _small_outputs(base, accumulator, _decode_small(planned['small_rhs'], base.q), stage, dt)
    labels = ('final',) if stage == 3 else ('input', 'accumulator')
    _require(isinstance(complete.get('outputs'), dict) and set(complete['outputs']) == set(labels),
             'completed RK4 output labels differ')
    for label in labels:
        binding = authority.at(directory/(label+'.json'))
        record = authority.read(binding)
        _require(complete['outputs'][label] == binding, 'changed physical/evaluation successor binding')
        _layout(record, seed)
        time = base.time+dt*((.5, .5, 1.)[stage] if label == 'input' else 1.)
        kind = ('mixed-rk4-evaluation-node' if label == 'input' else
                'mixed-rk4-physical-endpoint' if stage == 3 else 'mixed-rk4-quadrature-accumulator')
        provenance = dict(operation=STEP_SCHEMA, step=step_binding, stage=stage, label=label,
            base=base.binding, current=current.binding,
            accumulator=None if accumulator is None else accumulator.binding,
            prepared_stage=plan_binding, transitions=complete['transitions'], physical_endpoint=stage == 3)
        _require(record.get('time') == time and record.get('kind') == kind
                 and record.get('identity') == identity and record.get('provenance') == provenance
                 and record.get('legs') == {c: [t['successors'][label][c] for t in transitions] for c in COMPONENTS}
                 and record.get('small') == _encode_small(small[label]),
                 'successor is not the authenticated full RK4 physical/evaluation output')
        # Validate conventions, typed payload headers, and authentic retirement
        # sidecars for historical nodes, without reading retired source values.
        with MappedMixedState(binding):
            pass
    return complete['outputs']


def _endpoint(directory, owner_binding, request, saved, history, slots, initial, seed,
              state_binding, authority):
    steps = time_steps(history, request['time'], request['subdivisions'])
    paths = [directory/'steps'/f'{index:06d}'/'stage-3'/'final.json' for index in range(len(steps))]
    _require(state_binding['path'] in [str(p) for p in paths], 'state is not an original physical endpoint path')
    source = _native_source(request, saved, slots)
    work_identity = digest(dict(owner=owner_binding, dynamics=source, initial=initial,
        dependencies=request['source_bindings'], actual_frame_reuse=_frame_provenance(request, authority)))
    expected_sources = {key: _source_binding(request, HERE/name) for key, name in (
        ('producer', 'mixed_response_rk4.py'), ('transition', 'mixed_response_edge_transition.py'),
        ('storage', 'mixed_response_storage.py'), ('payload', 'mixed_payload_output.py'),
        ('dynamics', 'mixed_material_dynamics.py'))}
    predecessor = initial
    for index, step in enumerate(steps):
        step_directory = directory/'steps'/f'{index:06d}'
        step_binding = authority.at(step_directory/'request.json')
        _require(authority.read(step_binding) == dict(schema=STEP_SCHEMA, initial=predecessor,
            dt=step['dt'], identity=saved['model_identity'], work_identity=work_identity, sources=expected_sources),
            'original owner-bound physical step request differs')
        with ExitStack() as stack:
            base = stack.enter_context(MappedMixedState(predecessor))
            _layout(base.record, seed)
            _require(base.time == step['start'], 'broken contiguous physical predecessor chain')
            current, accumulator = base, None
            for stage in range(4):
                prepared_identity = dict(schema=NATIVE_SCHEMA, source=source, input=current.binding,
                    baseline=_baseline(request, history, step['interval'], current.time, authority))
                stage_directory = step_directory/f'stage-{stage}'
                plan = authority.read(authority.at(stage_directory/'plan.json'))
                prepared = authority.read(plan['prepared'])
                _require(prepared.get('schema') == NATIVE_SCHEMA and prepared.get('complete') is True
                         and prepared.get('identity') == prepared_identity,
                         'foreign Native prepared state/model/source/mean')
                outputs = _stage(stage_directory, step_binding=step_binding,
                    stage=stage, base=base, current=current, accumulator=accumulator, dt=step['dt'],
                    identity=saved['model_identity'], seed=seed, authority=authority)
                if stage < 3:
                    current = stack.enter_context(MappedMixedState(outputs['input']))
                    accumulator = stack.enter_context(MappedMixedState(outputs['accumulator']))
            predecessor = outputs['final']
        _require(authority.read(predecessor)['time'] == step['end'], 'physical endpoint time differs')
        if predecessor == state_binding:
            return
    raise ValueError('requested state is not in the original complete physical chain')


def admit_saved_boundary(producer_directory, owner_sha256, state_binding, *, source_archive=None):
    """Return authenticated original data, or raise on any missing/foreign input.

    This function is read-only. It does not acquire a persistent writer lock;
    callers must exclude the producer before positive-time admission/consumption.
    Missing or retired requested legs are refused even when headers survive.
    """
    directory = _plain_path(producer_directory)
    authority = _Authority()
    owner_binding = file_binding(directory/'owner.json')
    _require(owner_binding['sha256'] == owner_sha256, 'explicit original producer owner digest differs')
    owner = authority.read(owner_binding)
    if source_archive is not None:
        authority.recover_sources(source_archive, owner_binding, owner['request']['source_bindings'])
    state_binding = authority.bind(_binding(state_binding))
    state = authority.read(state_binding)
    _require(state.get('kind') in ('physical-local-pulse-raw-electronic-zero', 'mixed-rk4-physical-endpoint'),
             'evaluation node, quadrature accumulator or partial source state is not a physical boundary')
    request = owner['request']
    saved, original, history, slots = _original_authority(request, authority)
    initial, seed = _seed(directory, request, saved, original, slots, authority)
    _layout(state, seed)
    if state_binding == initial:
        boundary_kind = 'initial'
    else:
        _require(state.get('kind') == 'mixed-rk4-physical-endpoint' and state.get('time', 0) > 0,
                 'foreign initial state or nonphysical endpoint')
        _endpoint(directory, owner_binding, request, saved, history, slots, initial, seed,
                  state_binding, authority)
        boundary_kind = 'physical-endpoint'
    with MappedMixedState(state_binding) as mapped:
        for index in range(len(mapped._shapes['mixed'])):
            for component in COMPONENTS:
                entry = mapped.leg_binding(component, index)
                if payload_path(entry) is not None:
                    authority.bind(entry)
                mapped.leg(component, index)
            mapped.release_edge(index)
    return dict(producer_directory=str(directory), owner_binding=owner_binding, owner=owner,
                state_binding=state_binding, state=state, boundary_kind=boundary_kind,
                authority=list(authority.bindings.values()), source_archive=authority.source_archive)

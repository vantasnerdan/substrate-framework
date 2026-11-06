"""Recoverable, actual-QWork, one-step mixed RK4 producer.

Only a completed stage-3 checkpoint is a physical endpoint. Earlier current
states are RK4 evaluation nodes; accumulators carry a quadrature target time,
not an evolved physical state. Initial/base payloads are never retired. The
caller owns source-bound dynamics/work identities and the qualified QBudget;
this module neither resets resource counters nor estimates Hamiltonian costs.
"""
from contextlib import ExitStack
import inspect
import json
from pathlib import Path

import numpy as np

from substrate_framework.partitioned_q_storage import durable_publish
from substrate_framework.streamed_material_state import authenticate, file_binding
import mixed_response_edge_transition as edge_transition
import mixed_response_storage as storage
import mixed_payload_output as payload_output

SCHEMA = 'P254-mixed-rk4-step-v2'
STAGE_SCHEMA = 'P254-mixed-rk4-stage-plan-v2'
STAGE_COMPLETE_SCHEMA = 'P254-mixed-rk4-stage-complete-v2'


def _read(path):
    path = Path(path)
    if path.is_symlink():
        raise ValueError('RK4 journal cannot be a symlink')
    binding = file_binding(path)
    return binding, json.loads(authenticate(binding).read_text())


def _publish_exact(path, value, budget):
    if Path(path).exists():
        binding, previous = _read(path)
        if previous != value:
            raise ValueError('existing mixed RK4 journal/source/step differs')
        return binding
    return durable_publish(path, value, budget)


def _sources(dynamics):
    provider = inspect.getsourcefile(dynamics.prepare)
    if provider is None or not Path(provider).is_file():
        raise ValueError('dynamics.prepare requires an actual source-file binding')
    return dict(producer=file_binding(Path(__file__).resolve()),
                transition=file_binding(Path(edge_transition.__file__).resolve()),
                storage=file_binding(Path(storage.__file__).resolve()),
                payload=file_binding(Path(payload_output.__file__).resolve()),
                dynamics=file_binding(Path(provider).resolve()))


def _small(state):
    return {c: {name: getattr(state.small[c], name) for name in storage.SMALL_NAMES}
            for c in storage.COMPONENTS}


def _small_outputs(base, accumulator, rhs, stage, dt):
    outputs = {}
    for label in (('final',) if stage == 3 else ('input', 'accumulator')):
        scale = dt*(edge_transition.NODES[stage] if label == 'input'
                    else edge_transition.WEIGHTS[stage])
        previous = base if label == 'input' or accumulator is None else accumulator
        original = _small(previous)
        outputs[label] = {c: {name: original[c][name]+scale*rhs[c][name]
                             for name in storage.SMALL_NAMES} for c in storage.COMPONENTS}
        storage._small_arrays(outputs[label], base.q)
    return outputs


def _stage_plan(path, *, step, stage, current, accumulator, dynamics, work, budget):
    predecessors = dict(current=current.binding,
                        accumulator=None if accumulator is None else accumulator.binding)
    if Path(path).exists():
        binding, record = _read(path)
        if (record.get('schema') != STAGE_SCHEMA or record.get('step') != step
                or record.get('stage') != stage or record.get('predecessors') != predecessors):
            raise ValueError('prepared RK4 stage identity differs')
        authenticate(record['prepared'])
        prepared = dynamics.restore(record['prepared'], current, work)
        rhs = storage._decode_small(record['small_rhs'], current.q)
        actual = storage._small_arrays(prepared.small_rhs, current.q)
        if prepared.binding != record['prepared'] or any(
                not np.array_equal(rhs[c][name], actual[c][name])
                for c in storage.COMPONENTS for name in storage.SMALL_NAMES):
            raise ValueError('restored prepared RHS differs from saved complete derivative')
        return binding, prepared, rhs
    prepared = dynamics.prepare(current, work)
    authenticate(prepared.binding)
    rhs = storage._small_arrays(prepared.small_rhs, current.q)
    # Save all twelve small derivatives, including genuinely empty nonGamma
    # field directions, before any worker can retire even one source edge.
    record = dict(schema=STAGE_SCHEMA, step=step, stage=stage,
                  predecessors=predecessors, prepared=prepared.binding,
                  small_rhs=storage._encode_small(rhs))
    return durable_publish(path, record, budget), prepared, rhs


def _completed_stage(path, *, step, stage, current, accumulator, initial, dt,
                     identity, budget):
    binding, record = _read(path)
    predecessors = dict(current=current, accumulator=accumulator)
    if (record.get('schema') != STAGE_COMPLETE_SCHEMA or record.get('complete') is not True
            or record.get('step') != step or record.get('stage') != stage
            or record.get('predecessors') != predecessors):
        raise ValueError('completed RK4 stage source/step identity differs')
    labels = {'final'} if stage == 3 else {'input', 'accumulator'}
    if not isinstance(record.get('outputs'), dict) or set(record['outputs']) != labels:
        raise ValueError('completed RK4 stage checkpoint labels differ')
    plan = record['plan']
    planned = json.loads(authenticate(plan).read_text())
    if (planned.get('schema') != STAGE_SCHEMA or planned.get('step') != step
            or planned.get('stage') != stage or planned.get('predecessors') != predecessors):
        raise ValueError('completed stage prepared context identity differs')
    authenticate(planned['prepared'])
    with ExitStack() as stack:
        base_state = stack.enter_context(storage.MappedMixedState(initial, budget=budget))
        current_state = base_state if current == initial else stack.enter_context(
            storage.MappedMixedState(current, budget=budget))
        accumulator_state = None if accumulator is None else stack.enter_context(
            storage.MappedMixedState(accumulator, budget=budget))
        count = len(base_state._shapes['mixed'])
        if (len(record['work_receipts']) != count or len(record['transitions']) != count):
            raise ValueError('completed stage has an incomplete source barrier')
        transitions = []
        for index in range(count):
            expected = record['transitions'][index]
            authenticate(expected)
            transition = edge_transition.read_transition(expected['path'], plan=plan,
                stage=stage, index=index, base=base_state, current=current_state,
                accumulator=accumulator_state, require_retired=True, live_successors=False)
            if transition[0] != expected:
                raise ValueError('completed transition binding differs')
            receipt = json.loads(authenticate(record['work_receipts'][index]).read_text())
            if (receipt.get('complete') is not True or receipt.get('index') != index
                    or receipt.get('request') != record['work_request']
                    or set(receipt['arrays']) != {'transition_sha256'}):
                raise ValueError('completed stage actual QWork receipt differs')
            marker = np.load(authenticate(receipt['arrays']['transition_sha256']), allow_pickle=False)
            if (marker.dtype != np.uint8 or marker.shape != (32,)
                    or not np.array_equal(marker, edge_transition.transition_marker(expected))):
                raise ValueError('completed stage marker does not authenticate its transition')
            transitions.append(transition)
        for label, checkpoint in record['outputs'].items():
            if checkpoint['path'] != str(Path(path).parent/(label+'.json')):
                raise ValueError('completed successor checkpoint escaped its assigned stage')
            authenticate(checkpoint)
        rhs = storage._decode_small(planned['small_rhs'], base_state.q)
        small_outputs = _small_outputs(base_state, accumulator_state, rhs, stage, dt)
        outputs = _publish_outputs(Path(path).parent, step=step, stage=stage,
            base=base_state, current=current_state, accumulator=accumulator_state,
            small_outputs=small_outputs, transitions=transitions, plan=plan, dt=dt,
            identity=identity, budget=budget)
        if outputs != record['outputs']:
            raise ValueError('completed successor binding differs')
    return binding, record


def _publish_outputs(directory, *, step, stage, base, current, accumulator,
                     small_outputs, transitions, plan, dt, identity, budget):
    checkpoints = {}
    labels = ('final',) if stage == 3 else ('input', 'accumulator')
    for label in labels:
        legs = {c: [record['successors'][label][c] for _, record in transitions]
                for c in storage.COMPONENTS}
        if label == 'input':
            time = base.time+dt*edge_transition.NODES[stage]
            kind = 'mixed-rk4-evaluation-node'
        else:
            time = base.time+dt
            kind = 'mixed-rk4-physical-endpoint' if stage == 3 else 'mixed-rk4-quadrature-accumulator'
        provenance = dict(operation=SCHEMA, step=step, stage=stage, label=label,
            base=base.binding, current=current.binding,
            accumulator=None if accumulator is None else accumulator.binding,
            prepared_stage=plan, transitions=[binding for binding, _ in transitions],
            physical_endpoint=stage == 3)
        path = Path(directory)/(label+'.json')
        if path.exists():
            binding, record, shapes, _, _, _, arrays = storage._metadata(file_binding(path))
            if (record['legs'] != legs or shapes != base._shapes
                    or record['time'] != time or record['kind'] != kind
                    or record['identity'] != identity or record['provenance'] != provenance
                    or record['q'] != base.q.tolist()
                    or record['target_plus'] != list(base.target_plus)
                    or record['target_minus'] != list(base.target_minus)
                    or any(not np.array_equal(arrays[c][name], small_outputs[label][c][name])
                           for c in storage.COMPONENTS for name in storage.SMALL_NAMES)):
                raise ValueError('existing RK4 successor checkpoint differs')
        else:
            binding = storage.publish_checkpoint(path, legs=legs, small=small_outputs[label],
                q=base.q, target_plus=base.target_plus, target_minus=base.target_minus,
                time=time, identity=identity, kind=kind, provenance=provenance, budget=budget)
        checkpoints[label] = binding
    return checkpoints


def advance_step(initial_binding, dt, *, dynamics, work, directory, identity, budget):
    """Advance exactly one recoverable classical RK4 step through actual QWork.

    ``dynamics.prepare(state, work)`` and ``restore(binding, state, work)`` return
    a durable ``binding``, complete ``small_rhs``, and ``edge_rhs(index, state)``.
    All global jet reductions must finish inside prepare, before it returns;
    edge_rhs may then read only that index's plus/minus/mixed jet payloads.
    ``directory`` exclusively owns this step's new checkpoints and payloads.
    Reuse the same directory, initial binding, dt, identities and source bytes
    after interruption; changed inputs are refused rather than silently mixed.
    """
    if (isinstance(dt, (bool, np.bool_)) or not np.isscalar(dt) or np.iscomplexobj(dt)
            or not np.isfinite(dt) or dt <= 0 or not isinstance(identity, str) or not identity):
        raise ValueError('positive finite real dt and explicit nonempty step identity required')
    if work.budget is not budget:
        raise ValueError('producer and actual QWork must share the original QBudget')
    directory = Path(directory)
    if directory.is_symlink():
        raise ValueError('step owner cannot be a symlink')
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    budget.check()
    # Validate even a completed invocation against the currently supplied source
    # files/work source identity, not merely a mutable caller-selected filename.
    initial_binding = storage._metadata(initial_binding)[0]
    request = dict(schema=SCHEMA, initial=initial_binding, dt=float(dt), identity=identity,
                   work_identity=work.identity, sources=_sources(dynamics))
    step = _publish_exact(directory/'request.json', request, budget)
    current_binding, accumulator_binding = initial_binding, None
    for stage in range(4):
        stage_directory = directory/f'stage-{stage}'
        stage_directory.mkdir(exist_ok=True)
        complete = stage_directory/'complete.json'
        if complete.exists():
            _, record = _completed_stage(complete, step=step, stage=stage,
                current=current_binding, accumulator=accumulator_binding,
                initial=initial_binding, dt=float(dt), identity=identity, budget=budget)
            outputs = record['outputs']
        else:
            with ExitStack() as stack:
                # Open every state BEFORE fork, while no local worker can unlink
                # a predecessor inspected by another state's constructor.
                base = stack.enter_context(storage.MappedMixedState(initial_binding, budget=budget))
                current = base if current_binding == initial_binding else stack.enter_context(
                    storage.MappedMixedState(current_binding, budget=budget))
                accumulator = None if accumulator_binding is None else stack.enter_context(
                    storage.MappedMixedState(accumulator_binding, budget=budget))
                plan, prepared, rhs = _stage_plan(stage_directory/'plan.json', step=step,
                    stage=stage, current=current, accumulator=accumulator,
                    dynamics=dynamics, work=work, budget=budget)
                small_outputs = _small_outputs(base, accumulator, rhs, stage, float(dt))
                count = len(base._shapes['mixed'])
                # prepare can legitimately map every edge for global moments.
                # Children inherit metadata and EMPTY map caches, never global
                # parent mappings whose payloads another child may retire.
                for state in (base, current, accumulator):
                    if state is not None:
                        for index in range(count):
                            state.release_edge(index)
                def advance_edge(index):
                    return edge_transition.execute_edge(index, plan=plan, stage=stage,
                        dt=float(dt), prepared=prepared, base=base, current=current,
                        accumulator=accumulator, directory=stage_directory,
                        budget=budget, owner=directory)
                receipts = work.execute('mixed-rk4-edge-transition',
                    dict(step=step, stage=stage, plan=plan), count, advance_edge)
                key, work_directory, work_count = work.last
                if work_count != count or len(receipts) != count:
                    raise ValueError('incomplete actual RK4 QWork barrier')
                transitions = []
                for index, receipt in enumerate(receipts):
                    transition = edge_transition.read_transition(
                        stage_directory/f'edge-{index:06d}'/'transition.json', plan=plan,
                        stage=stage, index=index, base=base, current=current,
                        accumulator=accumulator, require_retired=True)
                    if set(receipt['arrays']) != {'transition_sha256'}:
                        raise ValueError('QWork must return the actual transition SHA marker only')
                    marker_binding = receipt['arrays']['transition_sha256']
                    authenticate(marker_binding)
                    marker = np.load(marker_binding['path'], allow_pickle=False)
                    if (marker.dtype != np.uint8 or marker.shape != (32,)
                            or not np.array_equal(marker, edge_transition.transition_marker(transition[0]))):
                        raise ValueError('actual QWork marker does not authenticate transition bytes')
                    transitions.append(transition)
                outputs = _publish_outputs(stage_directory, step=step, stage=stage,
                    base=base, current=current, accumulator=accumulator,
                    small_outputs=small_outputs, transitions=transitions, plan=plan,
                    dt=float(dt), identity=identity, budget=budget)
                # Keep observed original QWork receipts AND pre-receipt transition
                # measurements. Recovery never fabricates costs for skipped RHS.
                value = dict(schema=STAGE_COMPLETE_SCHEMA, complete=True, step=step,
                    stage=stage, predecessors=dict(current=current_binding,
                        accumulator=accumulator_binding), plan=plan, outputs=outputs,
                    work_request=key, work_receipts=[file_binding(work_directory/
                        ('item-%04d' % i)/'complete.json') for i in range(count)],
                    transitions=[binding for binding, _ in transitions], budget=budget.snapshot())
                durable_publish(complete, value, budget)
        if stage == 3:
            final = outputs['final']
            # A completed return authenticates every full-shape final value:
            # actual file bytes or an explicit exact-zero descriptor, never
            # checkpoint headers, compact markers or absent-file inference.
            with storage.MappedMixedState(final, budget=budget) as state:
                for index in range(len(state._shapes['mixed'])):
                    for component in storage.COMPONENTS:
                        state.leg(component, index)
                    state.release_edge(index)
            return final
        current_binding, accumulator_binding = outputs['input'], outputs['accumulator']
    raise RuntimeError('RK4 physical endpoint was not published')

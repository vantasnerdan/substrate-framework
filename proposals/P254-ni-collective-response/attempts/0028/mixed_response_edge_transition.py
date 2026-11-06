"""Assigned-edge RK4 writes, authenticated successor journals, and retirement.

A transition is published before any old payload is unlinked. Its SHA256 is the
actual QWork ndarray result. Replaying it never evaluates the RHS or opens old
jet payloads; original execution costs remain in the immutable transition even
when the subsequently recovered QWork receipt measures only recovery work.
"""
import json
import os
from pathlib import Path
import resource
from time import monotonic, process_time
from uuid import uuid4

import numpy as np

from substrate_framework.partitioned_q_storage import (
    durable_publish, retire_payloads, sync_directory)
from substrate_framework.streamed_material_state import authenticate, file_binding
from mixed_response_storage import (
    COMPONENTS, _barrier, _owned, _retirement, authenticate_payload,
    payload_binding, payload_path, publish_edge_barrier, retire_edge)
from mixed_payload_output import write_payload

SCHEMA = 'P254-mixed-rk4-edge-transition-v2'
WEIGHTS = (1/6, 1/3, 1/3, 1/6)
NODES = (.5, .5, 1.)


def transition_marker(binding):
    return np.frombuffer(bytes.fromhex(binding['sha256']), dtype=np.uint8).copy()


def _predecessors(base, current, accumulator):
    return dict(base=base.binding, current=current.binding,
                accumulator=None if accumulator is None else accumulator.binding)


def _old_states(base, current, accumulator):
    result = {}
    if current.binding != base.binding:
        result['current'] = current
    if accumulator is not None:
        if accumulator.binding in (base.binding, current.binding):
            raise ValueError('RK4 accumulator must be a distinct owned checkpoint')
        result['accumulator'] = accumulator
    return result


def read_transition(path, *, plan, stage, index, base, current, accumulator,
                    require_retired=False, live_successors=True):
    """Authenticate replacements, or retained journals after later retirement."""
    path = Path(path)
    if path.is_symlink():
        raise ValueError('transition cannot be a symlink')
    binding = file_binding(path)
    record = json.loads(authenticate(binding).read_text())
    if (record.get('schema') != SCHEMA or record.get('complete') is not True
            or record.get('plan') != plan or record.get('stage') != stage
            or type(record.get('index')) is not int or record['index'] != index
            or record.get('predecessors') != _predecessors(base, current, accumulator)):
        raise ValueError('transition source/step/stage/predecessors differ')
    labels = {'final'} if stage == 3 else {'input', 'accumulator'}
    successors = record.get('successors')
    if not isinstance(successors, dict) or set(successors) != labels:
        raise ValueError('transition replacement labels differ')
    seen = set()
    old_paths = {path for state in (base, current, accumulator) if state is not None
                 for component in COMPONENTS for entry in state.record['legs'][component]
                 if (path := payload_path(entry)) is not None}
    for triple in successors.values():
        if not isinstance(triple, dict) or set(triple) != set(COMPONENTS):
            raise ValueError('complete successor triple required')
        for component in COMPONENTS:
            entry = payload_binding(triple[component])
            entry_path = payload_path(entry)
            if entry_path is not None:
                if entry_path in seen or entry_path in old_paths:
                    raise ValueError('replacement aliases an original or another replacement')
                seen.add(entry_path)
            # Typed zeros retain full geometry even after later stages retire
            # the real successor files. Missing real files are never zeros.
            if live_successors or entry_path is None:
                authenticate_payload(entry, base._shapes[component][index])
    old = _old_states(base, current, accumulator)
    if not isinstance(record.get('barriers'), dict) or set(record['barriers']) != set(old):
        raise ValueError('separate old-checkpoint barriers required')
    for label, state in old.items():
        barrier = record['barriers'][label]
        _, value = _barrier(barrier, state.binding, state.record, state._shapes,
                            index, live_successors=live_successors)
        if value['successors'] != successors:
            raise ValueError('retirement barrier replacements differ from transition')
        if require_retired:
            _retirement(state.binding, state.record, state._shapes, index,
                        expected_barrier=barrier)
            if any(Path(path).exists() for c in COMPONENTS
                   if (path := payload_path(state.record['legs'][c][index])) is not None):
                raise ValueError('completed edge still has unretired owned payloads')
    measured = record.get('observed_execution')
    if not isinstance(measured, dict) or any(name not in measured for name in (
            'pid', 'cpu', 'affinity', 'hamiltonian_actions', 'action_columns',
            'cpu_seconds', 'elapsed_seconds', 'process_peak_rss_bytes')):
        raise ValueError('actual edge execution measurements required')
    return binding, record


def _finish_retirement(record, *, index, base, current, accumulator, budget, owner):
    for label, state in _old_states(base, current, accumulator).items():
        state.release_edge(index)
        retire_edge(state, index, record['barriers'][label], budget=budget, owner_root=owner)


def _discard_uncommitted(edge, plan, budget):
    """Recover only this edge's unpublished attempts via qualified accounting."""
    owner = edge.resolve()
    for attempt in sorted(edge.glob('attempt-*')):
        if attempt.is_symlink() or not attempt.is_dir() or attempt.parent.resolve() != owner:
            raise ValueError('uncommitted attempt escaped its assigned edge')
        receipt = attempt/'abandoned.json'
        paths = sorted(attempt.glob('*.npy'))
        if receipt.exists():
            _owned(receipt, owner)
            previous = json.loads(receipt.read_text())
            for entry in previous['payloads']:
                _owned(entry['path'], owner)
        for path in paths:
            _owned(path, owner)
        if paths or receipt.exists():
            retire_payloads(receipt, plan, paths, owner, budget)


def execute_edge(index, *, plan, stage, dt, prepared, base, current, accumulator,
                 directory, budget, owner):
    """Called by the actual pinned QWork worker for its exclusively owned edge."""
    edge = Path(directory)/f'edge-{index:06d}'
    if edge.is_symlink() or not edge.resolve().is_relative_to(Path(owner).resolve()):
        raise ValueError('assigned edge directory escaped the step owner')
    edge.mkdir(parents=True, exist_ok=True)
    transition = edge/'transition.json'
    if transition.exists():
        binding, record = read_transition(transition, plan=plan, stage=stage, index=index,
            base=base, current=current, accumulator=accumulator)
        _finish_retirement(record, index=index, base=base, current=current,
                           accumulator=accumulator, budget=budget, owner=owner)
        return {'transition_sha256': transition_marker(binding)}

    # No transition means these attempts were never published as successors.
    # Their compact journals/costs survive; credit only actually unlinked bytes.
    _discard_uncommitted(edge, plan, budget)
    attempt = edge/('attempt-'+uuid4().hex)
    attempt.mkdir()
    started, cpu_started = monotonic(), process_time()
    actions, columns = budget.local_actions, budget.local_columns
    rhs = prepared.edge_rhs(index, current)
    if not isinstance(rhs, dict) or set(rhs) != set(COMPONENTS):
        raise ValueError('edge RHS requires exactly plus/minus/mixed')
    labels = ('final',) if stage == 3 else ('input', 'accumulator')
    successors = {label: {} for label in labels}
    for component in COMPONENTS:
        derivative = np.asarray(rhs[component])
        if (derivative.shape != base._shapes[component][index]
                or derivative.dtype.kind not in 'fciub' or not np.isfinite(derivative).all()):
            raise ValueError('finite original full-row numerical RHS required')
        for label in labels:
            scale = dt*(NODES[stage] if label == 'input' else WEIGHTS[stage])
            previous = base if label == 'input' or accumulator is None else accumulator
            value = np.multiply(derivative, scale, dtype=np.complex128)
            np.add(value, previous.leg(component, index), out=value)
            entry = write_payload(attempt/(label+'-'+component+'.npy'), value, budget)
            successors[label][component] = entry
            del value
        del derivative
    del rhs
    sync_directory(attempt)
    barriers = {label: publish_edge_barrier(attempt/(label+'-barrier.json'),
        checkpoint=state.binding, index=index, successors=successors, budget=budget)
        for label, state in _old_states(base, current, accumulator).items()}
    affinity = sorted(os.sched_getaffinity(0))
    record = dict(schema=SCHEMA, complete=True, plan=plan, stage=stage, index=index,
        predecessors=_predecessors(base, current, accumulator), successors=successors,
        barriers=barriers, observed_execution=dict(pid=os.getpid(),
            cpu=affinity[0] if len(affinity) == 1 else None, affinity=affinity,
            hamiltonian_actions=budget.local_actions-actions,
            action_columns=budget.local_columns-columns,
            cpu_seconds=process_time()-cpu_started, elapsed_seconds=monotonic()-started,
            process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))
    binding = durable_publish(transition, record, budget)
    # RHS may have borrowed current views; all such local references are gone.
    _finish_retirement(record, index=index, base=base, current=current,
                       accumulator=accumulator, budget=budget, owner=owner)
    base.release_edge(index)
    return {'transition_sha256': transition_marker(binding)}

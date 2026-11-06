"""Bounded Native measurements of original and split full-PW edge equations."""
import os
import resource
from time import perf_counter, process_time

import numpy as np
from threadpoolctl import threadpool_info

from substrate_framework.partitioned_joint_budget import resident_bytes
from substrate_framework.partitioned_q_storage import bounded_array_identity
from substrate_framework.source_local_q_edge import dual_moments, dual_free
from substrate_framework.source_local_q_state import edge_state, small_state, terminal_edge
from substrate_framework.streamed_joint_response import NAMES


def _process_evidence():
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return dict(pid=os.getpid(), affinity=sorted(os.sched_getaffinity(0)),
        current_rss_bytes=resident_bytes(os.getpid()), peak_rss_bytes=int(usage.ru_maxrss)*1024,
        user_cpu_seconds=usage.ru_utime, system_cpu_seconds=usage.ru_stime,
        threadpools=threadpool_info())


def _array_evidence(value):
    result = bounded_array_identity(value)
    result.update(nbytes=int(value.nbytes), finite=bool(np.isfinite(value).all()))
    if not result['finite']:
        raise FloatingPointError('non-finite original full-PW edge measurement')
    return result


def _typed_array(value, shape):
    if not isinstance(value, np.ndarray) or value.shape != tuple(shape) or value.dtype != np.complex128:
        raise ValueError('original full-shape complex128 edge result required')
    return _array_evidence(value)


def _split_comparison(executor, index, time, state, reference):
    budget = executor.budget
    before = budget.snapshot()
    started = perf_counter()
    moments = dual_moments(executor, index, time, state)
    moments_done = perf_counter()
    middle = budget.snapshot()
    for name in ('n', 'g', 'ell'):
        if not np.array_equal(moments[name], reference[name]):
            raise ValueError('Native moment-only result differs: '+name)
    moment_evidence = {name: _array_evidence(value) for name, value in moments.items()}
    del moments
    free_started = perf_counter()
    free = dual_free(executor, index, time, state)
    free_done = perf_counter()
    after = budget.snapshot()
    for name, value in zip(('ket', 'bra'), free):
        if not np.array_equal(value, reference[name]):
            raise ValueError('Native free-only result differs: '+name)
    free_evidence = {name: _array_evidence(value) for name, value in zip(('ket', 'bra'), free)}
    if (middle['hamiltonian_actions'] != before['hamiltonian_actions']
            or after['hamiltonian_actions']-middle['hamiltonian_actions'] != 2):
        raise ValueError('split edge did not preserve its exact physical-action accounting')
    return dict(array_exact=True, signed_zero_byte_identity_not_claimed=True,
        moment_only=moment_evidence, free_only=free_evidence,
        moment_only_seconds=moments_done-started, free_only_seconds=free_done-free_started,
        moment_only_actions=middle['hamiltonian_actions']-before['hamiltonian_actions'],
        free_only_actions=after['hamiltonian_actions']-middle['hamiltonian_actions'],
        moment_only_columns=middle['action_columns']-before['action_columns'],
        free_only_columns=after['action_columns']-middle['action_columns'])


def measure_edge(runtime, *, index, time=0.0):
    """Measure one unreduced terminal/independent edge, never whole noise/RHS."""
    executor, probe = runtime.executor, runtime.probe
    response, budget = executor.r, executor.budget
    if type(index) is not int or not 0 <= index < len(response.edges):
        raise ValueError('selected original explicit edge index required')
    time = float(time)
    if not np.isfinite(time) or not 0 <= time <= executor.edge_stage.final.time:
        raise ValueError('selected probe time must be covered by actual saved history')
    probe.validate(response)
    budget.check_processes(os.getpid(), [])
    before, process_before = budget.snapshot(), _process_evidence()
    started, cpu_started = perf_counter(), process_time()
    executor.prepare_time(time)
    frames_done = perf_counter()
    pair = terminal_edge(executor, probe, time, index)
    terminal_done = perf_counter()
    small = small_state(response, {name: np.array(getattr(probe, name), copy=True) for name in NAMES})
    values = executor.independent_edge(index, time, edge_state(index, pair, small), True)
    independent_done = perf_counter()
    edge = response.edges[index]
    source, target = response.sectors[edge.source], response.sectors[edge.target]
    shapes = dict(ket=(target.shape[0], source.shape[1]),
                  bra=(target.shape[1], source.shape[0]),
                  n=response.density_shape, g=(response.ndof,), ell=(3,))
    if set(values) != set(shapes):
        raise ValueError('original independent-edge typed output keys changed')
    terminal = dict(ket=_typed_array(pair[0], shapes['ket']), bra=_typed_array(pair[1], shapes['bra']),
        small={name: _array_evidence(np.asarray(getattr(small, name))) for name in NAMES})
    output = {name: _typed_array(value, shapes[name]) for name, value in values.items()}
    split_comparison = _split_comparison(
        executor, index, time, edge_state(index, pair, small), values)
    budget.check_processes(os.getpid(), [])
    after, process_after = budget.snapshot(), _process_evidence()
    counters = ('hamiltonian_actions', 'action_columns', 'cumulative_write_bytes_reserved',
                'cumulative_allocation_bytes_reserved', 'mapped_overwrite_bytes_reserved',
                'retired_payload_bytes', 'live_bytes_reserved')
    return dict(schema='P254-source-local-material-one-edge-probe-v1', status='complete',
        measurement='actual original full independent edge and exact consumed split projections',
        executor_identity=executor.identity, baseline=runtime.provenance['baseline'],
        source_bindings=runtime.provenance['source_bindings'], stage=executor.edge_stage.stage_identity(time),
        time=time, index=index, source=edge.source, target=edge.target,
        source_slot=executor.edge_stage.initial.slots[edge.source].identity,
        target_slot=executor.edge_stage.initial.slots[edge.target].identity,
        source_full_shape=list(source.shape), target_full_shape=list(target.shape),
        original_sector_count=len(response.sectors), original_edge_count=len(response.edges),
        full_original_pw_and_columns=True, observable=probe.identity, dynamics=executor.dynamics,
        terminal=terminal, independent_edge=output, split_comparison=split_comparison,
        budget_before=before, budget_after=after,
        budget_delta={key: after[key]-before[key] for key in counters},
        process_before=process_before, process_after=process_after,
        timing=dict(frame_barrier_seconds=frames_done-started,
            terminal_probe_seconds=terminal_done-frames_done,
            independent_edge_seconds=independent_done-terminal_done,
            total_wall_seconds=perf_counter()-started, parent_cpu_seconds=process_time()-cpu_started),
        whole_rhs_executed=False, pair_executed=False, whole_noise_result=False,
        numerical_accuracy_certified=False,
        limitations='One full original edge and genuine terminal probe only; no complete global moments, feedback barrier, RK4 pullback or two-sided covariance. Interior-time frame preparation may charge all original 128 source exponentials; measured budget deltas include that work. RSS evidence is before/after process residency and OS lifetime peak, not a separately sampled per-operation peak.')

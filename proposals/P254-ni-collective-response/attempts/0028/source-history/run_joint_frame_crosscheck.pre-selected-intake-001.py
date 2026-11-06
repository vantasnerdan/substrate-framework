#!/usr/bin/env python3
"""Compare one actual completed joint frame with the independent SciPy action.

Consumes a genuine completed matched mean history, not a fabricated partial
boundary. Reconstructs its last frozen interval in the original full PW basis.
This also measures the Taylor method used for interior moving-history frames.
One sector is a numerical crosscheck, not a whole-material accuracy certificate.
"""
import argparse
import json
import os
from pathlib import Path
import resource
import sys
from time import perf_counter

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

import numpy as np
from run_joint_q_material import bindings, selected_initialization
from matched_joint_case import load_case
from substrate_framework.material_control_reference import bound_json
from substrate_framework.material_physical_operators import PhysicalFrozenAction
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.streamed_joint_reduction import JointStage
from substrate_framework.streamed_material_evolution import exponential_action
from substrate_framework.streamed_material_state import file_binding


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    for name in ('initialization', 'history', 'producer-model'):
        parser.add_argument('--'+name, type=Path, required=True)
        parser.add_argument('--'+name+'-sha256', required=True)
    parser.add_argument('--slot-index', type=int, required=True)
    parser.add_argument('--cpu', type=int, choices=range(6, 11), required=True)
    parser.add_argument('--max-actions', type=int, required=True)
    parser.add_argument('--max-seconds', type=float, required=True)
    parser.add_argument('--max-memory-bytes', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True)
    return parser.parse_args()


def crosscheck(args, budget):
    initialization, saved = selected_initialization(args)
    history_binding = file_binding(args.history)
    model_binding = file_binding(args.producer_model)
    if history_binding['sha256'] != args.history_sha256 or model_binding['sha256'] != args.producer_model_sha256:
        raise ValueError('explicit actual history or producer-model digest differs')
    producer = bound_json(model_binding)
    context = load_case(args.case)
    if (producer['identity'] != context.identity or saved['model_identity'] != context.identity
            or producer['specification'] != context.specification
            or producer['initialization'] != initialization):
        raise ValueError('history producer and initialization must use this exact matched scalar')
    if any(file_binding(item['path']) != item for item in producer['implementation']):
        raise ValueError('actual mean implementation changed before the independent crosscheck')
    limits = producer['column_limits']
    if not limits:
        raise ValueError('actual original-column producer and its numerical allocation required')
    allowance = float(limits['absolute_error_target'])
    history = context.make_history(history_binding, budget, block_size=producer['block_size'])
    if file_binding(history.initial.path) != saved['joint_boundary']:
        raise ValueError('history does not start at the selected literal initialization')
    if not 0 <= args.slot_index < len(history.final.slots):
        raise ValueError('explicit original slot index is outside the actual history')
    prior = history.priors[-1]
    record = history.records[-1]
    stage = JointStage.restore(record)
    slot = prior.slots[args.slot_index]
    target_slot = history.final.slots[args.slot_index]
    initial = prior.arrays(slot)[0]
    target = history.final.arrays(target_slot)[0]
    source = bindings()+[file_binding(Path(__file__).resolve()), model_binding]
    result = dict(case=args.case, source_bindings=source, initialization=initialization,
        producer_model=model_binding, history=history_binding, model_identity=context.identity,
        original_slot=slot.identity, original_shape=list(slot.shape), initial_frame=slot.frame,
        actual_target_frame=target_slot.frame, original_D=slot.coefficient,
        interval=record['dt_ry_inverse'], absolute_comparison_allowance=allowance,
        numerical_accuracy_certified=False)
    (args.output/'started.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('Authenticated actual frame; starting independent full-PW Taylor action.', flush=True)
    mean = context.make_mean(budget, block_size=producer['block_size'])
    before = budget.snapshot()
    started = perf_counter()
    with mean.factory(slot) as (coupling, h0, _):
        frozen = PhysicalFrozenAction(h0, coupling, potential=stage.potential[:, slot.component],
            direction=stage.u, shift=stage.shift)
        def counted(vectors):
            budget.check_processes(os.getpid(), [])
            budget.action(vectors)
            return frozen.apply(vectors)
        reference = exponential_action(counted, initial, record['dt_ry_inverse'])
    elapsed = perf_counter()-started
    error = float(np.linalg.norm(reference-target))
    scale = float(np.linalg.norm(reference))
    wrong_time_error = float(np.linalg.norm(reference-initial))
    if not np.isfinite([error, scale, wrong_time_error]).all():
        raise FloatingPointError('nonfinite actual/reference frame comparison')
    checks = dict(actual_frame_agrees_with_independent_action=error <= allowance,
        rejects_original_zero_interval_frame=wrong_time_error > allowance)
    result.update(status='passed' if all(checks.values()) else 'failed', checks=checks,
        absolute_frobenius_difference=error, relative_frobenius_difference=error/max(scale, np.finfo(float).tiny),
        zero_interval_mutation_difference=wrong_time_error, reference_method='existing SciPy full-PW scaling/Taylor exponential_action',
        reference_action_seconds=elapsed, budget_before_reference=before, budget_after_reference=budget.snapshot(),
        qualification='Floating independent-method corroboration at the producer declared target; neither method supplies an interval-certified material error bound. This is one actual sector and one frozen interval, not whole-history refinement or a physical coherence result.')
    return result


def main():
    args = parse_args()
    os.sched_setaffinity(0, {args.cpu})
    args.output.mkdir(parents=True, exist_ok=False)
    budget = QBudget(max_actions=args.max_actions, max_seconds=args.max_seconds,
        max_bytes=1048576, max_memory_bytes=args.max_memory_bytes,
        worker_memory_bytes=args.max_memory_bytes)
    budget.reserve(524288)
    started = perf_counter()
    result = dict(schema='P254-actual-joint-frame-crosscheck-v1', status='failed')
    try:
        result.update(crosscheck(args, budget))
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result.update(budget=budget.snapshot(), seconds=perf_counter()-started,
        process_peak_rss_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        cpu_affinity=sorted(os.sched_getaffinity(0)))
    output = args.output/'result.json'
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(output, flush=True)
    if result['status'] != 'passed':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()

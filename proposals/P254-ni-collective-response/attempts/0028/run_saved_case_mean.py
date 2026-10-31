#!/usr/bin/env python3
"""Execute actual original matched means without changing historical CLI bytes.

The existing source-qualified case, original-column midpoint and contraction
reuse APIs own the physics. Reduce/evolve/refine share one original immutable
allocation; changing CPU placement does not reset its wall clock or actions.
No new DFT, refilling, positive-time symmetry or clamped substitute is used.
"""
import argparse
import json
import math
import os
from pathlib import Path
import sys
from time import time
from uuid import uuid4

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

from mixed_material_driver_support import (
    checked_json, check_loaded_sources, owned_directory, source_closure)
from run_mixed_q_material import watchdog
from saved_case_mean_execution import execute, seed_sources
from saved_case_mean_owner import (
    SCHEMA, CAP_NAMES, execution_caps, original_owner, continuation_budget,
    terminal_signals, finish_receipt, isolated_output, reap_owned_workers)
from substrate_framework.streamed_material_state import file_binding


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument('--mode', choices=('reduce', 'evolve', 'refine'), required=True)
    result.add_argument('--case', required=True)
    for name in ('initialization', 'qualification'):
        result.add_argument('--'+name, type=Path, required=True)
        result.add_argument('--'+name+'-sha256', required=True)
    for name in ('oracle', 'column-oracle'):
        result.add_argument('--'+name, type=Path, required=True)
    result.add_argument('--reuse-initial-contractions', type=Path)
    result.add_argument('--reuse-initial-sha256')
    result.add_argument('--output', type=Path, required=True)
    result.add_argument('--resume', action='store_true')
    result.add_argument('--cpus', nargs='+', type=int, choices=range(6, 11), required=True)
    for name in CAP_NAMES:
        result.add_argument('--'+name.replace('_', '-'),
            type=float if name == 'max_seconds' else int, required=True)
    result.add_argument('--column-dimension', type=int, required=True)
    result.add_argument('--column-workspace-bytes', type=int, required=True)
    result.add_argument('--column-error-target', type=float, required=True)
    result.add_argument('--interval', type=float, default=.001)
    result.add_argument('--steps', type=int, default=1)
    result.add_argument('--max-steps', type=int, default=2)
    result.add_argument('--max-iterations', type=int, default=4)
    result.add_argument('--tolerance', type=float, default=1e-9)
    result.add_argument('--reality-tolerance', type=float, default=1e-10)
    result.add_argument('--block-size', type=int, default=8)
    return result


def inputs_and_request(args, sources, caps):
    if (not args.cpus or len(set(args.cpus)) != len(args.cpus)
            or not set(args.cpus) <= os.sched_getaffinity(0)):
        raise ValueError('unique explicitly available allocated CPUs required')
    if (not math.isfinite(args.interval) or args.interval <= 0
            or not 0 < args.steps <= args.max_steps or args.max_iterations < 1
            or not math.isfinite(args.tolerance) or args.tolerance <= 0 or args.block_size < 1):
        raise ValueError('positive original mean interval, step and solver limits required')
    if args.mode == 'refine' and args.max_steps < 2:
        raise ValueError('refine requires the declared two-step branch')
    if bool(args.reuse_initial_contractions) != bool(args.reuse_initial_sha256):
        raise ValueError('original contraction reuse requires its explicit path and digest')
    inputs = {name: checked_json(getattr(args, name), getattr(args, name+'_sha256'))[0]
              for name in ('initialization', 'qualification')}
    inputs.update(oracle=file_binding(args.oracle), column_oracle=file_binding(args.column_oracle),
        initial_reuse=(checked_json(args.reuse_initial_contractions, args.reuse_initial_sha256)[0]
                       if args.reuse_initial_contractions else None))
    numeric = {name: getattr(args, name) for name in (
        'interval', 'steps', 'max_steps', 'max_iterations', 'tolerance',
        'reality_tolerance', 'block_size', 'column_dimension',
        'column_workspace_bytes', 'column_error_target')}
    request = dict(schema=SCHEMA+'-request', case=args.case, inputs=inputs,
        numeric_plan=numeric, permitted_modes=['reduce', 'evolve', 'refine'],
        source_bindings=sources, caps=caps, blas_threads=1,
        phase_policy='Mode and CPU placement are invocation choices; original inputs/numerics/caps/clock remain fixed')
    return inputs, request


def run_owned(args, directory, invocation, started_unix, result):
    budget = None
    try:
        caps = execution_caps(args)
        sources = source_closure(seed_sources())
        inputs, request = inputs_and_request(args, sources, caps)
        owner, binding = original_owner(directory, request, resume=args.resume,
                                        started_unix=started_unix)
        result.update(owner=binding, original_started_unix=owner['started_unix'],
                      source_bindings=sources, explicit_inputs=inputs)
        budget = continuation_budget(directory, owner)
        os.sched_setaffinity(0, {args.cpus[0]})
        with watchdog(budget):
            execute(args, directory, inputs, sources, budget, result)
            check_loaded_sources(sources)
        result.update(status='complete', material_trajectory=bool(result.get('branches')))
    except BaseException as error:
        result.update(status='interrupted' if isinstance(error, (InterruptedError, KeyboardInterrupt)) else 'partial',
                      error=dict(type=type(error).__name__, message=str(error)[:16384]))
    finally:
        result['owned_worker_cleanup'] = reap_owned_workers()
        binding = finish_receipt(directory, invocation, result, budget)
    return binding


def main(arguments=None):
    started_unix = time()
    args = parser().parse_args(arguments)
    invocation = uuid4().hex
    result = dict(schema='P254-joint-partitioned-material-v1', status='partial', mode=args.mode,
        invocation=invocation, invocation_started_unix=started_unix, resume=args.resume,
        cpus=args.cpus, blas_threads=1, material_trajectory=False,
        requested_interval=args.interval if args.mode != 'reduce' else 0.,
        scope='All original slots, original physical mean; not whole local-q/noise/coherence or a coupled lattice ranking',
        numerical_accuracy_certified=False, automatic_restart=False)
    try:
        directory = isolated_output(args, ROOT)
        with terminal_signals(), owned_directory(directory, args.resume):
            binding = run_owned(args, directory, invocation, started_unix, result)
        print(json.dumps(dict(status=result['status'], receipt=binding)), flush=True)
    except BaseException as error:
        result.update(status='partial', receipt_location='stderr; owned output unavailable',
            receipt_error=dict(type=type(error).__name__, message=str(error)[:16384]))
        print(json.dumps(result, allow_nan=False), file=sys.stderr, flush=True)
    return 0 if result['status'] == 'complete' else 2


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        raise SystemExit(main())

#!/usr/bin/env python3
"""Finish the original material pair with an explicit cumulative wall allocation.

Reuse the original execute(), immutable owner, saved contexts, work identity and
append-only action ledger. Only the separately bound cumulative wall ceiling
changes; no physical propagation implementation or numerical plan is replaced.
"""
import argparse
import json
import os
from pathlib import Path
import sys
from time import time
from uuid import uuid4

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'

from mixed_material_driver_support import check_loaded_sources, exact_publish, owned_directory
from run_source_local_q_material import execute
from run_mixed_q_material import watchdog
from saved_case_mean_owner import reap_owned_workers
from source_local_continuation_allocation import SCHEMA, load_allocation
from source_local_material_owner import (
    context_bindings, encode, error_record, finish_receipt, material_budget,
    reserve_terminal_and_ledger, terminal_signals)


def continue_owned(allocation, grant, directory, owner, args, sources, result):
    budget, reserved = None, False
    with owned_directory(directory, True):
        try:
            budget = material_budget(directory, owner, grant['caps'])
            before = budget.snapshot()
            from substrate_framework.streamed_material_state import authenticate
            previous = json.loads(authenticate(grant['previous_terminal']).read_text())
            if any(before[name] < previous['budget'][name]
                   for name in ('hamiltonian_actions', 'action_columns')):
                raise ValueError('original charged actions and columns cannot disappear')
            reserve_terminal_and_ledger(budget)
            reserved = True
            copy = exact_publish(directory/'allocations'/allocation['sha256']/'allocation.json',
                                 grant, budget)
            contexts = context_bindings(args, directory)
            result.update(allocation=allocation, owned_allocation_copy=copy,
                previous_terminal=grant['previous_terminal'], owner=grant['original_owner'],
                original_started_unix=owner['started_unix'], cumulative_caps=grant['caps'],
                budget_before=before, source_bindings=sources,
                continuation_source_bindings=grant['source_bindings'],
                explicit_source_local_contexts=contexts,
                explicit_inputs=owner['request']['explicit_inputs'])
            executor = directory/'scientific-work'/'executor'
            if executor.exists():
                # Only an empty, released scratch container can be retired.
                # rmdir refuses retained scientific files or a nonempty executor.
                executor.rmdir()
                result['removed_empty_executor_directory'] = str(executor)
            os.sched_setaffinity(0, {args.cpus[0]})
            check_loaded_sources(grant['source_bindings'])
            with watchdog(budget):
                from source_local_material_runtime import material_inputs
                inputs = material_inputs(args)
                execute(args, inputs, sources, budget, directory, contexts, result)
                check_loaded_sources(grant['source_bindings'])
            result['status'] = 'complete'
        except BaseException as error:
            result.update(status='failed', error=error_record(error))
        finally:
            result['owned_worker_cleanup'] = reap_owned_workers()
            result['finished_unix'] = time()
            if not reserved:
                raise RuntimeError('terminal allocation unavailable; stderr only')
            receipt = finish_receipt(directory, result['invocation'], result, budget)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grant', type=Path, required=True)
    parser.add_argument('--grant-sha256', required=True)
    parser.add_argument('--cpus', nargs='+', type=int, required=True)
    cli = parser.parse_args()
    result = dict(schema=SCHEMA+'-result', status='failed', mode='pair',
        invocation=uuid4().hex, invocation_started_unix=time(), resume=True,
        cpus=cli.cpus, blas_threads=1, covariance_executed=False,
        original_owner_rewritten=False, whole_spatial_reconstruction=False,
        required_q_for_whole_local=64, numerical_accuracy_certified=False,
        heat_claim=False, coherence_claim=False, production_authorized=False,
        automatic_restart=False)
    try:
        allocation, grant, directory, owner, args, sources = load_allocation(
            cli.grant, cli.grant_sha256, cli.cpus)
        with terminal_signals():
            receipt = continue_owned(allocation, grant, directory, owner, args, sources, result)
        print(json.dumps(dict(status=result['status'], receipt=receipt)), flush=True)
    except BaseException as error:
        result.update(status='failed', receipt_location='stderr', receipt_error=error_record(error))
        print(json.dumps(encode(result), allow_nan=False), file=sys.stderr, flush=True)
    return 0 if result['status'] == 'complete' else 2


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        raise SystemExit(main())

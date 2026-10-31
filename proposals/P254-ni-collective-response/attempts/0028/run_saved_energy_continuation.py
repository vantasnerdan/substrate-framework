#!/usr/bin/env python3
"""Finish a saved energy readout on its original clock, ledger and scientific work.

An explicit bound cumulative allocation extends the exhausted original consumer;
its original owner/failed receipt are never rewritten. No state initialization,
propagation or global moment replay is performed. Existing QWork restores the
original prepared context in its own namespace and retains every charged action.
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
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'

from mixed_energy_block_work import EnergyBlockWork
from mixed_energy_saved_owner import producer_exclusion, terminal_receipt
from mixed_energy_saved_runtime import build_saved_runtime
from mixed_material_driver_support import (
    ContinuationBudget, checked_json, check_loaded_sources, exact_publish,
    owned_directory, owned_tree, source_closure)
from mixed_material_energy import mixed_energy_readout
from run_mixed_energy_readout import checked_energy_report
from run_mixed_q_material import watchdog
from saved_case_mean_owner import terminal_signals, reap_owned_workers
from substrate_framework.streamed_material_state import authenticate, file_binding

SCHEMA = 'P254-saved-energy-cumulative-continuation-v1'
TERMINAL_BYTES = 2*1024**2


def load_grant(path, sha, cpus):
    binding, grant = checked_json(path, sha)
    if grant.get('schema') != SCHEMA+'-allocation':
        raise ValueError('explicit original-clock cumulative energy allocation required')
    directory = Path(grant['consumer_directory']).resolve()
    owner_binding = file_binding(directory/'owner.json')
    owner = json.loads(authenticate(owner_binding).read_text())
    authority = json.loads((directory/'owner-authority.json').read_text())
    previous = json.loads(authenticate(grant['previous_terminal']).read_text())
    if (grant['original_consumer_owner'] != owner_binding
            or authority != dict(schema='P254-saved-native-mixed-energy-consumer-v1-owner-authority', owner=owner_binding)
            or previous.get('consumer_owner') != owner_binding
            or previous.get('complete') is not False or previous.get('status') != 'failed'
            or grant['original_started_unix'] != owner['started_unix']
            or previous['consumer_original_started_unix'] != owner['started_unix']):
        raise ValueError('unchanged original failed consumer/owner/start required')
    if (grant['admission'] != previous['admission'] or grant['state'] != previous['state']
            or dict(path=grant['mean_ledger']['path'], sha256=grant['mean_ledger']['sha256'])
                != owner['request']['inputs']['mean_energy_ledger']):
        raise ValueError('unchanged original admission, endpoint and saved mean required')
    if Path(grant['energy_subdirectory']).name != grant['energy_subdirectory']:
        raise ValueError('owned energy subdirectory must be one relative name')
    caps, old = grant['caps'], owner['request']['caps']
    if (not cpus or len(set(cpus)) != len(cpus)
            or cpus != grant['cpus'] or any(type(cpu) is not int for cpu in cpus)
            or not set(cpus) <= os.sched_getaffinity(0) or not set(cpus) <= {6, 7, 8}
            or type(caps['max_actions']) is not int or caps['max_actions'] <= old['max_actions']
            or not math.isfinite(caps['max_seconds']) or not old['max_seconds'] < caps['max_seconds']
            or any(caps[name] != old[name] for name in ('max_output_bytes', 'max_memory_bytes', 'worker_memory_bytes'))):
        raise ValueError('explicit cumulative action/time extension and unchanged memory/output ceilings required')
    proof = json.loads(authenticate(grant['native_block_comparison']).read_text())
    sample = proof.get('result', {})
    if (proof.get('status') != 'complete' or sample.get('status') != 'complete'
            or sample.get('comparison', {}).get('passed') is not True
            or sample.get('parent_model_unchanged') is not True
            or sample.get('new_readout_scope', {}).get('block_size') != grant['energy_column_block_size']
            or sample.get('inputs', {}).get('original_state') != grant['state']
            or sample.get('inputs', {}).get('original_admission') != grant['admission']):
        raise ValueError('matching actually completed original Native block comparison required')
    numerical = json.loads(authenticate(grant['sensitive_block_action_oracle']).read_text())
    if (numerical.get('schema') != 'P254-saved-energy-block-oracle-v1'
            or numerical.get('status') != 'passed'
            or numerical.get('original_prepared') != sample['inputs']['original_prepared']
            or int(numerical['cumulative_consumer_action_upper_bound_with_donors']) > caps['max_actions']):
        raise ValueError('source-bound exact block/action proof required')
    for source in numerical['source_bindings']:
        authenticate(source)
    if grant['source_bindings'] != source_closure([file_binding(__file__)]):
        raise ValueError('cumulative implementation differs from allocated source bytes')
    return binding, grant, directory, owner, previous


def execute(grant, directory, owner, cpus, budget, result):
    admission = json.loads(authenticate(grant['admission']).read_text())
    inputs = owner['request']['inputs']
    if (admission['producer_directory'] != inputs['producer_directory']
            or admission['owner_binding']['sha256'] != inputs['producer_owner_sha256']
            or admission['owner_binding'] != file_binding(Path(inputs['producer_directory'])/'owner.json')):
        raise ValueError('unchanged original producer authority required')
    if admission['state_binding'] != grant['state']:
        raise ValueError('unchanged genuine physical endpoint required')
    state = json.loads(authenticate(grant['state']).read_text())
    with producer_exclusion(admission['producer_directory'], state):
        dynamics, original_work, runtime = build_saved_runtime(admission, budget=budget,
            cpus=cpus, work_directory=directory/'work')
        dynamics.energy_mean_readout = grant['mean_ledger']
        work = EnergyBlockWork(original_work.directory, budget, cpus, original_work.identity,
            dynamics=dynamics, column_block_size=grant['energy_column_block_size'],
            prior_request=grant['original_bare_request'], prior_receipts=grant['retained_bare_receipts'])
        if dynamics.block_size != grant['registered_dynamics_block_size']:
            raise ValueError('registered mean/dynamics block size must remain unchanged')
        result.update(runtime_evidence=runtime, original_work_identity=work.identity,
            original_state=grant['state'], retained_scalar_sources=len(work.donors),
            energy_only_column_block_size=work.column_block_size)
        check_loaded_sources(grant['source_bindings'])
        report_binding = mixed_energy_readout(grant['state'], dynamics=dynamics, work=work,
            directory=directory/grant['energy_subdirectory'], budget=budget)
        report = checked_energy_report(report_binding, admission)
        if report.get('complete') is not True or report.get('status') != 'complete':
            raise RuntimeError('complete genuine eleven-term positive energy report required')
        check_loaded_sources(grant['source_bindings'])
        result.update(status='complete', complete=True, energy_report=report_binding,
                      time_ry_inverse=report['time_ry_inverse'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--grant', type=Path, required=True)
    parser.add_argument('--grant-sha256', required=True)
    parser.add_argument('--cpus', nargs='+', type=int, required=True)
    args = parser.parse_args()
    result = dict(schema=SCHEMA, status='failed', complete=False, invocation=uuid4().hex,
        invocation_started_unix=time(), propagation_performed=False,
        raw_initialization_performed=False, original_owner_rewritten=False,
        scope='Actual complex mixed eleven-term response-energy coefficient at the saved physical endpoint; not finite-dose energy, heat, noise, spatial closure or a conserved-energy certificate')
    budget, reserved = None, False
    try:
        allocation, grant, directory, owner, previous = load_grant(args.grant, args.grant_sha256, args.cpus)
        with terminal_signals(), owned_directory(directory, True):
            try:
                caps = grant['caps']
                budget = ContinuationBudget(directory, started_at=owner['started_unix'],
                    max_actions=caps['max_actions'], max_seconds=caps['max_seconds'],
                    max_bytes=caps['max_output_bytes'], max_memory_bytes=caps['max_memory_bytes'],
                    worker_memory_bytes=caps['worker_memory_bytes'], existing_bytes=sum(owned_tree(directory).values()))
                before = budget.snapshot()
                if before['hamiltonian_actions'] < previous['consumer_budget']['hamiltonian_actions']:
                    raise ValueError('original charged actions cannot disappear')
                remaining = max(0, budget.max_actions*budget.RECORD.size-budget.ledger.stat().st_size)
                budget.reserve(TERMINAL_BYTES+remaining)
                reserved = True
                extension = exact_publish(directory/'allocations'/allocation['sha256']/'allocation.json', grant, budget)
                result.update(allocation=allocation, owned_allocation_copy=extension,
                    original_consumer_owner=grant['original_consumer_owner'], previous_terminal=grant['previous_terminal'],
                    original_started_unix=owner['started_unix'], cumulative_caps=caps,
                    budget_before=before, cpus=args.cpus, blas_threads=1,
                    diagnostic_cost_scope='Native comparison is separately charged in its original diagnostic ledger; no costs erased')
                os.sched_setaffinity(0, {args.cpus[0]})
                with watchdog(budget):
                    execute(grant, directory, owner, args.cpus, budget, result)
            except BaseException as error:
                result.update(status='failed', complete=False,
                              error=type(error).__name__+': '+str(error)[:16384])
            finally:
                result['owned_worker_cleanup'] = reap_owned_workers()
                result['finished_unix'] = time()
                result['budget'] = budget.snapshot() if budget is not None else None
                if not reserved:
                    raise RuntimeError('terminal allocation unavailable; stderr only')
                receipt = terminal_receipt(directory, result['invocation'], result)
        print(json.dumps(dict(status=result['status'], complete=result['complete'], receipt=receipt)), flush=True)
    except BaseException as error:
        result.update(status='failed', complete=False, receipt_location='stderr',
                      receipt_error=type(error).__name__+': '+str(error)[:16384])
        print(json.dumps(result, allow_nan=False), file=sys.stderr, flush=True)
    return 0 if result['complete'] else 2


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        raise SystemExit(main())

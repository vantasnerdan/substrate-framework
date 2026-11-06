#!/usr/bin/env python3
"""Read the actual eleven-term mixed energy of one authentic saved boundary.

--output is a NEW isolated consumer directory, never the producer's directory.
--resume retains the exact request, caps, source bytes, independent consumer
identity, original consumer wall deadline and cumulative pre-action admissions.
The producer keeps its separate original identity, start/deadline and costs.
An explicit --producer-source-archive/--producer-source-archive-sha256 pair
authenticates retained historical source bytes after a source cutover. It never
selects code to execute; current imported code remains separately source-bound.

No raw initialization, RHS, propagation, interpolation or producer-work mutation
is performed. An initial-boundary readout is an engineering result at t0, not a
positive-time propagated material result. Only an actual complete eleven-term
mixed_energy_readout report permits a successful exit.
"""
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import resource
import sys
from time import perf_counter, process_time, time
from uuid import uuid4

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'

from mixed_energy_saved_boundary import admit_saved_boundary
from mixed_energy_saved_runtime import build_saved_runtime
from mixed_energy_saved_owner import (
    SCHEMA, consumer_budget, consumer_owner, energy_attempt_evidence, execution_caps, execution_costs,
    execution_snapshot, isolated_output, producer_cost_snapshot, producer_exclusion,
    reserve_terminal_and_ledger, terminal_receipt, terminal_signals)
from mixed_material_driver_support import (
    check_loaded_sources, exact_publish, nested_bindings, owned_directory,
    owned_tree, source_closure)
from mixed_material_energy import SCHEMA as ENERGY_SCHEMA, TERMS, mixed_energy_readout
from run_mixed_q_material import checked_json, watchdog
from substrate_framework.streamed_material_state import authenticate, file_binding


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument('--producer', type=Path, required=True)
    result.add_argument('--producer-owner-sha256', required=True)
    result.add_argument('--producer-source-archive', type=Path)
    result.add_argument('--producer-source-archive-sha256')
    result.add_argument('--mean-energy-ledger', type=Path, required=True)
    result.add_argument('--mean-energy-ledger-sha256', required=True)
    result.add_argument('--state', type=Path, required=True)
    result.add_argument('--state-sha256', required=True)
    result.add_argument('--output', type=Path, required=True,
                        help='new exclusively owned consumer directory, outside the repository')
    result.add_argument('--cpus', nargs='+', type=int, choices=range(6, 11), required=True)
    result.add_argument('--max-actions', type=int, required=True)
    result.add_argument('--max-seconds', type=float, required=True)
    result.add_argument('--max-output-bytes', type=int, required=True)
    result.add_argument('--max-memory-bytes', type=int, required=True)
    result.add_argument('--worker-memory-bytes', type=int, required=True)
    result.add_argument('--resume', action='store_true')
    return result


def requested_inputs(args):
    archive = (None if args.producer_source_archive is None
               and args.producer_source_archive_sha256 is None else dict(
                   path=str(args.producer_source_archive.resolve()) if args.producer_source_archive else None,
                   sha256=args.producer_source_archive_sha256))
    return dict(producer_directory=str(args.producer.resolve()),
        producer_owner_sha256=args.producer_owner_sha256,
        mean_energy_ledger=dict(path=str(args.mean_energy_ledger.resolve()),
                               sha256=args.mean_energy_ledger_sha256),
        producer_source_archive=archive,
        state=dict(path=str(args.state.resolve()), sha256=args.state_sha256))


def state_binding(args):
    for value in (args.producer_owner_sha256, args.state_sha256):
        if len(value) != 64 or any(character not in '0123456789abcdef' for character in value):
            raise ValueError('explicit lowercase SHA256 owner and state digests required')
    binding = file_binding(args.state.resolve())
    if binding['sha256'] != args.state_sha256:
        raise ValueError('explicit requested state SHA256 differs from original bytes')
    return binding


def source_archive_binding(args):
    path, sha = args.producer_source_archive, args.producer_source_archive_sha256
    if path is None and sha is None:
        return None
    if (path is None or sha is None or len(sha) != 64
            or any(character not in '0123456789abcdef' for character in sha)):
        raise ValueError('source archive requires both its path and explicit lowercase SHA256')
    binding = file_binding(path.resolve())
    if binding['sha256'] != sha:
        raise ValueError('explicit source archive SHA256 differs from retained bytes')
    return binding


def admit_inputs(args, output, budget):
    admission = admit_saved_boundary(args.producer.resolve(),
        args.producer_owner_sha256, state_binding(args), source_archive=source_archive_binding(args))
    # The admission is immutable consumer metadata, not a rewritten Native state.
    for binding in nested_bindings(admission):
        if Path(binding['path']).resolve().is_relative_to(output):
            raise ValueError('borrowed producer/input authority overlaps the consumer namespace')
    binding = exact_publish(output/'admission.json', admission, budget)
    return admission, binding


def checked_energy_report(binding, admission):
    """A successful invocation must bind the genuine complete eleven-term report."""
    report = json.loads(authenticate(binding).read_text())
    if (report.get('schema') != ENERGY_SCHEMA
            or report.get('input') != admission['state_binding']):
        raise ValueError('actual energy report does not bind the authenticated requested state')
    if report.get('status') != 'complete':
        return report
    if (report.get('complete') is not True
            or report.get('time_ry_inverse') != admission['state']['time']
            or report.get('model_identity') != admission['state']['identity']
            or set(report.get('components_global_ry', {})) != set(TERMS)):
        raise ValueError('success requires the actual complete eleven-term physical-boundary energy')
    for value in [*report['components_global_ry'].values(), report.get('total_global_ry')]:
        if (not isinstance(value, list) or len(value) != 2
                or any(type(part) not in (int, float) or not math.isfinite(part) for part in value)):
            raise ValueError('actual complete report requires finite complex eleven-term values and total')
    request = json.loads(authenticate(report['request']).read_text())
    if request.get('input') != admission['state_binding']:
        raise ValueError('energy request authority differs from the admitted raw boundary')
    return report


def execute(args, directory, result, budget, sources, stack):
    """Admit first, then assemble the exact isolated runtime and call the real API."""
    preliminary_state = json.loads(authenticate(state_binding(args)).read_text())
    result['producer_writer_policy'] = stack.enter_context(
        producer_exclusion(args.producer.resolve(), preliminary_state))
    mean_ledger, _ = checked_json(args.mean_energy_ledger, args.mean_energy_ledger_sha256)
    admission, admission_binding = admit_inputs(args, directory, budget)
    result.update(admission=admission_binding, state=admission['state_binding'],
        boundary_kind=admission['boundary_kind'], time_ry_inverse=admission['state']['time'],
        producer_owner=admission['owner_binding'], boundary_authority=admission['authority'])
    result['historical_producer_costs'] = producer_cost_snapshot(admission)
    result['scope']['positive_time_material_result'] = admission['boundary_kind'] == 'physical-endpoint'
    budget.check_processes(os.getpid(), [])
    dynamics, work, evidence = build_saved_runtime(admission, budget=budget,
        cpus=args.cpus, work_directory=directory/'work')
    dynamics.energy_mean_readout = mean_ledger
    result.update(runtime_evidence=evidence, qwork_identity=work.identity)
    if (not Path(work.directory).resolve().is_relative_to(directory/'work')
            or work.budget is not budget or dynamics.budget is not budget):
        raise ValueError('runtime must use this isolated consumer work namespace and actual budget')
    check_loaded_sources(sources)
    result['budget_before_energy'] = budget.snapshot()
    binding = mixed_energy_readout(admission['state_binding'], dynamics=dynamics,
        work=work, directory=directory/'energy', budget=budget)
    result['energy_report'] = binding
    report = checked_energy_report(binding, admission)
    result['energy_report_status'] = report.get('status')
    check_loaded_sources(sources)
    result.update(status='complete' if report.get('status') == 'complete' else 'unavailable',
        complete=report.get('status') == 'complete')
    if not result['complete']:
        result['reason'] = report.get('reason', 'real energy API returned no complete eleven-term result')


def owned_invocation(args, directory, result, started_unix):
    """Keep all admitted exits and resource failures inside the immutable receipt."""
    budget, before = None, None
    sources, work_before = [], execution_snapshot(directory/'work')
    receipt = None
    # Exclusion is released only after the final consumer receipt is durable.
    with ExitStack() as stack:
        try:
            caps = execution_caps(args)
            os.sched_setaffinity(0, {args.cpus[0]})
            sources = source_closure([file_binding(HERE/'run_mixed_energy_readout.py')])
            request = dict(schema=SCHEMA, inputs=requested_inputs(args), cpus=args.cpus,
                caps=caps, source_bindings=sources, output_directory=str(directory),
                operation='propagation-free actual eleven-term mixed_energy_readout')
            owner, binding, identity = consumer_owner(directory, request,
                resume=args.resume, started_unix=started_unix)
            result.update(consumer_owner=binding, consumer_identity=identity,
                consumer_original_started_unix=owner['started_unix'],
                consumer_deadline_unix=owner['started_unix']+caps['max_seconds'],
                source_bindings=sources, execution_caps=caps)
            budget = consumer_budget(directory, owner, caps)
            before = budget.snapshot()
            result['budget_at_invocation_start'] = before
            reserve_terminal_and_ledger(budget)
            result['terminal_and_ledger_space_reserved'] = True
            with watchdog(budget):
                execute(args, directory, result, budget, sources, stack)
        except BaseException as error:
            result.update(status='failed', complete=False,
                          error=type(error).__name__+': '+str(error))
        finally:
            result['consumer_budget'] = budget.snapshot() if budget is not None else None
            if budget is not None and before is not None:
                result['current_invocation_admitted_actions'] = budget.actions-before['hamiltonian_actions']
                result['current_invocation_admitted_columns'] = budget.columns-before['action_columns']
            try:
                result['actual_energy_api_attempts'] = energy_attempt_evidence(directory/'energy')
                result['consumer_execution_costs'] = execution_costs(directory/'work', work_before)
                result['owned_bytes_before_terminal_receipt'] = sum(owned_tree(directory).values())
            except BaseException as error:
                result.update(status='failed', complete=False,
                    cost_accounting_error=type(error).__name__+': '+str(error))
            result.update(terminal_at=datetime.now(timezone.utc).isoformat(),
                invocation_wall_seconds=perf_counter()-result.pop('_wall_start'),
                invocation_parent_cpu_seconds=process_time()-result.pop('_cpu_start'),
                parent_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                recovery='Only identical --resume reopens this consumer owner, request, caps, '
                         'sources, original deadline and durable action ledger. Completed isolated '
                         'consumer work may be adopted; producer caches/costs are not relabelled.')
            receipt = terminal_receipt(directory, result['invocation'], result)
    return receipt


def main(argv=None):
    args = parser().parse_args(argv)
    started_unix = time()
    result = dict(schema=SCHEMA, status='failed', complete=False, invocation=uuid4().hex,
        requested_inputs=requested_inputs(args), cpus=args.cpus,
        _wall_start=perf_counter(), _cpu_start=process_time(),
        terminal_and_ledger_space_reserved=False,
        scope=dict(quantity='complex mixed d2E/(dlambda deta), physical eta=lambda*, no half',
            propagation_performed=False, raw_initialization_performed=False,
            whole_spatial_reconstruction=False, numerical_accuracy_certified=False,
            production_authorization=False, positive_time_material_result=False,
            not_claimed=['finite-dose energy', 'energy transfer/destination', 'heat/noise power',
                         'conserved energy certificate', 'full spatial E/D/C closure']),
        budget_scope='Independent consumer original-start wall deadline and cumulative durable '
                     'pre-action admissions. Original producer identity/clock/caps/costs stay '
                     'historical and separate. Only isolated consumer output is charged.',
        terminal_receipt_scope='Every catchable owned execution exit is recorded; SIGKILL or '
                               'host loss cannot write a terminal receipt, but the owner and '
                               'pre-action ledger remain durable for identical resume.')
    receipt = None
    try:
        args.output = isolated_output(args, ROOT)
        with terminal_signals(), owned_directory(args.output, args.resume) as directory:
            receipt = owned_invocation(args, directory, result, started_unix)
    except BaseException as error:
        result.update(status='failed', complete=False,
                      terminal_error=type(error).__name__+': '+str(error))
    summary = dict(schema=SCHEMA, status=result['status'], complete=result['complete'],
        output_directory=str(args.output), invocation_receipt=receipt,
        energy_report=result.get('energy_report'), error=result.get('error'),
        terminal_error=result.get('terminal_error'))
    if receipt is None:
        # Never write through an unacquired owner lock or into an unsafe namespace.
        summary['unpersisted_failure'] = {key: value for key, value in result.items()
                                          if not key.startswith('_')}
    print(json.dumps(summary, allow_nan=False, sort_keys=True))
    return 0 if result['complete'] else 2


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        raise SystemExit(main())

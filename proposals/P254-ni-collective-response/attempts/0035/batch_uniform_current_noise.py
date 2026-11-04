#!/usr/bin/env python3
"""Serial, explicitly pinned and finite A35 current-noise production.

Reuses the A34 real-subprocess worker unchanged (including its A34-named
operational outcome schema); only the serial loop and A35 science admission are
owned here. Resume admits committed cases, never unfinished action arrays. One
attempt per pending record per invocation, no automatic retry after failure.
Exit: requested subset complete=0, valid incomplete subset=2, invalid/failed=1.
"""
import argparse
import fcntl
import json
import math
import os
from pathlib import Path
import shutil
import signal
import sys
import threading
import time
import uuid

from aggregate_uniform_current_noise import (
    A34, BATCH_SCHEMA, SCOPE, add_input_arguments, authenticate, binding, load_case,
    load_context, load_json, save_json)
from batch_mode_noise import BLAS_ENV, execute_case, parse_indices, positive_seconds, utc_now
from checkpoint_directories import make_directory

MEASUREMENTS = ('elapsed_seconds', 'user_cpu_seconds', 'system_cpu_seconds', 'peak_rss_bytes')
OUTCOME_SCHEMA = 'P254-A34-mode-noise-subprocess-outcome-v1'


def read_exit(generation, record):
    """Authenticate actual successful durable exit and streams before science."""
    intent = load_json(generation / 'intent.json')
    outcome = load_json(generation / 'outcome.json')
    if (outcome['schema'] != OUTCOME_SCHEMA or type(outcome['original_record']) is not int
            or outcome['original_record'] != record or intent['original_record'] != record
            or type(outcome['exit_code']) is not int or outcome['exit_code'] != 0
            or outcome['termination'] is not None or outcome['error'] is not None
            or outcome['argv'] != intent['argv'] or outcome['cpu'] != intent['cpu']):
        raise ValueError('generation has no successful durable real-subprocess exit')
    authenticate(outcome['stdout'], generation / 'stdout.log')
    authenticate(outcome['stderr'], generation / 'stderr.log')
    return intent, outcome


def reusable(root, record, current, context, rejected):
    directory = root / 'records' / ('record-%03d' % record)
    if not directory.exists():
        return None
    compatible = []
    for generation in sorted(directory.glob('generation-*')):
        done = generation / 'case' / 'done.json'
        if not done.exists():
            continue
        try:
            intent, _ = read_exit(generation, record)
            for key in ('inputs', 'contract', 'scientific_compatibility_identity'):
                if intent['binding_contract'][key] != current[key]:
                    raise ValueError('generation original inputs/plan/current scientific identity differs')
            case = load_case(done, context, record)
            compatible.append(dict(original_record=record, done=case['done'],
                                   result=case['result'], generation=str(generation),
                                   reused=True, work=case['payload']['work']))
        except Exception as error:
            rejected.append(dict(original_record=record, generation=str(generation), error=repr(error)))
    if len(compatible) > 1:
        raise ValueError('multiple compatible committed generations for record %d' % record)
    return compatible[0] if compatible else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_input_arguments(parser)
    parser.add_argument('--python', type=Path, default=Path(sys.executable))
    parser.add_argument('--output', type=Path, required=True,
                        help='preserved batch root; each invocation creates a new run generation')
    parser.add_argument('--records', type=lambda s: parse_indices(s, 128), default=list(range(128)))
    parser.add_argument('--cpu', type=int, required=True,
                        help='one explicitly granted allowed logical CPU; exactly one worker')
    parser.add_argument('--wall-seconds', type=positive_seconds, required=True)
    parser.add_argument('--case-timeout-seconds', type=positive_seconds, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    make_directory(root, exist_ok=True)
    run = root / 'runs' / ('run-' + uuid.uuid4().hex)
    make_directory(run)
    started = time.monotonic()
    deadline = started + args.wall_seconds
    stop = threading.Event()
    completed, failed, rejected, outcomes = [], [], [], []
    launched, interrupted = set(), []
    error, context = None, None
    old_handlers = {}
    events = (run / 'progress.jsonl').open('xb')
    lock = (root / '.batch.lock').open('a+b')

    def cost():
        def total(field):
            values = [x[field] for x in outcomes]
            return None if any(x is None for x in values) else sum(values)

        rss = [x['peak_rss_bytes'] for x in outcomes]
        return dict(subprocess_elapsed_seconds=total('elapsed_seconds'),
                    user_cpu_seconds=total('user_cpu_seconds'),
                    system_cpu_seconds=total('system_cpu_seconds'),
                    peak_case_rss_bytes=None if any(x is None for x in rss) else max(rss, default=0),
                    measured_subprocesses=sum(not x['missing_measurements'] for x in outcomes),
                    unmeasured_subprocesses=sum(bool(x['missing_measurements']) for x in outcomes),
                    new_successful_current_action_evaluations=sum(
                        x['work']['operator_work']['current_action_evaluations']
                        for x in completed if not x['reused']),
                    new_successful_adjoint_reuses_by_hermiticity=sum(
                        x['work']['operator_work']['adjoint_reuses_by_hermiticity']
                        for x in completed if not x['reused']),
                    failed_generation_current_action_count='unknown' if failed else 0,
                    reused_generation_costs_included=False)

    def emit(event, **fields):
        value = dict(event=event, utc=utc_now(), elapsed_seconds=time.monotonic() - started,
                     completed_original_records=sorted(x['original_record'] for x in completed),
                     cost=cost(), **fields)
        raw = (json.dumps(value, allow_nan=False) + '\n').encode()
        events.write(raw)
        events.flush()
        os.fsync(events.fileno())
        try:
            print(raw.decode().rstrip('\n'), flush=True)
        except BrokenPipeError:
            pass

    def consume(record, generation, worker_error):
        def reject(exception):
            item = dict(original_record=record, generation=str(generation), error=repr(exception))
            failed.append(item)
            save_json(generation / 'rejected.json', item)
            emit('case-failed', **item)

        try:
            outcome = load_json(generation / 'outcome.json')
            if outcome['original_record'] != record:
                raise ValueError('durable outcome belongs to another record')
        except BaseException as exception:
            outcomes.append(dict(original_record=record, exit_code=None,
                                 elapsed_seconds=None, user_cpu_seconds=None,
                                 system_cpu_seconds=None, peak_rss_bytes=None,
                                 missing_measurements=list(MEASUREMENTS),
                                 measurement_error=repr(exception)))
            reject(exception)
            raise
        for field in MEASUREMENTS:
            value = outcome.get(field)
            if (type(value) not in (int, float) or not math.isfinite(value) or value < 0):
                outcome[field] = None
        outcome['missing_measurements'] = [f for f in MEASUREMENTS if outcome[f] is None]
        outcomes.append(outcome)
        try:
            intent, _ = read_exit(generation, record)
            if intent['binding_contract'] != current:
                raise ValueError('durable launch binding contract differs')
            if worker_error is not None:
                raise worker_error
            # Durable actual exit and both streams are admitted before result readout.
            case = load_case(generation / 'case' / 'done.json', context, record)
            item = dict(original_record=record, done=case['done'], result=case['result'],
                        generation=str(generation), reused=False, work=case['payload']['work'])
        except BaseException as exception:
            reject(exception)
            raise
        completed.append(item)
        save_json(generation / 'accepted.json', item)
        emit('case-complete', original_record=record, done=case['done'])

    def on_signal(signum, _frame):
        interrupted.append(signum)
        stop.set()

    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if args.cpu < 0 or args.cpu not in os.sched_getaffinity(0):
            raise ValueError('explicit CPU outside inherited actual CPU affinity')
        taskset = shutil.which('taskset')
        if taskset is None:
            raise ValueError('Linux taskset required for explicit child CPU affinity')
        for sig in (signal.SIGINT, signal.SIGTERM):
            old_handlers[sig] = signal.signal(sig, on_signal)
        context = load_context(args.initial_result, args.coefficient_manifest, args.contract)
        current = dict(inputs=context['inputs'], contract=context['contract_binding'],
                       scientific_compatibility_identity=context['scientific_compatibility_identity'],
                       producer=binding(Path(__file__).with_name('export_uniform_current_noise.py')),
                       python=binding(args.python), python_command=str(args.python.absolute()),
                       runner=binding(__file__),
                       aggregator=binding(Path(__file__).with_name('aggregate_uniform_current_noise.py')),
                       reused_subprocess_worker=binding(A34 / 'batch_mode_noise.py'))
        save_json(run / 'intent.json', dict(schema='P254-A35-uniform-current-noise-batch-intent-v1',
                                           binding_contract=current, records=args.records,
                                           cpu=args.cpu, workers=1, wall_seconds=args.wall_seconds,
                                           case_timeout_seconds=args.case_timeout_seconds,
                                           environment_overrides=BLAS_ENV, argv=sys.argv,
                                           started_utc=utc_now(),
                                           operational_outcome_schema=OUTCOME_SCHEMA))
        emit('batch-start', requested_original_records=args.records, run=str(run))
        pending = []
        for record in args.records:
            case = reusable(root, record, current, context, rejected)
            if case is None:
                pending.append(record)
            else:
                completed.append(case)
                emit('case-reused', original_record=record, done=case['done'])
        for record in pending:
            if stop.is_set() or time.monotonic() >= deadline:
                break
            generation = root / 'records' / ('record-%03d' % record) / ('generation-' + uuid.uuid4().hex)
            argv = [taskset, '--cpu-list', str(args.cpu), current['python_command'],
                    current['producer']['path'], '--initial-result',
                    current['inputs']['initial_result']['path'], '--coefficient-manifest',
                    current['inputs']['coefficient_manifest']['path'], '--original-record',
                    str(record), '--output', str(generation / 'case')]
            intent = dict(binding_contract=current, original_record=record, started_utc=utc_now())
            emit('case-start', original_record=record, cpu=args.cpu, generation=str(generation), argv=argv)
            launched.add(record)
            worker_error = None
            try:
                # Deliberately synchronous: no executor, parallel scheduler or second child.
                execute_case(record, args.cpu, generation, argv, intent, stop, deadline,
                             args.case_timeout_seconds)
            except BaseException as exception:
                worker_error = exception
            consume(record, generation, worker_error)
        if interrupted:
            error = 'interrupted by signal(s): ' + repr(interrupted)
    except BaseException as exception:
        stop.set()
        error = repr(exception)
    finally:
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
        done_records = {item['original_record'] for item in completed}
        valid = not failed and error is None
        full = valid and done_records == set(range(128))
        subset_complete = valid and done_records == set(args.records)
        code = 0 if subset_complete else (2 if valid else 1)
        summary = dict(schema=BATCH_SCHEMA,
                       status='requested-subset-complete' if subset_complete else (
                           'incomplete' if valid else 'invalid'),
                       full_original128_coverage=full,
                       requested_subset_complete=subset_complete,
                       requested_original_records=args.records, launched_original_records=sorted(launched),
                       incomplete_original_records=sorted(set(args.records) - done_records),
                       missing_original_records=sorted(set(range(128)) - done_records),
                       completed=sorted(completed, key=lambda x: x['original_record']),
                       failed=failed, unstarted_original_records=sorted(set(args.records) - done_records - launched),
                       rejected_reuse_candidates=rejected, error=error, interrupted_signals=interrupted,
                       elapsed_seconds=time.monotonic() - started, cost=cost(),
                       scientific_compatibility_identity=(context['scientific_compatibility_identity']
                                                          if context is not None else None),
                       subprocess_outcomes=[{key: x.get(key) for key in
                           ('original_record', 'exit_code', *MEASUREMENTS,
                            'missing_measurements', 'measurement_error')} for x in outcomes],
                       exit_code=code, scope=SCOPE)
        try:
            receipt = save_json(run / 'summary.json', summary)
            emit('batch-exit', status=summary['status'], exit_code=code, summary=receipt, error=error)
        finally:
            events.close()
            lock.close()
    return code


if __name__ == '__main__':
    raise SystemExit(main())

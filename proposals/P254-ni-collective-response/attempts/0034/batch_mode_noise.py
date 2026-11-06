#!/usr/bin/env python3
"""Finite, explicitly pinned A34 original-sector force-noise batch.

No work launches on import. Supply a CPU list and actual wall/case grants after
the separately authorized pilot. Each producer gets a fresh nonexistent case
directory inside a preserved generation container. Resume skips only committed,
binding-compatible sectors; it never resumes an unfinished in-case history.
"""
import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
import fcntl
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import threading
import time
import uuid

from aggregate_mode_noise import (SCOPE, add_input_arguments, authenticate, binding,
                                  load_case, load_context, load_json,
                                  same_binding, save_json)
from checkpoint_directories import make_directory

BLAS_ENV = {key: '1' for key in (
    'OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
    'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS',
    'GOTO_NUM_THREADS', 'OMP_THREAD_LIMIT')}


def parse_indices(text, upper=None):
    values = []
    for field in text.split(','):
        if not field:
            raise argparse.ArgumentTypeError('empty index field')
        try:
            if '-' in field:
                left, right = field.split('-', 1)
                start, stop = int(left), int(right)
                if stop < start:
                    raise ValueError('descending range')
                values.extend(range(start, stop + 1))
            else:
                values.append(int(field))
        except ValueError as error:
            raise argparse.ArgumentTypeError('invalid index list: ' + text) from error
    if (not values or len(values) != len(set(values)) or min(values) < 0
            or (upper is not None and max(values) >= upper)):
        raise argparse.ArgumentTypeError('indices must be unique and in range')
    return sorted(values)


def positive_seconds(text):
    try:
        value = float(text)
    except ValueError as error:
        raise argparse.ArgumentTypeError('positive finite seconds required') from error
    if not math.isfinite(value) or value <= 0:
        raise argparse.ArgumentTypeError('positive finite seconds required')
    return value


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def execute_case(record, cpu, generation, argv, intent, stop, deadline, timeout):
    """Persist launch intent, raw streams and wait4 exit/cost before readout."""
    make_directory(generation)
    save_json(generation / 'intent.json', dict(intent, argv=argv, cpu=cpu,
                                             created_utc=utc_now(),
                                             environment_overrides=BLAS_ENV))
    started = time.monotonic()
    process = None
    exit_code = None
    usage = None
    termination = None
    failure = None
    try:
        with (generation / 'stdout.log').open('xb') as stdout, \
                (generation / 'stderr.log').open('xb') as stderr:
            try:
                if stop.is_set() or time.monotonic() >= deadline:
                    termination = 'batch stopped before spawn'
                else:
                    process = subprocess.Popen(argv, stdout=stdout, stderr=stderr,
                                               env=dict(os.environ, **BLAS_ENV),
                                               start_new_session=True)
                    save_json(generation / 'spawn.json', dict(pid=process.pid,
                                                             spawned_utc=utc_now()))
                    while True:
                        pid, status, current_usage = os.wait4(process.pid, os.WNOHANG)
                        if pid:
                            exit_code = os.waitstatus_to_exitcode(status)
                            process.returncode = exit_code
                            usage = current_usage
                            break
                        now = time.monotonic()
                        if termination is None:
                            if stop.is_set():
                                termination = 'batch interrupted'
                            elif now >= deadline:
                                termination = 'batch wall grant exhausted'
                            elif now - started >= timeout:
                                termination = 'case wall grant exhausted'
                            if termination is not None:
                                try:
                                    os.killpg(process.pid, signal.SIGKILL)
                                except ProcessLookupError:
                                    pass
                        time.sleep(0.05)
            except BaseException as error:
                failure = repr(error)
                if process is not None and process.returncode is None:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    _, status, usage = os.wait4(process.pid, 0)
                    exit_code = os.waitstatus_to_exitcode(status)
                    process.returncode = exit_code
            finally:
                stdout.flush()
                stderr.flush()
                os.fsync(stdout.fileno())
                os.fsync(stderr.fileno())
    except BaseException as error:
        failure = repr(error)
    outcome = dict(schema='P254-A34-mode-noise-subprocess-outcome-v1',
                   original_record=record, argv=argv, cpu=cpu,
                   started_utc=intent['started_utc'], finished_utc=utc_now(),
                   exit_code=exit_code, termination=termination, error=failure,
                   elapsed_seconds=time.monotonic() - started,
                   user_cpu_seconds=usage.ru_utime if usage is not None else None,
                   system_cpu_seconds=usage.ru_stime if usage is not None else None,
                   peak_rss_bytes=usage.ru_maxrss * 1024 if usage is not None else None,
                   stdout=binding(generation / 'stdout.log')
                          if (generation / 'stdout.log').exists() else None,
                   stderr=binding(generation / 'stderr.log')
                          if (generation / 'stderr.log').exists() else None)
    # Consumer sees only this durable receipt, never an in-memory process outcome.
    save_json(generation / 'outcome.json', outcome)
    return generation


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
            intent = load_json(generation / 'intent.json')
            # Implementation byte bindings are provenance, not scientific compatibility.
            # load_case admits the science identity or an exact authorized legacy result.
            for key in ('inputs', 'contract', 'source0_action_packet',
                        'source0_mode1_packet'):
                if intent['binding_contract'][key] != current[key]:
                    raise ValueError('generation fixed input/plan/paid-packet contract differs')
            if intent['original_record'] != record:
                raise ValueError('generation original record differs')
            outcome = load_json(generation / 'outcome.json')
            if (outcome['original_record'] != record
                    or type(outcome['exit_code']) is not int or outcome['exit_code'] != 0
                    or outcome['termination'] is not None or outcome['error'] is not None
                    or outcome['argv'] != intent['argv']):
                raise ValueError('generation has no successful durable subprocess exit')
            authenticate(outcome['stdout'], generation / 'stdout.log')
            authenticate(outcome['stderr'], generation / 'stderr.log')
            case = load_case(done, context, record)
            compatible.append(dict(original_record=record, done=case['done'],
                                   result=case['result'], generation=str(generation),
                                   reused=True, work=case['payload']['work']))
        except Exception as error:
            rejected.append(dict(original_record=record, generation=str(generation),
                                 error=repr(error)))
    if len(compatible) > 1:
        raise ValueError('multiple compatible committed generations for record %d; '
                         'select an unambiguous batch directory' % record)
    return compatible[0] if compatible else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    add_input_arguments(parser)
    parser.add_argument('--reuse-source0-action-packet', type=Path, required=True)
    parser.add_argument('--reuse-source0-mode1-packet', type=Path, required=True)
    parser.add_argument('--producer', type=Path,
                        default=Path(__file__).with_name('export_full_source_modes.py'))
    parser.add_argument('--python', type=Path, default=Path(sys.executable))
    parser.add_argument('--output', type=Path, required=True,
                        help='batch root, preserved across invocations for sector-level resume')
    parser.add_argument('--records', type=lambda s: parse_indices(s, 128),
                        default=list(range(128)), help='original records, e.g. 0-127 or 0,7,64')
    parser.add_argument('--cpus', type=parse_indices, required=True,
                        help='explicit allowed logical CPU list; one supplied CPU per active case')
    parser.add_argument('--workers', type=int, default=1,
                        help='bounded concurrent independent cases, default 1, at most supplied CPUs')
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
    drain_errors = []
    launched = set()
    error = None
    interrupted = []
    events = (run / 'progress.jsonl').open('xb')
    lock = (root / '.batch.lock').open('a+b')
    old_handlers = {}

    def cost():
        def total(field):
            values = [x[field] for x in outcomes]
            return None if any(x is None for x in values) else sum(values)

        rss = [x['peak_rss_bytes'] for x in outcomes]
        return dict(subprocess_elapsed_seconds=total('elapsed_seconds'),
                    user_cpu_seconds=total('user_cpu_seconds'),
                    system_cpu_seconds=total('system_cpu_seconds'),
                    peak_case_rss_bytes=None if any(x is None for x in rss)
                                        else max(rss, default=0),
                    measured_subprocesses=sum(not x['missing_measurements'] for x in outcomes),
                    unmeasured_subprocesses=sum(bool(x['missing_measurements']) for x in outcomes),
                    new_successful_physical_force_forward_actions=sum(
                        x['work']['physical_force_forward_actions'] for x in completed if not x['reused']),
                    new_successful_physical_force_adjoint_actions=sum(
                        x['work']['physical_force_adjoint_actions'] for x in completed if not x['reused']),
                    failed_generation_physical_action_count='unknown' if failed else 0)

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

    def consume(future, record, generation):
        """Reap one owned worker and account it before fallible publication."""
        worker_error = None
        try:
            future.result()
        except BaseException as exception:
            worker_error = exception

        def reject(exception):
            item = dict(original_record=record, generation=str(generation),
                        error=repr(exception))
            # Keep accounting even if either subsequent publication fails.
            failed.append(item)
            save_json(generation / 'rejected.json', item)
            emit('case-failed', **item)

        measurement_fields = ('elapsed_seconds', 'user_cpu_seconds',
                              'system_cpu_seconds', 'peak_rss_bytes')
        try:
            outcome = load_json(generation / 'outcome.json')
            if outcome['original_record'] != record:
                raise ValueError('durable outcome original record differs')
        except BaseException as exception:
            # A submitted future must never disappear from cost accounting.
            outcomes.append(dict(original_record=record, exit_code=None,
                                 elapsed_seconds=None, user_cpu_seconds=None,
                                 system_cpu_seconds=None, peak_rss_bytes=None,
                                 missing_measurements=list(measurement_fields),
                                 measurement_error=repr(exception)))
            reject(exception)
            raise
        for field in measurement_fields:
            value = outcome.get(field)
            if (value is None or isinstance(value, bool)
                    or not isinstance(value, (int, float))
                    or not math.isfinite(value) or value < 0):
                outcome[field] = None
        outcome['missing_measurements'] = [
            field for field in measurement_fields if outcome[field] is None]
        outcomes.append(outcome)
        try:
            successful = (type(outcome['exit_code']) is int and outcome['exit_code'] == 0
                          and outcome['termination'] is None and outcome['error'] is None)
        except BaseException as exception:
            outcome['exit_code'] = None
            reject(exception)
            raise
        if not successful:
            reject(ValueError('producer did not exit successfully: ' + repr(outcome)))
            if worker_error is not None:
                raise worker_error
            return
        try:
            intent = load_json(generation / 'intent.json')
            if outcome['argv'] != intent['argv']:
                raise ValueError('durable outcome argv differs from launch intent')
            authenticate(outcome['stdout'], generation / 'stdout.log')
            authenticate(outcome['stderr'], generation / 'stderr.log')
            # Only a successful durable exit/stream receipt permits case admission.
            case = load_case(generation / 'case' / 'done.json', context, record)
            item = dict(original_record=record, done=case['done'],
                        result=case['result'], generation=str(generation),
                        reused=False, work=case['payload']['work'])
        except BaseException as exception:
            reject(exception)
            raise
        # Publication errors must not erase a scientifically admitted successful case.
        completed.append(item)
        save_json(generation / 'accepted.json', item)
        emit('case-complete', original_record=record, done=case['done'])
        if worker_error is not None:
            raise worker_error

    def on_signal(signum, _frame):
        interrupted.append(signum)
        stop.set()

    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if not 1 <= args.workers <= min(len(args.cpus), 128):
            raise ValueError('workers must be between1 and min(number of supplied CPUs,128)')
        if not set(args.cpus).issubset(os.sched_getaffinity(0)):
            raise ValueError('supplied CPUs are outside the actual inherited CPU affinity')
        taskset = shutil.which('taskset')
        if taskset is None:
            raise ValueError('Linux taskset is required for explicit child CPU affinity')
        context = load_context(args.initial_result, args.coefficient_manifest,
                               args.channel_archive, args.contract)
        mode0 = binding(args.reuse_source0_action_packet)
        mode1 = binding(args.reuse_source0_mode1_packet)
        if (not same_binding(mode0, context['contract']['mode0_packet'])
                or not same_binding(mode1, context['contract']['mode1_packet'])):
            raise ValueError('source0 paid reuse packets differ from A34 contract')
        current = dict(inputs=context['inputs'], contract=context['contract_binding'],
                       source0_action_packet=mode0, source0_mode1_packet=mode1,
                       scientific_compatibility_identity=context['scientific_compatibility_identity'],
                       producer=binding(args.producer), python=binding(args.python),
                       python_command=str(args.python.absolute()),
                       runner=binding(__file__),
                       aggregator=binding(Path(__file__).with_name('aggregate_mode_noise.py')))
        save_json(run / 'intent.json', dict(schema='P254-A34-mode-noise-batch-intent-v1',
                                           binding_contract=current, records=args.records,
                                           cpus=args.cpus, workers=args.workers,
                                           wall_seconds=args.wall_seconds,
                                           case_timeout_seconds=args.case_timeout_seconds,
                                           environment_overrides=BLAS_ENV,
                                           argv=sys.argv, started_utc=utc_now()))
        for sig in (signal.SIGINT, signal.SIGTERM):
            old_handlers[sig] = signal.signal(sig, on_signal)
        emit('batch-start', requested_original_records=args.records, run=str(run))
        pending = []
        for record in args.records:
            case = reusable(root, record, current, context, rejected)
            if case is None:
                pending.append(record)
            else:
                completed.append(case)
                emit('case-reused', original_record=record, done=case['done'])
        active = {}
        available = list(args.cpus[:args.workers])
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            try:
                while pending or active:
                    while (pending and available and not stop.is_set()
                           and time.monotonic() < deadline):
                        record = pending.pop(0)
                        cpu = available.pop(0)
                        generation = root / 'records' / ('record-%03d' % record) / (
                            'generation-' + uuid.uuid4().hex)
                        argv = [taskset, '--cpu-list', str(cpu), current['python_command'],
                                current['producer']['path'], '--initial-result',
                                current['inputs']['initial_result']['path'], '--coefficient-manifest',
                                current['inputs']['coefficient_manifest']['path'], '--channel-archive',
                                current['inputs']['channel_archive']['path'], '--original-record',
                                str(record), '--output', str(generation / 'case'),
                                '--reuse-source0-action-packet', mode0['path'],
                                '--reuse-source0-mode1-packet', mode1['path']]
                        intent = dict(binding_contract=current, original_record=record,
                                      started_utc=utc_now())
                        future = executor.submit(execute_case, record, cpu, generation, argv,
                                                 intent, stop, deadline, args.case_timeout_seconds)
                        active[future] = (record, cpu, generation)
                        launched.add(record)
                        emit('case-start', original_record=record, cpu=cpu,
                             generation=str(generation), argv=argv)
                    if not active:
                        break
                    settled, _ = wait(active, timeout=0.25, return_when=FIRST_COMPLETED)
                    for future in settled:
                        record, cpu, generation = active.pop(future)
                        available.append(cpu)
                        consume(future, record, generation)
            except BaseException:
                # This must precede ThreadPoolExecutor.__exit__'s blocking shutdown.
                stop.set()
                for future, (record, _cpu, generation) in list(active.items()):
                    try:
                        consume(future, record, generation)
                    except BaseException as exception:
                        # One broken receipt/progress write cannot strand other workers.
                        drain_errors.append(dict(original_record=record,
                                                 generation=str(generation),
                                                 error=repr(exception)))
                    finally:
                        del active[future]
                raise
        if interrupted:
            error = 'interrupted by signal(s): ' + repr(interrupted)
    except BaseException as exception:
        stop.set()
        error = repr(exception)
    finally:
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
        done_records = {item['original_record'] for item in completed}
        unstarted = sorted(set(args.records) - done_records - launched)
        complete = done_records == set(args.records) and not failed and error is None
        code = 0 if complete else (130 if interrupted else 1)
        summary = dict(schema='P254-A34-mode-noise-batch-summary-v1',
                       status='requested-subset-complete' if complete else 'incomplete',
                       full_original128_coverage=complete and done_records == set(range(128)),
                       requested_original_records=args.records,
                       launched_original_records=sorted(launched),
                       incomplete_original_records=sorted(set(args.records) - done_records),
                       completed=sorted(completed, key=lambda x: x['original_record']),
                       failed=failed, unstarted_original_records=unstarted,
                       rejected_reuse_candidates=rejected, error=error,
                       drain_errors=drain_errors,
                       elapsed_seconds=time.monotonic() - started, cost=cost(),
                       subprocess_outcomes=[dict(original_record=x['original_record'],
                                                 exit_code=x['exit_code'],
                                                 elapsed_seconds=x['elapsed_seconds'],
                                                 user_cpu_seconds=x['user_cpu_seconds'],
                                                 system_cpu_seconds=x['system_cpu_seconds'],
                                                 peak_rss_bytes=x['peak_rss_bytes'],
                                                 missing_measurements=x['missing_measurements'],
                                                 measurement_error=x.get('measurement_error'))
                                            for x in outcomes],
                       exit_code=code, scope=SCOPE)
        try:
            receipt = save_json(run / 'summary.json', summary)
            emit('batch-exit', status=summary['status'], exit_code=code, summary=receipt,
                 unstarted_original_records=unstarted, error=error)
        finally:
            events.close()
            lock.close()
    return code


if __name__ == '__main__':
    raise SystemExit(main())

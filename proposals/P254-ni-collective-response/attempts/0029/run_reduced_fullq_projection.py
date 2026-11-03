"""Produce explicit full-q shards from a SETTLED complete128-sector Gamma bank.

Examples (replace paths with the actual case's authenticated original inputs):

  python run_reduced_fullq_projection.py export --case loaded-ni31h6 \
    --boundary /absolute/original/boundary.json --cache /var/tmp/fullq-loaded-q01 \
    --gamma /var/tmp/p254-activeband-loaded-contrast-projection-001/manifest.json \
    --producer-history /absolute/frozen-producer-snapshots \
    --q 1 --sources all --kernels --seconds 7200 --maxrss-mib 3072

  python run_reduced_fullq_projection.py export --case loaded-ni31h6 \
    --boundary /absolute/original/boundary.json --cache /var/tmp/fullq-loaded-q02a \
    --gamma /absolute/settled-Gamma/manifest.json \
    --producer-history /absolute/frozen-producer-snapshots \
    --q 2 --sources 0-63 --kernels --seconds 7200 --maxrss-mib 3072

  python run_reduced_fullq_projection.py union --manifest /absolute/shard-a/manifest.json \
    --manifest /absolute/shard-b/manifest.json --output /var/tmp/fullq-union/manifest.json

The source-split companion for q2 must use sources64-127 WITHOUT --kernels;
otherwise union refuses duplicate kernel ownership. Enumerate every required
finite-q edge and exactly one owner per kernel across the real shard inputs.
The union example rejects incomplete coverage until ALL64q/all128-slot products
are present; Gamma is shared, not fabricated. Resume by repeating the SAME
explicit command/assignment and --gamma input. The destination's explicit
pilot seed is always that settled Gamma bank, whether originally fresh or
pilot-derived. No automatic retries, workers, launchers or watchers are created.
One writer owns each cache; repeating a command does not reset its deadline.

Scientific scope: approved reduced fixed12-band prepared electronic model with
nested8/12 and7/27 readout, original64q/128 sectors, actual populations, original
full paired ionic support and independently evaluated complex reverse maps.
This entrypoint changes storage/ownership only; MaterialExporter supplies all
physics, and missing admitted Native/Gram/H0 preparations are hard failures.
No witness mode is exposed because it would rebuild admitted full factors.
"""
from time import perf_counter
COMMAND_STARTED = perf_counter()

import argparse
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import json
import math
from pathlib import Path
import resource
import signal
import sys
from time import time
import traceback
from uuid import uuid4

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0028'), str(HERE)]
from reduced_projection_cache import atomic_json
from reduced_material_projection import MaterialExporter, selected
from reduced_fullq_cache import adopt_complete_gamma, complete_keys, FullQCache, union_manifests
from reduced_preparation_storage import MappedPreparations


class CommandInterrupted(BaseException):
    def __init__(self, signum):
        self.signum = signum
        super().__init__('projection command interrupted by signal %d' % signum)


class PreparedExporter(MaterialExporter):
    def check_budget(self):
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
        if time() >= self.command_deadline_unix:
            raise TimeoutError('original command-start wall deadline reached')
        if rss > self.args.maxrss_mib:
            raise MemoryError('projection command RSS budget reached')

    def _admit_basis(self, index):
        raise ValueError('missing admitted original preparation; raw basis/Gram rebuild forbidden: %03d' % index)

    def _export_slot(self, index):
        raise ValueError('missing admitted original preparation; Native/H0 rebuild forbidden: %03d' % index)


def _utc():
    return datetime.now(timezone.utc).isoformat()


@contextmanager
def _writer_lock(directory):
    directory = Path(directory).resolve()
    if directory == ROOT or ROOT in directory.parents:
        raise ValueError('numerical outputs must remain outside Git')
    directory.mkdir(parents=True, exist_ok=True)
    # flock is released by every process exit, including SIGKILL. The persistent
    # empty lock file is not an item receipt and never indicates completion.
    with (directory/'.fullq-writer.lock').open('a') as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError('another writer owns this full-q cache') from exc
        yield


@contextmanager
def _deadline(seconds, started):
    if not hasattr(signal, 'setitimer'):
        raise RuntimeError('POSIX original-start wall deadline required')
    if signal.getitimer(signal.ITIMER_REAL) != (0.0, 0.0):
        raise RuntimeError('cannot replace an existing process alarm')
    previous = {number: signal.getsignal(number)
                for number in (signal.SIGALRM, signal.SIGINT, signal.SIGTERM)}

    def interrupt(signum, _frame):
        if signum == signal.SIGALRM:
            raise TimeoutError('original command-start wall deadline reached')
        raise CommandInterrupted(signum)

    try:
        for number in previous:
            signal.signal(number, interrupt)
        remaining = seconds-(perf_counter()-started)
        if remaining <= 0:
            raise TimeoutError('original command-start wall deadline already reached')
        signal.setitimer(signal.ITIMER_REAL, remaining)
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        for number, handler in previous.items():
            signal.signal(number, handler)


def _release_q(exporter):
    exporter.last_bundle = None
    exporter.last_index = None
    exporter.local_cache.clear()
    exporter.prepared_maps.recent.clear()


def _save_command(path, run, exporter=None):
    run['elapsed_seconds'] = perf_counter()-run['_started']
    run['maxrss_mib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
    payload = _public_run(run)
    atomic_json(path, payload)
    if exporter is not None and hasattr(exporter, 'cache'):
        manifest = exporter.cache.manifest
        commands = manifest.setdefault('fullq_commands', [])
        for index, command in enumerate(commands):
            if command['id'] == run['id']:
                commands[index] = payload
                break
        else:
            commands.append(payload)
        exporter.cache.publish()


def _public_run(run):
    return {name: value for name, value in run.items() if not name.startswith('_')}


def _original_allocation(directory, args, started):
    """Persist one original cache/union allocation, including failed attempts."""
    path = directory/'original-allocation.json'
    identity = vars(args).copy()
    if path.exists():
        allocation = json.loads(path.read_bytes())
        if allocation['arguments'] != identity:
            raise ValueError('original allocation arguments/budget differ on resume')
    else:
        original = time()-(perf_counter()-started)
        allocation = dict(arguments=identity, original_started_unix=original,
                          wall_deadline_unix=original+args.seconds)
        atomic_json(path, allocation)
    return allocation


def _execute(args, started):
    destination = Path(args.cache if args.command == 'export' else args.output).resolve()
    directory = destination if args.command == 'export' else destination.parent
    exporter = None
    run = dict(id=uuid4().hex, command=args.command, status='running',
               original_started_utc=datetime.fromtimestamp(
                   time()-(perf_counter()-started), timezone.utc).isoformat(),
               wall_deadline_unix=time()+args.seconds-(perf_counter()-started),
               seconds=args.seconds, maxrss_limit_mib=args.maxrss_mib,
               arguments=vars(args).copy(), progress=dict(stage='initializing'), _started=started)
    command_path = directory/('command-'+run['id']+'.json')
    code = 1
    with _writer_lock(directory):
        allocation = _original_allocation(directory, args, started)
        run['wall_deadline_unix'] = allocation['wall_deadline_unix']
        run['original_started_utc'] = datetime.fromtimestamp(
            allocation['original_started_unix'], timezone.utc).isoformat()
        effective_seconds = perf_counter()-started+run['wall_deadline_unix']-time()
        _save_command(command_path, run)
        try:
            with _deadline(effective_seconds, started):
                if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024 > args.maxrss_mib:
                    raise MemoryError('projection command RSS budget reached')
                if args.command == 'union':
                    run['progress'] = dict(stage='authenticating-and-uniting', inputs=len(args.manifest))
                    _save_command(command_path, run)
                    result = union_manifests(args.manifest, args.output)
                    if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024 > args.maxrss_mib:
                        raise MemoryError('union command RSS budget reached')
                    run['progress'] = dict(stage='complete', manifest=str(result))
                else:
                    qs = selected(args.q, 64)
                    sources = selected(args.sources, 128)
                    # Retain a partial object if constructor fails after it has
                    # opened its cache: record that exit rather than losing it.
                    exporter = PreparedExporter.__new__(PreparedExporter)
                    exporter.command_deadline_unix = run['wall_deadline_unix']
                    exporter.__init__(args)
                    # Exclusive writer ownership proves that an older running
                    # receipt is abandoned, not still being produced. Preserve
                    # its unknown exit honestly; never claim a caught failure.
                    for previous in exporter.cache.manifest.get('fullq_commands', []):
                        if previous.get('status') == 'running':
                            previous['status'] = 'abandoned'
                            previous['exit_code'] = None
                            previous['observed_abandoned_by_command'] = run['id']
                    exporter.started = started
                    run['progress'] = dict(stage='adopting-complete-Gamma')
                    _save_command(command_path, run, exporter)
                    admitted = adopt_complete_gamma(exporter, args.gamma, args.producer_history)
                    exporter.cache = FullQCache(exporter.cache, qs, sources, args.kernels)
                    exporter.prepared_maps = MappedPreparations(exporter.cache,
                        receipt_bindings=admitted['prepared_map_receipts'])
                    exporter.cache.manifest['fullq_worker_status'] = 'running'
                    _save_command(command_path, run, exporter)
                    print(json.dumps(dict(event='fullq-start', case=args.case,
                        q=qs, assigned_sources=len(sources),
                        original_wall_deadline_unix=run['wall_deadline_unix'],
                        command_receipt=str(command_path))), flush=True)
                    for q in qs:
                        _release_q(exporter)
                        exporter.check_budget()
                        for offset, source in enumerate(sources):
                            run['progress'] = dict(stage='edge', q=q, source=source)
                            # The edge itself checkpoints the complete item map.
                            # Do not rewrite that growing map twice per edge.
                            _save_command(command_path, run)
                            exporter.edge(q, source, witness=False)
                            exporter.check_budget()
                            if (offset+1) % 8 == 0 or offset+1 == len(sources):
                                print(json.dumps(dict(event='fullq-sources-complete',
                                    case=args.case, q=q, source_count=offset+1,
                                    assigned_sources=len(sources))), flush=True)
                        if args.kernels:
                            # Release prepared columns before potentially large
                            # physical kernel work; keep no previous-q gradient.
                            _release_q(exporter)
                            run['progress'] = dict(stage='kernel', q=q)
                            _save_command(command_path, run, exporter)
                            exporter.kernel(q)
                            exporter.check_budget()
                        _release_q(exporter)
                        run['progress'] = dict(stage='q-complete', q=q)
                        _save_command(command_path, run, exporter)
                    # Both passive splits/chunks are already admitted. This
                    # current-interface finish therefore cannot reconstruct them.
                    exporter.finish()
                    exporter.check_budget()
                    exporter.cache.manifest['status'] = (
                        'complete' if complete_keys().issubset(exporter.cache.manifest['items'])
                        else 'gamma-complete')
                    exporter.cache.manifest['fullq_worker_status'] = 'settled'
                    run['progress'] = dict(stage='complete', manifest=str(exporter.cache.path),
                                           completed_items=len(exporter.cache.manifest['items']))
                run['status'] = 'settled'
                code = 0
        except BaseException as exc:
            run['status'] = ('budget-exit' if isinstance(exc, (TimeoutError, MemoryError))
                             else 'interrupted' if isinstance(exc, (CommandInterrupted, KeyboardInterrupt))
                             else 'failed')
            run['failure'] = dict(kind=type(exc).__name__, message=str(exc),
                                  traceback=traceback.format_exc())
            code = (3 if run['status'] == 'budget-exit' else
                    128+exc.signum if isinstance(exc, CommandInterrupted) else
                    130 if isinstance(exc, KeyboardInterrupt) else 1)
            if exporter is not None and hasattr(exporter, 'cache'):
                exporter.cache.manifest['fullq_worker_status'] = run['status']
                exporter.cache.manifest.setdefault('failures', []).append(dict(
                    command_id=run['id'], **run['failure'],
                    elapsed_seconds=perf_counter()-started))
            print(str(exc), file=sys.stderr)
        finally:
            if exporter is not None and hasattr(exporter, 'local_cache'):
                _release_q(exporter)
            run['finished_utc'] = _utc()
            run['exit_code'] = code
            # This is outside the alarm context: preserve the failure receipt
            # rather than repeatedly interrupting atomic publication at deadline.
            _save_command(command_path, run, exporter)
    print(json.dumps(_public_run(run), sort_keys=True))
    return code


def main(argv=None):
    started = COMMAND_STARTED
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)
    export = sub.add_parser('export', help='one explicit prepared-only q/source shard')
    export.add_argument('--case', required=True)
    export.add_argument('--boundary', required=True)
    export.add_argument('--cache', required=True)
    export.add_argument('--gamma', required=True, help='settled complete128-sector Gamma manifest')
    export.add_argument('--producer-history', required=True)
    export.add_argument('--q', required=True, help='explicit all or comma/ranges, q indices0..63')
    export.add_argument('--sources', required=True, help='explicit all or comma/ranges, slots0..127')
    export.add_argument('--kernels', action='store_true', help='own selected finite-q kernels too')
    export.set_defaults(bands=12, features=27)
    union = sub.add_parser('union', help='refuse anything short of a complete settled full-q union')
    union.add_argument('--manifest', action='append', required=True, help='repeat for each independent shard')
    union.add_argument('--output', required=True, help='new manifest.json outside Git')
    for command in (export, union):
        command.add_argument('--seconds', type=float, default=7200.)
        command.add_argument('--maxrss-mib', type=float, default=3072.)
    export.set_defaults(reuse_gamma_pilot=None)
    args = parser.parse_args(argv)
    if (not math.isfinite(args.seconds) or args.seconds <= 0 or
            not math.isfinite(args.maxrss_mib) or args.maxrss_mib <= 0):
        parser.error('finite positive wall/RSS budgets required')
    if args.command == 'export':
        args.reuse_gamma_pilot = args.gamma
    return _execute(args, started)


if __name__ == '__main__':
    raise SystemExit(main())

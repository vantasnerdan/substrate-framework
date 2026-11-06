"""Exercise F1 at the actual summary rotation boundary; never recompute gPsi.

Reuse the isolated fork/rename-boundary pattern from attempt0021's
run_state_recovery_smoke.py:218-265. Original material files are read-only;
only private summary copies are written. No interrupted-operation, Native
history, adversarial path-replacement or power-loss recovery is claimed.
"""
import time
BOOT = time.monotonic()
import argparse
import json
import os
from pathlib import Path
import select
import signal
import sys
import tempfile
import traceback
sys.dont_write_bytecode = True
import material_pilot_io as checkpoints
from substrate_framework.qe_saved_native import read_bound_bytes


def recover_all(original, target):
    records = [original['coordinates'], *original['blocks'], *original['saved_band_matrices']]
    observed = []
    for expected in records:
        key = Path(expected['archive']['path']).stem
        actual = checkpoints.recover_block(target, key, original['recipe_sha256'], expected['metadata'])
        if actual != expected:
            raise ValueError('completed archive/receipt differs: '+key)
        observed.append({'archive': actual['archive'],
                         'receipt': checkpoints.binding(target/(key+'.json'))})
    return observed


def crash_summary(scratch, original):
    summary = scratch/'summary.json'
    checkpoints.save_json(summary, original)
    old = checkpoints.binding(summary)
    pipe_read, pipe_write = os.pipe()
    child = os.fork()
    if child == 0:
        os.close(pipe_read)
        signal.signal(signal.SIGALRM, lambda *_: os._exit(124))
        signal.setitimer(signal.ITIMER_REAL, 10)
        real_replace = checkpoints.os.replace

        def pause_after_rotation(source, destination):
            real_replace(source, destination)
            if Path(source) == summary and Path(destination) == summary.with_name('summary.json.prev'):
                os.write(pipe_write, b'rotated')
                while True:
                    signal.pause()

        checkpoints.os.replace = pause_after_rotation
        try:
            checkpoints.save_json(summary, dict(original, isolated_summary_probe=True), rotate=True)
        finally:
            os._exit(2)
    os.close(pipe_write)
    killed_pid = child
    try:
        ready, _, _ = select.select([pipe_read], [], [], 5)
        if not ready or os.read(pipe_read, 64) != b'rotated':
            raise RuntimeError('actual writer did not reach the rotation boundary')
        os.kill(child, signal.SIGKILL)
        _, status = os.waitpid(child, 0)
        child = None
        if not os.WIFSIGNALED(status) or os.WTERMSIG(status) != signal.SIGKILL:
            raise RuntimeError('owned summary writer did not die by SIGKILL')
        previous = summary.with_name('summary.json.prev')
        temporary = list(scratch.glob('summary.json.tmp-*'))
        if summary.exists() or not previous.exists() or len(temporary) != 1:
            raise RuntimeError('unexpected actual summary crash boundary')
        previous_bytes = previous.read_bytes()
        if checkpoints.binding(previous, previous_bytes)['sha256'] != old['sha256']:
            raise ValueError('previous committed summary changed')
        recovered_summary = json.loads(previous_bytes)
        if recovered_summary != original:
            raise ValueError('previous summary cannot recover the original recipe/receipts')
        return recovered_summary, {'child_pid': killed_pid, 'wait_status': status,
            'signal': 'SIGKILL', 'current_summary_absent': True,
            'previous_summary': checkpoints.binding(previous, previous_bytes),
            'interrupted_temporary': checkpoints.binding(temporary[0]),
            'boundary': 'after real summary-to-prev rename; before temporary-to-summary rename',
            'instrumentation': 'Only the child rename boundary pauses; the real writer and real rename run.'}
    finally:
        if child is not None:
            try:
                os.kill(child, signal.SIGKILL)
            except ProcessLookupError:
                pass
            os.waitpid(child, 0)
        os.close(pipe_read)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--admission', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    report = {'passed': False, 'scope': __doc__}
    try:
        def timed_out(*_):
            raise TimeoutError('isolated crash-check allocation exhausted')
        signal.signal(signal.SIGALRM, timed_out)
        signal.setitimer(signal.ITIMER_REAL, max(0.01, 27-(time.monotonic()-BOOT)))
        raw = args.admission.read_bytes()
        grant = json.loads(raw)
        report['admission'] = checkpoints.binding(args.admission, raw)
        for source in grant['sources']:
            read_bound_bytes(source)
        terminal_raw, terminal_receipt = read_bound_bytes(grant['original_terminal'])
        original = json.loads(terminal_raw)
        target = Path(grant['original_terminal']['path']).parent
        if args.output.resolve().is_relative_to(target.resolve()):
            raise ValueError('crash receipt must be outside original material outputs')
        report['implementation_sources'] = checkpoints.implementation()
        report['original_terminal'] = dict(terminal_receipt)
        before = recover_all(original, target)
        report['original_pairs_before'] = before
        with tempfile.TemporaryDirectory(prefix='material-summary-crash-', dir=args.output.parent) as directory:
            recovered, report['crash'] = crash_summary(Path(directory), original)
            after = recover_all(recovered, target)
            if after != before:
                raise ValueError('original archive/receipt identities changed after the crash')
            report['original_pairs_after'] = after
        report['scratch_removed'] = True
        report['exact_recovered_pairs'] = len(after)
        report['original_terminal_unchanged'] = checkpoints.binding(Path(terminal_receipt['path'])) == dict(terminal_receipt)
        if not report['original_terminal_unchanged']:
            raise ValueError('original terminal bytes changed')
        report['passed'] = True
    except BaseException as error:
        report['failure'] = {'type': type(error).__name__, 'message': str(error),
                             'traceback': traceback.format_exc()}
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        report['elapsed_seconds'] = time.monotonic()-BOOT
        checkpoints.save_json(args.output, report)
        print(json.dumps({k: report.get(k) for k in
            ('passed', 'failure', 'crash', 'exact_recovered_pairs', 'scratch_removed',
             'original_terminal_unchanged', 'elapsed_seconds')}), flush=True)
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())

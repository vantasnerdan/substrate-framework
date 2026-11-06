#!/usr/bin/env python3
"""Root-only completed-cache/IO smoke; never an unfinished-trajectory restart.

Run with the original pilot's Python/PYTHONPATH. --receipt-dir must not exist.
Only this smoke's private scratch is mutable; no original CLI --resume, input
loading, model reconstruction, RHS, Native, or scientific production is invoked.
The retained source copy and receipt are outside scratch, which is removed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import select
import shutil
import signal
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ORIGINAL = HERE / 'pilot-001'
FINGERPRINT = '49ef49241b4aef1fd45633bf6b6fb4da1bb6726dd7c65d2957f6b8ca81e955f5'
SOURCE_HASHES = {
    HERE / 'time_pilot_checkpoints.py': '999ad8ba9242d77346c456cb6e86ca1a6aa9869a541b48f230e21d0ed563af7d',
    HERE / 'run_local_kick.py': '9cc0d74ee70b7308ea6e0337c12ab27a048b239796aec98a11370c516f586807',
    HERE.parent / '0020/pilot_inputs.py': '68b4a221e7073e422b7eab539c30f966c3b6c851bdc22c08c3a78029bb1153fc',
    HERE.parent / '0020/pilot_checkpoints.py': '69ac9bd90cd2b0267cd3674c8c2b54874609159a1dc60703b1a97c155dc75f7e',
}
FINDINGS = [
    {
        'id': 'SR-01', 'kind': 'bounded_recovery_scope',
        'finding': 'Committed receipts index full states independently of result.json and rotating state pointers. This smoke exercises completed-cache recovery and exact saved-array reconstruction only, not numerical continuation.',
        'evidence': ['time_pilot_checkpoints.py:28-72', 'run_local_kick.py:116-120'],
        'impact_on_completed_result': 'None identified by this storage assessment.',
    },
    {
        'id': 'SR-02', 'kind': 'operational_bug',
        'finding': 'save_json rotates result.json to result.json.prev before publishing the temporary summary. SIGKILL in that gap leaves committed samples recoverable but CLI resume unconditionally reads missing result.json before admission and fails.',
        'evidence': ['../0020/pilot_checkpoints.py:13-18', 'run_local_kick.py:81-91'],
        'minimal_fix_for_root': 'At the resume reference read, on FileNotFoundError only, read result.json.prev instead. Keep recover_samples authoritative. If both summaries are absent, explicitly decide how advisory reference metadata is recovered; do not discard committed states or silently reset the generation. No original source is changed here.',
        'impact_on_completed_result': 'No evidence this crash occurred; result.json exists and reports complete.',
    },
    {
        'id': 'SR-03', 'kind': 'unproven_numerical_continuation',
        'finding': 'Electronic state coordinates are target/source eigenbasis amplitudes. Unfinished resume rebuilds eigh frames and checks only saved vector shape. Eigenvectors/gauge and solver-library identity are not persisted with state; identical source/input fingerprints alone do not establish identical coordinates after rebuilding (phase/degenerate-subspace freedom). No actual gauge mismatch was observed or produced.',
        'evidence': ['run_local_kick.py:109-137', '../../../../src/substrate_framework/qe_pw_response_frames.py:180-194,265-270,293-295,306-337', '../../../../src/substrate_framework/pw_electronic_dynamics.py:10-13,63-79,117-125', 'time_pilot_checkpoints.py:35-43'],
        'required_before_claim': 'Demonstrate coordinate-compatible reconstructed frames or preserve/transform the actual basis. Any future basis identity check is justified by preventing reinterpretation of saved amplitudes, not by adding an arbitrary frame signature. Not exercised by this smoke.',
        'impact_on_completed_result': 'No restart was needed for the banked completed generation; this is a future-continuation limitation, not a rejection of its saved first moments.',
    },
    {
        'id': 'SR-04', 'kind': 'declared_loss_boundary',
        'finding': 'Adaptive DOP853 history and in-flight stage work are not saved. Only completed quarter-tau states commit. Work since the latest receipt, potentially the entire current interval, can be lost; restarted integration need not reproduce the interrupted numerical path bitwise.',
        'evidence': ['run_local_kick.py:61,132-166', 'time_pilot_checkpoints.py:28-46'],
        'impact_on_completed_result': 'All ten declared output states are already committed; no missing interval is inferred.',
    },
    {
        'id': 'SR-05', 'kind': 'operational_bug',
        'finding': 'RHS action recovery parses every JSONL line unconditionally. A torn final append can raise JSONDecodeError and block unfinished continuation despite valid committed states. The counter is reserved before derivative evaluation, so complete records may conservatively include an action whose derivative never completed.',
        'evidence': ['run_local_kick.py:128-131,145-154'],
        'minimal_fix_for_root': 'Recognize and report a malformed final ledger record distinctly, retain the raw ledger, and require an explicit accounting decision before continuation. A bounded repair must conservatively charge a possibly reserved action rather than silently ignoring it or resetting the cap; malformed interior records must still refuse. Do not add automatic retries.',
        'impact_on_completed_result': 'This is a crash-path limitation, not evidence of a corrupt original ledger.',
    },
    {
        'id': 'SR-06', 'kind': 'intentional_fixed_allocation',
        'finding': 'Wall allocations are absolute time.time timestamps, not active CPU time: downtime consumes the original total and sector windows. Expired unfinished resume refuses before frame construction; completed-cache reuse skips those budget checks but still loads original inputs and rewrites result/trajectory in the actual CLI. Clock adjustments also affect these wall checks. Caps are cooperative checks, not preemption of an in-flight RHS/frame/observable operation.',
        'evidence': ['run_local_kick.py:54-55,72-79,85,93-115,124,147,169-173', 'pilot-001/allocation.json'],
        'impact_on_completed_result': 'No allocation extension is requested or exercised; the smoke has a separate IO-only budget.',
    },
    {
        'id': 'SR-07', 'kind': 'durability_and_binding_boundary',
        'finding': 'NPZ publication precedes the immutable receipt commit; orphan NPZ/tmp files are not accepted results. Absolute archive paths are part of binding, so copied original receipts are not portable. Writes/renames do not fsync files or directories: process-SIGKILL behavior is not host/power-loss durability.',
        'evidence': ['time_pilot_checkpoints.py:28-46,50-67', '../0020/pilot_checkpoints.py:13-18', '../0020/pilot_inputs.py:18-20'],
        'impact_on_completed_result': 'Original receipts remain untouched. All scratch rebinding is explicit fixture construction and is not claimed as portable original-generation recovery.',
    },
]


def bind(path):
    raw = path.read_bytes()
    return {'path': str(path), 'bytes': str(len(raw)), 'sha256': hashlib.sha256(raw).hexdigest()}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def original_inventory():
    paths = sorted(ORIGINAL.rglob('*'))
    require(all(not p.is_symlink() for p in paths), 'original tree contains symlinks')
    return [bind(p) for p in paths if p.is_file()]


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    os.replace(temporary, path)


def timeout_handler(signum, frame):
    raise TimeoutError('bounded IO smoke wall allowance exhausted')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt-dir', required=True, type=Path)
    parser.add_argument('--wall-seconds', type=int, default=90, choices=range(10, 121), metavar='10..120')
    args = parser.parse_args()
    output = args.receipt_dir.resolve()
    require(output != ORIGINAL and ORIGINAL not in output.parents,
            'receipt directory must be outside original data')
    output.mkdir(parents=False, exist_ok=False)
    source_copy = output / Path(__file__).name
    source_copy.write_bytes(Path(__file__).read_bytes())
    scratch = output / 'owned-scratch'
    scratch.mkdir()
    start = time.monotonic()
    report = {
        'schema': 'p254-complete0021-state-recovery-smoke-v1',
        'status': 'running', 'executed_by': 'Root invocation only',
        'command': [sys.executable, *sys.argv],
        'original_directory': str(ORIGINAL), 'original_fingerprint': FINGERPRINT,
        'scope': 'Completed-cache/IO only; no midtrajectory numerical continuation claim.',
        'scientific_boundary': 'Held nonstationary two-atom/two-cell first moments only; no physical loss, domain extent, equilibrium noise, or total actual electronic energy claim.',
        'wall_budget_seconds': args.wall_seconds,
        'cleanup_allowance_seconds': 10,
        'preserved_smoke_source': bind(source_copy),
        'original_smoke_source': bind(Path(__file__).resolve()),
        'findings': FINDINGS, 'checks': [], 'fixture_rebindings': [],
        'scratch_removed': False,
    }
    child = None
    pipe_read = None
    before = None
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.setitimer(signal.ITIMER_REAL, args.wall_seconds)
    try:
        before = original_inventory()
        report['original_bindings_before'] = before
        sources = [bind(path) for path in SOURCE_HASHES]
        for record in sources:
            require(record['sha256'] == SOURCE_HASHES[Path(record['path'])],
                    'original recovery implementation changed: ' + record['path'])
        report['original_source_bindings'] = sources
        identity = json.loads((ORIGINAL / 'identity.json').read_bytes())
        require(identity['fingerprint'] == FINGERPRINT, 'unexpected original generation')
        identity_payload = {key: value for key, value in identity.items() if key != 'fingerprint'}
        require(hashlib.sha256(json.dumps(identity_payload, sort_keys=True).encode()).hexdigest()
                == FINGERPRINT, 'original identity fingerprint does not close')
        report['original_identity'] = identity
        report['original_allocation'] = json.loads((ORIGINAL / 'allocation.json').read_bytes())
        for key in ('registration', 'input_inventory'):
            require(bind(Path(identity[key]['path'])) == identity[key], 'changed ' + key)
        original_report = json.loads((ORIGINAL / 'result.json').read_bytes())
        require(original_report['status'] == 'complete', 'original result is not complete')
        require(bind(ORIGINAL / 'two-cell-trajectory.npz') == original_report['trajectory_archive'],
                'original trajectory binding differs')

        # Prevent imports from creating __pycache__ in the original source tree.
        # Importing these helpers does not call load_inputs or construct models.
        sys.dont_write_bytecode = True
        import numpy as np
        import time_pilot_checkpoints as checkpoints
        import pilot_checkpoints as summary_writer

        with np.load(ORIGINAL / 'two-cell-trajectory.npz', allow_pickle=False) as saved:
            times = saved['tau_ry'].copy()
            masses = saved['masses_ry'].copy()
        require(times.tolist() == [0., .25, .5, .75, 1.], 'unexpected saved output times')
        require(np.array_equal(masses, np.asarray(original_report['masses_ry'])), 'mass metadata differs')
        records = checkpoints.recover_samples(ORIGINAL, FINGERPRINT, times)
        expected = {f'{point}-t{index:02d}' for point in ('gamma', 'half') for index in range(5)}
        require(set(records) == expected, 'not all ten original committed states recovered')
        report['checks'].append({'name': 'original_read_only_recover_samples', 'status': 'pass',
                                 'samples': sorted(records)})

        rebuilt = scratch / 'rebuilt'
        rebuilt.mkdir()
        rebuilt_binding = checkpoints.cell_trajectory(rebuilt, records, masses, times)
        comparisons = []
        with np.load(ORIGINAL / 'two-cell-trajectory.npz', allow_pickle=False) as old, \
                np.load(rebuilt_binding['path'], allow_pickle=False) as new:
            require(set(old.files) == set(new.files), 'trajectory member sets differ')
            for name in sorted(old.files):
                a, b = old[name], new[name]
                require(a.dtype == b.dtype and a.shape == b.shape and
                        a.tobytes(order='C') == b.tobytes(order='C'),
                        'trajectory array is not byte-exact: ' + name)
                comparisons.append({'array': name, 'dtype': str(a.dtype), 'shape': list(a.shape),
                                    'logical_C_bytes_sha256': hashlib.sha256(a.tobytes(order='C')).hexdigest()})
        report['checks'].append({'name': 'rebuilt_saved_trajectory_arrays_exact', 'status': 'pass',
                                 'arrays': comparisons,
                                 'note': 'Exact decoded array bytes, not necessarily identical ZIP container metadata.'})

        # A copied archive gets a NEW fixture receipt binding its private path.
        # The fingerprint is retained solely to exercise the same reader; the
        # original generation is neither moved nor re-pinned by this fixture.
        def fixture(directory, names):
            directory.mkdir()
            for name in names:
                source = records[name]
                archive = directory / (name + '.npz')
                shutil.copyfile(source['archive']['path'], archive)
                receipt = dict(source)
                receipt['archive'] = checkpoints.binding(archive)
                require(receipt['archive']['sha256'] == source['archive']['sha256'],
                        'fixture copy bytes differ')
                write_json(directory / (name + '.json'), receipt)
                report['fixture_rebindings'].append({
                    'sample': name, 'original_archive_binding': source['archive'],
                    'fixture_archive_binding': receipt['archive'],
                    'classification': 'isolated storage fixture, NOT original-generation portability',
                })

        killed = scratch / 'summary-kill'
        fixture(killed, sorted(records))
        require(len(checkpoints.recover_samples(killed, FINGERPRINT, times)) == 10,
                'summary-kill fixture invalid before fault injection')
        summary = killed / 'result.json'
        summary_writer.save_json(summary, {'fixture': True, 'completed_samples': list(records.values())})
        pipe_read, pipe_write = os.pipe()
        child = os.fork()
        if child == 0:
            os.close(pipe_read)
            signal.setitimer(signal.ITIMER_REAL, 0)
            real_replace = summary_writer.os.replace

            def stop_after_summary_rotation(src, dst):
                real_replace(src, dst)
                if Path(src) == summary and Path(dst) == summary.with_suffix('.json.prev'):
                    os.write(pipe_write, b'rotated')
                    while True:
                        signal.pause()

            summary_writer.os.replace = stop_after_summary_rotation
            try:
                # This is the REAL original writer; only the rename boundary is
                # instrumented so the parent kills deterministically in its gap.
                summary_writer.save_json(summary, {'fixture': True, 'new_summary': True})
            finally:
                os._exit(2)
        os.close(pipe_write)
        ready, _, _ = select.select([pipe_read], [], [], min(10, args.wall_seconds))
        require(bool(ready) and os.read(pipe_read, 64) == b'rotated',
                'child did not reach original summary rotation boundary')
        os.kill(child, signal.SIGKILL)
        _, status = os.waitpid(child, 0)
        child = None
        require(os.WIFSIGNALED(status) and os.WTERMSIG(status) == signal.SIGKILL,
                'summary writer did not die by SIGKILL')
        require(not summary.exists() and summary.with_suffix('.json.prev').exists()
                and summary.with_suffix('.json.tmp').exists(), 'wrong injected summary crash boundary')
        recovered = checkpoints.recover_samples(killed, FINGERPRINT, times)
        require(set(recovered) == expected, 'committed states lost after summary SIGKILL')
        for name, record in recovered.items():
            require(record['archive']['sha256'] == records[name]['archive']['sha256'],
                    'committed fixture archive changed after kill')
        report['checks'].append({'name': 'real_child_summary_writer_SIGKILL', 'status': 'pass',
                                 'signal': 'SIGKILL', 'recovered_samples': len(recovered),
                                 'result_json_absent': True, 'previous_and_temporary_present': True,
                                 'actual_CLI_resume_tested': False,
                                 'boundary': 'after original summary -> .prev; before .tmp -> summary'})

        # Independent negative cases: wrong generation uses intact originals;
        # corruption uses the CORRECT generation and a separate private fixture.
        try:
            checkpoints.recover_samples(ORIGINAL, 'wrong-generation-' + FINGERPRINT, times)
        except ValueError as error:
            require('source/bytes differ' in str(error), 'wrong-generation refusal was unrelated')
            report['checks'].append({'name': 'wrong_generation_refused', 'status': 'pass',
                                     'exception': str(error)})
        else:
            raise AssertionError('wrong generation was accepted')
        corrupt = scratch / 'corrupt-archive'
        fixture(corrupt, ['gamma-t00'])
        require(len(checkpoints.recover_samples(corrupt, FINGERPRINT, times)) == 1,
                'corrupt fixture invalid before corruption')
        archive = corrupt / 'gamma-t00.npz'
        with archive.open('r+b') as stream:
            stream.seek(-1, os.SEEK_END)
            last = stream.read(1)
            stream.seek(-1, os.SEEK_END)
            stream.write(bytes([last[0] ^ 1]))
        try:
            checkpoints.recover_samples(corrupt, FINGERPRINT, times)
        except ValueError as error:
            require('source/bytes differ' in str(error), 'corrupt archive refusal was unrelated')
            report['checks'].append({'name': 'corrupt_archive_refused_correct_generation', 'status': 'pass',
                                     'exception': str(error), 'corrupt_binding': bind(archive)})
        else:
            raise AssertionError('corrupt archive was accepted')
        report['status'] = 'pass'
    except BaseException as error:
        report['status'] = 'fail'
        report['failure'] = {'type': type(error).__name__, 'message': str(error),
                             'traceback': traceback.format_exc()}
    finally:
        # Separate, finite shutdown allowance; no retries or numerical work.
        signal.setitimer(signal.ITIMER_REAL, 10)
        try:
            if child:
                try:
                    os.kill(child, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                os.waitpid(child, 0)
            if pipe_read is not None:
                os.close(pipe_read)
            if before is not None:
                after = original_inventory()
                report['original_bindings_after'] = after
                report['original_data_byte_immutable'] = before == after
                if before != after:
                    report['status'] = 'fail'
            else:
                report['original_data_byte_immutable'] = None
        except BaseException as error:
            report['status'] = 'fail'
            report['cleanup_failure'] = {'type': type(error).__name__, 'message': str(error)}
        finally:
            try:
                shutil.rmtree(scratch)
                report['scratch_removed'] = not scratch.exists()
            except BaseException as error:
                report['status'] = 'fail'
                report['scratch_removal_failure'] = {'type': type(error).__name__, 'message': str(error)}
            signal.setitimer(signal.ITIMER_REAL, 0)
            report['elapsed_seconds'] = time.monotonic() - start
            write_json(output / 'receipt.json', report)
    print(json.dumps({'status': report['status'], 'receipt': str(output / 'receipt.json'),
                      'scratch_removed': report['scratch_removed']}), flush=True)
    return 0 if report['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())

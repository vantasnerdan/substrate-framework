#!/usr/bin/env python3
"""Durable canonical oracle: immutable receipts and actual signal/checkpoint recovery."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SRC = ROOT / 'src/substrate_framework'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exclusive_json(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())


def source_paths():
    local = list(HERE.glob('*.py'))+[HERE/'registration.json']
    reusable = [SRC/f'finite_fermion_{name}.py' for name in
                ('state', 'response', 'thermodynamics', 'bath')]
    existing = [SRC/'material_matter_field.py', SRC/'verification.py',
                HERE.parent/'0001/derive_blocks.py',
                HERE.parents[1]/'sources/giustino2017-electron-phonon.txt']
    return sorted(local+reusable+existing)


def prepare(manifest_path):
    value = {'schema_version': 1, 'created_ns': time.time_ns(),
             'purpose': 'immutable source registration before any oracle execution',
             'files': {str(p.relative_to(ROOT)): digest(p) for p in source_paths()}}
    exclusive_json(manifest_path, value)
    print('REGISTERED', digest(manifest_path), flush=True)


class Tee:
    def __init__(self, file, terminal):
        self.file, self.terminal = file, terminal
    def write(self, text):
        if not self.file.closed:
            self.file.write(text)
            self.file.flush()
        self.terminal.write(text)
        self.terminal.flush()
    def flush(self):
        if not self.file.closed:
            self.file.flush()
        self.terminal.flush()


def block_sources(name):
    dependencies = {
        'state': ('state',),
        'response': ('state', 'response'),
        'thermodynamics': ('state', 'response', 'thermodynamics'),
        'bath': ('state', 'response', 'thermodynamics', 'bath')}
    paths = [HERE/'fixtures.py', HERE/'registration.json', HERE/f'check_{name}.py',
             SRC/'verification.py']
    paths += [SRC/f'finite_fermion_{part}.py' for part in dependencies[name]]
    if name in ('thermodynamics', 'bath'):
        paths.append(HERE/'free_energy_oracle.py')
    if name == 'bath':
        paths.append(SRC/'material_matter_field.py')
    return {str(p.relative_to(ROOT)): digest(p) for p in paths}


def cached_block(name, function, ledger, manifest_hash):
    sources = block_sources(name)
    key = hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest()
    path = HERE/f'{name}-checkpoint-{key}.json'
    if path.exists():
        saved = json.loads(path.read_text())
        if digest(path) != path.with_suffix('.sha256').read_text().strip():
            raise RuntimeError(f'{name} checkpoint hash mismatch')
        ledger.passed.extend(saved['checks'])
        print('RECOVERED_BLOCK', name, digest(path), flush=True)
        return saved['result']
    offset = len(ledger.passed)
    result = function(ledger)
    exclusive_json(path, {'manifest_sha256': manifest_hash, 'source_hashes': sources,
                         'result': result, 'checks': ledger.passed[offset:]})
    with path.with_suffix('.sha256').open('x') as stream:
        stream.write(digest(path)+'\n')
        stream.flush()
        os.fsync(stream.fileno())
    return result


def execute(interrupt, manifest_path):
    started = time.monotonic()
    run_id = f'{time.time_ns()}-{os.getpid()}'
    manifest_hash = digest(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    if any(digest(ROOT/p) != h for p, h in manifest['files'].items()):
        raise RuntimeError('source manifest mismatch; new immutable generation required')
    with (HERE/f'invocation-{run_id}.log').open('x') as log:
        sys.stdout = Tee(log, sys.stdout)
        sys.stderr = Tee(log, sys.stderr)
        evidence = {'run_id': run_id, 'pid': os.getpid(), 'argv': sys.argv,
                    'affinity': sorted(os.sched_getaffinity(0)), 'manifest_sha256': manifest_hash,
                    'environment': {k: os.environ.get(k) for k in
                        ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                         'NUMEXPR_NUM_THREADS', 'PYTHONPATH')},
                    'proc_status': Path('/proc/self/status').read_text(),
                    'cgroup': Path('/proc/self/cgroup').read_text()}
        print('ACTUAL_ORACLE_PID', os.getpid(), 'AFFINITY', evidence['affinity'], flush=True)
        exclusive_json(HERE/f'invocation-{run_id}-start.json', evidence)
        try:
            if evidence['affinity'] != [5] or any(evidence['environment'][k] != '1'
                    for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')):
                raise RuntimeError('singleton CPU5 / BLAS1 contract violated')
            from substrate_framework.verification import CheckLedger
            from check_state import run as state
            from check_response import run as response
            from check_thermodynamics import run as thermodynamics
            from check_bath import run as bath
            ledger = CheckLedger('P254-attempt0007-finite-fermion')
            checkpoint = HERE/'state-checkpoint.json'
            if checkpoint.exists():
                expected = (HERE/'state-checkpoint.sha256').read_text().strip()
                if digest(checkpoint) != expected:
                    raise RuntimeError('checkpoint content hash mismatch')
                saved = json.loads(checkpoint.read_text())
                if saved['source_hashes'] != block_sources('state'):
                    raise RuntimeError('state checkpoint scientific source mismatch')
                ledger.passed.extend(saved['checks'])
                state_result = saved['state']
                print('RECOVERED_STATE_CHECKPOINT', expected, 'STATE_BLOCK_NOT_RECOMPUTED', flush=True)
                recovery = {'reused': True, 'checkpoint_sha256': expected,
                            'original_pid': saved['pid'], 'current_pid': os.getpid()}
            else:
                state_result = state(ledger)
                exclusive_json(checkpoint, {'manifest_sha256': manifest_hash,
                    'source_hashes': block_sources('state'),
                    'checks': ledger.passed, 'state': state_result, 'pid': os.getpid()})
                with (HERE/'state-checkpoint.sha256').open('x') as stream:
                    stream.write(digest(checkpoint)+'\n')
                    stream.flush()
                    os.fsync(stream.fileno())
                recovery = {'reused': False, 'checkpoint_sha256': digest(checkpoint)}
                print('SAFE_STATE_CHECKPOINT_READY', digest(checkpoint), flush=True)
            if interrupt:
                def terminate(signum, frame):
                    exclusive_json(HERE/f'invocation-{run_id}-signal.json',
                        {'signal': signum, 'checkpoint_sha256': digest(checkpoint),
                         'pid': os.getpid(), 'manifest_sha256': manifest_hash})
                    print('ACTUAL_SIGNAL_CHECKPOINT_PRESERVED', signum, flush=True)
                    raise SystemExit(128+signum)
                signal.signal(signal.SIGTERM, terminate)
                print('WAITING_FOR_ACTUAL_SIGTERM', flush=True)
                while True:
                    signal.pause()
            results = {'state': state_result}
            for name, function in (('response', response), ('thermodynamics', thermodynamics), ('bath', bath)):
                results[name] = cached_block(name, function, ledger, manifest_hash)
            ledger.finish()
            registration = json.loads((HERE/'registration.json').read_text())
            receipt = {'schema_version': 1, 'status': 'completed',
                'route_verdict': 'established_as_stated_for_encoded_finite_fixtures',
                'evidence_scope': ['exact finite symbolic algebra', 'independent numerical crosschecks'],
                'source_manifest_sha256': manifest_hash, 'registration_sha256': digest(HERE/'registration.json'),
                'registration': registration, 'run': evidence, 'recovery': recovery,
                'checks': ledger.passed, 'check_count': len(ledger.passed), 'results': results,
                'elapsed_seconds': time.monotonic()-started,
                'max_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                'limitations': registration['does_not_license']}
            exclusive_json(HERE/'receipt.json', receipt)
            print('COMPLETE_RECEIPT_SHA256', digest(HERE/'receipt.json'), flush=True)
        except SystemExit:
            raise
        except BaseException:
            failure = {'run_id': run_id, 'manifest_sha256': manifest_hash,
                       'traceback': traceback.format_exc(), 'elapsed_seconds': time.monotonic()-started}
            exclusive_json(HERE/f'invocation-{run_id}-failure.json', failure)
            traceback.print_exc()
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--interrupt-after-state', action='store_true')
    parser.add_argument('--manifest', type=Path, default=HERE/'source-manifest.json')
    args = parser.parse_args()
    prepare(args.manifest) if args.prepare else execute(args.interrupt_after_state, args.manifest)

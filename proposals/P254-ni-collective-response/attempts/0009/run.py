#!/usr/bin/env python3
"""Run only the preregistered force-print enclosure correction, never DFT."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import tempfile

import numpy as np
import scipy
import sympy

import check_enclosure

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / '0008'


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def atomic_json(path, value):
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix='.' + path.name)
    with os.fdopen(descriptor, 'w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    directory = os.open(path.parent, os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def runtime():
    group = next(line.split(':', 2)[2] for line in Path('/proc/self/cgroup').read_text().splitlines()
                 if line.startswith('0::'))
    root = Path('/sys/fs/cgroup' + group)
    return {'pid': os.getpid(), 'affinity': sorted(os.sched_getaffinity(0)),
            'cgroup': group, 'cgroup_observed': {
                name: (root / name).read_text().strip()
                for name in ('memory.max', 'memory.current', 'memory.peak', 'memory.events', 'cpu.stat')},
            'environment': {name: os.environ.get(name) for name in
                            ('PYTHONPATH', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                             'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')},
            'python': platform.python_version(), 'numpy': np.__version__,
            'scipy': scipy.__version__, 'sympy': sympy.__version__}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest_path = HERE / 'source-manifest.json'
    manifest = json.loads(manifest_path.read_text())
    manifest_sha = digest(manifest_path)
    for path, expected in manifest['files'].items():
        if digest(path) != expected:
            raise ValueError('Frozen source/input mismatch: ' + path)
    args.output.mkdir(parents=True, exist_ok=False)
    observed_start = runtime()
    assert observed_start['affinity'] == [5]
    assert observed_start['cgroup_observed']['memory.max'] == '2147483648'
    assert all(observed_start['environment'][name] == '1' for name in
               ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'))
    atomic_json(args.output / 'runtime-start.json', observed_start)
    primary = json.loads((PRIOR / 'receipt.json').read_text())
    factor = primary['blocks']['units']['result']['force_factor']
    prior_force = primary['blocks']['forces']['result']
    registration = json.loads((HERE / 'registration.json').read_text())
    stdout = next(path for path in registration['inputs'] if path.endswith('/pw.stdout'))
    parser_spec = importlib.util.spec_from_file_location('frozen_primary_native_parse', PRIOR / 'native_parse.py')
    native_parser = importlib.util.module_from_spec(parser_spec)
    parser_spec.loader.exec_module(native_parser)
    rows = native_parser.parse_stdout(Path(stdout).read_text())['forces'] * factor
    assert rows.shape == (37, 3)
    radius = 5e-9 * factor
    gate = prior_force['gate_ev_per_angstrom']
    ledger, blocks = [], {}
    tasks = [('component_box', lambda: check_enclosure.theory(ledger)),
             ('native_enclosure', lambda: check_enclosure.native(
                 ledger, rows, radius, gate, prior_force['max_rounding_upper_ev_per_angstrom'])),
             ('exposing_witnesses', lambda: check_enclosure.mutations(ledger, radius, gate))]
    for name, work in tasks:
        blocks[name] = work()
        atomic_json(args.output / (name + '.json'), {
            'source_manifest_sha256': manifest_sha, 'result': blocks[name]})
        print('FORCE_BOX_BLOCK_COMPLETE', name, flush=True)
    for path, expected in manifest['files'].items():
        if digest(path) != expected:
            raise ValueError('Source/input changed during execution: ' + path)
    assert digest(manifest_path) == manifest_sha
    receipt = {'schema_version': 1, 'status': 'complete', 'source_manifest_sha256': manifest_sha,
               'checks': ledger, 'check_count': len(ledger), 'blocks': blocks,
               'component_half_quantum_ev_per_angstrom': radius,
               'runtime_start': observed_start, 'runtime_end': runtime(),
               'scope': 'Bounded source-print enclosure correction only; original native results/receipts unchanged. No tested checkpoint restoration or physical convergence license.'}
    atomic_json(args.output / 'receipt.json', receipt)
    print('FORCE_BOX_CORRECTION_COMPLETE', len(ledger), digest(args.output / 'receipt.json'), flush=True)


if __name__ == '__main__':
    main()

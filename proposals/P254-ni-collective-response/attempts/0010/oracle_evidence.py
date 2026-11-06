"""Check receipts, atomic publication and directly observed process resources."""
import hashlib
import json
import os
from pathlib import Path
import platform
import tempfile

import numpy as np
import sympy as sp


def require(ledger, name, condition):
    if not condition:
        raise AssertionError(name)
    ledger.append(name)


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
                name: (root / name).read_text().strip() for name in
                ('memory.max', 'memory.current', 'memory.peak', 'memory.events',
                 'memory.swap.max', 'memory.swap.current', 'cpu.stat')},
            'environment': {name: os.environ.get(name) for name in
                            ('PYTHONPATH', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                             'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')},
            'python': platform.python_version(), 'numpy': np.__version__, 'sympy': sp.__version__}


def source_guard(files):
    for path, expected in files.items():
        if digest(path) != expected:
            raise ValueError('Frozen source/input mismatch: ' + path)

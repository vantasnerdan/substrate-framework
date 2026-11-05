"""Literal source ownership and actual runtime for paid long-history columns."""
from pathlib import Path
import json
import os
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/n) for n in
    ('0060','0052','0051','0046','0044','0045','0041','0030','0029')]
from moving_reference_inputs import atomic_json, file_binding

CHANNELS = ('greater','lesser','retarded')
BUNDLE = Path('/var/tmp/p254-a52-requested-reference-001/producer/bundle.json')


def bound_items(tree):
    if isinstance(tree, dict):
        if {'path','bytes','sha256'} <= tree.keys():
            yield {key: tree[key] for key in ('path','bytes','sha256')}
        for value in tree.values():
            yield from bound_items(value)
    elif isinstance(tree, list):
        for value in tree:
            yield from bound_items(value)


def authenticate(item):
    if file_binding(item['path']) != item:
        raise ValueError('paid input bytes changed: '+item['path'])
    return Path(item['path'])


def begin(output, cpu):
    output = Path(output)
    output.mkdir(exist_ok=False)
    group = Path('/proc/self/cgroup').read_text().strip().split('::')[-1]
    runtime = dict(affinity=sorted(os.sched_getaffinity(0)),
        threads={key:os.environ.get(key) for key in
            ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')},
        limits={key:(Path('/sys/fs/cgroup')/group.lstrip('/')/key).read_text().strip()
                for key in ('memory.max','memory.swap.max')})
    atomic_json(output/'actual-runtime.json', runtime)
    if runtime['affinity'] != [cpu] or set(runtime['threads'].values()) != {'1'} or runtime['limits'] != {'memory.max':'3221225472','memory.swap.max':'0'}:
        raise ValueError('actual assigned one-thread3GiB/noSwap runtime differs')
    spec = json.loads((HERE/'preregistration.json').read_bytes())
    return output, spec

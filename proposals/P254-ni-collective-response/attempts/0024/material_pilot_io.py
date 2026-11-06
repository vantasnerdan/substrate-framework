"""Source closure, own context and immutable exact-readback material blocks.

Atomic archive/receipt protocol extends 0023 noise_pilot_io.py, original
4785 bytes SHA256 ae9d2571d40951ec06a886ac17132b3f91ff36acce1eccc3496da481b4e62280.
No notification owner is defined here; Root wraps the established notifier.
Admission schema is documented in validate_admission; no implicit resource grant.
"""
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import uuid
import zipfile
import numpy as np
import numpy.lib.format
import numpy.fft
import numpy.linalg
import xml.etree.ElementTree
from material_pilot_inputs import HERE, DIRECTIONS, bound_json

# Exercise ONLY serialization dispatch, not scientific computation, before
# source inventory. This forces NPZ writer/reader lazy consumers to load.
_buffer = io.BytesIO()
np.savez(_buffer, preload=np.empty(0, dtype=np.float64))
_buffer.seek(0)
with np.load(_buffer, allow_pickle=False) as _loaded:
    _loaded['preload']
del _buffer, _loaded


def binding(path, raw=None):
    path = Path(path).resolve()
    raw = path.read_bytes() if raw is None else raw
    return {'path': str(path), 'bytes': str(len(raw)), 'sha256': hashlib.sha256(raw).hexdigest()}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def implementation():
    paths = {HERE/name for name in ('material_pilot_inputs.py', 'material_pilot_io.py',
                                     'run_material_ionic_pilot.py', 'run_material_ionic_oracle.py')}
    for module in tuple(sys.modules.values()):
        filename = getattr(module, '__file__', None)
        if filename:
            path = Path(filename).resolve()
            if path.is_file() and path.suffix in ('.py', '.so', '.pyc'):
                paths.add(path)
    return [binding(p) for p in sorted(paths)]


def runtime_context():
    relative = next(line.split(':', 2)[2] for line in Path('/proc/self/cgroup').read_text().splitlines()
                    if line.startswith('0::'))
    base = Path('/sys/fs/cgroup')/relative.lstrip('/')
    resources = {}
    for key in ('memory.max', 'memory.swap.max', 'cpuset.cpus.effective'):
        try:
            resources[key] = {'availability': 'available', 'value': (base/key).read_text().strip()}
        except OSError as error:
            resources[key] = {'availability': 'unavailable', 'reason': str(error)}
    return {'pid': os.getpid(), 'uid': os.getuid(), 'gid': os.getgid(), 'groups': os.getgroups(),
            'cpu_affinity': sorted(os.sched_getaffinity(0)), 'cgroup': relative,
            'thread_variables': {k: os.environ.get(k) for k in
                ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')},
            'INVOCATION_ID': os.environ.get('INVOCATION_ID'), 'resources': resources,
            'kernel_limits_text': Path('/proc/self/limits').read_text(), 'python': sys.version,
            'numpy': np.__version__}


def validate_admission(path, sources, context, context_only=False):
    """JSON REQUIRED fields:

    registration, source_access: {path,bytes,sha256}; implementation_sources:
    list of same bindings for the preloaded implementation() result (superset
    allowed, every entry verified); remaining_seconds:0<number<=300;
    first_block_seconds:0<number<=90; output_max_bytes:integer<=1073741824;
    directions:exact DIRECTIONS; block_bands:integer1..16;
    resume:bool; resume_identity:null for fresh or prior summary recipe_sha256;
    context:{uid:int,gid:int,groups:[int],cgroup:str,cpu_affinity:[int]} exact.
    Context-only emits observed sources/context even when admission fails, but
    never reads material or builds radial tables. A corrected explicit grant
    is required for producer; remaining_seconds is NOT reset on recovery.
    """
    raw = Path(path).read_bytes()
    grant = json.loads(raw)
    registration, _ = bound_json(grant['registration'])
    inventory, _ = bound_json(grant['source_access'])
    if grant['registration']['path'] != str(HERE/'implementation-registration.json'):
        raise ValueError('wrong registration path')
    if grant['source_access'] != registration['access']:
        raise ValueError('source-access differs from registered authority')
    expected = {str(Path(x['path']).resolve()): x for x in grant['implementation_sources']}
    if len(expected) != len(grant['implementation_sources']):
        raise ValueError('duplicate implementation source')
    for record in expected.values():
        if binding(record['path']) != record:
            raise ValueError('implementation byte binding mismatch: '+record['path'])
    for record in sources:
        if expected.get(record['path']) != record:
            raise ValueError('unbound preloaded consumer: '+record['path'])
    for key, cap in (('remaining_seconds', 300), ('first_block_seconds', 90)):
        if type(grant[key]) not in (int, float) or not 0 < grant[key] <= cap:
            raise ValueError('explicit bounded '+key+' required')
    if type(grant['output_max_bytes']) is not int or not 0 < grant['output_max_bytes'] <= 1073741824:
        raise ValueError('output cap invalid')
    if type(grant['block_bands']) is not int or not 1 <= grant['block_bands'] <= 16:
        raise ValueError('band block cap invalid')
    if grant['directions'] != DIRECTIONS or type(grant['resume']) is not bool:
        raise ValueError('directions/recovery admission invalid')
    if (grant['resume'] and not isinstance(grant['resume_identity'], str)) or (not grant['resume'] and grant['resume_identity'] is not None):
        raise ValueError('explicit recovery identity required')
    for key in ('uid', 'gid', 'groups', 'cgroup', 'cpu_affinity'):
        if grant['context'][key] != context[key]:
            raise ValueError('actual manager context differs: '+key)
    if not context_only:
        if not context['INVOCATION_ID'] or '/user.slice/' not in context['cgroup']:
            raise ValueError('actual user-manager invocation required')
        if not 1 <= len(context['cpu_affinity']) <= 6:
            raise ValueError('CPU affinity exceeds six')
        for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
            if context['thread_variables'][key] != '1':
                raise ValueError('single-thread BLAS/OpenMP required')
        for key, limit in (('memory.max', 4294967296), ('memory.swap.max', 0)):
            entry = context['resources'][key]
            if entry['availability'] != 'available' or not entry['value'].isdigit() or int(entry['value']) > limit:
                raise ValueError('unadmitted cgroup resource: '+key)
    recipe = {'registration': grant['registration'], 'source_access': grant['source_access'],
              'implementation_sources': grant['implementation_sources'], 'directions': DIRECTIONS,
              'block_bands': grant['block_bands'], 'inputs': inventory['material_inputs'],
              'algorithm': 'Gamma-bare-local-plus-both-NC-derivatives-v1',
              'units': 'Ry/bohr', 'bands': 346, 'spins': [1, 2], 'npw': 18613, 'nat': 37}
    identity = hashlib.sha256(canonical(recipe)).hexdigest()
    if grant['resume'] and grant['resume_identity'] != identity:
        raise ValueError('recovery source/input/recipe identity mismatch')
    return grant, inventory, recipe, identity, binding(path, raw)


def save_json(path, value, rotate=False):
    path = Path(path)
    temporary = path.with_name(path.name+'.tmp-'+uuid.uuid4().hex)
    with temporary.open('xb') as stream:
        stream.write(canonical(value)+b'\n')
        stream.flush()
        os.fsync(stream.fileno())
    if rotate and path.exists():
        os.replace(path, path.with_name(path.name+'.prev'))
    os.replace(temporary, path)


def array_record(array):
    contiguous = np.ascontiguousarray(array)
    return {'dtype': array.dtype.str, 'shape': list(array.shape),
            'sha256_C_bytes': hashlib.sha256(memoryview(contiguous).cast('B')).hexdigest()}


def recover_block(target, key, identity, metadata):
    receipt_path, archive = target/(key+'.json'), target/(key+'.npz')
    if not receipt_path.exists() and not archive.exists():
        return None
    if not receipt_path.exists() or not archive.exists():
        raise ValueError('incomplete atomic block; no interrupted-operation continuation: '+key)
    receipt = json.loads(receipt_path.read_bytes())
    if receipt['recipe_sha256'] != identity or receipt['metadata'] != metadata or binding(archive) != receipt['archive']:
        raise ValueError('completed block identity mismatch: '+key)
    with np.load(archive, allow_pickle=False) as loaded:
        if set(loaded.files) != set(receipt['arrays']):
            raise ValueError('recovery array membership differs')
        for name in loaded.files:
            if array_record(loaded[name]) != receipt['arrays'][name]:
                raise ValueError('recovery exact array digest differs: '+name)
    return receipt


def commit_block(target, key, arrays, metadata, identity, cap):
    """Archive is not visible until exact dtype/shape/value readback succeeds."""
    path = target/(key+'.npz')
    receipt_path = target/(key+'.json')
    if path.exists() or receipt_path.exists():
        raise FileExistsError('immutable block already exists: '+key)
    used = sum(p.stat().st_size for p in target.iterdir() if p.is_file())
    if used + sum(a.nbytes for a in arrays.values()) + 65536 > cap:
        raise RuntimeError('output byte cap would be exceeded')
    temp = target/(key+'.tmp-'+uuid.uuid4().hex)
    try:
        with temp.open('xb') as stream:
            np.savez(stream, **arrays)
            stream.flush()
            os.fsync(stream.fileno())
        with np.load(temp, allow_pickle=False) as loaded:
            if set(loaded.files) != set(arrays):
                raise ValueError('output readback membership differs')
            for name, expected in arrays.items():
                actual = loaded[name]
                if actual.dtype != expected.dtype or actual.shape != expected.shape or not np.array_equal(actual, expected):
                    raise ValueError('exact saved array readback failed: '+name)
        receipt = {'recipe_sha256': identity, 'archive': binding(path, temp.read_bytes()),
                   'arrays': {k: array_record(v) for k, v in arrays.items()}, 'metadata': metadata,
                   'readback': 'exact dtype/shape/value; allow_pickle=False'}
        # Hard-link is an atomic no-overwrite commit on this same filesystem.
        os.link(temp, path)
        with receipt_path.open('xb') as stream:
            stream.write(canonical(receipt)+b'\n')
            stream.flush()
            os.fsync(stream.fileno())
        return receipt
    finally:
        if temp.exists():
            temp.unlink()

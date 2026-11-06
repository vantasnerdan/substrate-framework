"""Source-bound spectral recipes and exact saved-array readback for P2540023."""
from dataclasses import fields
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
# Preload the chosen saved-output reader before the implementation inventory.
import numpy.lib.format
import zipfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'0020'))
from pilot_inputs import binding, load_inputs
from pilot_checkpoints import save_json
from substrate_framework.qe_frozen_curvature import verify_buffer


def bound_registration():
    raw = (HERE/'registration.json').read_bytes()
    registration = json.loads(raw)
    record = registration['source_inventory']
    inventory_raw = Path(record['path']).read_bytes()
    verify_buffer(inventory_raw, record)
    for source in json.loads(inventory_raw)['files']:
        verify_buffer(Path(source['path']).read_bytes(), source)
    return registration, binding(HERE/'registration.json', raw)


def implementation():
    paths = {Path(__file__).resolve(), HERE/'run_prepared_noise.py'}
    for name, module in list(sys.modules.items()):
        filename = getattr(module, '__file__', None)
        if filename and (name.startswith(('substrate_framework', 'catdft', 'numpy.lib')) or name in ('pilot_inputs','pilot_checkpoints')):
            path = Path(filename).resolve()
            if path.suffix == '.py':
                paths.add(path)
    return [binding(path) for path in sorted(paths)]


def runtime_context():
    """Observe this actual non-Native consumer under the user service manager."""
    relative = next(line.split(':',2)[2] for line in Path('/proc/self/cgroup').read_text().splitlines()
                    if line.startswith('0::'))
    cgroup = Path('/sys/fs/cgroup')/relative.lstrip('/')
    resources = {}
    for name in ('memory.max','memory.swap.max','cpuset.cpus.effective'):
        try:
            resources[name] = {'availability':'available','value':(cgroup/name).read_text().strip()}
        except FileNotFoundError:
            resources[name] = {'availability':'unavailable','reason':'not exposed in own cgroup'}
    return {'pid':str(os.getpid()),'uid':str(os.getuid()),'gid':str(os.getgid()),
            'groups':[str(value) for value in os.getgroups()],
            'cpu_affinity':[str(cpu) for cpu in sorted(os.sched_getaffinity(0))],
            'thread_variables':{key:os.environ.get(key) for key in
                ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')},
            'invocation':os.environ.get('INVOCATION_ID'),'cgroup':relative,'resources':resources,
            'kernel_limits_text':Path('/proc/self/limits').read_text(),
            'python':sys.version,'numpy':np.__version__}


def coordinate_bindings(model, cache):
    """Identify actual spectral coordinates without claiming amplitude interchange."""
    def array(value):
        key = id(value)
        if key not in cache:
            contiguous = np.ascontiguousarray(value)
            cache[key] = {'dtype':value.dtype.str,'shape':list(value.shape),
                          'sha256_C_bytes':hashlib.sha256(memoryview(contiguous).cast('B')).hexdigest()}
        return cache[key]
    return {'q':array(model.q_cart_2pi_alat),
            'channels':[{field.name:(array(value) if isinstance(value,np.ndarray) else value)
                         for field in fields(channel) for value in [getattr(channel,field.name)]}
                        for channel in model.channels],
            'scope':'Actual basis identity only; no map to prior delta-psi amplitudes or coupled quantum state is implied.'}


def save_arrays(target, key, arrays, metadata):
    path = target/(key+'.npz')
    with path.open('xb') as stream:
        np.savez_compressed(stream, **arrays)
    with np.load(path, allow_pickle=False) as loaded:
        if set(loaded.files) != set(arrays):
            raise ValueError('saved spectral recipe membership differs')
        for name, expected in arrays.items():
            expected = np.asarray(expected)
            actual = loaded[name]
            if actual.dtype != expected.dtype or actual.shape != expected.shape or not np.array_equal(actual, expected):
                raise ValueError('saved spectral recipe changed: '+key+'/'+name)
    receipt = {'archive':binding(path),'arrays':{name:{'dtype':np.asarray(value).dtype.str,
                 'shape':list(np.shape(value))} for name,value in arrays.items()},
               'metadata':metadata,'readback':'every array exact dtype/shape/value; allow_pickle=False'}
    save_json(target/(key+'.json'),receipt)
    return receipt

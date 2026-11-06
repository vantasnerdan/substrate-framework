"""Once-global joint on an actual saved provider, using cancellation-free Sigma.

API: produce_joint(provider, output, *, resume=False). Resume only with an
identical request; bitmap-committed pairs are read/validated, never recomputed.
CLI: python produce_joint.py --provider G1/provider.json --output NEW [--resume]
"""
import argparse
from contextlib import contextmanager
import fcntl
import json
import hashlib
import os
from pathlib import Path
import resource
import time

from joint_blocks import HERE, joint_pair
import numpy as np
from history_banks import ORDERS, authenticate, binding, save
from history_sources import open_reference
from joint_material_inputs import load_inputs
from joint_gamma_adapter import vertices
from electron_ion_functional import evaluate_components

PAIRS = tuple((a,b) for a in range(31) for b in range(31))


def atomic_save(path, record):
    """Owning atomic JSON writer, plus directory durability after rename."""
    path = Path(path)
    save(path, record)
    fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def contract():
    return json.loads((HERE/'preregistration.json').read_bytes())


@contextmanager
def output_lock(output, resume):
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=resume)
    if resume and not (out/'checkpoint.json').is_file():
        raise ValueError('resume requires an existing durable checkpoint')
    with (out/'writer.lock').open('a+b') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield out


def material_reference(provider):
    """Original coefficient ownership; actual electronic/Gamma input is G1."""
    c, baseline, old, coefficients = load_inputs()
    old.close()
    del old, baseline
    ref = open_reference(c, coefficients['manifest'], provider=provider)
    if ((c.ns,c.m,c.d,c.f,c.n) != (128,12,111,3,64) or
            not np.array_equal(ref.times, np.arange(31))):
        ref.close()
        raise ValueError('full31/all128/full111+3 original axes required')
    return c, ref, coefficients


def validate_pair(bank, a, b):
    for values in bank.values():
        if not np.isfinite(values).all():
            raise ValueError('nonfinite persisted scientific pair')
    expected = bank['greater']-bank['lesser'] if a >= b else np.zeros_like(bank['retarded'])
    scale = max(1., float(np.max(np.abs(expected))))
    if float(np.max(np.abs(bank['retarded']-expected))) > 2e-12*scale:
        raise ValueError('persisted raw theta(0)=1 causal identity differs')


def _pair_digest(values):
    digest = hashlib.sha256()
    for name in sorted(values):
        digest.update(name.encode('ascii'))
        digest.update(memoryview(np.ascontiguousarray(values[name])).cast('B'))
    return np.frombuffer(digest.digest(), dtype=np.uint8)


def retain_dependencies(record, dependencies, stores):
    """Never relabel committed pairs with a newly selected implementation."""
    saved = record.get('implementation_bindings')
    if saved is None:
        if any(store.seen.any() for store in stores):
            raise ValueError('committed pairs lack original dependency bindings')
        record['implementation_bindings'] = dependencies
    else:
        for item in saved:
            authenticate(item)
        if saved != dependencies:
            raise ValueError('selected implementation differs from saved dependencies')


class PairStore:
    """Numeric fsync -> pair bitmap fsync -> caller's atomic checkpoint.

    Each store commits all declared arrays for a pair together. On resume the
    durable bitmap is authoritative even after a crash before JSON publication.
    Uncommitted slices may contain interrupted writes, never scientific coverage.
    """
    def __init__(self, directory, tails, *, resume):
        self.directory = Path(directory)
        self.directory.mkdir(exist_ok=resume)
        self.arrays, self.fds, self.seen, self.digests = {}, {}, None, None
        try:
            for name, tail in tails.items():
                path = self.directory/(name+'.npy')
                shape = (31,31,*tail)
                arr = (np.load(path, mmap_mode='r+', allow_pickle=False) if resume else
                       np.lib.format.open_memmap(path, mode='w+', dtype=np.complex128, shape=shape))
                self.arrays[name] = arr
                if arr.shape != shape or arr.dtype != np.complex128:
                    raise ValueError('resumed numeric axes/dtype differ: '+name)
                self.fds[name] = os.open(path, os.O_RDWR)
                if not resume:
                    arr.flush()
                    os.fsync(self.fds[name])
            path = self.directory/'pair-digests.npy'
            self.digests = (np.load(path, mmap_mode='r+', allow_pickle=False) if resume else
                            np.lib.format.open_memmap(path, mode='w+', dtype=np.uint8, shape=(31,31,32)))
            if self.digests.shape != (31,31,32) or self.digests.dtype != np.uint8:
                raise ValueError('invalid persisted pair digests')
            self.fds['digests'] = os.open(path, os.O_RDWR)
            if not resume:
                self.digests[:] = 0
                self.digests.flush()
                os.fsync(self.fds['digests'])
            path = self.directory/'pair-coverage.npy'
            self.seen = (np.load(path, mmap_mode='r+', allow_pickle=False) if resume else
                         np.lib.format.open_memmap(path, mode='w+', dtype=bool, shape=(31,31)))
            if self.seen.shape != (31,31) or self.seen.dtype != np.bool_:
                raise ValueError('invalid persisted pair bitmap')
            self.fds['coverage'] = os.open(path, os.O_RDWR)
            if not resume:
                self.seen[:] = False
                self.seen.flush()
                os.fsync(self.fds['coverage'])
            fd = os.open(self.directory, os.O_RDONLY | os.O_DIRECTORY)
            try:
                os.fsync(fd)
            finally:
                os.close(fd)
            self.validate()
        except BaseException:
            self.close()
            raise

    def completed(self):
        return [[a,b] for a,b in PAIRS if self.seen[a,b]]

    def validate(self):
        prefixes = {name.rsplit('-',1)[0] for name in self.arrays}
        for a,b in PAIRS:
            if self.seen[a,b]:
                actual = _pair_digest({name:arr[a,b] for name,arr in self.arrays.items()})
                if not np.array_equal(actual, self.digests[a,b]):
                    raise ValueError('committed scientific pair bytes changed')
                for prefix in prefixes:
                    validate_pair({o:self.arrays[prefix+'-'+o][a,b] for o in ORDERS}, a,b)

    def commit(self, a, b, values):
        if self.seen[a,b] or set(values) != set(self.arrays):
            raise ValueError('duplicate pair or incomplete numerical commit')
        for name, arr in self.arrays.items():
            value = np.asarray(values[name])
            if value.shape != arr.shape[2:] or not np.isfinite(value).all():
                raise ValueError('invalid numerical pair: '+name)
            arr[a,b] = value
        for name, arr in self.arrays.items():
            arr.flush()
            os.fsync(self.fds[name])
        # Validate the persisted mappings before claiming coverage.
        for prefix in {name.rsplit('-',1)[0] for name in self.arrays}:
            validate_pair({o:self.arrays[prefix+'-'+o][a,b] for o in ORDERS}, a,b)
        self.digests[a,b] = _pair_digest({name:arr[a,b] for name,arr in self.arrays.items()})
        self.digests.flush()
        os.fsync(self.fds['digests'])
        self.seen[a,b] = True
        self.seen.flush()
        os.fsync(self.fds['coverage'])

    def descriptors(self, prefix):
        return {o:dict(binding(self.directory/(prefix+'-'+o+'.npy')),
                       shape=list(self.arrays[prefix+'-'+o].shape), dtype='complex128')
                for o in ORDERS}

    def bitmap(self):
        return dict(binding(self.directory/'pair-coverage.npy'), shape=[31,31], dtype='bool')

    def close(self):
        for arr in self.arrays.values():
            arr._mmap.close()
        self.arrays.clear()
        if self.seen is not None:
            self.seen._mmap.close()
            self.seen = None
        if self.digests is not None:
            self.digests._mmap.close()
            self.digests = None
        for fd in self.fds.values():
            os.close(fd)
        self.fds.clear()


def produce_joint(provider, output, *, resume=False):
    started = time.monotonic()
    source = binding(provider)
    prereg = contract()
    request = dict(input_dataset=source, preparation_origin=prereg['original_preparation'],
                   numerical_reference=prereg['original_numerical_reference'],
                   implementation=[binding(HERE/n) for n in ('joint_blocks.py','produce_joint.py')])
    for name in ('original_preparation','original_numerical_reference'):
        authenticate(prereg[name])
    with output_lock(output, resume) as out:
        if resume:
            record = json.loads((out/'checkpoint.json').read_bytes())
            if record.get('request') != request or record.get('component') != 'joint':
                raise ValueError('resume request/source/implementation differs')
            for item in record.get('implementation_bindings', ()):
                authenticate(item)
        else:
            record = dict(schema='p254-additive-history-component-v1', component='joint',
                q_additive=False, status='partial-global-joint-component', full_mesh=False,
                full_solve_eligible=False, full_time_square=False, request=request,
                input_dataset=source, preparation_origin=prereg['original_preparation'],
                numerical_reference=prereg['original_numerical_reference'],
                times_inverse_ry=list(range(31)), slots=128, bands=12, spins=2,
                ionic_directions=111, gamma_coordinates=114,
                axes=dict(sigma=['time','time','slot','band','band'],
                          pi=['time','time','coordinate','coordinate']),
                sigma_convention='joint_minus_ionic_Gamma', pi_convention='full_joint_Gamma_once',
                sigma_evaluation='same A41 monomial: ionic-field + field-ionic + field-field; no II subtraction',
                normalization='A50 ionic/sqrt(N), normalized field vertices; A41 n=1',
                coverage=dict(owned_pairs=[list(p) for p in PAIRS], completed_pairs=[],
                    external_targets=list(range(128)), all114_coordinates=True,
                    checkpoint='six numeric slices flush/fsync, bitmap flush/fsync, atomic checkpoint+directory fsync'),
                scope='Once-global joint correction and full Pi; not a q component or full collision map',
                timings=[], boundary_loss=dict(incomplete_pair=None, static_rebuild='original coefficients and vertices rebuilt; committed pairs never recomputed'))
        store = ref = None
        active = None
        try:
            store = PairStore(out/'numeric', {kind+'-'+o:tail for kind,tail in
                (('sigma',(128,12,12)),('pi',(114,114))) for o in ORDERS}, resume=resume)
            record['coverage']['completed_pairs'] = store.completed()
            atomic_save(out/'checkpoint.json', record)
            c, ref, coefficients = material_reference(provider)
            canonical = vertices(c)
            dependencies = [binding(HERE/'produce_joint.py'),
                binding(HERE/'joint_blocks.py'), binding(Path(vertices.__code__.co_filename)),
                binding(Path(evaluate_components.__code__.co_filename)),
                binding(Path(open_reference.__code__.co_filename)),
                binding(Path(ref.electron.__func__.__code__.co_filename)),
                binding(Path(load_inputs.__code__.co_filename))]
            retain_dependencies(record, dependencies, (store,))
            record.update(source_bindings=ref.bindings, coefficient_bindings=coefficients,
                canonical_routes=dict(source=list(range(128)), target=np.asarray(c.target[0]).tolist(),
                    slot_spin=np.asarray(c.slot_spin).tolist(), slot_k=np.asarray(c.slot_k).tolist()))
            atomic_save(out/'checkpoint.json',record)
            for a,b in PAIRS:
                if store.seen[a,b]:
                    continue
                active = [a,b]
                tick = time.monotonic()
                values = joint_pair(c,ref,a,b,canonical)
                store.commit(a,b,{key.replace('_','-',1):value for key,value in values.items()})
                record['coverage']['completed_pairs'] = store.completed()
                record['timings'].append(dict(pair=[a,b],seconds=time.monotonic()-tick))
                record['elapsed_seconds_this_invocation'] = time.monotonic()-started
                atomic_save(out/'checkpoint.json',record)
                active = None
                print(json.dumps(dict(event='actual-dressed-joint-pair-committed',pair=[a,b],
                    completed_pairs=len(record['coverage']['completed_pairs']))),flush=True)
            record['failure'] = None
        except BaseException as exc:
            record['failure'] = dict(type=type(exc).__name__,message=str(exc))
            record['boundary_loss']['incomplete_pair'] = active
            raise
        finally:
            if ref is not None:
                ref.close()
            if store is not None:
                try:
                    completed = store.completed()
                    complete = len(completed) == 961
                    record.update(status='complete-global-joint-component' if complete else 'partial-global-joint-component',
                        full_time_square=complete, elapsed_seconds_this_invocation=time.monotonic()-started,
                        native_cpu_seconds=time.process_time(),
                        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
                    record['coverage'].update(completed_pairs=completed,array=store.bitmap())
                    record['sigma'],record['pi'] = store.descriptors('sigma'),store.descriptors('pi')
                    atomic_save(out/'checkpoint.json',record)
                    atomic_save(out/'catalog.json',record)
                finally:
                    store.close()
        # Exercise the existing saved-output reader on the actual published
        # payload, not freshly computed temporary values.
        from importlib.util import spec_from_file_location, module_from_spec
        spec = spec_from_file_location('a87_joint_reader', HERE.parent/'0075/component_reader.py')
        reader = module_from_spec(spec)
        spec.loader.exec_module(reader)
        report = reader.consume(out/'catalog.json')
        atomic_save(out/'validation.json',report)
        return out/'catalog.json'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--provider',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--resume',action='store_true')
    args = parser.parse_args()
    print(produce_joint(args.provider,args.output,resume=args.resume))

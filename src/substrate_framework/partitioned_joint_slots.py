"""Independent pinned slot contractions/propagation, immutable receipts and barriers.

No worker evaluates K or feeds partial n/G/ell back into a Hamiltonian. Each
receipt belongs to an exact state, frozen stage, interval and source identity.
Failed attempt files are retained and never confused with completed receipts.
"""
from dataclasses import asdict, fields
from hashlib import sha256
import json
import multiprocessing as mp
import os
from pathlib import Path
from queue import Empty
import resource
from time import monotonic
from uuid import uuid4

import numpy as np

from .material_physical_operators import PhysicalFrozenAction
from .qe_material_hamiltonian import MaterialHamiltonian
from .streamed_material_state import authenticate, file_binding


def digest(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def array_identity(value):
    a = np.ascontiguousarray(value)
    return dict(shape=list(a.shape), dtype=a.dtype.str, sha256=sha256(memoryview(a).cast('B')).hexdigest())


def stage_identity(stage):
    if stage is None:
        return None
    return {f.name: array_identity(getattr(stage, f.name)) if isinstance(getattr(stage, f.name), np.ndarray)
            else getattr(stage, f.name) for f in fields(stage)}


def publish(path, value, budget):
    """Unique staging name survives interrupted earlier publications without reuse."""
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    temporary = path.with_name(path.name+'.'+uuid4().hex+'.incomplete')
    budget.json(temporary, value)
    # Only the parent or the uniquely assigned slot worker can publish this name.
    temporary.rename(path)
    return file_binding(path)


class SlotExecutor:
    def __init__(self, model, directory, *, cpus, source_identity):
        cpus = tuple(cpus)
        if (not cpus or len(set(cpus)) != len(cpus) or any(type(c) is not int or c not in range(6,11) for c in cpus)
                or not source_identity):
            raise ValueError('explicit unique allocated CPUs6–10 and source identity required')
        self.model, self.directory = model, Path(directory).resolve()
        self.cpus, self.source_identity = cpus, source_identity
        self.directory.mkdir(parents=True, exist_ok=True)
        self._verified = {}

    def verify(self, binding):
        path = Path(binding['path'])
        stat = path.stat()
        key = (str(path), binding['sha256'], binding['bytes'])
        signature = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
        if self._verified.get(key) != signature:
            authenticate(binding)
            self._verified[key] = signature

    def request(self, state, *, frozen_stage=None, interval=0., coordinates=None, operation='reduce'):
        self.model.validate(state)
        if operation not in ('reduce', 'frames') or not np.isfinite(interval) or interval < 0:
            raise ValueError('known operation and nonnegative finite interval required')
        if (frozen_stage is None) != (interval == 0):
            raise ValueError('a positive propagation interval requires its exact frozen stage')
        for slot in state.slots:
            self.verify(slot.frame)
            self.verify(slot.coefficient)
        coords = (state.u, state.p, state.coordinate, state.momentum) if coordinates is None else coordinates
        m, r = self.model, self.model.reference
        parameters = dict(reference={name: array_identity(getattr(r,name)) for name in
            ('nbar','masses','curvature','ionic_force','support_force','held_density','held_gradient')},
            coefficients={name: array_identity(getattr(m.coefficients,name)) for name in
                ('delta_curvature','contact','mixed')}, charges=array_identity(m.charges),
            polarizations=array_identity(m.field.polarizations), field_mass=m.field.mass,
            field_volume=m.field.volume_bohr3, measure=m.measure, ncells=r.ncells,
            reality_tolerance=m.reality_tolerance)
        request = dict(schema='partitioned-joint-slot-request-v1', operation=operation,
            source_identity=self.source_identity, model_identity=state.identity, model_parameters=parameters,
            boundary=file_binding(state.path), slots=[asdict(s) for s in state.slots],
            time=state.time, initial_coordinates=[array_identity(x) for x in
                (state.u, state.p, state.coordinate, state.momentum)],
            frozen_stage=stage_identity(frozen_stage), interval=float(interval),
            reduction_coordinates=[array_identity(x) for x in coords])
        key = digest(request)
        directory = self.directory/key
        directory.mkdir(exist_ok=True)
        manifest = directory/'request.json'
        if manifest.exists():
            if json.loads(manifest.read_text()) != json.loads(json.dumps(request)):
                raise ValueError('stale or modified slot request')
        else:
            publish(manifest, request, self.model.budget)
        return key, directory, coords

    def receipt(self, directory, key, index, slot):
        path = directory/('slot-%04d' % index)/'complete.json'
        if not path.exists():
            return None
        record = json.loads(path.read_text())
        if (record.get('complete') is not True or record.get('request') != key
                or record.get('index') != index or record.get('slot') != slot.identity):
            raise ValueError('stale/mutated completed slot receipt')
        for binding in record['arrays'].values():
            self.verify(binding)
        return record

    def collect(self, directory, key, state):
        """Strict complete barrier; never silently renormalize or reduce a prefix."""
        records = [self.receipt(directory, key, i, slot) for i, slot in enumerate(state.slots)]
        if any(record is None for record in records):
            raise ValueError('incomplete all-slot barrier: no scalar/stage feedback permitted')
        return records

    def _worker(self, cpu, queue, state, stage, interval, operation, directory, key):
        os.sched_setaffinity(0, {cpu})
        from threadpoolctl import threadpool_limits
        with threadpool_limits(limits=1):
            while True:
                try:
                    index = queue.get_nowait()
                except Empty:
                    return
                slot = state.slots[index]
                target = directory/('slot-%04d' % index)
                target.mkdir(exist_ok=True)
                attempt = target/('attempt-'+uuid4().hex)
                attempt.mkdir()
                started = monotonic()
                action_start = self.model.budget.local_actions
                column_start = self.model.budget.local_columns
                self.model.budget.check()
                psi, d = state.arrays(slot)
                with self.model.operators(slot) as (coupling, h0, vertex):
                    if stage is not None:
                        psi = self.model.propagate(coupling, h0, vertex, slot, stage, psi, interval)
                    if operation == 'frames':
                        if psi.shape != slot.shape or psi.dtype != np.complex128 or not np.isfinite(psi).all():
                            raise ValueError('propagated raw frame changed shape/type or is nonfinite')
                        values = dict(frame=psi)
                    else:
                        n, g = (coupling.prepared_moments(psi, d) if hasattr(coupling, 'prepared_moments')
                                else coupling.moments(psi, psi, d))
                        # First L moments only: never recompute contacts/curvatures.
                        ell = vertex.moments(psi, d)['field_gradient']
                        bare = PhysicalFrozenAction(h0).apply if isinstance(h0, MaterialHamiltonian) else h0
                        e = np.einsum('ij,ji->', d, psi.conj().T @ self.model.counted(bare)(psi))
                        values = dict(density=n, gradient=g, ell=ell, one_body=np.asarray(e),
                                      number=np.asarray(self.model.measure*np.sum(n)))
                    arrays = {name: self.model.budget.array(attempt/(name+'.npy'), value)
                              for name, value in values.items()}
                publish(target/'complete.json', dict(complete=True, request=key, index=index,
                    slot=slot.identity, cpu=cpu, blas_threads=1, arrays=arrays,
                    hamiltonian_actions=self.model.budget.local_actions-action_start,
                    action_columns=self.model.budget.local_columns-column_start,
                    array_output_bytes=sum(b['bytes'] for b in arrays.values()),
                    process_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                    elapsed_seconds=monotonic()-started), self.model.budget)
                # Drop this slot's action closures/factors and full frame before
                # constructing the next slot; context exit alone keeps locals.
                del psi, d, coupling, h0, vertex, values, arrays
                if operation == 'reduce':
                    del n, g, ell, bare, e

    def execute(self, state, *, frozen_stage=None, interval=0., coordinates=None, operation='reduce'):
        key, directory, coords = self.request(state, frozen_stage=frozen_stage, interval=interval,
                                             coordinates=coordinates, operation=operation)
        missing = [i for i,s in enumerate(state.slots) if self.receipt(directory,key,i,s) is None]
        if missing:
            ctx = mp.get_context('fork')
            next_index = ctx.Value('q', 0)
            class Tasks:
                def get_nowait(self):
                    with next_index.get_lock():
                        i = next_index.value
                        if i >= len(missing):
                            raise Empty
                        next_index.value += 1
                    return missing[i]
            allocated_cpus = self.cpus[:len(missing)]
            workers = [ctx.Process(target=self._worker, args=(cpu, Tasks(), state, frozen_stage,
                       interval, operation, directory, key)) for cpu in allocated_cpus]
            before = self.model.budget.snapshot()
            journal = directory/('execution-'+uuid4().hex+'.json')
            self.model.budget.reserve(16384)
            publish(journal, dict(status='allocated', cpus=list(allocated_cpus),
                slots=missing, allocation=before), self.model.budget)
            failure = None
            try:
                self.model.budget.check_processes(os.getpid(), [])
                for process in workers:
                    process.start()
                while any(process.is_alive() for process in workers):
                    self.model.budget.check_processes(os.getpid(), [p.pid for p in workers if p.is_alive()])
                    if any(p.exitcode not in (None, 0) for p in workers):
                        raise RuntimeError('slot worker failed; completed receipts retained')
                    for process in workers:
                        process.join(timeout=.1/len(workers))
                if any(process.exitcode != 0 for process in workers):
                    raise RuntimeError('slot worker failed; completed receipts retained')
            except BaseException as error:
                failure = type(error).__name__+': '+str(error)
                raise
            finally:
                for process in workers:
                    if process.is_alive():
                        process.terminate()
                for process in workers:
                    if process.pid is not None:
                        process.join()
                # Terminal accounting must survive an exhausted action/time cap.
                terminal = dict(status='failed' if failure else 'workers-finished', error=failure,
                    allocation=file_binding(journal), budget=self.model.budget.snapshot(),
                    exitcodes=[p.exitcode for p in workers])
                encoded = (json.dumps(terminal, indent=2)+'\n').encode()
                if len(encoded) > 16384:
                    raise RuntimeError('execution terminal journal exceeds reserved space')
                with journal.with_suffix('.terminal.json').open('xb') as stream:
                    stream.write(encoded)
        return key, directory, coords, self.collect(directory, key, state)

    def stage(self, directory, key, state, coords, records):
        m, r = self.model, self.model.reference
        density = np.zeros_like(r.nbar, dtype=complex)
        gradient, ell = np.zeros(r.ndof,complex), np.zeros(3,complex)
        one_body, number = 0j, 0j
        for slot, record in zip(state.slots, records):
            values = {name: np.load(binding['path'], mmap_mode='r', allow_pickle=False)
                      for name,binding in record['arrays'].items()}
            if (values['density'].shape != (r.nbar.shape[0],) or values['gradient'].shape != (r.ndof,)
                    or values['ell'].shape != (3,) or values['one_body'].shape != () or values['number'].shape != ()):
                raise ValueError('completed contraction dimensions differ')
            density[:,slot.component] += slot.weight*values['density']
            gradient += slot.weight*values['gradient']
            ell += slot.weight*values['ell']
            one_body += slot.weight*values['one_body'].item()
            number += slot.weight*values['number'].item()
        m.reductions += 1
        return m.stage(density, gradient, ell, one_body, number, *coords)

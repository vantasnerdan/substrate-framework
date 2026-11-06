"""Bounded Linux process maps for independent sectors, never global stages.

Fork isolates lazy physical caches without copying the material background up
front. Pools are stage-local so closures inherit the current frozen stage.
Only small moments/bindings cross IPC; full successor frames are saved in place.
"""
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from itertools import islice
import json
import multiprocessing as mp
import os
from pathlib import Path
import threading

from .streamed_material_state import DiskMean, authenticate, save_array


COUNTERS = ('operator_constructions', 'h0_calls', 'h0_vector_columns',
            'moment_reductions', 'kernel_calls')
_FORK_FUNCTION = None


class WorkCounts:
    """One shared lock covers the cap test AND reservation before each action."""
    def __init__(self):
        self.values = mp.get_context('fork').Array('q', len(COUNTERS), lock=True)

    def add(self, name, amount=1):
        with self.values.get_lock():
            self.values[COUNTERS.index(name)] += amount

    def reserve_h0(self, columns, cap):
        with self.values.get_lock():
            index = COUNTERS.index('h0_calls')
            if cap is not None and self.values[index] >= cap:
                raise RuntimeError('full-PW H0 action cap exhausted; no successor published')
            self.values[index] += 1
            self.values[COUNTERS.index('h0_vector_columns')] += columns

    def snapshot(self):
        with self.values.get_lock():
            return dict(zip(COUNTERS, self.values[:]))


def _pin_worker(cpus, next_cpu):
    with next_cpu.get_lock():
        index = next_cpu.value
        next_cpu.value += 1
    os.sched_setaffinity(0, {cpus[index]})


def _invoke(value):
    return _FORK_FUNCTION(value)


def _publish_boundary(directory, slots, identity, time, u, p, evidence):
    """Publish a complete JSON inode atomically without overwriting a state."""
    staging = directory/'.boundary-publication'
    staging.mkdir()
    DiskMean._publish(staging, slots, identity, time, u, p, evidence)
    pending = staging/'boundary.json'
    with pending.open('rb') as stream:
        os.fsync(stream.fileno())
    os.link(pending, directory/'boundary.json')
    pending.unlink()
    staging.rmdir()


class SectorExecutor:
    """Explicit CPU allocation, ordered bounded batches, fail-closed publication.

    At most len(cpus) task results reside in the pool/parent. A batch barrier
    trades some tail utilization for a simple auditable bound. No executor is
    retained between stages; no updated potential is consumed before reduction.
    Native thread limits are process-wide and require an otherwise idle caller.
    """
    def __init__(self, cpus):
        self.cpus = tuple(cpus)
        if (not self.cpus or len(set(self.cpus)) != len(self.cpus)
                or any(type(cpu) is not int or cpu not in range(6, 11) for cpu in self.cpus)
                or not set(self.cpus) <= os.sched_getaffinity(0)):
            raise ValueError('explicit distinct available unclaimed CPUs in 6..10 required')
        self._active = False

    def map(self, function, values):
        global _FORK_FUNCTION
        from threadpoolctl import threadpool_limits
        if self._active or _FORK_FUNCTION is not None or threading.active_count() != 1:
            raise RuntimeError('sector fork requires a single-threaded, non-nested caller')
        context = mp.get_context('fork')
        names = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
        prior = {name: os.environ.get(name) for name in names}
        self._active, _FORK_FUNCTION = True, function
        try:
            for name in names:
                os.environ[name] = '1'
            with threadpool_limits(limits=1):
                with ProcessPoolExecutor(len(self.cpus), mp_context=context,
                        initializer=_pin_worker, initargs=(self.cpus, context.Value('i', 0))) as pool:
                    iterator = iter(values)
                    while True:
                        batch = tuple(islice(iterator, len(self.cpus)))
                        if not batch:
                            break
                        # map returns in input order regardless of finish order.
                        results = list(pool.map(_invoke, batch, chunksize=1))
                        yield from results
                        del results
        finally:
            self._active, _FORK_FUNCTION = False, None
            for name, value in prior.items():
                if value is None:
                    os.environ.pop(name, None)
                else:
                    os.environ[name] = value

    def successor(self, state, directory, frame, *, time, u, p, evidence, counts):
        """Workers write unique frames; only parent publishes a complete state."""
        import numpy as np
        directory = Path(directory).resolve()
        directory.mkdir(parents=True, exist_ok=False)
        def save(item):
            index, slot = item
            value = np.asarray(frame(slot))
            if (value.shape != slot.shape or value.dtype != np.complex128
                    or not np.isfinite(value).all()):
                raise ValueError('successor changed full-frame shape or contains nonfinite values')
            return dict(slot.metadata(), coefficient=slot.coefficient,
                        frame=save_array(directory/('frame-%04d.npy' % index), value))
        slots = list(self.map(save, enumerate(state.slots)))
        evidence['work_counts'] = counts()
        _publish_boundary(directory, slots, state.identity, time, u, p, evidence)
        return DiskMean(directory/'boundary.json', expected_identity=state.identity, verify=False)


PHYSICAL_KEYS = ('aggregate', 'manifest', 'preparation', 'initial_ions', 'scalar')
KNOWN_T0_SHA256 = '550391812d873a1cf10ed02576025c68c23d23600da5b9bd3c4ba2e649832869'
KNOWN_T0_IDENTITY = '9c79fe0f79d68c624727356e143c13e5339a80294ee2205aadaedd1392399257'


def physical_specification(specification):
    return {key: specification[key] for key in PHYSICAL_KEYS}


def physical_identity(specification):
    return sha256(json.dumps(physical_specification(specification), sort_keys=True).encode()).hexdigest()


def adopt_initial_boundary(boundary_binding, directory, *, specification, u, p,
                           expected_slots, execution, authority_binding=None):
    """Copy-free narrow t=0 migration. Production authority is fixed, not CLI-set.

    authority_binding is for independent finite-oracle fixtures; production
    callers MUST omit it. It must be a separately bound old model authority
    containing the exact initial boundary binding and its old specification.
    """
    import numpy as np
    boundary = authenticate(boundary_binding)
    state = DiskMean(boundary)
    if authority_binding is None:
        if (boundary_binding['sha256'] != KNOWN_T0_SHA256
                or state.identity != KNOWN_T0_IDENTITY):
            raise ValueError('only the known original t0 boundary may cross execution identities')
        old_specification = state.evidence['model']
    else:
        authority = json.loads(authenticate(authority_binding).read_text())
        if (authority['initial_boundary'] != boundary_binding
                or authority['identity'] != state.identity):
            raise ValueError('initial boundary differs from independently bound authority')
        old_specification = authority['specification']
    old_identity = sha256(json.dumps(old_specification, sort_keys=True).encode()).hexdigest()
    if (old_identity != state.identity or state.time != 0.
            or physical_specification(old_specification) != physical_specification(specification)
            or not np.array_equal(state.u, u) or not np.array_equal(state.p, p)
            or tuple(s.identity for s in state.slots) != tuple(expected_slots)
            or not state.evidence.get('original_raw_frames_unchanged')):
        raise ValueError('changed physical initial state cannot be adopted')
    # The exact bound boundary authenticates source/basis/gauge, frame and D.
    # No rewritten slot metadata, re-gauging, frame copying or D recomputation.
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=False)
    slots = [dict(s.metadata(), frame=s.frame, coefficient=s.coefficient) for s in state.slots]
    evidence = dict(state.evidence, model=physical_specification(specification),
                    adoption=dict(original_boundary=boundary_binding,
                                  original_model_identity=state.identity,
                                  method='exact authenticated t0 bindings; no array copies'),
                    execution=execution)
    identity = physical_identity(specification)
    _publish_boundary(directory, slots, identity, 0., state.u, state.p, evidence)
    return DiskMean(directory/'boundary.json', expected_identity=identity, verify=False)

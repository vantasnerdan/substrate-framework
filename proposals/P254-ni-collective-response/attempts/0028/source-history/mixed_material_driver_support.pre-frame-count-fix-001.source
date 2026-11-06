"""Source/input-bound local-pulse seeds, metadata accounting and read-only frames.

These helpers never infer raw jets from old D-weighted response factors. The
initial electronic variation is exact full-shape complex128 zero metadata.
Resource forecasts describe logical retention and aliases, not a storage-fit,
accuracy, whole-q, energy, or production admission. No helper launches physics.
"""
import ast
from contextlib import contextmanager
import fcntl
import json
import math
import os
from pathlib import Path
import struct
from uuid import uuid4
import sys
from time import monotonic, time

import numpy as np

from substrate_framework.partitioned_joint_slots import digest
from substrate_framework.partitioned_q_storage import QBudget, durable_publish, sync_directory
from substrate_framework.partitioned_q_work import QWork
from substrate_framework.streamed_material_state import authenticate, file_binding
from mixed_response_storage import (
    COMPONENTS, MappedMixedState, publish_checkpoint, publish_edge_barrier, retire_edge)
from local_pulse_energy import fourier_pulse
from mixed_payload_output import exact_zero_payload

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SEARCH = (HERE, HERE.parent/'0025', HERE.parent/'0024', ROOT/'src')


def canonical(value):
    return json.loads(json.dumps(value, sort_keys=True, allow_nan=False))


def checked_json(path, expected_sha256):
    binding = file_binding(Path(path).resolve())
    if binding['sha256'] != expected_sha256:
        raise ValueError('explicit input digest differs: '+str(path))
    return binding, json.loads(authenticate(binding).read_text())


def exact_publish(path, value, budget):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        binding = file_binding(path)
        if json.loads(authenticate(binding).read_text()) != canonical(value):
            raise ValueError('changed immutable owned request: '+str(path))
        return binding
    if budget is not None:
        return durable_publish(path, value, budget)
    temporary = path.with_name(path.name+'.'+uuid4().hex+'.incomplete')
    with temporary.open('xb') as stream:
        stream.write((json.dumps(value, indent=2, allow_nan=False)+'\n').encode())
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    sync_directory(path.parent)
    return file_binding(path)


def source_closure(seeds):
    """Recursively bind repository imports, including function-local imports.

    This is a static dependency walk, not a claim that loaded modules constitute
    the closure. The driver also refuses an unbound repository module at use.
    Explicit seeds retain the existing source-registration/admission convention.
    """
    bindings, pending = {}, [Path(b['path']).resolve() for b in seeds]
    def locate(name):
        if not name:
            return None
        for base in SEARCH:
            relative = Path(*name.split('.'))
            for candidate in (base/relative.with_suffix('.py'), base/relative/'__init__.py'):
                if candidate.is_file():
                    return candidate.resolve()
        return None
    while pending:
        path = pending.pop()
        if str(path) in bindings:
            continue
        bindings[str(path)] = file_binding(path)
        if path.suffix != '.py':
            continue
        tree = ast.parse(path.read_text(), filename=str(path))
        relative = path.relative_to(ROOT/'src') if path.is_relative_to(ROOT/'src') else None
        package = list(relative.parts[:-1]) if relative else []
        names = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names.extend(n.name for n in node.names)
            elif isinstance(node, ast.ImportFrom):
                prefix = package[:len(package)-node.level+1] if node.level else []
                module = '.'.join(prefix+([node.module] if node.module else []))
                names.append(module)
                names.extend('.'.join(filter(None, (module, n.name))) for n in node.names if n.name != '*')
            elif isinstance(node, ast.Call) and node.args:
                function = node.func
                dynamic = ((isinstance(function, ast.Name) and function.id == '__import__')
                           or (isinstance(function, ast.Attribute) and function.attr == 'import_module'))
                if dynamic and isinstance(node.args[0], ast.Constant) and isinstance(node.args[0].value, str):
                    names.append(node.args[0].value)
        for name in names:
            candidate = locate(name)
            if candidate is not None and str(candidate) not in bindings:
                pending.append(candidate)
            # Importing a submodule also executes each package initializer.
            bits = name.split('.')
            for end in range(1, len(bits)):
                candidate = locate('.'.join(bits[:end]))
                if candidate is not None and str(candidate) not in bindings:
                    pending.append(candidate)
    return [bindings[key] for key in sorted(bindings)]


def check_loaded_sources(bindings):
    expected = {b['path']: b for b in bindings}
    for module in tuple(sys.modules.values()):
        value = getattr(module, '__file__', None)
        if not value:
            continue
        path = Path(value).resolve()
        if path.suffix == '.py' and path.is_relative_to(ROOT):
            if str(path) not in expected:
                raise ValueError('repository runtime dependency absent from static closure: '+str(path))
    for binding in bindings:
        authenticate(binding)


def nested_bindings(*values):
    found = {}
    def visit(value):
        if isinstance(value, dict):
            if {'path', 'bytes', 'sha256'} <= value.keys():
                binding = dict(path=value['path'], bytes=int(value['bytes']), sha256=value['sha256'])
                old = found.get(binding['path'])
                if old is not None and old != binding:
                    raise ValueError('conflicting input file identity: '+binding['path'])
                found[binding['path']] = binding
            else:
                for item in value.values():
                    visit(item)
        elif isinstance(value, (list, tuple)):
            for item in value:
                visit(item)
    for value in values:
        visit(value)
    return [found[k] for k in sorted(found)]


def owned_tree(directory):
    """Logical path accounting: hardlinked aliases each remain charged."""
    result = {}
    if not Path(directory).exists():
        return result
    for path in Path(directory).rglob('*'):
        if path.is_symlink():
            raise ValueError('symlink in owned work tree: '+str(path))
        if path.is_file():
            result[str(path.resolve())] = path.stat().st_size
    return result


@contextmanager
def owned_directory(path, resume):
    path = Path(path)
    for part in (path, *path.parents):
        if part.is_symlink():
            raise ValueError('owned path cannot traverse a symlink')
    path = path.resolve()
    path.mkdir(parents=True, exist_ok=resume)
    if resume and not (path/'owner.json').is_file():
        raise ValueError('resume requires the original immutable owner request')
    with (path/'owner.lock').open('a+b') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield path


class ContinuationBudget(QBudget):
    """Original deadline and durable action admissions survive process loss.

    The append-only ledger is fsynced BEFORE a Hamiltonian action is admitted;
    a killed action remains charged, as in the original pre-action QBudget.
    Per-item native/QWork observed costs remain separate immutable receipts.
    Storage on restart is remeasured, never credited for retained aliases.
    """
    RECORD = struct.Struct('<QQ')

    def __init__(self, directory, *, started_at, **kwargs):
        self.ledger = Path(directory)/'action-admissions.bin'
        self.ledger_lock = Path(directory)/'action-admissions.lock'
        super().__init__(**kwargs)
        self.started = monotonic()-max(0., time()-started_at)
        self.started_at = started_at
        with self.ledger_lock.open('a+b') as lock:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
            with self.ledger.open('a+b') as stream:
                size = stream.seek(0, os.SEEK_END)
                if size % self.RECORD.size:
                    raise ValueError('incomplete durable action ledger; refusing an uncharged restart')
                if size:
                    stream.seek(-self.RECORD.size, os.SEEK_END)
                    actions, columns = self.RECORD.unpack(stream.read(self.RECORD.size))
                    self._counts[0], self._counts[1] = actions, columns

    def action(self, vectors):
        self.check()
        columns = 1 if vectors.ndim == 1 else vectors.shape[1]
        with self.ledger_lock.open('a+b') as lock:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
            with self.ledger.open('r+b') as stream:
                size = stream.seek(0, os.SEEK_END)
                if size % self.RECORD.size:
                    raise ValueError('damaged durable action ledger')
                previous = (0, 0)
                if size:
                    stream.seek(-self.RECORD.size, os.SEEK_END)
                    previous = self.RECORD.unpack(stream.read(self.RECORD.size))
                if previous[0] >= self.max_actions:
                    raise RuntimeError('original cumulative Hamiltonian action cap exhausted')
                counts = previous[0]+1, previous[1]+columns
                stream.seek(0, os.SEEK_END)
                stream.write(self.RECORD.pack(*counts))
                stream.flush()
                os.fsync(stream.fileno())
                with self._counts.get_lock():
                    self._counts[0], self._counts[1] = counts
        self.local_actions += 1
        self.local_columns += columns

    def snapshot(self):
        value = super().snapshot()
        value.update(original_started_unix=self.started_at,
                     deadline_unix=self.started_at+self.max_seconds,
                     action_accounting='durable pre-action admissions; killed/incomplete actions stay charged',
                     historical_worker_costs='original immutable QWork/transition receipts, not current execution')
        return value


def metadata_geometry(initial, q):
    slots = initial['slots']
    lookup = {s['identity']: i for i, s in enumerate(slots)}
    if len(lookup) != len(slots):
        raise ValueError('duplicate original explicit source slots')
    step = np.rint(4*np.asarray(q)).astype(int)
    targets = {}
    for component, sign in (('plus', 1), ('minus', -1)):
        target = []
        for slot in slots:
            spin, indices = slot['identity'].split(':')
            mesh = np.asarray([int(i) for i in indices.split(',')])
            key = spin+':'+','.join(str(int(i)) for i in (mesh+sign*step) % 4)
            if key not in lookup:
                raise ValueError('missing original physical target slot: '+key)
            target.append(lookup[key])
        targets[component] = target
    shapes = {c: [] for c in COMPONENTS}
    for i, slot in enumerate(slots):
        for c in ('plus', 'minus'):
            shapes[c].append([slots[targets[c][i]]['shape'][0], slot['shape'][1]])
        shapes['mixed'].append(list(slot['shape']))
    return shapes, targets['plus'], targets['minus']


def time_steps(history, final_time, subdivisions):
    if (type(subdivisions) is not int or subdivisions < 1 or not math.isfinite(final_time)
            or not 0 < final_time <= history['time_ry_inverse']):
        raise ValueError('positive selected time and subdivisions within completed actual history required')
    result, last = [], 0.
    for index, record in enumerate(history['history']):
        lo, span = record['start_time_ry_inverse'], record['dt_ry_inverse']
        if record['kind'] != 'actual-frozen-midpoint' or lo != last or span <= 0:
            raise ValueError('broken actual midpoint history coverage')
        last = lo+span
        hi = min(last, final_time)
        if hi > lo:
            grid = np.linspace(lo, hi, subdivisions+1)
            result.extend(dict(interval=index, start=float(a), end=float(b), dt=float(b-a))
                          for a, b in zip(grid[:-1], grid[1:]))
        if hi == final_time:
            break
    if not result or result[-1]['end'] != final_time:
        raise ValueError('selected time is not covered by the actual history')
    return result




def frame_requirements(history, steps, *, block_size, admissions):
    """Exact q-independent interior frame inputs and their last step consumers."""
    admitted = {digest(proof['expected_request']['input']) for _, proof in admissions}
    required = {}
    for index, step in enumerate(steps):
        interval = step['interval']
        record = history['history'][interval]
        lo = record['start_time_ry_inverse']
        hi = lo+record['dt_ry_inverse']
        for node in (step['start'], step['start']+step['dt']/2, step['end']):
            if node == lo or node == hi:
                continue
            identity = canonical(dict(time=float(node), interval=interval, stage=record,
                baseline_boundary=nested_bindings(record['prior_boundary_binding'])[0],
                model=history['model_identity'], origin=history['original_electronic_boundary'],
                block_size=block_size))
            key = digest(identity)
            if key not in required:
                required[key] = dict(input=identity, first_step=index, last_step=index,
                    consumer_steps=[], retained_admission=key in admitted)
            required[key]['last_step'] = index
            if index not in required[key]['consumer_steps']:
                required[key]['consumer_steps'].append(index)
    return required


def initial_small(initial, *, q, atom, axis, amplitude, cell, ncells=64):
    ndof, nfield = len(initial['p_canonical']), len(initial['field_coordinate'])
    if (not 1 <= atom <= ndof//3 or axis not in 'xyz' or not math.isfinite(amplitude)
            or len(cell) != 3 or any(type(i) is not int for i in cell) or ncells != 64):
        raise ValueError('original atom/axis, finite pulse and original four-cube cell required')
    gamma = bool(np.all(np.asarray(q) == 0))
    small = {c: dict(u=np.zeros(ndof, complex), p=np.zeros(ndof, complex),
                    coordinate=np.zeros(nfield if gamma or c == 'mixed' else 0, complex),
                    momentum=np.zeros(nfield if gamma or c == 'mixed' else 0, complex))
             for c in COMPONENTS}
    # The canonical executed function returns orthonormal collective DFT
    # coefficients. Its documented inverse contributes the second 1/sqrt(N),
    # giving the FIXED physical per-cell J*exp(-2*pi*i*q.cell)/N convention.
    mesh_index = tuple(int(i) for i in np.rint(4*np.asarray(q)).astype(int) % 4)
    row = next(row for row in fourier_pulse((4, 4, 4), cell, amplitude)
               if tuple(row['mesh_index']) == mesh_index)
    kick = row['coefficient']/math.sqrt(ncells)
    index = 3*(atom-1)+'xyz'.index(axis)
    small['plus']['p'][index] = kick
    small['minus']['p'][index] = kick.conjugate()
    return small


def initialize_raw_zero(directory, *, shapes, small, q, target_plus, target_minus,
                        identity, provenance, budget):
    """Recoverable exact initial variation without artificial zero payload files."""
    directory = Path(directory)
    directory.mkdir(exist_ok=True)
    request = dict(shapes=shapes, small={c: {k: [[z.real, z.imag] for z in a]
                    for k, a in values.items()} for c, values in small.items()},
                   q=list(q), target_plus=list(target_plus), target_minus=list(target_minus),
                   identity=identity, provenance=provenance)
    exact_publish(directory/'request.json', request, budget)
    legs = {c: [exact_zero_payload(tuple(shape)) for shape in shapes[c]] for c in COMPONENTS}
    checkpoint = directory/'checkpoint.json'
    if checkpoint.exists():
        binding = file_binding(checkpoint)
        with MappedMixedState(binding, budget=budget) as state:
            record = state.record
            if (record.get('identity') != identity
                    or record.get('provenance') != canonical(provenance)
                    or record.get('small') != canonical(request['small']) or record.get('time') != 0.
                    or record.get('shapes') != canonical(shapes) or record.get('q') != list(q)
                    or record.get('target_plus') != list(target_plus)
                    or record.get('target_minus') != list(target_minus)
                    or record.get('legs') != legs):
                raise ValueError('owned initial checkpoint differs from exact physical pulse request')
        return binding
    return publish_checkpoint(checkpoint, legs=legs, small=small, q=q,
        target_plus=target_plus, target_minus=target_minus, time=0., identity=identity,
        kind='physical-local-pulse-raw-electronic-zero', provenance=provenance, budget=budget)


def retire_previous_endpoint(directory, *, previous, successor, successor_receipt, budget):
    """Retire one obsolete full raw base after a durable actual next endpoint.

    Every edge uses the qualified storage successor barrier and write-ahead
    retirement journal. Only old payloads are unlinked; all checkpoints, small
    states, source identities, action/worker costs and journals survive. Replay
    never opens retired bases, and already-retired successors need not be live.
    The original exact-zero initial checkpoint is excluded by the caller.
    """
    directory = Path(directory).resolve()
    completed = json.loads(authenticate(successor_receipt).read_text())
    if (completed.get('final') != successor or completed.get('predecessor') != previous
            or completed.get('complete') is not True):
        raise ValueError('durable actual successor/readout commit required before cross-step retirement')
    for binding in completed.get('readouts', {}).values():
        if json.loads(authenticate(binding).read_text()).get('status') != 'complete':
            raise ValueError('requested readout is not complete; previous raw base must stay live')
    old = json.loads(authenticate(previous).read_text())
    new = json.loads(authenticate(successor).read_text())
    if (old.get('identity') != new.get('identity') or not 0 < old['time'] < new['time']
            or old['shapes'] != new['shapes'] or old['q'] != new['q']
            or old['target_plus'] != new['target_plus'] or old['target_minus'] != new['target_minus']):
        raise ValueError('cross-step retirement requires unchanged raw geometry and a later actual endpoint')
    root = directory/'endpoint-retirements'/previous['sha256']
    request = exact_publish(root/'request.json', dict(previous=previous, successor=successor,
        successor_receipt=successor_receipt, owner=str(directory)), budget)
    receipts = []
    for index in range(len(old['shapes']['mixed'])):
        path = root/f'edge-{index:06d}.json'
        successors = {'final': {c: new['legs'][c][index] for c in COMPONENTS}}
        if path.exists():
            barrier = file_binding(path)
            value = json.loads(authenticate(barrier).read_text())
            if (value.get('checkpoint') != previous or value.get('index') != index
                    or value.get('successors') != successors):
                raise ValueError('changed owned cross-step successor barrier')
        else:
            barrier = publish_edge_barrier(path, checkpoint=previous, index=index,
                                          successors=successors, budget=budget)
        receipts.append(retire_edge(previous, index, barrier, budget=budget, owner_root=directory))
    return exact_publish(root/'complete.json', dict(schema='P254-mixed-endpoint-retirement-v1',
        complete=True, request=request, previous=previous, successor=successor,
        successor_receipt=successor_receipt, edge_retirements=receipts), budget)



class RetainedFrameWork:
    """Read-only admitted frames plus genuine missing-only owned QWork frames.

    Original namespaces/payloads/costs never change. Missing interior times use
    the original history's real exponential callback, with the same budget and
    full source-bound owned identity. Owned frames retire after their final
    step/readout consumer commits; compact actual provenance/costs survive.
    """
    verify = QWork.verify

    def __init__(self, admissions, *, history_binding, budget):
        from actual_frame_cache_provenance import _complete_cache
        self.budget, self._verified, self.entries, self.last = budget, {}, {}, None
        self.provenance, self.borrowed = [], []
        self.generated, self.requirements = None, {}
        self.execution_log = []
        for proof_binding, proof in admissions:
            if (proof.get('schema') != 'P254-real-frame-consumer-admission-v1'
                    or proof.get('complete') is not True):
                raise ValueError('completed real-QWork full-PW frame admission required')
            for source in proof['implementation_sources']:
                authenticate(source)
            descriptor, items = _complete_cache(Path(proof['cache']), proof['owner'])
            original = descriptor['original_request']
            request = proof['expected_request']
            if (descriptor['history'] != history_binding or request['count'] != original['count']
                    or request['input'] != original['input'] or proof['original_request'] != original):
                raise ValueError('retained frames differ from explicit original history/input/count/producer')
            # The original producer and the admitted consumer identities are both
            # kept, even when the latter adopted the same physical frame bytes.
            frame_work = object.__new__(QWork)
            frame_work.directory = Path(proof['target']).resolve()
            frame_work.budget, frame_work.cpus = budget, ()
            frame_work.identity, frame_work._verified, frame_work.last = request['source_identity'], {}, None
            key = digest(request)
            root = frame_work.directory/key
            generations = [p for p in root.glob('generation-*') if not (p/'retire.json').exists()]
            if len(generations) != 1:
                raise ValueError('one live immutable admitted frame request generation required: '+str(root))
            generation = generations[0]
            if json.loads((generation/'request.json').read_text()) != request:
                raise ValueError('admitted actual QWork request bytes changed')
            records = frame_work.collect(generation, key, request['count'])
            for i, (record, (cached, origin)) in enumerate(zip(records, items)):
                a, b = record['arrays']['frame'], cached['frame']
                if (a['sha256'] != b['sha256'] or a['bytes'] != b['bytes']
                        or not os.path.samefile(a['path'], b['path'])):
                    raise ValueError('admitted frame is not the original immutable cached payload')
                shape = tuple(descriptor['slots'][i]['shape'])
                value = np.load(a['path'], mmap_mode='r', allow_pickle=False)
                try:
                    if value.shape != shape or value.dtype != np.complex128:
                        raise ValueError('retained original full-PW frame header differs')
                finally:
                    value._mmap.close()
                self.borrowed.extend(nested_bindings(cached, origin, record))
            identity = digest(request['input'])
            if identity in self.entries:
                raise ValueError('duplicate retained frame stage input')
            self.entries[identity] = (request, frame_work, generation, records)
            self.provenance.append(dict(proof=proof_binding, cache=descriptor,
                consumer_request=request, original_request=original,
                cache_complete=file_binding(Path(proof['cache'])/'cache-complete.json')))

    def configure_owned(self, work, requirements, source_count):
        if work.budget is not self.budget:
            raise ValueError('missing-frame QWork must share the original cumulative budget')
        self.generated, self.requirements = work, canonical(requirements)
        if type(source_count) is not int or source_count < 1:
            raise ValueError('original complete source count required for owned baseline frames')
        self.source_count = source_count
        for key, value in self.requirements.items():
            if (key != digest(value['input'])
                    or value['retained_admission'] != (key in self.entries)):
                raise ValueError('planned actual frame input/admission differs')

    def retire_after_commit(self, index, commit_binding, shared_records):
        """Retire ONLY owned times whose last requested consumer just committed."""
        result = []
        if self.generated is None:
            raise ValueError('source-bound owned missing-frame QWork was not configured')
        for value in self.requirements.values():
            if value['retained_admission'] or value['last_step'] != index:
                continue
            request = dict(schema='partitioned-q-work-v2', source_identity=self.generated.identity,
                operation='actual-baseline-frames', input=value['input'], count=self.source_count)
            root = self.generated.directory/digest(request)
            generations = list(root.glob('generation-*'))
            if not generations:
                raise ValueError('committed step lacks its actual missing-time frame provenance')
            for generation in generations:
                if json.loads((generation/'request.json').read_text()) != request:
                    raise ValueError('owned baseline-frame request changed before retirement')
                if not (generation/'retire.json').exists():
                    self.generated.collect(generation, digest(request), self.source_count)
                result.append(self.generated.retire(generation, commit_binding))
            identity = value['input']
            shared_records.pop((identity['interval'], identity['time']), None)
        return result

    def execute(self, operation, identity, count, function):
        self.budget.check()
        entry = self.entries.get(digest(identity))
        if operation != 'actual-baseline-frames':
            raise ValueError('baseline frame dispatcher only accepts genuine actual-baseline-frames')
        if entry is None:
            key = digest(identity)
            planned = self.requirements.get(key)
            if (self.generated is None or planned is None
                    or planned['input'] != canonical(identity) or count != self.source_count
                    or not callable(function)):
                raise ValueError('exact supported missing-time input and original exponential callback required')
            request_key, generation = self.generated.request(operation, identity, count)
            log = dict(request=file_binding(generation/'request.json'), input=canonical(identity),
                count=count, execution='original qualified frozen-mean exponential; missing retained time only',
                status='allocated', historical_costs='immutable original QWork worker/attempt receipts')
            self.execution_log.append(log)
            try:
                records = self.generated.execute(operation, identity, count, function)
            except BaseException as error:
                log.update(status='failed', error=type(error).__name__+': '+str(error))
                raise
            self.last = self.generated.last
            log.update(status='complete', worker_receipts=[file_binding(
                generation/('item-%04d' % i)/'complete.json') for i in range(count)])
            return records
        request, frame_work, generation, _ = entry
        if request['input'] != canonical(identity) or request['count'] != count:
            raise ValueError('retained frame exact input/count mismatch')
        records = frame_work.collect(generation, digest(request), count)
        self.last = digest(request), generation, count
        return records


def resource_plan(shapes, steps, borrowed, *, existing, workers, caps, history,
                  gamma, frame_requirements, retain_endpoints=False):
    """Exact raw-state/transition counts plus disclosed non-certified workspace."""
    edge = [sum(16*math.prod(shapes[c][i])+512 for c in COMPONENTS)
            for i in range(len(shapes['mixed']))]
    state = sum(edge)
    # The zero seed remains. Default rolling retention keeps one nonzero base;
    # current+accumulator and per-worker successor triples are separate. Before
    # unlinking a base, its next endpoint/readouts and edge journals are durable.
    transition = 2*sum(sorted(edge, reverse=True)[:workers])
    retained = len(steps)+1 if retain_endpoints else min(2, len(steps)+1)
    raw_peak = ((len(steps)+2) if retain_endpoints else min(4, len(steps)+2))*state+transition
    missing = [value for value in frame_requirements.values() if not value['retained_admission']]
    frame_set = sum(16*math.prod(shape)+512 for shape in shapes['mixed'])
    live_sets = len(missing) if retain_endpoints else max((sum(
        value['first_step'] <= index <= value['last_step'] for value in missing)
        for index in range(len(steps))), default=0)
    new_frame_peak = live_sets*frame_set
    density_binding = history['history'][0]['arrays']['density']
    density = np.load(density_binding['path'], mmap_mode='r', allow_pickle=False)
    try:
        density_shape, ngrid = tuple(density.shape), math.prod(density.shape)
    finally:
        density._mmap.close()
    ndof, nfield = len(history['p_canonical']), len(history['field_coordinate'])
    moment_per_source = 3*(16*(ngrid+ndof+3)+3*512)
    moment_peak = len(edge)*moment_per_source
    context = 16*(6*ngrid+9*ndof+9+2*nfield*(3 if gamma else 1))+24*512
    contexts = 4*len(steps)*context
    paths = dict(existing)
    for binding in borrowed:
        paths[binding['path']] = binding['bytes']
    logical = sum(paths.values())
    inodes, allocated, missing = set(), 0, []
    for path in paths:
        try:
            stat = Path(path).stat()
        except FileNotFoundError:
            # Retired producer aliases remain immutable origin provenance;
            # the admitted cache payload is authenticated by its real consumer.
            missing.append(path)
            continue
        inode = stat.st_dev, stat.st_ino
        if inode not in inodes:
            inodes.add(inode)
            allocated += stat.st_blocks*512
    return dict(schema='P254-mixed-material-resource-plan-v1', metadata_only=True,
        source_count=len(edge), raw_components=list(COMPONENTS), raw_shapes=shapes,
        raw_state_logical_bytes_with_header_allowance=state,
        simultaneous_base_current_accumulator_bytes=3*state,
        per_edge_triple_bytes_with_header_allowance=edge,
        concurrent_two_successor_transition_bytes=transition,
        retain_all_physical_endpoints=retain_endpoints,
        retained_physical_endpoints_including_initial=retained,
        cross_step_retirement='qualified per-edge barriers after actual next endpoint and requested readouts; '
            'checkpoints/small/source/cost metadata survive; resume uses latest actual endpoint',
        raw_payload_peak_upper_bound_bytes=raw_peak,
        density_response_shape=density_shape,
        dense_all_source_moment_barrier_peak_bytes=moment_peak,
        prepared_context_bytes_per_rk_stage=context,
        retained_prepared_context_bytes=contexts,
        baseline_frame_requirements=frame_requirements,
        new_missing_baseline_time_count=len(missing),
        new_baseline_frame_set_logical_bytes=frame_set,
        peak_owned_baseline_frame_sets=live_sets,
        peak_owned_baseline_frame_bytes=new_frame_peak,
        cumulative_new_baseline_frame_payload_bytes=len(missing)*frame_set,
        owned_frame_retention='explicit all-endpoint retention' if retain_endpoints else
            'last step/readout consumer only; actual QWork retirement journals preserve all observed costs',
        declared_readonly_and_existing_paths=len(paths),
        declared_readonly_and_existing_logical_alias_bytes=logical,
        current_readonly_and_existing_unique_inodes=len(inodes),
        current_readonly_and_existing_allocated_bytes=allocated,
        retired_declared_alias_paths=missing,
        existing_owned_output_logical_bytes=sum(existing.values()),
        owned_payload_context_and_moment_forecast_bytes=sum(existing.values())+raw_peak+contexts+moment_peak+new_frame_peak,
        selected_steps=steps, completed_baseline_time=history['time_ry_inverse'], caps=caps,
        owned_output_forecast_within_live_cap=sum(existing.values())+raw_peak+contexts+moment_peak+new_frame_peak <= caps['max_output_bytes'],
        workspace_accounting=dict(native_fft_projector_kernel_and_optional_energy='not certified by metadata',
            density_classical_stages='all four required; no density/PW/rank/precision reduction',
            rss='actual parent/worker watchdog enforces explicit caps; not inferred from raw payload'),
        storage_semantics='only owned live logical output is charged; readonly aliases and unique physical blocks are disclosed separately; '
            'initial holes do not imply dense successor storage fit; compact metadata/receipts additionally charged; '
            'moment/context upper bounds assume all dense nonzero responses',
        scope='resource accounting only; no storage-fit/accuracy/scientific gate or production authorization')

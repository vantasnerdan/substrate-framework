"""Serial, immutable additive-history accumulation; no physics evaluation.

accumulate(root, catalog_path) admits one explicitly q-owned component.
publish(root, output_catalog) emits a partial receipt (no six arrays) until the
literal original 31 x 31 x 64 ownership is complete. Complete publication uses
A62's Sigma/Pi convention. Storage zeros are never evidence of ownership.

Each transaction copies only its changed family. The atomic current.json holds
current and previous immutable state bindings; those two states own all retained
family objects. An interrupted precommit is unreferenced and collected on reopen.
Callers must not modify source files during admission. flock serializes writers
and publishers; no daemon, evaluator, interpolation, or hidden launch is used.
"""
from contextlib import contextmanager
from pathlib import Path
import fcntl
import hashlib
import json
import os
import shutil
import sys
import uuid

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / '0062'))
from history_banks import (ORDERS, SCHEMA, REFERENCE, PI_CONVENTION,
                           binding, authenticate, save)

STATE_SCHEMA = 'p254-serial-additive-bank-v1'
COMPONENTS = ('ee', 'mixed_ei', 'cubic', 'joint')
SIGMA_AXES = ['time', 'time', 'slot', 'band', 'band']
PI_AXES = ['time', 'time', 'coordinate', 'coordinate']
TAILS = {'sigma': (128, 12, 12), 'cubic': (111, 111), 'pi': (114, 114)}


def _sync_dir(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _save(path, record):
    save(path, record)
    _sync_dir(Path(path).parent)


@contextmanager
def _locked(root):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    with (root / 'serial.lock').open('a+b') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        yield root


def _reference():
    prereg = json.loads((HERE / 'preregistration.json').read_bytes())
    descriptor = prereg['reference']
    record = json.loads(authenticate(descriptor).read_bytes())
    if record['times_inverse_ry'] != list(range(31)):
        raise ValueError('original literal 31-node grid differs')
    # This original manifest fixes CAR slot/band bases and all coordinate routes.
    authenticate(record['coefficient_manifest'])
    authenticate(record['U'])
    return descriptor, record


def _load(root):
    path = root / 'current.json'
    if not path.exists():
        return None, None
    pointer = json.loads(path.read_bytes())
    state = json.loads(authenticate(pointer['current']).read_bytes())
    if state['schema'] != STATE_SCHEMA:
        raise ValueError('accumulator schema differs')
    return pointer, state


def reopen(root, *, require_complete=False):
    """Read committed state; optionally require every physical ownership cell."""
    with _locked(root) as directory:
        pointer, state = _load(directory)
        if state is None:
            raise ValueError('no committed accumulation')
        if require_complete and not all(item['complete'] for item in _summary(state).values()):
            raise ValueError('partial accumulated bank is ineligible for a complete solve')
        return dict(pointer=pointer, state=state)


def _collect(root, pointer):
    states, objects = set(), set()
    for key in ('current', 'previous'):
        descriptor = (pointer or {}).get(key)
        if descriptor:
            state_path = authenticate(descriptor)
            states.add(state_path.resolve())
            state = json.loads(state_path.read_bytes())
            objects.update(Path(f['directory']).resolve() for f in state['families'].values())
    for directory in (root / 'objects', root / 'states'):
        directory.mkdir(exist_ok=True)
    for path in (root / 'states').iterdir():
        if path.resolve() not in states:
            path.unlink()
    for path in (root / 'objects').iterdir():
        if path.resolve() not in objects:
            shutil.rmtree(path)


def _bindings(value, found):
    if isinstance(value, dict):
        # A paid archive is the input; its metadata records transitive producer
        # history, not current operands. Freeze that history via catalog SHA,
        # without reinterpreting obsolete live producer paths as source files.
        if {'archive', 'metadata'} <= value.keys():
            _bindings(value['archive'], found)
            return
        if {'archived_snapshot', 'original_binding'} <= value.keys():
            original, snapshot = value['original_binding'], value['archived_snapshot']
            if any(original[key] != snapshot[key] for key in ('bytes', 'sha256')):
                raise ValueError('historical snapshot does not preserve original bytes')
            _bindings(snapshot, found)
            return
        if {'path', 'bytes', 'sha256'} <= value.keys():
            path = value['path']
            old = found.get(path)
            descriptor = {key: value[key] for key in ('path', 'bytes', 'sha256')}
            if old is not None and old != descriptor:
                raise ValueError('conflicting source descriptors: ' + path)
            found[path] = descriptor
        for child in value.values():
            _bindings(child, found)
    elif isinstance(value, list):
        for child in value:
            _bindings(child, found)


def _indices(value, count, name):
    if not isinstance(value, list) or not value or any(type(x) is not int or not 0 <= x < count for x in value):
        raise ValueError('explicit original ' + name + ' indices required')
    if len(set(value)) != len(value):
        raise ValueError('duplicate ' + name)
    return value


def _normalize(catalog, reference):
    if catalog.get('schema') != 'p254-additive-history-component-v1':
        raise ValueError('component information contract differs')
    component = catalog['component']
    if component not in COMPONENTS:
        raise ValueError('unknown additive family')
    if catalog['input_dataset'] != reference:
        raise ValueError('changed original input_dataset identity')
    times = catalog['times_inverse_ry']
    if (not isinstance(times, list) or not times or
            any(type(t) not in (int, float) or not np.isfinite(t) or t not in range(31) for t in times) or
            any(b <= a for a, b in zip(times, times[1:]))):
        raise ValueError('literal increasing source-bank time subset required')
    grid = [int(t) for t in times]
    nt = len(times)
    coverage = catalog['coverage']
    pairs = coverage['completed_pairs']
    if pairs == 'all':
        pairs = [[a, b] for a in range(nt) for b in range(nt)]
    if not isinstance(pairs, list) or not pairs:
        raise ValueError('no persisted numeric ownership')
    if any(not isinstance(p, list) or len(p) != 2 or any(type(i) is not int or not 0 <= i < nt for i in p) for p in pairs):
        raise ValueError('invalid completed source-time pair')
    if len(set(map(tuple, pairs))) != len(pairs):
        raise ValueError('duplicate pair coverage')
    owned = coverage.get('owned_pairs', catalog.get('owned_pairs'))
    if owned != 'all' and (not isinstance(owned, list) or not set(map(tuple, pairs)) <= set(map(tuple, owned))):
        raise ValueError('completion outside declared ownership')
    if component == 'joint':
        if (catalog.get('q_additive') is not False or catalog.get('owned_q') or catalog.get('completed_q') or
                catalog.get('sigma_convention') != 'joint_minus_ionic_Gamma' or
                catalog.get('pi_convention') != 'full_joint_Gamma_once'):
            raise ValueError('joint must be once-only global A62 correction')
        qs = [0]  # storage singleton, never a physical momentum partition
        if catalog['axes'] != dict(sigma=SIGMA_AXES, pi=PI_AXES):
            raise ValueError('joint axes differ')
        if coverage.get('all114_coordinates') is not True:
            raise ValueError('explicit full114 joint coordinate ownership required')
        slots = _indices(coverage['external_targets'], 128, 'slot')
        selectors = [('sigma', s, s) for s in slots] + [('pi', s, 128+s) for s in range(114)]
        kinds = {'sigma': catalog['sigma'], 'pi': catalog['pi']}
    else:
        qs = _indices(catalog['owned_q'], 64, 'q')
        # One component array is one q contribution, not an unidentified sum.
        if len(qs) != 1 or catalog.get('completed_q', qs) != qs:
            raise ValueError('single explicit persisted original q contribution required')
        if component == 'cubic':
            if catalog['axes'] not in (['time', 'time', 'ionic_coordinate', 'ionic_coordinate'],
                                        ['time', 'time', 'ionic', 'ionic']):
                raise ValueError('cubic111 coordinate order differs')
            # Cubic is full111 by its declared canonical matrix shape. Historical
            # emitters use external_targets=[0] only as a pair-commit marker.
            slots = list(range(111))
            selectors = [('cubic', s, s) for s in slots]
            kinds = {'cubic': catalog['arrays']}
        else:
            if catalog['axes'] != SIGMA_AXES:
                raise ValueError('CAR slot/band axes differ')
            if catalog.get('slots', 128) != 128 or catalog.get('bands', 12) != 12 or catalog.get('spins', 2) != 2:
                raise ValueError('original128/12/bothspins required')
            slots = coverage.get('external_targets', coverage.get('completed_slots', catalog.get('external_slots')))
            if slots is None:
                raise ValueError('explicit completed external slots required')
            slots = _indices(slots, 128, 'slot')
            selectors = [('sigma', s, s) for s in slots]
            kinds = {'sigma': catalog['sigma']}
    for kind, descriptors in kinds.items():
        if set(descriptors) != set(ORDERS):
            raise ValueError('three contour orders required')
        for descriptor in descriptors.values():
            if descriptor['shape'] != [nt, nt, *TAILS[kind]] or descriptor['dtype'] != 'complex128':
                raise ValueError('canonical component shape/dtype differs')
    return component, qs[0], grid, pairs, selectors, kinds


def _open(descriptor, shape=None, dtype=None, mode='r'):
    array = np.load(authenticate(descriptor), mmap_mode=mode, allow_pickle=False)
    expected_shape = shape if shape is not None else descriptor['shape']
    expected_dtype = dtype if dtype is not None else descriptor['dtype']
    if tuple(array.shape) != tuple(expected_shape) or array.dtype != np.dtype(expected_dtype):
        array._mmap.close()
        raise ValueError('stored NPY axes/dtype differ')
    return array


def _close(maps):
    for array in maps:
        array._mmap.close()


def _digest(values):
    digest = hashlib.sha256()
    for value in values:
        if not np.isfinite(value).all():
            raise ValueError('nonfinite admitted component')
        digest.update(np.ascontiguousarray(value).view(np.uint8))
    return digest.digest()


def _family_shape(component):
    return (1 if component == 'joint' else 64, 31, 31,
            111 if component == 'cubic' else 242 if component == 'joint' else 128)


def _new_map(path, shape, dtype):
    # New/truncated regular files have zero-filled unwritten payloads. These
    # are storage-only zeros; ownership is an independent authenticated map.
    array = np.lib.format.open_memmap(path, mode='w+', dtype=dtype, shape=shape)
    array.flush()
    array._mmap.close()


def _persist(path):
    with Path(path).open('rb') as handle:
        os.fsync(handle.fileno())


def _summary(state):
    result = {}
    for component in COMPONENTS:
        family = state['families'].get(component)
        count = int(np.prod(_family_shape(component)))
        if family is None:
            result[component] = dict(covered_cells=0, required_cells=count, complete=False, completed_q=[])
            continue
        owned = _open(family['coverage'], _family_shape(component), 'uint8')
        try:
            qdone, covered = [], 0
            for q in range(owned.shape[0]):
                block = owned[q]
                if not np.isin(block, [0, 1]).all():
                    raise ValueError('noncanonical coverage map')
                covered += int(np.count_nonzero(block))
                if bool(block.all()):
                    qdone.append(q)
            result[component] = dict(covered_cells=covered, required_cells=count,
                                     complete=covered == count, completed_q=[] if component == 'joint' else qdone)
        finally:
            owned._mmap.close()
    return result


def accumulate(root, catalog_path):
    """Atomically admit actual component cells; identical overlap adds nothing.

    Changed numeric overlap or an already-bound source path changing bytes is
    rejected before commit. Catalog/array/source bindings remain in provenance;
    cell digests cover all three contour orders, individually per external slot.
    """
    with _locked(root) as root:
        reference, reference_record = _reference()
        catalog_binding = binding(catalog_path)
        catalog = json.loads(Path(catalog_path).read_bytes())
        component, q, grid, pairs, selectors, kinds = _normalize(catalog, reference)
        pointer, state = _load(root)
        if state is None:
            state = dict(schema=STATE_SCHEMA, generation=0, input_dataset=reference,
                         original_B_identity=reference_record['coefficient_manifest'],
                         times_inverse_ry=reference_record['times_inverse_ry'], families={}, sources={}, imports=[])
        if state['input_dataset'] != reference or state['original_B_identity'] != reference_record['coefficient_manifest']:
            raise ValueError('accumulator original source/B identity differs')
        _collect(root, pointer)
        sources = {catalog_binding['path']: catalog_binding}
        _bindings(catalog, sources)
        # Rehash actual bound operands/archives, deduplicated by absolute path.
        # Transitive archive metadata stays unchanged in the bound input catalog.
        for path, descriptor in sources.items():
            if path in state['sources'] and descriptor != state['sources'][path]:
                raise ValueError('changed previously admitted source: ' + path)
            authenticate(descriptor)
        source_maps, old_maps = {}, []
        family = state['families'].get(component)
        coverage = digests = None
        additions, duplicates = [], 0
        try:
            for kind, descriptors in kinds.items():
                source_maps[kind] = {}
                for order in ORDERS:
                    source_maps[kind][order] = _open(descriptors[order])
            if family:
                for bank in family['arrays'].values():
                    for descriptor in bank.values():
                        authenticate(descriptor)
                coverage = _open(family['coverage'], _family_shape(component), 'uint8')
                old_maps.append(coverage)
                digests = _open(family['digests'], _family_shape(component), 'V32')
                old_maps.append(digests)
            for a, b in pairs:
                i, j = grid[a], grid[b]
                for kind, slot, external in selectors:
                    digest = _digest([source_maps[kind][order][a, b, slot] for order in ORDERS])
                    key = (q, i, j, external)
                    if coverage is not None and coverage[key]:
                        if digests[key].tobytes() != digest:
                            raise ValueError('changed overlapping numeric cell: ' + str((component, key)))
                        duplicates += 1
                    else:
                        additions.append((a, b, kind, slot, external, digest))
            if not additions:
                return dict(status='idempotent', component=component, added_cells=0,
                            identical_overlap_cells=duplicates, pointer=pointer, catalog=catalog_binding)
            token = uuid.uuid4().hex
            target = root / 'objects' / (component + '-' + token)
            target.mkdir()
            if family:
                for path in Path(family['directory']).iterdir():
                    shutil.copyfile(path, target / path.name)
            else:
                _new_map(target / 'coverage.npy', _family_shape(component), 'uint8')
                _new_map(target / 'digests.npy', _family_shape(component), 'V32')
                for kind in kinds:
                    for order in ORDERS:
                        _new_map(target / (kind + '-' + order + '.npy'), (31, 31, *TAILS[kind]), 'complex128')
            _close(old_maps)
            old_maps.clear()
            # Bound mappings to one input family and one new family only.
            for kind in kinds:
                for order in ORDERS:
                    path = target / (kind + '-' + order + '.npy')
                    output = np.load(path, mmap_mode='r+', allow_pickle=False)
                    try:
                        for a, b, cell_kind, slot, external, digest in additions:
                            if cell_kind == kind:
                                output[grid[a], grid[b], slot] += source_maps[kind][order][a, b, slot]
                        output.flush()
                    finally:
                        output._mmap.close()
                    _persist(path)
            for name in ('coverage', 'digests'):
                path = target / (name + '.npy')
                output = np.load(path, mmap_mode='r+', allow_pickle=False)
                try:
                    for a, b, kind, slot, external, digest in additions:
                        output[q, grid[a], grid[b], external] = 1 if name == 'coverage' else np.void(digest)
                    output.flush()
                finally:
                    output._mmap.close()
                _persist(path)
            new_family = dict(directory=str(target), arrays={})
            for kind in kinds:
                new_family['arrays'][kind] = {}
                for order in ORDERS:
                    path = target / (kind + '-' + order + '.npy')
                    new_family['arrays'][kind][order] = dict(binding(path), shape=[31, 31, *TAILS[kind]], dtype='complex128')
            for name, dtype in (('coverage', 'uint8'), ('digests', 'V32')):
                new_family[name] = dict(binding(target / (name + '.npy')), shape=list(_family_shape(component)), dtype=dtype)
            _sync_dir(target)
            _sync_dir(target.parent)
            state['families'][component] = new_family
            state['generation'] += 1
            state['sources'].update(sources)
            state['imports'].append(dict(catalog=catalog_binding, component=component, q=None if component == 'joint' else q,
                added_cells=len(additions), identical_overlap_cells=duplicates))
            state_path = root / 'states' / (token + '.json')
            _save(state_path, state)
            new_pointer = dict(current=binding(state_path), previous=pointer['current'] if pointer else None)
            _save(root / 'current.json', new_pointer)
            _collect(root, new_pointer)
            return dict(status='committed', component=component, added_cells=len(additions),
                        identical_overlap_cells=duplicates, pointer=new_pointer, coverage=_summary(state), catalog=catalog_binding)
        finally:
            _close(old_maps)
            _close([array for bank in source_maps.values() for array in bank.values()])


def publish(root, output_catalog):
    """Publish only a partial receipt until every actual cell exists.

    A complete publication streams one contour order/family at a time, with the
    literal A62 combination and no normalization or represented-control change.
    The output catalog is atomically committed only after all six arrays persist.
    """
    with _locked(root) as root:
        pointer, state = _load(root)
        if state is None:
            raise ValueError('no committed accumulation')
        reference, reference_record = _reference()
        if state['input_dataset'] != reference or state['original_B_identity'] != reference_record['coefficient_manifest']:
            raise ValueError('changed source identity at publication')
        summary = _summary(state)
        complete = all(item['complete'] for item in summary.values())
        path = Path(output_catalog).resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        record = dict(schema='p254-accumulated-history-partial-v1', status='partial',
                      full_solve_eligible=False, input_dataset=state['input_dataset'],
                      original_B_identity=state['original_B_identity'],
                      times_inverse_ry=state['times_inverse_ry'], coverage=summary,
                      committed_pointer=pointer, families=state['families'],
                      missing_families=[k for k, v in summary.items() if not v['complete']],
                      storage_zeros='not physical values; not complete-solve eligible')
        if not complete:
            _save(path, record)
            return record
        destination = path.parent / (path.stem + '-arrays-' + uuid.uuid4().hex)
        destination.mkdir()
        arrays = {}
        for kind, tail in (('sigma', TAILS['sigma']), ('pi', TAILS['pi'])):
            arrays[kind] = {}
            for order in ORDERS:
                output_path = destination / (kind + '-' + order + '.npy')
                _new_map(output_path, (31, 31, *tail), 'complex128')
                output = np.load(output_path, mmap_mode='r+', allow_pickle=False)
                try:
                    terms = [('ee', 'sigma'), ('mixed_ei', 'sigma'), ('joint', 'sigma')] if kind == 'sigma' else [('joint', 'pi'), ('cubic', 'cubic')]
                    for component, source_kind in terms:
                        source = _open(state['families'][component]['arrays'][source_kind][order])
                        try:
                            for a in range(31):
                                for b in range(31):
                                    if source_kind == 'cubic':
                                        output[a, b, :111, :111] += source[a, b]
                                    else:
                                        output[a, b] += source[a, b]
                        finally:
                            source._mmap.close()
                    output.flush()
                finally:
                    output._mmap.close()
                _persist(output_path)
                arrays[kind][order] = dict(binding(output_path), shape=[31, 31, *tail], dtype='complex128')
        _sync_dir(destination)
        _sync_dir(destination.parent)
        coverage = {name: dict(owned_pairs='all', completed_pairs='all') for name in COMPONENTS}
        for name in ('ee', 'mixed_ei', 'cubic'):
            coverage[name].update(owned_q=list(range(64)), completed_q=list(range(64)))
        record.update(schema=SCHEMA, status='complete', full_solve_eligible=True,
                      reference_source_kind=REFERENCE, pi_convention=PI_CONVENTION,
                      coverage=coverage, exact_cell_coverage=summary, missing_families=[],
                      slots=128, bands=12, gamma_coordinates=114, global_q=64,
                      axes=dict(sigma=SIGMA_AXES, pi=PI_AXES), **arrays)
        _save(path, record)
        return record

"""Immutable COW storage for complete dressed functional contributions.

Only storage arithmetic lives here. Numeric files are fsynced before the state
and current/previous pointer; a failed transaction never owns scientific cells.
The A77 locking, authentication and low-level NPY primitives remain unchanged.
"""
from contextlib import contextmanager
from pathlib import Path
import hashlib
import json
import shutil
import sys
import uuid

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / n) for n in ('0077', '0062')]
from component_accumulator import (_locked as _a77_locked, _save, _open, _new_map,
                                   _persist, _sync_dir, _collect, _bindings, _digest)
from history_banks import ORDERS, binding, authenticate

SCHEMA = 'p254-dressed-functional-bank-v1'
COMPONENTS = ('ee', 'mixed_ei', 'joint', 'cubic')
TAILS = {'sigma': (128, 12, 12), 'pi': (114, 114), 'cubic': (111, 111)}


def load(root):
    path = Path(root) / 'current.json'
    if not path.exists():
        return None, None
    pointer = json.loads(path.read_bytes())
    state = json.loads(authenticate(pointer['current']).read_bytes())
    if state['schema'] != SCHEMA:
        raise ValueError('not a fresh dressed functional bank; old banks are immutable')
    identity_path = Path(root) / 'identity.json'
    if json.loads(identity_path.read_bytes()) != state['identity']:
        raise ValueError('bank identity marker differs')
    return pointer, state


@contextmanager
def _locked(root):
    directory = Path(root).resolve()
    # Refuse old/unrelated storage before A77 opens or creates its lock file.
    if directory.exists():
        if not directory.is_dir():
            raise ValueError('bank root must be a fresh directory')
        if (directory / 'current.json').exists():
            load(directory)
        elif not (directory / 'identity.json').exists() and any(
                p.name != 'serial.lock' for p in directory.iterdir()):
            raise ValueError('refusing any write to old or unrelated storage')
    with _a77_locked(directory) as locked:
        yield locked


def initialize(root, identity):
    marker = root / 'identity.json'
    if marker.exists():
        if json.loads(marker.read_bytes()) != identity:
            raise ValueError('fresh bank is already bound to another actual source')
    else:
        if any(p.name != 'serial.lock' for p in root.iterdir()):
            raise ValueError('refusing to initialize over an existing non-bank directory')
        _save(marker, identity)
    return dict(schema=SCHEMA, identity=identity, generation=0, families={},
                sources={}, imports=[])


def summary(state):
    result = {}
    for component in COMPONENTS:
        family = state['families'].get(component, {})
        required = 64 if component in ('ee', 'mixed_ei') else 1
        entries = family.get('entries', {})
        expected = {str(q) for q in range(required)}
        if not set(entries) <= expected:
            raise ValueError('invalid committed component ownership')
        complete = set(entries) == expected
        result[component] = dict(complete=complete,
            completed_contributions=len(entries), required_contributions=required,
            completed_pairs=961 if entries else 0,
            completed_q=sorted(map(int, entries)) if required == 64 else
                (list(range(64)) if component == 'cubic' and complete else []),
            all128_slots=bool(entries) if component != 'cubic' else None,
            all111_coordinates=bool(entries) if component == 'cubic' else None,
            all114_coordinates=bool(entries) if component == 'joint' else None)
    return result


def reopen(root, *, require_complete=False):
    with _locked(root) as directory:
        pointer, state = load(directory)
        if state is None:
            raise ValueError('no committed dressed accumulation')
        coverage = summary(state)
        if require_complete and not all(v['complete'] for v in coverage.values()):
            raise ValueError('partial dressed bank is ineligible for a complete solve')
        return dict(pointer=pointer, state=state, coverage=coverage)


def numeric_digest(maps):
    digest = hashlib.sha256()
    for kind in sorted(maps):
        digest.update(kind.encode('ascii'))
        for a in range(31):
            for b in range(31):
                digest.update(_digest([maps[kind][order][a, b] for order in ORDERS]))
    return digest.hexdigest()


def commit(root, catalog_path, catalog, identity, component, ownership, kinds):
    """Add one complete q, one global joint, or one complete cubic family.

    Whole-contribution numeric digests are sufficient because admission permits
    no partial q/pair/slot scopes. Duplicate numeric ownership is never re-added.
    """
    with _locked(root) as root:
        pointer, state = load(root)
        if state is None:
            state = initialize(root, identity)
        if state['identity'] != identity:
            raise ValueError('actual source/original preparation identity differs')
        catalog_binding = binding(catalog_path)
        if json.loads(authenticate(catalog_binding).read_bytes()) != catalog:
            raise ValueError('component catalog changed during admission')
        sources = {catalog_binding['path']: catalog_binding}
        _bindings(catalog, sources)
        for path, descriptor in sources.items():
            if path in state['sources'] and state['sources'][path] != descriptor:
                raise ValueError('changed previously admitted source: ' + path)
            authenticate(descriptor)
        _collect(root, pointer)
        maps = {}
        try:
            for kind, descriptors in kinds.items():
                maps[kind] = {}
                for order in ORDERS:
                    maps[kind][order] = _open(descriptors[order],
                        (31, 31, *TAILS[kind]), 'complex128')
            digest = numeric_digest(maps)
            family = state['families'].get(component)
            key = str(ownership)
            if family and key in family['entries']:
                if family['entries'][key]['numeric_sha256'] != digest:
                    raise ValueError('changed overlapping complete contribution: ' + component + '/' + key)
                return dict(status='idempotent', component=component,
                    added_contributions=0, catalog=catalog_binding, pointer=pointer,
                    coverage=summary(state))
            token = uuid.uuid4().hex
            target = root / 'objects' / (component + '-' + token)
            target.mkdir()
            if family:
                for descriptors in family['arrays'].values():
                    for descriptor in descriptors.values():
                        authenticate(descriptor)
                for path in Path(family['directory']).iterdir():
                    shutil.copyfile(path, target / path.name)
            arrays = {}
            for kind in kinds:
                arrays[kind] = {}
                for order in ORDERS:
                    path = target / (kind + '-' + order + '.npy')
                    if family is None:
                        _new_map(path, (31, 31, *TAILS[kind]), 'complex128')
                    output = np.load(path, mmap_mode='r+', allow_pickle=False)
                    try:
                        for a in range(31):
                            for b in range(31):
                                output[a, b] += maps[kind][order][a, b]
                        output.flush()
                    finally:
                        output._mmap.close()
                    _persist(path)
                    arrays[kind][order] = dict(binding(path),
                        shape=[31, 31, *TAILS[kind]], dtype='complex128')
            entries = dict(family['entries']) if family else {}
            entries[key] = dict(numeric_sha256=digest, catalog=catalog_binding)
            state['families'][component] = dict(directory=str(target), arrays=arrays,
                                                 entries=entries)
            state['generation'] += 1
            state['sources'].update(sources)
            state['imports'].append(dict(catalog=catalog_binding, component=component,
                owned_q=ownership if component in ('ee', 'mixed_ei') else None,
                numeric_sha256=digest))
            _sync_dir(target)
            _sync_dir(target.parent)
            state_path = root / 'states' / (token + '.json')
            _save(state_path, state)
            new_pointer = dict(current=binding(state_path),
                               previous=pointer['current'] if pointer else None)
            _save(root / 'current.json', new_pointer)
            _collect(root, new_pointer)
            return dict(status='committed', component=component, added_contributions=1,
                catalog=catalog_binding, pointer=new_pointer, coverage=summary(state))
        finally:
            for bank in maps.values():
                for array in bank.values():
                    array._mmap.close()

"""Compose complete actual-G1 contractions without evaluating physical terms.

The EE and raw EI families each own 64 complete q contributions. Joint is a
once-global Sigma correction/full Pi; cubic is one inherited/replaced all64
family. Partial state is a receipt, never a six-array solve input.
"""
from pathlib import Path
import argparse
import json
import sys
import uuid

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / n) for n in ('0076', '0075', '0062')]
from history_banks import (ORDERS, SCHEMA, REFERENCE, PI_CONVENTION,
                           binding, authenticate)
from component_consumer import consume_component
from component_reader import JointComponent
from bank_storage import (TAILS, COMPONENTS, _locked, _save, _open, _new_map,
                          _persist, _sync_dir, load, summary, commit, reopen)

SIGMA_AXES = ['time', 'time', 'slot', 'band', 'band']
PI_AXES = ['time', 'time', 'coordinate', 'coordinate']
ALL_PAIRS = {(a, b) for a in range(31) for b in range(31)}


def complete_pairs(value):
    if value == 'all':
        return True
    return (isinstance(value, list) and len(value) == 961 and
            all(isinstance(p, list) and len(p) == 2 and
                all(type(i) is int for i in p) for p in value) and
            set(map(tuple, value)) == ALL_PAIRS)


def source_identity(provider=None, preparation=None):
    settled = json.loads((HERE / 'preregistration.json').read_bytes())
    actual = binding(provider) if provider is not None else settled['actual_G1']
    origin = binding(preparation) if preparation is not None else settled['original_preparation']
    record = json.loads(authenticate(actual).read_bytes())
    if (record['schema'] != 'p254-dressed-history-provider-v1' or
            record['status'] != 'complete' or
            record['times_inverse_ry'] != list(range(31))):
        raise ValueError('actual complete full31 dressed provider required')
    numerical = record['input_identity']['numerical_reference']
    original = json.loads(authenticate(numerical).read_bytes())
    bare = json.loads(authenticate(origin).read_bytes())
    if (original['schema'] != 'p254-a52-requested-time-reference-v1' or
            original['times_inverse_ry'] != list(range(31)) or
            bare['schema'] != 'p254-prescribed-bare-history-source-v1' or
            bare['status'] != 'complete' or bare['reference_bundle'] != numerical or
            bare['times_inverse_ry'] != list(range(31))):
        raise ValueError('original numerical reference/preparation differs')
    cursor, visited = actual, set()
    while cursor != origin:
        if cursor['path'] in visited:
            raise ValueError('cyclic provider ancestry')
        visited.add(cursor['path'])
        node = json.loads(authenticate(cursor).read_bytes())
        if (node['schema'] != 'p254-dressed-history-provider-v1' or
                node['status'] != 'complete' or
                node['input_identity']['numerical_reference'] != numerical):
            raise ValueError('actual provider does not descend from original preparation')
        cursor = node['input_previous']
    authenticate(cursor)
    authenticate(original['coefficient_manifest'])
    authenticate(original['U'])
    return dict(input_dataset=actual, numerical_reference=numerical,
                preparation_origin=origin,
                original_B_identity=original['coefficient_manifest'],
                times_inverse_ry=list(range(31)))


def _orders(descriptors, tail):
    if set(descriptors) != set(ORDERS):
        raise ValueError('all three contour orders required')
    for descriptor in descriptors.values():
        if descriptor['shape'] != [31, 31, *tail] or descriptor['dtype'] != 'complex128':
            raise ValueError('complete original shape/precision required')


def _admit(catalog_path, catalog, identity):
    if (catalog['input_dataset'] != identity['input_dataset'] or
            catalog['times_inverse_ry'] != identity['times_inverse_ry']):
        raise ValueError('actual dressed input/grid identity differs')
    for key in ('numerical_reference', 'preparation_origin'):
        if key in catalog and catalog[key] != identity[key]:
            raise ValueError('separate original identity differs: ' + key)
    component = catalog.get('component')
    coverage = catalog['coverage']
    if catalog.get('schema') == 'p254-dressed-cubic-family-v1':
        if (catalog['status'] != 'complete' or
                catalog.get('preparation_origin') != identity['preparation_origin'] or
                catalog.get('numerical_reference') != identity['numerical_reference'] or
                catalog['inherited_q'] != list(range(1, 64)) or
                catalog['replaced_q'] != [0] or
                coverage['owned_q'] != list(range(64)) or
                coverage['completed_q'] != list(range(64)) or
                coverage.get('all111_coordinates') is not True or
                not complete_pairs(coverage['owned_pairs']) or
                not complete_pairs(coverage['completed_pairs'])):
            raise ValueError('complete truthful inherited63/replaced0 cubic family required')
        if catalog.get('owned_q') == [0]:
            raise ValueError('all64 cubic family cannot masquerade as a single-q catalog')
        provenance = catalog['provenance']
        if not {'base', 'old_q0', 'new_q0'} <= set(provenance):
            raise ValueError('explicit retained-base/old0/new0 provenance required')
        for name in ('base', 'old_q0', 'new_q0'):
            manifest = json.loads(authenticate(provenance[name]).read_bytes())
            if (manifest['status'] != 'complete' or
                    manifest['preparation_origin'] != identity['preparation_origin'] or
                    not complete_pairs(manifest['coverage']['owned_pairs']) or
                    not complete_pairs(manifest['coverage']['completed_pairs'])):
                raise ValueError('incomplete or changed cubic provenance: ' + name)
            if name == 'base':
                if (manifest['schema'] != 'p254-retained-cubic-base-v1' or
                        manifest['original_input_dataset'] != identity['numerical_reference'] or
                        manifest['coverage']['owned_q'] != list(range(64)) or
                        manifest['coverage']['completed_q'] != list(range(64))):
                    raise ValueError('retained base must be the original full64 family')
                authenticate(manifest['base_catalog'])
            else:
                expected = identity['numerical_reference'] if name == 'old_q0' else identity['input_dataset']
                if (manifest['schema'] != 'p254-cubic-q0-cache-v1' or
                        manifest['input_dataset'] != expected or
                        manifest['numerical_reference'] != identity['numerical_reference'] or
                        manifest['owned_q'] != [0] or manifest['completed_q'] != [0] or
                        manifest['coverage']['owned_q'] != [0] or
                        manifest['coverage']['completed_q'] != [0] or
                        any(manifest['route'][key] != 0
                            for key in ('external_q', 'internal_q', 'partner_q'))):
                    raise ValueError('truthful complete original/dressed q0 provenance required')
            _orders(manifest['arrays'], TAILS['cubic'])
            for descriptor in manifest['arrays'].values():
                authenticate(descriptor)
        _orders(catalog['arrays'], TAILS['cubic'])
        return 'cubic', 0, {'cubic': catalog['arrays']}
    if catalog.get('schema') != 'p254-additive-history-component-v1':
        raise ValueError('unknown component schema')
    if component in ('ee', 'mixed_ei'):
        q, = catalog['owned_q']
        if (type(q) is not int or q not in range(64) or
                catalog['status'] != 'complete-owned-component-q%d' % q or
                catalog['completed_q'] != [q] or catalog['full_mesh'] is not False or
                catalog['full_time_square'] is not True or
                catalog['axes'] != SIGMA_AXES or
                (catalog['slots'], catalog['bands'], catalog['spins']) != (128, 12, 2) or
                not complete_pairs(coverage['owned_pairs']) or
                not complete_pairs(coverage['completed_pairs']) or
                coverage['completed_slot_cells'] != 961 * 128):
            raise ValueError('only complete full31/all128 single-q Sigma is admissible')
        _orders(catalog['sigma'], TAILS['sigma'])
        # The unchanged public consumer authenticates routes, pair/slot bitmap,
        # inputs, contour identity and every covered numeric value.
        consume_component(catalog_path)
        return component, q, {'sigma': catalog['sigma']}
    if component == 'joint':
        if (catalog.get('owned_q') or catalog.get('completed_q') or
                coverage.get('all114_coordinates') is not True or
                coverage.get('external_targets') != list(range(128)) or
                not complete_pairs(coverage['owned_pairs']) or
                not complete_pairs(coverage['completed_pairs'])):
            raise ValueError('once-global all961/all128/full114 joint ownership required')
        with JointComponent(catalog_path, require_complete=True):
            pass
        return 'joint', 0, {'sigma': catalog['sigma'], 'pi': catalog['pi']}
    raise ValueError('single-q cubic or unknown component is not a dressed family')


def accumulate(root, catalog_path, *, provider=None, preparation=None):
    identity = source_identity(provider, preparation)
    catalog_path = Path(catalog_path).resolve()
    catalog = json.loads(catalog_path.read_bytes())
    component, ownership, kinds = _admit(catalog_path, catalog, identity)
    return commit(root, catalog_path, catalog, identity, component, ownership, kinds)


def publish(root, output_catalog):
    with _locked(root) as root:
        pointer, state = load(root)
        if state is None:
            raise ValueError('no committed dressed contributions')
        identity = state['identity']
        current_identity = source_identity(identity['input_dataset']['path'],
                                           identity['preparation_origin']['path'])
        if identity != current_identity:
            raise ValueError('changed actual source or original numerical reference')
        coverage_summary = summary(state)
        complete = all(item['complete'] for item in coverage_summary.values())
        path = Path(output_catalog).resolve()
        if path.is_relative_to(root / 'objects') or path.is_relative_to(root / 'states'):
            raise ValueError('publication may not overwrite internal immutable storage')
        if path.parent == root and path.name in ('current.json', 'identity.json', 'serial.lock'):
            raise ValueError('publication may not overwrite bank control files')
        if str(path) in state['sources']:
            raise ValueError('publication may not overwrite an admitted input')
        if path.exists():
            old = json.loads(path.read_bytes())
            if (old.get('schema') not in ('p254-dressed-functional-partial-v1', SCHEMA) or
                    old.get('input_dataset') != identity['input_dataset'] or
                    old.get('preparation_origin') != identity['preparation_origin'] or
                    'committed_pointer' not in old):
                raise ValueError('publication may only replace this dressed bank publication')
        path.parent.mkdir(parents=True, exist_ok=True)
        record = dict(schema='p254-dressed-functional-partial-v1', status='partial',
            full_solve_eligible=False, **identity, coverage=coverage_summary,
            committed_pointer=pointer, families=state['families'],
            missing_families=[k for k, v in coverage_summary.items() if not v['complete']],
            storage_zeros='not physical values; partial state is not solve eligible',
            source_bindings=list(state['sources'].values()),
            composition='Sigma=EE+rawEI+jointMinusIonicGamma; Pi=joint+embed111cubic')
        if not complete:
            _save(path, record)
            return record
        destination = path.parent / (path.stem + '-arrays-' + uuid.uuid4().hex)
        destination.mkdir()
        arrays = {}
        for kind in ('sigma', 'pi'):
            tail = TAILS[kind]
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
                    for a in range(31):
                        if not np.isfinite(output[a]).all():
                            raise ValueError('nonfinite composed canonical kernels')
                    output.flush()
                finally:
                    output._mmap.close()
                _persist(output_path)
                arrays[kind][order] = dict(binding(output_path), shape=[31, 31, *tail],
                                           dtype='complex128')
        _sync_dir(destination)
        _sync_dir(destination.parent)
        coverage = {name: dict(owned_pairs='all', completed_pairs='all') for name in COMPONENTS}
        for name in ('ee', 'mixed_ei', 'cubic'):
            coverage[name].update(owned_q=list(range(64)), completed_q=list(range(64)))
        coverage['joint'].update(q_additive=False, all114_coordinates=True,
                                  external_targets=list(range(128)))
        record.update(schema=SCHEMA, status='complete', full_solve_eligible=True,
            reference_source_kind=REFERENCE, pi_convention=PI_CONVENTION,
            coverage=coverage, exact_contribution_coverage=coverage_summary,
            missing_families=[], slots=128, bands=12, gamma_coordinates=114,
            global_q=64, axes=dict(sigma=SIGMA_AXES, pi=PI_AXES), **arrays)
        _save(path, record)
        return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    add = sub.add_parser('accumulate')
    add.add_argument('--bank', type=Path, required=True)
    add.add_argument('--catalog', type=Path, required=True)
    add.add_argument('--provider', type=Path)
    add.add_argument('--preparation', type=Path)
    pub = sub.add_parser('publish')
    pub.add_argument('--bank', type=Path, required=True)
    pub.add_argument('--output', type=Path, required=True)
    show = sub.add_parser('reopen')
    show.add_argument('--bank', type=Path, required=True)
    show.add_argument('--require-complete', action='store_true')
    args = parser.parse_args()
    if args.command == 'accumulate':
        result = accumulate(args.bank, args.catalog, provider=args.provider,
                            preparation=args.preparation)
    elif args.command == 'publish':
        result = publish(args.bank, args.output)
    else:
        result = reopen(args.bank, require_complete=args.require_complete)
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()

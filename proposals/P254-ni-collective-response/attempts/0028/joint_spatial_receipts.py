"""Authenticate existing selected-q receipts and their original case geometry.

The explicit manifest binds both result and request bytes: result-v2 omits
atom/axis and some seed metadata carried by request.json. Never infer them
from filenames, and never assemble missing q from conjugate symmetry.
"""
from hashlib import sha256
from pathlib import Path
import re

import numpy as np

from substrate_framework.material_control_reference import bound_json, bound_arrays
from substrate_framework.joint_spatial_readout import ORIGINS, complete_grid
from joint_spatial_provenance import (
    COMMITTED_SCHEMA, authenticate, binding_key, canonical, same_binding, source_key,
    validate_committed_split,
)
from joint_spatial_source_local import MATERIAL_SCHEMA, validate_source_local
from joint_spatial_dual_coherent import SCHEMA as DUAL_SCHEMA


def number(value):
    if isinstance(value, dict):
        if set(value) != {'real', 'imag'}:
            raise ValueError('expected encoded complex scalar')
        z = complex(value['real'], value['imag'])
    else:
        z = complex(value)
    if not np.isfinite(z):
        raise ValueError('nonfinite readout')
    return z




def read_entries(entries):
    """Return original forward/pair objects, never a synthesized producer receipt.

    All-mode uses the same original objects for both roles. Split entries bind
    four independent files. Both forms can coexist on a complete original grid.
    """
    records, requests, pairs, pair_requests = [], [], [], []
    for entry in entries:
        if set(entry) == {'result', 'request'}:
            r, request = bound_json(entry['result']), bound_json(entry['request'])
            if request.get('mode') != 'all':
                raise ValueError('single-receipt entry requires genuine all-mode output')
            pair, pair_request = r, request
        elif set(entry) == {'coherent', 'pair'}:
            if set(entry['coherent']) != {'result', 'seed'} or set(entry['pair']) != {'result', 'owner'}:
                raise ValueError('native dual coherent entry requires exact result/seed and pair/owner')
            r, request = bound_json(entry['coherent']['result']), bound_json(entry['coherent']['seed'])
            pair, pair_request = bound_json(entry['pair']['result']), bound_json(entry['pair']['owner'])
            if (r.get('schema') != DUAL_SCHEMA or not same_binding(r['seed'], entry['coherent']['seed'])
                    or not same_binding(r['owner'], entry['pair']['owner'])
                    or not same_binding(r['context'], pair['covariance']['left_context'])
                    or not same_binding(pair['owner'], entry['pair']['owner'])):
                raise ValueError('native coherent original owner/context/seed bindings differ')
        elif set(entry) == {'forward', 'pair'}:
            if set(entry['forward']) != {'result', 'request'}:
                raise ValueError('forward side must bind its original result and request')
            r, request = bound_json(entry['forward']['result']), bound_json(entry['forward']['request'])
            if set(entry['pair']) == {'result', 'owner'}:
                pair = bound_json(entry['pair']['result'])
                pair_request = bound_json(entry['pair']['owner'])
                if pair.get('schema') != MATERIAL_SCHEMA or not same_binding(entry['pair']['owner'], pair['owner']):
                    raise ValueError('native source-local pair must bind its own original owner')
            elif set(entry['pair']) == {'result', 'request'}:
                pair, pair_request = bound_json(entry['pair']['result']), bound_json(entry['pair']['request'])
                if pair_request.get('mode') != 'pair':
                    raise ValueError('legacy split pair requires genuine pair mode')
            else:
                raise ValueError('pair side must bind result/request or native result/owner')
            if request.get('mode') != 'forward':
                raise ValueError('split entry requires genuine forward mode')
            if (r.get('schema') == COMMITTED_SCHEMA
                    and not same_binding(entry['forward']['request'],
                                         r['historical_producer_request']['binding'])):
                raise ValueError('native forward entry must bind its original producer request')
        else:
            raise ValueError('entry must bind either all-mode or separate forward/pair outputs')
        records.append(r)
        requests.append(request)
        pairs.append(pair)
        pair_requests.append(pair_request)
    return records, requests, pairs, pair_requests


def validate_split(record, request, pair, pair_request):
    """Compare scientific identities BEFORE joining coherent and quantum values."""
    if record.get('schema') == COMMITTED_SCHEMA:
        return validate_committed_split(record, request, pair, pair_request)
    mode = request.get('mode')
    if mode == 'all':
        if pair_request != request or pair != record:
            raise ValueError('all-mode roles must refer to the same original receipt/request')
        return record['coherent']
    if (mode != 'forward' or pair_request.get('mode') != 'pair'
            or pair.get('status') != 'complete' or 'error' in pair
            or 'quantum' in record or 'coherent' in pair):
        raise ValueError('successful genuine split forward/pair receipts required')
    fields = ('schema', 'status', 'required_q_for_whole_local', 'fixed_original_D',
              'all_explicit_edges', 'case', 'model_identity', 'scalar', 'q', 'dynamics',
              'time_grid', 'second_time_grid', 'observable')
    if any(record[key] != pair[key] for key in fields):
        raise ValueError('split case/scalar/q/dynamics/time/observable differs')
    if any(not same_binding(record[key], pair[key])
           for key in ('initialization', 'baseline', 'matched_aggregate')):
        raise ValueError('split original initialization/baseline/aggregate differs')
    if source_key(record['source_bindings']) != source_key(pair['source_bindings']):
        raise ValueError('split producer sources differ')
    ignored = {'mode', 'left_checkpoint', 'left_checkpoint_sha256',
               'right_checkpoint', 'right_checkpoint_sha256'}
    if ({key: value for key, value in request.items() if key not in ignored}
            != {key: value for key, value in pair_request.items() if key not in ignored}):
        raise ValueError('split scientific request or original coherent seed differs')
    return record['coherent']


def validate_records(records, requests, saved, initialization, mesh, *,
                     pair_records=None, pair_requests=None, source_snapshots=()):
    """Join original all-mode or matched split receipts after byte authentication."""
    if not records or len(records) != len(requests):
        raise ValueError('one result and request per q required')
    if pair_records is None and pair_requests is None:
        pair_records, pair_requests = records, requests
    if (pair_records is None or pair_requests is None or len(pair_records) != len(records)
            or len(pair_requests) != len(records)):
        raise ValueError('one independently authenticated pair side per q required')
    complete_grid([r['q'] for r in records], mesh)
    n = int(np.prod(mesh))
    common, ab, ba, coherent, source_families = None, [], [], [], {}
    for r, request, pair, pair_request in zip(records, requests, pair_records, pair_requests):
        family = pair.get('schema')
        sources = pair['source_bindings']
        selected_sources = source_key(sources)
        if family in source_families and source_families[family] != selected_sources:
            raise ValueError('mixed selected-q producer sources within one algorithm')
        source_families[family] = selected_sources
        if family == MATERIAL_SCHEMA:
            selected, seed, rows = validate_source_local(r, request, pair, pair_request,
                saved, initialization, mesh, source_snapshots)
            if common is not None and common != selected:
                raise ValueError('mixed selected-q scientific request')
            common = selected
            coherent.append(number(seed['readout']))
            ab.append(rows[0])
            ba.append(rows[1])
            continue
        seed = validate_split(r, request, pair, pair_request)
        r, request = pair, pair_request
        if (r.get('schema') != 'P254-joint-q-material-v2' or r.get('status') != 'complete'
                or request.get('mode') not in ('all', 'pair') or 'error' in r
                or r.get('required_q_for_whole_local') != n
                or r.get('fixed_original_D') is not True
                or r.get('all_explicit_edges') != 2*n):
            raise ValueError('successful complete original selected-q receipt required')
        if (r['case'] != saved['case'] or request['case'] != saved['case']
                or r['model_identity'] != saved['model_identity']
                or r['scalar'] != saved['specification']['scalar']
                or not same_binding(r['matched_aggregate'], saved['specification']['aggregate'])
                or not same_binding(r['initialization'], initialization)
                or not same_binding(request['initialization'], initialization)
                or not same_binding(request['baseline'], r['baseline'])
                or request['q'] != r['q']):
            raise ValueError('mixed original case/scalar/initialization/baseline/q')
        quantum = pair['quantum']
        expected_identity = sha256(canonical(dict(model=r['model_identity'], q=r['q'],
            dynamics=r['dynamics'], sources=r['source_bindings'])).encode()).hexdigest()
        if (request['identity'] != expected_identity or quantum['propagation_identity'] != expected_identity
                or r['dynamics'] not in ('screened-joint', 'bare-electronic')
                or quantum['dynamics'] != r['dynamics']
                or quantum['preparation_identity'] != 'original-matched-product-'+saved['model_identity']):
            raise ValueError('mixed propagation/dynamics/original preparation')
        for times, key, endpoint in ((r['time_grid'], 'times', 'time'),
                                     (r['second_time_grid'], 'second_times', 'second_time')):
            t = np.asarray(times, float)
            if (t.ndim != 1 or not len(t) or not np.isfinite(t).all() or t[0] != 0
                    or np.any(np.diff(t) <= 0) or request[key] != times
                    or quantum[endpoint] != times[-1]):
                raise ValueError('mixed or invalid original two-sided time grid')
        observable = re.fullmatch(r'original-([A-Z][a-z]?)([1-9][0-9]*)([xyz])-total-force', r['observable'])
        if (observable is None or int(observable[2]) != request['atom'] or observable[3] != request['axis']
                or seed['amplitude'] != request['amplitude'] or seed['cell'] != request['cell']
                or not np.isfinite(request['amplitude']) or request.get('checkpoint') != 'None'
                or request.get('checkpoint_sha256') is not None):
            raise ValueError('mixed observable/coherent original kick or unsupported forward restart')
        cell = np.asarray(request['cell'])
        if cell.shape != (3,) or not np.array_equal(cell, np.rint(cell)):
            raise ValueError('original integer kick cell required')
        selected = dict(case=r['case'], baseline=binding_key(r['baseline']),
            model_identity=r['model_identity'], initialization=binding_key(initialization),
            preparation=quantum['preparation_identity'], observable=r['observable'],
            dynamics=r['dynamics'], times=r['time_grid'], second_times=r['second_time_grid'],
            atom=request['atom'], axis=request['axis'],
            seed_atom=request['atom'], seed_axis=request['axis'],
            amplitude=request['amplitude'], cell=request['cell'], subdivisions=request['subdivisions'],
            block_size=request['block_size'])
        if common is not None and common != selected:
            raise ValueError('mixed selected-q scientific request')
        common = selected
        for ordering, destination in (('AB', ab), ('BA', ba)):
            row = [number(quantum[ordering][key]) for key in ORIGINS]
            if abs(sum(row[:3])-row[3]) > 1e-11*max(1., sum(abs(x) for x in row)):
                raise ValueError('inconsistent initial-origin total')
            destination.append(row)
        if (abs(number(quantum['symmetrized'])-(ab[-1][3]+ba[-1][3])/2) >
                1e-11*max(1., abs(ab[-1][3]), abs(ba[-1][3]))
                or abs(number(quantum['commutator'])-(ab[-1][3]-ba[-1][3])) >
                1e-11*max(1., abs(ab[-1][3]), abs(ba[-1][3]))):
            raise ValueError('inconsistent symmetrized/commutator receipt')
        coherent.append(number(seed['readout']))
    return common, np.asarray(coherent), np.asarray(ab), np.asarray(ba)


def load_material(manifest_binding):
    """Read actual authenticated initialization, geometry and ALL64 q receipts."""
    manifest = bound_json(manifest_binding)
    if manifest.get('schema') != 'P254-joint-spatial-input-v1':
        raise ValueError('explicit spatial input manifest required')
    ib = manifest['initialization']
    initial = bound_json(ib)
    if (initial.get('schema') != 'P254-matched-joint-case-v1' or initial.get('status') != 'complete'
            or initial.get('mode') != 'initialize' or initial['case'] != initial['saved']['case']):
        raise ValueError('complete actual matched initialization required')
    saved = initial['saved']
    if saved['model_identity'] != sha256(canonical(saved['specification']).encode()).hexdigest():
        raise ValueError('initialization model/specification identity differs')
    origin = bound_json(saved['joint_boundary'])
    if (origin.get('schema') != 'streamed-joint-mean-v1' or origin.get('complete') is not True
            or origin['time_ry_inverse'] != 0 or origin['model_identity'] != saved['model_identity']
            or not same_binding(origin['original_electronic_boundary'], saved['original_boundary'])):
        raise ValueError('original joint boundary/preparation differs')
    authenticate(saved['original_boundary'])
    authenticate(saved['preparation'])
    slot_keys = [tuple(map(int, s['identity'].replace(':', ',').split(','))) for s in origin['slots']]
    mesh = np.max(np.asarray(slot_keys)[:, 1:], axis=0)+1
    expected = {(spin, *index) for spin in (1, 2)
                for index in np.indices(tuple(mesh)).reshape(3, -1).T}
    if len(slot_keys) != 128 or len(set(slot_keys)) != 128 or set(slot_keys) != expected or int(np.prod(mesh)) != 64:
        raise ValueError('complete original 128-spin-slot/64-cell material grid required')
    aggregate_binding = saved['specification']['aggregate']
    aggregate = bound_json(aggregate_binding)
    if aggregate.get('status') != 'complete':
        raise ValueError('complete matched aggregate required')
    arrays = bound_arrays(aggregate['archive'])
    cell_bohr = arrays['cell_bohr']
    if tuple(arrays['common_force_lattice_ry_bohr2'].shape[:3]) != tuple(mesh):
        raise ValueError('original lattice and original slot mesh differ')
    entries = manifest['selected_q']
    records, requests, pairs, pair_requests = read_entries(entries)
    common, coherent, ab, ba = validate_records(records, requests, saved, ib, mesh,
                                               pair_records=pairs, pair_requests=pair_requests,
                                               source_snapshots=manifest.get('source_snapshots', ()))
    atom, seed_atom = common['atom'], common['seed_atom']
    if (not 1 <= atom <= len(arrays['atom_names'])
            or not 1 <= seed_atom <= len(arrays['atom_names'])
            or common['observable'] != 'original-%s%d%s-total-force' % (arrays['atom_names'][atom-1], atom, common['axis'])
            or not np.array_equal(arrays['atom_indices'], np.arange(1, len(arrays['atom_names'])+1))):
        raise ValueError('selected force is not in original atom order')
    authenticated = set()
    for record in records+pairs:
        for source in record['source_bindings']:
            key = binding_key(source)
            if key not in authenticated:
                authenticate(source, manifest.get('source_snapshots', ()))
                authenticated.add(key)
    for pair in pairs:
        producer = Path(__file__).resolve().with_name(
            'run_source_local_q_material.py' if pair.get('schema') == MATERIAL_SCHEMA
            else 'run_joint_q_material.py')
        if not any(Path(source['path']).resolve() == producer for source in pair['source_bindings']):
            raise ValueError('selected-q receipt does not bind its original material producer')
    baseline_binding = records[0]['baseline']
    baseline = bound_json(baseline_binding)
    if (baseline.get('schema') != 'streamed-joint-mean-v1' or baseline.get('complete') is not True
            or baseline['model_identity'] != saved['model_identity']
            or not same_binding(baseline['original_electronic_boundary'], saved['original_boundary'])
            or baseline['time_ry_inverse'] < max(common['times'][-1], common['second_times'][-1])):
        raise ValueError('completed exact moving baseline must cover both requested times')
    if baseline['time_ry_inverse'] == 0:
        if not same_binding(baseline_binding, saved['joint_boundary']):
            raise ValueError('t0 baseline differs from original initialization')
    elif not same_binding(baseline['history'][0]['prior_boundary_binding'], saved['joint_boundary']):
        raise ValueError('moving history does not start at original matched boundary')
    mean = None
    if 'mean_diagnostics' in manifest:
        receipt = bound_json(manifest['mean_diagnostics'])
        if (receipt.get('schema') != 'P254-joint-partitioned-material-v1'
                or receipt.get('status') != 'complete' or receipt.get('case') != common['case']
                or not same_binding(receipt['initialization'], ib)):
            raise ValueError('actual matched mean diagnostic receipt required')
        branches = [branch for branch in receipt['branches'] if same_binding(branch['boundary'], baseline_binding)]
        if len(branches) != 1:
            raise ValueError('mean diagnostics must match exact q baseline boundary, not merely time')
        mean = dict(receipt=manifest['mean_diagnostics'], boundary=baseline_binding,
                    time_ry_inverse=baseline['time_ry_inverse'], diagnostics=branches[0]['diagnostics'],
                    scope='homogeneous baseline mean diagnostics only, not coherent-kick or local energy flow')
    geometry = dict(aggregate=aggregate_binding, archive=aggregate['archive'], cell_bohr=cell_bohr,
                    atom_one_based=atom, atom_name=str(arrays['atom_names'][atom-1]),
                    atom_position_bohr=arrays['positions_bohr'][:, atom-1],
                    seed_atom_one_based=seed_atom, seed_atom_name=str(arrays['atom_names'][seed_atom-1]),
                    seed_atom_position_bohr=arrays['positions_bohr'][:, seed_atom-1],
                    convention='cell lattice vectors are columns; coherent distance uses probe minus seed intra-cell position; same-probe connected offsets cancel')
    return dict(manifest=manifest, common=common, saved=saved, mesh=mesh, geometry=geometry,
                receipt_assembly_by_q=['all' if 'result' in entry else
                    ('dual-coherent+source-local-pair' if record.get('schema') == DUAL_SCHEMA else
                     'committed-forward+source-local-pair' if pair.get('schema') == MATERIAL_SCHEMA else
                     'committed-forward+pair' if record.get('schema') == COMMITTED_SCHEMA else 'forward+pair')
                    for entry, record, pair in zip(entries, records, pairs)],
                q=np.asarray([r['q'] for r in records]), coherent=coherent, ab=ab, ba=ba,
                mean_diagnostics=mean, baseline=baseline_binding,
                baseline_time=baseline['time_ry_inverse'],
                producer_diagnostics=[dict(q=record['q'], numerical_accuracy_certified=r.get('numerical_accuracy_certified'),
                    execution_scheme=r.get('covariance', r.get('quantum', {})).get('execution_scheme'),
                    raw_frame_compatibility=dict(
                        edge_count=len(r.get('covariance', r.get('quantum', {}))['raw_frame_compatibility']),
                        left_max=max(row['left'] for row in r.get('covariance', r.get('quantum', {}))['raw_frame_compatibility']),
                        right_max=max(row['right'] for row in r.get('covariance', r.get('quantum', {}))['raw_frame_compatibility'])))
                    for record, r in zip(records, pairs)])

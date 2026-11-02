"""Read native source-local publications without translating producer schemas.

Execution identities authenticate each algorithm independently. The saved model
specification, not equality of algorithm wrappers, is the physical authority.
"""
from hashlib import sha256
import io

import numpy as np

from substrate_framework.material_control_reference import bound_json, read_bound_bytes
from substrate_framework.partitioned_joint_slots import digest
from substrate_framework.joint_spatial_readout import ORIGINS
from joint_spatial_provenance import (
    COMMITTED_SCHEMA, authenticate, binding_key, canonical, committed_seed,
    same_binding, source_key,
)
from joint_spatial_dual_coherent import SCHEMA as DUAL_SCHEMA, validate_readout
from joint_spatial_dual_authority import native_proof, physical_sources

MATERIAL_SCHEMA = 'P254-source-local-q-material-v1-result'
OWNER_SCHEMA = 'P254-source-local-q-material-v1-owner'
CONVENTION = '<O_A(q,t) O_B(q,s)^dagger> connected; origins refer to original t0'


def scalar(value):
    """Native pair files use [real, imag], material receipts use {real, imag}."""
    if isinstance(value, dict) and set(value) == {'real', 'imag'}:
        value = [value['real'], value['imag']]
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ValueError('native encoded complex scalar required')
    result = complex(*value)
    if not np.isfinite(result):
        raise ValueError('nonfinite native covariance')
    return result


def _grid(grid, endpoint):
    values = np.asarray(grid, float)
    if (values.ndim != 1 or not len(values) or not np.isfinite(values).all()
            or values[0] != 0 or np.any(np.diff(values) <= 0) or values[-1] != endpoint):
        raise ValueError('native original two-sided time grid differs')



def _source_reduction(covariance, edge_count):
    reduction = covariance['source_contraction']
    request = bound_json(reduction['request'])
    if (reduction['count'] != edge_count or len(reduction['receipts']) != edge_count
            or request['schema'] != 'partitioned-q-work-v2'
            or request['operation'] != 'source-local-initial-contraction'
            or request['source_identity'] != covariance['propagation_identity']
            or request['count'] != edge_count
            or request['input']['preparation'] != covariance['preparation_identity']
            or not same_binding(request['input']['left'], covariance['left_context'])
            or not same_binding(request['input']['right'], covariance['right_context'])):
        raise ValueError('complete original electronic source reduction differs')
    ordered, defects = np.zeros(2, complex), []
    for index, binding in enumerate(reduction['receipts']):
        item = bound_json(binding)
        if item.get('complete') is not True or item['index'] != index or item['request'] != digest(request):
            raise ValueError('native original source order/identity differs')
        values = {}
        for name in ('ordered', 'compatibility'):
            raw, _ = read_bound_bytes(item['arrays'][name], max_bytes=4096)
            value = np.load(io.BytesIO(raw), allow_pickle=False)
            if value.shape != (2,) or not np.isfinite(value).all():
                raise ValueError('native scalar source payload differs')
            values[name] = value
        ordered += values['ordered']
        defects.append(dict(edge=index, left=float(values['compatibility'][0]),
                            right=float(values['compatibility'][1])))
    if defects != covariance['raw_frame_compatibility']:
        raise ValueError('native raw frame defects differ from source receipts')
    if any(abs(ordered[i]-scalar(covariance[o]['electronic'])) > 1e-11*max(1., abs(ordered[i]))
           for i, o in enumerate(('AB', 'BA'))):
        raise ValueError('native electronic AB/BA differs from original source reduction')


def validate_source_local(record, request, pair, owner, saved, initialization, mesh,
                          snapshots=()):
    """Join an actual material result+owner to a native committed forward scalar."""
    n = int(np.prod(mesh))
    dual = record.get('schema') == DUAL_SCHEMA
    if dual:
        seed = validate_readout(record, request)
        if (not same_binding(record['owner'], pair['owner'])
                or not same_binding(record['context'], pair['covariance']['left_context'])):
            raise ValueError('dual coherent owner/context differs from this noise pair')
    elif record.get('schema') == COMMITTED_SCHEMA:
        seed = committed_seed(record, request)
    else:
        raise ValueError('source-local pair requires native committed or dual coherent readout')
    cell = np.asarray(request['cell'])
    if (cell.shape != (3,) or not np.isfinite(cell).all() or not np.array_equal(cell, np.rint(cell))
            or not np.isfinite(request['amplitude'])
            or (not dual and (request.get('checkpoint') != 'None' or request.get('checkpoint_sha256') is not None))
            or saved['model_identity'] != sha256(canonical(saved['specification']).encode()).hexdigest()):
        raise ValueError('original coherent kick and authenticated physical model required')
    if (pair.get('schema') != MATERIAL_SCHEMA or pair.get('status') != 'complete'
            or pair.get('mode') != 'pair' or pair.get('covariance_executed') is not True
            or 'error' in pair or pair.get('required_q_for_whole_local') != n
            or owner.get('schema') != OWNER_SCHEMA):
        raise ValueError('executed complete native material pair and owner required')
    owned = bound_json(pair['owner'])
    if owned != owner:
        raise ValueError('source-local result does not bind this original owner')
    scientific = owner['request']
    if (scientific.get('schema') != 'P254-source-local-q-material-v1-request'
            or scientific.get('original_time_origin') != 0.):
        raise ValueError('original source-local scientific request required')
    identity = digest(dict(schema='P254-source-local-q-material-v1', owner=pair['owner']))
    if pair['material_identity'] != identity:
        raise ValueError('source-local material owner identity differs')
    sources = pair['source_bindings']
    if source_key(sources) != source_key(scientific['source_bindings']):
        raise ValueError('source-local owner/producer source bindings differ')
    for source in sources:
        authenticate(source, snapshots)
    for source in record['source_bindings']:
        authenticate(source, snapshots)
    p = pair['provenance']
    selectors = scientific['selectors']
    edge_count = len(bound_json(saved['joint_boundary'])['slots'])
    runtime = p['request']
    q = record['q']
    if (selectors['mode'] != 'pair' or selectors['q'] != q
            or runtime['selectors'] != {k: v for k, v in selectors.items() if k not in ('mode', 'edge_index')}
            or p['case'] != saved['case'] or selectors['case'] != saved['case']
            or record['case'] != saved['case']
            or p['model_identity'] != saved['model_identity'] or record['model_identity'] != saved['model_identity']
            or (record['physical_specification'] if dual else record['request']['operator']['specification']) != saved['specification']
            or p['scalar'] != saved['specification']['scalar'] or record['scalar'] != p['scalar']
            or p['observable'] != record['observable'] or selectors['observable'] != 'total-force'
            or selectors['field_axis'] is not None or p['observable_metadata'] is not None
            or selectors['atom'] != record['atom'] or selectors['axis'] != record['axis']
            or selectors['dynamics'] != record['dynamics'] or p['fixed_original_D'] is not True
            or p['all_explicit_edges'] != edge_count or record['all_explicit_edges'] != edge_count
            or (not dual and (request['atom'] != record['atom'] or request['axis'] != record['axis']))):
        raise ValueError('native forward/pair exact physical case/scalar/probe/geometry differs')
    for key, expected in (('initialization', initialization), ('baseline', record['baseline']),
                          ('matched_aggregate', saved['specification']['aggregate'])):
        if not same_binding(p[key], expected):
            raise ValueError('native original input differs: '+key)
    if (not same_binding(record['matched_aggregate'], saved['specification']['aggregate'])
            or not same_binding(request['initialization'], initialization)
            or (not dual and not same_binding(request['baseline'], p['baseline']))):
        raise ValueError('native historical forward physical inputs differ')
    if not same_binding(record['initialization'], initialization):
        raise ValueError('native coherent initialization differs')
    for key, runtime_key in (('initialization', 'initialization'), ('history', 'baseline'),
                             ('physical_qualification', 'physical_qualification'),
                             ('source_local_oracle', 'source_local_oracle')):
        b = runtime[runtime_key]
        explicit = scientific['explicit_inputs'][key]
        if explicit['path'] != b['path'] or explicit['sha256'] != b['sha256']:
            raise ValueError('native explicit owner inputs differ')
    if (pair['explicit_inputs'] != scientific['explicit_inputs']
            or not same_binding(runtime['initialization'], initialization)
            or not same_binding(runtime['baseline'], p['baseline'])
            or not same_binding(runtime['physical_qualification'], p['physical_qualification'])
            or not same_binding(runtime['source_local_oracle'], p['source_local_numerical_proof'])
            or source_key(p['source_bindings']) != source_key(sources)):
        raise ValueError('native runtime provenance differs')
    physical = native_proof(p['physical_qualification'], 'P254-joint-partitioned-q-oracle-v1',
                            sources, snapshots, 'source_bindings')
    native_proof(p['source_local_numerical_proof'], 'P254-source-local-q-oracle-v1',
                 sources, snapshots, 'sources')
    physical_sources(saved['specification'], sources, physical, snapshots)
    forward_physical = physical if dual else bound_json(record['physical_operator_oracle'])
    physical_sources(saved['specification'], record['source_bindings'], forward_physical, snapshots)
    covariance = pair['covariance']
    publication = bound_json(covariance['context_pair_binding'])
    if (publication.get('schema') != 'source-local-q-pair-v1' or publication.get('complete') is not True):
        raise ValueError('complete native context-pair publication required')
    original = publication['covariance']
    if set(covariance) != set(original) | {'context_pair_binding'}:
        raise ValueError('native pair publication fields differ')
    for key in original:
        if key in ('AB', 'BA'):
            if set(original[key]) != set(ORIGINS) or set(covariance[key]) != set(ORIGINS):
                raise ValueError('all native independent origins required')
            equal = all(scalar(original[key][o]) == scalar(covariance[key][o]) for o in ORIGINS)
        elif key in ('symmetrized', 'commutator'):
            equal = scalar(original[key]) == scalar(covariance[key])
        else:
            equal = original[key] == covariance[key]
        if not equal:
            raise ValueError('native pair result differs from original publication: '+key)
    if (covariance['preparation_identity'] != 'original-matched-product-'+saved['model_identity']
            or covariance['convention'] != CONVENTION or covariance['dynamics'] != record['dynamics']
            or covariance['execution_scheme'] != 'source-local-discrete-rk4-v1'):
        raise ValueError('native preparation/ordering/normalization convention differs')
    propagation = sha256(canonical(dict(model=saved['model_identity'], q=q,
        dynamics=selectors['dynamics'], sources=sources, observable_metadata=runtime)).encode()).hexdigest()
    if covariance['propagation_identity'] != propagation:
        raise ValueError('native source-local execution identity differs')
    for side, time_key, grid_key in (('left', 'time', 'time_grid'), ('right', 'second_time', 'second_time_grid')):
        grid = pair[grid_key]
        _grid(grid, selectors[time_key])
        if runtime['times' if side == 'left' else 'second_times'] != grid or covariance[time_key] != grid[-1]:
            raise ValueError('native runtime two-sided times differ')
        context = bound_json(covariance[side+'_context'])
        context_request = bound_json(context['request'])
        if (context.get('schema') != 'source-local-q-pullback-v1' or context.get('complete') is not True
                or context['time'] != 0 or context['terminal_time'] != grid[-1]
                or context_request['schema'] != 'source-local-q-request-v1'
                or context_request['identity'] != propagation or context_request['q'] != q
                or context_request['dynamics'] != record['dynamics'] or context_request['times'] != grid
                or context_request['provenance'] != p
                or context_request['probe']['identity'] != p['observable']
                or context_request['probe']['units'] != 'Ry/bohr'):
            raise ValueError('native terminal probe/context provenance differs')
    if ((record['time_grid'] if dual else request['times']) != pair['time_grid']
            or record['time_ry_inverse'] != selectors['time']):
        raise ValueError('native coherent/pair left-time grid differs')
    _source_reduction(covariance, edge_count)
    rows = [np.asarray([scalar(covariance[order][o]) for o in ORIGINS]) for order in ('AB', 'BA')]
    for row in rows:
        if abs(row[:3].sum()-row[3]) > 1e-11*max(1., sum(abs(row))):
            raise ValueError('native origin sum differs')
    if (abs(scalar(covariance['symmetrized'])-(rows[0][3]+rows[1][3])/2) > 1e-11*max(1., abs(rows[0][3]), abs(rows[1][3]))
            or abs(scalar(covariance['commutator'])-(rows[0][3]-rows[1][3])) > 1e-11*max(1., abs(rows[0][3]), abs(rows[1][3]))):
        raise ValueError('native ordered symmetrized/commutator differs')
    common = dict(case=saved['case'], baseline=binding_key(p['baseline']),
        model_identity=saved['model_identity'], initialization=binding_key(initialization),
        preparation=covariance['preparation_identity'], observable=p['observable'],
        dynamics=record['dynamics'], times=pair['time_grid'], second_times=pair['second_time_grid'],
        atom=record['atom'], axis=record['axis'], seed_atom=request['atom'], seed_axis=request['axis'],
        amplitude=request['amplitude'], cell=request['cell'],
        subdivisions=selectors['subdivisions'] if dual else request['subdivisions'],
        block_size=selectors['block_size'] if dual else request['block_size'])
    return common, seed, rows

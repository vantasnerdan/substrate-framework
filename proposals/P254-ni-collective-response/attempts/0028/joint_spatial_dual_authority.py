"""Original material authority for a complete one-sided force pullback.

A coherent canonical impulse consumes this pullback, not a later electronic
covariance reduction. Combined spatial/noise joins still authenticate that pair.
"""
from hashlib import sha256
from pathlib import Path
import re

import numpy as np

from substrate_framework.streamed_material_state import file_binding
from substrate_framework.material_control_reference import bound_json
from joint_spatial_provenance import authenticate, canonical, same_binding, source_key


def native_proof(binding, schema, sources, snapshots, source_field):
    proof = bound_json(binding)
    if proof.get('schema') != schema or proof.get('status') != 'passed':
        raise ValueError('passed bound numerical/physical authority required')
    expected = {str(Path(b['path']).resolve()):
                (int(b['bytes']), b['sha256']) for b in sources}
    bindings = proof[source_field]
    source_key(bindings)
    for source in bindings:
        authenticate(source, snapshots)
        if expected.get(str(Path(source['path']).resolve())) != (
                int(source['bytes']), source['sha256']):
            raise ValueError('native proof sources differ from its own producer')
    return proof


def physical_sources(specification, sources, qualification, snapshots):
    """Preserve the existing, specifically proven history-owner continuation."""
    actual = {str(Path(b['path']).resolve()): b for b in sources}
    continuation = qualification.get('source_compatibility')
    for original in specification['implementation']:
        authenticate(original, snapshots)
        selected = actual.get(str(Path(original['path']).resolve()))
        if selected is not None and same_binding(original, selected):
            continue
        if (not isinstance(continuation, dict)
                or continuation.get('kind') != 'same-physical-scalar-new-execution'
                or Path(original['path']).name != 'material_joint_response.py'
                or not same_binding(continuation['before'], original)
                or selected is None or not same_binding(continuation['after'], selected)):
            raise ValueError('native physical constructor/preparation authority differs')


def _original_probe(scientific, context, request):
    p = request['provenance']
    runtime, selectors = p['request'], scientific['selectors']
    mode = selectors.get('mode')
    if (mode not in ('pair', 'pullback') or selectors.get('observable') != 'total-force'
            or selectors.get('field_axis') is not None or p.get('observable_metadata') is not None
            or runtime['selectors'] != {k: v for k, v in selectors.items()
                                        if k not in ('mode', 'edge_index')}
            or request['q'] != selectors['q'] or request['dynamics'] != selectors['dynamics']
            or request['probe']['identity'] != p['observable']
            or not re.fullmatch(r'original-[A-Za-z]+%d%s-total-force' % (
                selectors['atom'], selectors['axis']), p['observable'])
            or request['probe']['units'] != 'Ry/bohr'):
        raise ValueError('original force probe, q or dynamics differs')
    grid = request['times']
    values = np.asarray(grid, float)
    choices = [(runtime['times'], selectors['time'])]
    if mode == 'pair':
        choices.append((runtime['second_times'], selectors['second_time']))
    elif ('second_time' in selectors or 'second_times' in runtime
            or 'second_time' in runtime['selectors']):
        raise ValueError('one-sided owner cannot claim a covariance second leg')
    if (values.ndim != 1 or not len(values) or not np.isfinite(values).all()
            or values[0] != 0 or np.any(np.diff(values) <= 0)
            or context['terminal_time'] != values[-1]
            or not any(grid == expected and values[-1] == end for expected, end in choices)):
        raise ValueError('complete original pullback time grid differs')
    return selectors


def _source_authority(scientific, request, snapshots):
    p = request['provenance']
    runtime, sources = p['request'], p['source_bindings']
    if source_key(sources) != source_key(scientific['source_bindings']):
        raise ValueError('original material and pullback source bindings differ')
    for source in sources:
        authenticate(source, snapshots)
    for name, runtime_name in (('initialization', 'initialization'), ('history', 'baseline'),
                               ('physical_qualification', 'physical_qualification'),
                               ('source_local_oracle', 'source_local_oracle')):
        explicit, binding = scientific['explicit_inputs'][name], runtime[runtime_name]
        if explicit['path'] != binding['path'] or explicit['sha256'] != binding['sha256']:
            raise ValueError('original owner explicit input differs: '+name)
    if (not same_binding(runtime['initialization'], p['initialization'])
            or not same_binding(runtime['baseline'], p['baseline'])
            or not same_binding(runtime['physical_qualification'], p['physical_qualification'])
            or not same_binding(runtime['source_local_oracle'], p['source_local_numerical_proof'])):
        raise ValueError('original pullback runtime provenance differs')
    propagation = sha256(canonical(dict(model=p['model_identity'], q=request['q'],
        dynamics=request['dynamics'], sources=sources, observable_metadata=runtime)).encode()).hexdigest()
    if request['identity'] != propagation:
        raise ValueError('original source-local propagation identity differs')
    return sources


def _original_preparation(p, selectors, snapshots):
    initial = bound_json(p['initialization'])
    saved = initial['saved']
    if (initial.get('status') != 'complete' or saved['case'] != p['case']
            or p['case'] != selectors['case'] or saved['model_identity'] != p['model_identity']
            or saved['specification']['scalar'] != p['scalar']
            or sha256(canonical(saved['specification']).encode()).hexdigest() != p['model_identity']
            or not same_binding(saved['specification']['aggregate'], p['matched_aggregate'])
            or p.get('fixed_original_D') is not True):
        raise ValueError('original saved physical specification or preparation differs')
    for binding in (saved['original_boundary'], saved['preparation']):
        authenticate(binding, snapshots)
    return saved


def _source_geometry(saved, p):
    origin = bound_json(saved['joint_boundary'])
    keys = [tuple(int(i) for i in s['identity'].replace(':', ',').split(',')) for s in origin['slots']]
    if not keys or any(len(key) != 4 or min(key[1:]) < 0 for key in keys):
        raise ValueError('literal original spin/source geometry required')
    mesh = np.max(np.asarray(keys)[:, 1:], axis=0)+1
    expected = {(spin, *index) for spin in {key[0] for key in keys}
                for index in np.indices(tuple(mesh)).reshape(3, -1).T}
    if (len(keys) != len(set(keys)) or set(keys) != expected
            or origin['time_ry_inverse'] != 0 or origin['model_identity'] != p['model_identity']
            or p['all_explicit_edges'] != len(keys)):
        raise ValueError('complete original spin/source geometry differs')
    return mesh, len(keys)


def complete_pullback(owner_binding, context_binding, snapshots=()):
    """Authenticate original owner, preparation, sources and one completed dual."""
    owner, context = bound_json(owner_binding), bound_json(context_binding)
    if (owner.get('schema') != 'P254-source-local-q-material-v1-owner'
            or context.get('schema') != 'source-local-q-pullback-v1'
            or context.get('complete') is not True or context.get('time') != 0):
        raise ValueError('original material owner and complete original-time pullback required')
    scientific = owner['request']
    request = bound_json(context['request'])
    if (scientific.get('schema') != 'P254-source-local-q-material-v1-request'
            or scientific.get('original_time_origin') != 0.
            or request.get('schema') != 'source-local-q-request-v1'):
        raise ValueError('original source-local material and pullback requests required')
    p = request['provenance']
    selectors = _original_probe(scientific, context, request)
    sources = _source_authority(scientific, request, snapshots)
    saved = _original_preparation(p, selectors, snapshots)
    mesh, slots = _source_geometry(saved, p)
    point = np.asarray(request['q'], float)
    if (point.shape != (3,) or not np.isfinite(point).all()
            or not np.allclose(point*mesh, np.rint(point*mesh), atol=1e-12, rtol=0)
            or len(request['sectors']) != slots):
        raise ValueError('original commensurate q or complete source slots differ')
    if selectors['mode'] == 'pullback':
        namespace = Path(owner_binding['path']).resolve().parent/'scientific-work'/'pullback'
        if (scientific.get('scientific_namespace') != 'scientific-work'
                or Path(context_binding['path']).resolve() != namespace/'complete.json'
                or Path(context['request']['path']).resolve() != namespace/'request.json'):
            raise ValueError('one-sided context must belong to the original material owner namespace')
        authority = bound_json(file_binding(namespace.parent.parent/'owner-authority.json'))
        if (authority != dict(schema='P254-source-local-q-material-v1-owner-authority',
                              owner=owner_binding)
                or type(owner.get('started_unix')) not in (int, float)
                or not np.isfinite(owner['started_unix'])
                or not isinstance(owner.get('material_id'), str)
                or not scientific.get('caps')):
            raise ValueError('original one-sided owner/start/cap authority differs')
    physical = native_proof(p['physical_qualification'], 'P254-joint-partitioned-q-oracle-v1',
                            sources, snapshots, 'source_bindings')
    native_proof(p['source_local_numerical_proof'], 'P254-source-local-q-oracle-v1',
                 sources, snapshots, 'sources')
    physical_sources(saved['specification'], sources, physical, snapshots)
    return context, request, p, saved, mesh, slots

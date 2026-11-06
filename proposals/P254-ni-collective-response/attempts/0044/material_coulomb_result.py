"""Consume a complete material result and its authenticated q chunks, no kernel replay."""
import json
from pathlib import Path
import numpy as np
from coulomb_material_inputs import HERE, file_binding
from substrate_framework.streamed_material_state import authenticate
from reduced_projection_cache import ProjectionCache


def consume_complete_result(directory, manifest, externals, matching, bounds):
    directory = Path(directory)
    path = directory/'result.json'
    if not path.exists():
        return None
    result = json.loads(path.read_bytes())
    if result.get('status') != 'complete':
        return None
    if (result['schema'] != 'p254-a44-original-bare-coulomb-material-v1'
            or result['external_slots'] != externals or result['times_inverse_ry'] != [0.,1.,10.,30.]
            or result['matching'] != matching or result['source_bounds'] != bounds
            or result['implementation'] != [file_binding(p) for p in sorted(HERE.glob('*.py'))]
            or result['input_bindings']['manifest'] != file_binding(manifest)):
        raise ValueError('complete material source/state/model/time binding differs')
    inputs = result['input_bindings']
    for key in ('manifest','original_boundary','background_receipt','background_archive'):
        authenticate(inputs[key])
    for binding in inputs['authenticated_archives']+inputs['producer_implementation']:
        authenticate(binding)
    specification = {k:result[k] for k in ('schema','matching','source_bounds','input_bindings',
                                         'implementation','times_inverse_ry','external_slots')}
    cache = ProjectionCache(directory/'transfer-q-chunks',specification)
    if cache.manifest['status'] != 'complete':
        raise ValueError('complete result lacks complete immutable q chunks')
    keys = [o+'_'+k for o in ('greater','lesser','retarded') for k in ('direct','exchange','total')]
    if [item['external'] for item in result['items']] != externals:
        raise ValueError('complete final external inventory differs')
    observed = []
    for item in result['items']:
        external = item['external']
        total = {k:np.zeros((4,4,12,12),complex) for k in keys}
        for q in range(64):
            key = f'slot/{external:03d}/q/{q:02d}'
            arrays = cache.get(key)
            meta = cache.manifest['items'][key]['metadata']
            if (meta['external'] != external or meta['outer_transfer_q'] != q or meta['quartets'] != 192
                    or set(arrays) != set(keys) or any(x.shape != (4,4,12,12) for x in arrays.values())):
                raise ValueError('complete q chunk scope/layout differs')
            for name in keys:
                total[name] += arrays[name]
        with np.load(authenticate(item['archive']),allow_pickle=False) as saved:
            if set(saved.files) != set(keys) or not all(np.array_equal(saved[k],total[k]) for k in keys):
                raise ValueError('saved final result differs from authenticated original-order q sum')
            norms = {k:[[float(np.linalg.norm(saved[k][a,b])) for b in range(4)] for a in range(4)] for k in keys}
        observed.append(dict(external=external,archive=item['archive'],norms=norms,
                             authenticated_q_chunks=64,quartets=12288))
    return dict(result=file_binding(path),complete_result_reused=True,kernel_replays=0,observed=observed)

"""Authenticate persisted component descriptors and consume exactly covered cells."""
from pathlib import Path
import json
import numpy as np
from component_emitter import ORDERS, authenticate, binding


def consume_component(catalog):
    record = json.loads(Path(catalog).read_bytes())
    if record['schema'] != 'p254-additive-history-component-v1' or record['full_mesh'] is not False:
        raise ValueError('explicit additive component required')
    q, = record['owned_q']
    if type(q) is not int or not 0 <= q < 64:
        raise ValueError('original q required')
    authenticate(record['input_dataset'])
    inputs = json.loads(authenticate(record['input_bindings']).read_bytes())
    routes = record['canonical_routes']
    if (routes['source'] != list(range(128)) or sorted(routes['target']) != list(range(128))
            or any(routes['slot_spin'][s] != routes['slot_spin'][t] for s,t in zip(routes['source'],routes['target']))):
        raise ValueError('full canonical spin-preserving route required')
    times = np.asarray(record['times_inverse_ry'],float)
    nt = len(times)
    cubic = record['component'] == 'cubic'
    shape = (nt,nt,111,111) if cubic else (nt,nt,128,12,12)
    descriptors = record['arrays' if cubic else 'sigma']
    if set(descriptors) != set(ORDERS):
        raise ValueError('three contour arrays required')
    owned = [tuple(p) for p in record['coverage']['owned_pairs']]
    complete = [tuple(p) for p in record['coverage']['completed_pairs']]
    if len(set(owned)) != len(owned) or len(set(complete)) != len(complete) or not set(complete) <= set(owned):
        raise ValueError('duplicate/out-of-scope coverage')
    arrays = {}
    seen = np.load(authenticate(record['coverage']['array']),mmap_mode='r',allow_pickle=False)
    try:
        if seen.shape != (nt,nt,1 if cubic else 128) or seen.dtype != bool:
            raise ValueError('slot coverage shape/dtype mismatch')
        if [p for p in owned if seen[p].all()] != complete:
            raise ValueError('pair and slot coverage disagree')
        if record['coverage']['completed_slot_cells'] != int(seen.sum()):
            raise ValueError('slot count mismatch')
        for a in range(nt):
            for b in range(nt):
                if (a,b) not in owned and seen[a,b].any():
                    raise ValueError('unowned slot claimed')
        done = len(complete) == len(owned)
        if (record['status'] != ('complete-owned-component-q%d'%q if done else 'partial-owned-component-q%d'%q)
                or record['completed_q'] != ([q] if done else [])):
            raise ValueError('owned-scope status mismatch')
        norms = {key:0. for key in ORDERS}
        maxima = {key:0. for key in ORDERS}
        causal = 0.
        count = 0
        for order,item in descriptors.items():
            arrays[order] = np.load(authenticate(item),mmap_mode='r',allow_pickle=False)
            arr = arrays[order]
            if arr.shape != shape or arr.dtype != np.complex128 or item['shape'] != list(shape) or item['dtype'] != str(arr.dtype):
                raise ValueError('canonical numeric axes/dtype differ')
        for a,b in owned:
            slots = np.flatnonzero(seen[a,b])
            if not len(slots):
                continue
            values = {key:arr[a,b] if cubic else arr[a,b,slots] for key,arr in arrays.items()}
            count += values['greater'].size
            for order,value in values.items():
                if not np.isfinite(value).all():
                    raise ValueError('nonfinite saved component value')
                norms[order] += float(np.vdot(value,value).real)
                maxima[order] = max(maxima[order],float(np.max(np.abs(value))))
            expected = values['greater']-values['lesser'] if times[a]>=times[b] else np.zeros_like(values['retarded'])
            causal = max(causal,float(np.max(np.abs(expected-values['retarded']))))
        if causal != 0.:
            raise ValueError('saved retarded identity differs')
        return dict(schema='p254-a76-saved-component-consumer-v1',catalog=binding(catalog),
            component=record['component'],owned_q=[q],completed_pairs=len(complete),
            completed_slot_cells=int(seen.sum()),scalar_elements_consumed_per_order=count,
            all_descriptors_authenticated=True,canonical_routes_consumed=routes,
            original_input_binding=record['input_dataset'],inputs=binding(authenticate(record['input_bindings'])),
            norms={key:float(np.sqrt(value)) for key,value in norms.items()},max_abs=maxima,
            causal_max_abs_defect=causal,full_mesh=False,full_time_square=record['full_time_square'],
            actual_times=times.tolist(),source_input_sections=sorted(inputs))
    finally:
        seen._mmap.close()
        for arr in arrays.values():
            arr._mmap.close()

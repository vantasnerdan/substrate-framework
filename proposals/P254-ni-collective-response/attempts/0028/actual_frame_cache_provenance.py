"""Exact actual-history/request and original worker provenance for frame reuse."""
import math
import os
from pathlib import Path

from substrate_framework.partitioned_joint_slots import digest
from substrate_framework.streamed_joint_state import DiskJointMean
from substrate_framework.streamed_material_state import DiskMean, authenticate, file_binding
from actual_frame_cache_io import (
    _bound_json, _checked_binding, _directory, _frame, _json, _safe, _same, _snapshot)

SCHEMA = 'P254-actual-frame-cache-v1'


def _identity(value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError('explicit nonempty identity/owner required')
    return value


def _request(request, identity):
    if (set(request) != {'schema', 'source_identity', 'operation', 'input', 'count'}
            or request['schema'] != 'partitioned-q-work-v2'
            or request['operation'] != 'actual-baseline-frames'
            or request['source_identity'] != _identity(identity)
            or type(request['count']) is not int or request['count'] <= 0):
        raise ValueError('exact actual-baseline-frames QWork v2 request required')
    digest(request)
    return request


def _history_slots(identity, count, history_binding):
    """Validate q-independent stage/history and original per-slot shape bindings."""
    if set(identity) != {'time', 'interval', 'stage', 'baseline_boundary', 'model',
                         'origin', 'block_size'}:
        raise ValueError('unknown actual frame input identity')
    if (type(identity['interval']) is not int or type(identity['block_size']) is not int
            or identity['block_size'] <= 0 or type(identity['time']) not in (int, float)
            or not math.isfinite(identity['time'])):
        raise ValueError('invalid actual time/interval/block size')
    _bound_json(history_binding)
    final = DiskJointMean(history_binding['path'], expected_identity=identity['model'], verify=False)
    if not final.history or final.time <= 0 or not _same(identity['origin'], final.origin):
        raise ValueError('completed moving actual history/original source required')
    original = DiskMean(authenticate(_checked_binding(final.origin)), verify=False)
    priors, last = [], 0.
    for index, stage in enumerate(final.history):
        _bound_json(stage['prior_boundary_binding'])
        prior = DiskJointMean(stage['prior_boundary_binding']['path'],
                              expected_identity=identity['model'], verify=False)
        dt = stage['dt_ry_inverse']
        if (stage['kind'] != 'actual-frozen-midpoint' or not math.isfinite(dt) or dt <= 0
                or prior.time != last or stage['start_time_ry_inverse'] != last
                or stage['time_ry_inverse'] != last+dt/2
                or not _same(prior.history, final.history[:index])
                or not _same(prior.origin, final.origin)
                or len(prior.slots) != len(final.slots)):
            raise ValueError('broken completed actual history chain')
        for left, right in zip(prior.slots, final.slots):
            if (left.metadata() != right.metadata()
                    or not _same(left.coefficient, right.coefficient)):
                raise ValueError('actual history changed original slot/basis/gauge/D')
        priors.append(prior)
        last += dt
    if last != final.time or len(final.slots) != count or original.time != 0:
        raise ValueError('history coverage/count/original time differs')
    first = priors[0]
    if (len(original.slots) != count or any(a.metadata() != b.metadata()
            or not _same(a.frame, b.frame) or not _same(a.coefficient, b.coefficient)
            for a, b in zip(original.slots, first.slots))):
        raise ValueError('history is not bound to the literal original preparation')
    interval = identity['interval']
    if not 0 <= interval < len(priors):
        raise ValueError('actual frame interval not covered')
    stage, prior = final.history[interval], priors[interval]
    if (not prior.time < identity['time'] < prior.time+stage['dt_ry_inverse']
            or not _same(identity['stage'], stage)
            or not _same(identity['baseline_boundary'], file_binding(prior.path))):
        raise ValueError('frame input is not this actual interior stage/boundary')
    for binding in stage['arrays'].values():
        authenticate(_checked_binding(binding))
    return [dict(slot.metadata(), shape=list(slot.shape)) for slot in first.slots]


def _source_alive(directory, binding):
    if os.path.lexists(directory/'retire.json'):
        raise ValueError('retired/retiring source request is not pinnable')
    _bound_json(binding)


def _origin(directory, request_key, index):
    path = _safe(directory/('item-%04d' % index)/'complete.json', missing=True)
    if not os.path.lexists(path):
        return None
    record, completion_binding = _snapshot(path)
    if (record.get('complete') is not True or record.get('request') != request_key
            or type(record.get('index')) is not int or record['index'] != index
            or set(record.get('arrays', {})) != {'frame'}
            or not {'cpu', 'affinity', 'pid', 'blas_threads', 'cpu_seconds',
                    'elapsed_seconds', 'hamiltonian_actions',
                    'process_peak_rss_bytes'}.issubset(record)):
        raise ValueError('original completed QWork worker frame receipt required')
    binding = record['arrays']['frame']
    _checked_binding(binding)
    payload = Path(binding['path'])
    if (payload.name != 'frame.npy' or not payload.parent.name.startswith('attempt-')
            or payload.parent.parent != path.parent):
        raise ValueError('frame payload escaped its original item/attempt binding')
    return dict(request_key=request_key, index=index, completion=completion_binding,
                worker_receipt=record, worker_metadata='historical, not current execution')


def _cached_item(cache, descriptor, index):
    item = _safe(cache/('item-%04d' % index))
    receipt = _json(item/'complete.json')
    origin = _bound_json(receipt['origin'])
    if (receipt.get('schema') != SCHEMA or receipt.get('complete') is not True
            or receipt.get('index') != index
            or receipt['origin']['path'] != str(item/'origin.json')
            or origin['request_key'] != digest(descriptor['original_request'])
            or origin['index'] != index):
        raise ValueError('modified cached frame completion/provenance')
    old = origin['worker_receipt']
    if (old.get('complete') is not True or old.get('request') != origin['request_key']
            or old.get('index') != index or set(old.get('arrays', {})) != {'frame'}
            or receipt['frame']['path'] != str(item/'frame.npy')):
        raise ValueError('cached frame lost its original per-item binding')
    _frame(receipt['frame'], old['arrays']['frame'], descriptor['slots'][index]['shape'])
    return receipt, origin


def _complete_cache(cache, owner):
    cache = _directory(cache)
    descriptor = _json(cache/'cache.json')
    if (descriptor.get('schema') != SCHEMA or descriptor.get('owner') != owner
            or descriptor.get('cache') != str(cache)):
        raise ValueError('explicit cache owner/path differs')
    request = _request(descriptor['original_request'], descriptor['original_request']['source_identity'])
    slots = _history_slots(request['input'], request['count'], descriptor['history'])
    if not _same(slots, descriptor['slots']):
        raise ValueError('cache original slot/history identity differs')
    full = _json(cache/'cache-complete.json')
    if (full.get('schema') != SCHEMA or full.get('complete') is not True
            or full.get('count') != request['count'] or full.get('missing') != []
            or full.get('available') != list(range(request['count']))
            or not _same(full.get('cache'), file_binding(cache/'cache.json'))
            or len(full.get('items', [])) != request['count']):
        raise ValueError('complete immutable cache barrier required for adoption')
    items = []
    for index, binding in enumerate(full['items']):
        if binding['path'] != str(cache/('item-%04d' % index)/'complete.json'):
            raise ValueError('cache barrier item binding escaped owner')
        _bound_json(binding)
        items.append(_cached_item(cache, descriptor, index))
    return descriptor, items

"""Restore adopted-frame metadata and admit its real full-shape QWork payloads."""
import os
from pathlib import Path

import numpy as np

from substrate_framework.partitioned_joint_slots import digest
from substrate_framework.streamed_material_state import file_binding
from actual_frame_cache_io import (
    _bound_json, _checked_binding, _directory, _json, _safe, _same, _snapshot)
from actual_frame_cache_provenance import _origin


def restore_adoption_metadata(cache, target, owner, expected, descriptor, items, budget):
    """Bind one complete target barrier to its original cached worker receipts."""
    key = digest(expected)
    source_binding = descriptor['original_request_binding']
    source_directory = Path(source_binding['path']).parent
    old = descriptor['original_request']
    namespace = _json(target/'adoption.json')
    adopted, adopted_binding = _snapshot(target/'adopted.json')
    if (namespace.get('schema') != 'P254-actual-frame-adoption-v1'
            or namespace.get('owner') != owner or namespace.get('target') != str(target)
            or namespace.get('producer_identity') != expected['source_identity']
            or namespace.get('independent_consumer_authentication_required') is not True
            or not _same({name: adopted.get(name) for name in namespace}, namespace)
            or adopted.get('complete') is not True or len(adopted.get('requests', [])) != 1):
        raise ValueError('complete owned single-request adoption required')
    barrier = adopted['requests'][0]
    generation = _safe(barrier['request']['path']).parent
    root = _directory(target/key)
    if (generation.parent != root or not generation.name.startswith('generation-')
            or Path(barrier['request']['path']).name != 'request.json'
            or set(root.iterdir()) != {generation}
            or set(target.iterdir()) != {target/'adoption.json', target/'adopted.json', root}
            or barrier.get('complete') is not True or barrier.get('count') != 128
            or len(barrier.get('receipts', [])) != 128):
        raise ValueError('adopted request/generation/full barrier differs')
    _directory(generation)
    if set(path.name for path in generation.iterdir()) != {
            'request.json', 'adoption.json', *('item-%04d' % i for i in range(128))}:
        raise ValueError('adoption generation has missing/extra/retired work')
    actual = _bound_json(barrier['request'])
    if not _same(actual, expected):
        raise ValueError('explicit expected scientific request differs')
    provenance, provenance_binding = _snapshot(generation/'adoption.json')
    if (provenance.get('schema') != 'P254-actual-frame-request-adoption-v1'
            or provenance.get('producer_identity') != expected['source_identity']
            or not _same(provenance.get('new_request'), barrier['request'])
            or not _same(provenance.get('original_request'), old)
            or not _same(provenance.get('original_request_binding'), source_binding)
            or not _same(provenance.get('cache_complete'), file_binding(cache/'cache-complete.json'))):
        raise ValueError('adoption lost original cache/request/history provenance')
    _bound_json(provenance['cache_complete'])

    # Authenticate all completion metadata before QWork.request can see it.
    receipts = []
    for index, (cached, origin) in enumerate(items):
        budget.check_processes(os.getpid(), ())
        item = _directory(generation/('item-%04d' % index))
        binding = barrier['receipts'][index]
        if binding['path'] != str(item/'complete.json'):
            raise ValueError('adoption barrier receipt escaped its slot')
        receipt = _bound_json(binding)
        reuse = receipt.get('reuse_provenance', {})
        if (receipt.get('complete') is not True or receipt.get('request') != key
                or type(receipt.get('index')) is not int or receipt['index'] != index
                or set(receipt.get('arrays', {})) != {'frame'}
                or receipt.get('execution') != 'adopted-existing-frame; no current worker/cost/PID'
                or not _same(reuse.get('adoption'), provenance_binding)
                or not _same(reuse.get('cache_item'), file_binding(cache/('item-%04d' % index)/'complete.json'))
                or not _same(reuse.get('original_worker_provenance'), origin)
                or not _same(_origin(source_directory, digest(old), index), origin)):
            raise ValueError('adopted/original receipt or historical provenance differs')
        frame = _checked_binding(receipt['arrays']['frame'])
        original = _checked_binding(origin['worker_receipt']['arrays']['frame'])
        attempt = _directory(Path(frame['path']).parent)
        if (attempt.parent != item or not attempt.name.startswith('attempt-reuse-')
                or Path(frame['path']).name != 'frame.npy'
                or set(item.iterdir()) != {attempt, item/'complete.json'}
                or set(attempt.iterdir()) != {Path(frame['path'])}
                or frame['sha256'] != original['sha256']
                or int(frame['bytes']) != int(original['bytes'])
                or not os.path.samefile(frame['path'], cached['frame']['path'])
                or not os.path.samefile(frame['path'], original['path'])):
            raise ValueError('target/cache/original full-frame inode or SHA/size differs')
        receipts.append(receipt)
    return dict(request=actual, generation=generation, provenance=provenance_binding,
                cache_complete=provenance['cache_complete'],
                adopted=adopted_binding, receipts=receipts)


def collect_original_frames(work, metadata, derived_key, items, descriptor, budget):
    """Use actual QWork and authenticate every original/target full-frame alias."""
    derived_directory = metadata['generation']
    receipts = metadata['receipts']
    records = work.collect(derived_directory, derived_key, 128)
    admitted = []
    for index, (record, receipt, (cached, origin)) in enumerate(zip(records, receipts, items)):
        budget.check_processes(os.getpid(), ())
        if not _same(record, receipt):
            raise ValueError("QWork completion metadata changed during admission")
        original = origin['worker_receipt']['arrays']['frame']
        work.verify(original)
        value = np.load(original['path'], mmap_mode='r', allow_pickle=False)
        try:
            if (not isinstance(value, np.memmap) or value.flags.writeable
                    or value.dtype != np.dtype(np.complex128)
                    or list(value.shape) != descriptor['slots'][index]['shape']):
                raise ValueError('original readonly full complex128 slot header differs')
            shape, dtype = list(value.shape), value.dtype.str
        finally:
            if isinstance(value, np.memmap):
                value._mmap.close()
        info = _safe(original['path']).stat()
        if (not os.path.samefile(record['arrays']['frame']['path'], original['path'])
                or not os.path.samefile(cached['frame']['path'], original['path'])):
            raise ValueError('original inode linkage changed during QWork verification')
        admitted.append(dict(index=index, shape=shape, dtype=dtype,
            target_frame=record['arrays']['frame'], cache_frame=cached['frame'],
            original_frame=original, device=info.st_dev, inode=info.st_ino,
            link_count=info.st_nlink, original_worker_provenance=origin))
    return admitted

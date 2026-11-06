#!/usr/bin/env python3
"""Pin completed actual-baseline-frames, then adopt them into NEW q work.

No propagation, array copy, source deletion, polling, or producer admission is
performed. Each cache owns one exact source request. Repeating `pin` retains
already authenticated items and adds newly completed items; only a full set gets
cache-complete.json. `adopt` requires complete caches and explicit expected NEW
QWork requests. Its destination is the future scientific-work/frames directory,
which must not exist, not an active producer's namespace.
The physical producer identity may be unchanged: forward/pair modes at the same
q/observable share that identity. NEW always describes the target namespace.

API: pin(source_request, cache, *, owner, source_identity, expected_input,
         count, history_binding)
     adopt(caches, target, *, owner, producer_identity, expected_requests)

Inputs/expected requests are ordinary JSON objects, not namespace-derived
identities. history_binding authenticates the original completed actual mean
boundary covering the requested time. CLI JSON files provide these same objects.
Payloads must share a filesystem with their destinations: hardlink failure is
fatal, never a reason to copy. Links pin physical storage until ALL aliases are
removed, even after the original producer retires its own link. Do not count
producer retirement as physical disk reclamation while a cache/adoption survives.

Caches are private, explicitly owned, advisory-locked directories. All traversed
paths reject symlinks; existing unowned caches and ALL existing adoption targets
are refused. Interrupted owned work is retained, never cleaned by this tool.
Completion metadata is published only after target bytes match the ORIGINAL
expected SHA/size and full complex128 shape. Original request/completion bindings
and complete worker receipts survive in provenance; worker costs/PIDs are always
historical, never attributed to this invocation. The actual QWork consumer must
independently hash/admit adopted arrays. No scientific result is certified here.
"""
import argparse
from contextlib import ExitStack
import json
import os
from pathlib import Path
import sys
from uuid import uuid4

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'src'))

from substrate_framework.partitioned_joint_slots import digest
from substrate_framework.streamed_material_state import file_binding
from actual_frame_cache_io import (
    STORAGE, _bound_json, _directory, _disjoint, _json, _link, _locked, _mkdir,
    _publish, _safe, _same, _snapshot, _storage)
from actual_frame_cache_provenance import (
    SCHEMA, _cached_item, _complete_cache, _history_slots, _identity, _origin,
    _request, _source_alive)


def pin(source_request, cache, *, owner, source_identity, expected_input, count, history_binding):
    """Pin currently completed items; missing items remain an explicit partial set."""
    _identity(owner)
    source_request = _safe(source_request)
    request, source_binding = _snapshot(source_request)
    _request(request, source_identity)
    if (not _same(request['input'], expected_input) or type(count) is not int
            or request['count'] != count):
        raise ValueError('source actual input/count differs from explicit expectation')
    directory, key = source_request.parent, digest(request)
    if (source_request.name != 'request.json' or directory.parent.name != key
            or not directory.name.startswith('generation-')):
        raise ValueError('source request is not its content-addressed QWork generation')
    _source_alive(directory, source_binding)
    slots = _history_slots(expected_input, count, history_binding)
    cache = _safe(cache, missing=True)
    _disjoint(cache, directory.parent.parent)
    descriptor = dict(schema=SCHEMA, owner=owner, cache=str(cache),
        original_request=request, original_request_binding=source_binding,
        history=history_binding, slots=slots, storage=STORAGE)
    if not cache.exists():
        _directory(cache, new=True)
    with _locked(cache):
        marker = cache/'cache.json'
        if marker.exists():
            if not _same(_json(marker), descriptor):
                raise ValueError('cache owner/original request/history binding differs')
        else:
            if any(cache.iterdir()):
                raise ValueError('unowned nonempty cache refused')
            _publish(marker, descriptor)
        available, missing, bindings = [], [], []
        for index in range(count):
            _source_alive(directory, source_binding)
            item = cache/('item-%04d' % index)
            if os.path.lexists(item/'complete.json'):
                _cached_item(cache, descriptor, index)
            else:
                origin = _origin(directory, key, index)
                if origin is None:
                    missing.append(index)
                    continue
                _mkdir(item)
                origin_path = item/'origin.json'
                if os.path.lexists(origin_path):
                    if not _same(_json(origin_path), origin):
                        raise ValueError('original completion changed during interrupted pin')
                else:
                    _publish(origin_path, origin)
                expected = origin['worker_receipt']['arrays']['frame']
                frame = _link(expected['path'], item/'frame.npy', expected, slots[index]['shape'])
                _source_alive(directory, source_binding)
                _bound_json(origin['completion'])
                _publish(item/'complete.json', dict(schema=SCHEMA, complete=True, index=index,
                    origin=file_binding(origin_path), frame=frame))
            available.append(index)
            bindings.append(file_binding(item/'complete.json'))
        _source_alive(directory, source_binding)
        summary = dict(schema=SCHEMA, complete=not missing, cache=file_binding(marker),
                       count=count, available=available, missing=missing, items=bindings,
                       pinned_storage=_storage(cache/('item-%04d' % i)/'frame.npy' for i in available),
                       storage=STORAGE)
        completion = cache/'cache-complete.json'
        if not missing:
            if os.path.lexists(completion):
                if not _same(_json(completion), summary):
                    raise ValueError('modified full cache barrier')
            else:
                _publish(completion, summary)
        elif os.path.lexists(completion):
            raise ValueError('claimed complete cache has missing items')
        return summary


def adopt(caches, target, *, owner, producer_identity, expected_requests):
    """Adopt complete pinned sets into an absent, explicitly NEW frames namespace."""
    _identity(owner)
    _identity(producer_identity)
    caches = [_safe(path) for path in caches]
    expected_requests = list(expected_requests)
    if not caches or len(caches) != len(set(caches)) or len(caches) != len(expected_requests):
        raise ValueError('one distinct complete cache per explicit expected request required')
    target = _safe(target, missing=True)
    if os.path.lexists(target):
        raise FileExistsError('adoption target MUST be a new namespace: '+str(target))
    with ExitStack() as locks:
        for cache in sorted(caches):
            _disjoint(cache, target)
            locks.enter_context(_locked(cache))
        admitted, keys = [], set()
        # Refuse incompatible inputs/partial caches before creating any target.
        for cache, request in zip(caches, expected_requests):
            _request(request, producer_identity)
            descriptor, items = _complete_cache(cache, owner)
            old = descriptor['original_request']
            _disjoint(Path(descriptor['original_request_binding']['path']).parent.parent.parent, target)
            if not _same(request['input'], old['input']) or request['count'] != old['count']:
                raise ValueError('byte-exact actual input/count equality required')
            key = digest(request)
            if key in keys:
                raise ValueError('duplicate target scientific request')
            keys.add(key)
            admitted.append((cache, descriptor, items, request, key))
        _directory(target, new=True)
        with _locked(target):
            namespace = dict(schema='P254-actual-frame-adoption-v1', owner=owner,
                target=str(target), producer_identity=producer_identity,
                execution='hardlink adoption; no workers launched, no propagation',
                historical_worker_costs='retained only in original_worker_provenance',
                independent_consumer_authentication_required=True, storage=STORAGE)
            _publish(target/'adoption.json', namespace)
            completed = []
            for cache, descriptor, items, request, key in admitted:
                root = _mkdir(target/key)
                generation = _mkdir(root/('generation-'+uuid4().hex))
                request_binding = _publish(generation/'request.json', request)
                provenance = dict(schema='P254-actual-frame-request-adoption-v1',
                    new_request=request_binding, producer_identity=producer_identity,
                    original_request=descriptor['original_request'],
                    original_request_binding=descriptor['original_request_binding'],
                    cache_complete=file_binding(cache/'cache-complete.json'),
                    storage=STORAGE, execution='reuse only; original worker costs are historical')
                provenance_binding = _publish(generation/'adoption.json', provenance)
                receipts = []
                for index, (cached, origin) in enumerate(items):
                    # Reauthenticate after preflight; publication binds ORIGINAL bytes.
                    _cached_item(cache, descriptor, index)
                    item = _mkdir(generation/('item-%04d' % index))
                    attempt = _mkdir(item/('attempt-reuse-'+uuid4().hex))
                    expected = origin['worker_receipt']['arrays']['frame']
                    frame = _link(cached['frame']['path'], attempt/'frame.npy', expected,
                                  descriptor['slots'][index]['shape'])
                    receipt = dict(complete=True, request=key, index=index, arrays=dict(frame=frame),
                        execution='adopted-existing-frame; no current worker/cost/PID',
                        reuse_provenance=dict(adoption=provenance_binding,
                            cache_item=file_binding(cache/('item-%04d' % index)/'complete.json'),
                            original_worker_provenance=origin))
                    receipts.append(_publish(item/'complete.json', receipt))
                completed.append(dict(request=request_binding, complete=True, count=len(receipts),
                                      receipts=receipts))
            result = dict(namespace, complete=True, requests=completed,
                pinned_storage=_storage(cached['frame']['path']
                    for _, _, items, _, _ in admitted for cached, _ in items))
            _publish(target/'adopted.json', result)
            return result


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest='command', required=True)
    pin_parser = commands.add_parser('pin', help='Pin available completed items; safely repeat to fill missing items.')
    pin_parser.add_argument('--source-request', type=Path, required=True)
    pin_parser.add_argument('--source-identity', required=True, help='Explicit ORIGINAL producer identity.')
    pin_parser.add_argument('--expected-input', type=Path, required=True, help='JSON q-independent actual stage_identity object.')
    pin_parser.add_argument('--count', type=int, required=True)
    pin_parser.add_argument('--history', type=Path, required=True, help='Original completed actual mean boundary.')
    pin_parser.add_argument('--history-sha256', required=True)
    pin_parser.add_argument('--cache', type=Path, required=True, help='Distinct owned cache; parent must exist.')
    pin_parser.add_argument('--owner', required=True, help='Explicit cache owner token, reused during adoption.')
    adopt_parser = commands.add_parser('adopt', help='Adopt COMPLETE caches into an absent future frames namespace.')
    adopt_parser.add_argument('--cache', type=Path, action='append', required=True)
    adopt_parser.add_argument('--expected-request', type=Path, action='append', required=True,
        help='Full NEW QWork v2 JSON request, paired in order with --cache; repeat both as needed.')
    adopt_parser.add_argument('--producer-identity', required=True,
        help='Explicit physical producer identity; may equal the original identity.')
    adopt_parser.add_argument('--target', type=Path, required=True, help='NEW scientific-work/frames directory; parent must exist.')
    adopt_parser.add_argument('--owner', required=True)
    args = parser.parse_args()
    try:
        if args.command == 'pin':
            history = file_binding(_safe(args.history))
            if history['sha256'] != args.history_sha256:
                raise ValueError('completed actual history digest differs')
            result = pin(args.source_request, args.cache, owner=args.owner,
                source_identity=args.source_identity, expected_input=_json(args.expected_input),
                count=args.count, history_binding=history)
        else:
            result = adopt(args.cache, args.target, owner=args.owner,
                producer_identity=args.producer_identity,
                expected_requests=[_json(path) for path in args.expected_request])
        print(json.dumps(result, indent=2, allow_nan=False), flush=True)
        return 0
    except Exception as error:
        print(json.dumps(dict(status='refused-or-interrupted', error=type(error).__name__+': '+str(error),
            retained_owned_state=True, source_deleted=False)), file=sys.stderr, flush=True)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())

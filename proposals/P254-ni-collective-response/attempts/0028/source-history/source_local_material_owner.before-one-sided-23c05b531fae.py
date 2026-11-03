"""Immutable ownership and durable budgets for actual source-local q work.

Only this owner's isolated directory is writable. Original physical inputs and
source-local pullback contexts retain their own authenticated scientific meaning.
"""
from contextlib import contextmanager
import json
import math
import os
from pathlib import Path
import signal
from time import time
from uuid import uuid4

import numpy as np

from mixed_material_driver_support import (
    ContinuationBudget, canonical, exact_publish, owned_tree)
from substrate_framework.partitioned_joint_slots import digest
from substrate_framework.streamed_material_state import authenticate, file_binding

SCHEMA = 'P254-source-local-q-material-v1'
TERMINAL_BYTES = 2*1024**2
CAP_NAMES = ('max_actions', 'max_seconds', 'max_output_bytes',
             'max_memory_bytes', 'worker_memory_bytes')


def encode(value):
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(item) for item in value]
    if isinstance(value, np.ndarray):
        return encode(value.tolist())
    if isinstance(value, np.generic):
        return encode(value.item())
    if isinstance(value, complex):
        return [value.real, value.imag]
    if isinstance(value, Path):
        return str(value.resolve())
    return value


def execution_caps(args):
    if (not args.cpus or len(set(args.cpus)) != len(args.cpus)
            or not set(args.cpus) <= set(range(6, 11))
            or not set(args.cpus) <= os.sched_getaffinity(0)):
        raise ValueError('unique explicitly available CPUs6-10 required')
    if (args.max_actions < 1 or not math.isfinite(args.max_seconds)
            or args.max_seconds <= 0 or args.max_output_bytes < 1
            or not 0 < args.worker_memory_bytes <= args.max_memory_bytes):
        raise ValueError('positive explicit action/wall/output/aggregate/worker caps required')
    return {name: getattr(args, name) for name in CAP_NAMES}


def isolated_output(args, repository):
    output = args.output.absolute()
    if any(part.is_symlink() for part in (output, *output.parents)):
        raise ValueError('owned output cannot traverse a symlink')
    output = output.resolve()
    repository = repository.resolve()
    if (output.is_relative_to(repository) or repository.is_relative_to(output)):
        raise ValueError('source-local output must be outside and disjoint from the repository')
    borrowed = [args.initialization, args.history, args.physical_qualification,
                args.source_local_oracle]
    borrowed.extend(Path(path) for path, _ in args.frame_admission)
    if any(path.resolve().is_relative_to(output) for path in borrowed):
        raise ValueError('owned output cannot contain borrowed physical inputs or frame admissions')
    return output


def material_owner(directory, request, *, resume, started_unix):
    """Authenticate the first request and start; resuming cannot refresh either."""
    path = Path(directory)/'owner.json'
    authority_path = Path(directory)/'owner-authority.json'
    request = canonical(encode(request))
    if resume:
        binding = file_binding(path)
        owner = json.loads(authenticate(binding).read_text())
        authority = json.loads(authority_path.read_text())
        if authority != dict(schema=SCHEMA+'-owner-authority', owner=binding):
            raise ValueError('source-local owner bytes changed; original clock cannot be replaced')
        if (owner.get('schema') != SCHEMA+'-owner'
                or owner.get('request') != request
                or type(owner.get('started_unix')) not in (int, float)
                or not math.isfinite(owner['started_unix'])
                or owner['started_unix'] > time()
                or not isinstance(owner.get('material_id'), str)):
            raise ValueError('resume requires identical scientific request, sources and caps')
    else:
        owner = dict(schema=SCHEMA+'-owner', material_id=uuid4().hex,
                     request=request, started_unix=started_unix)
        binding = exact_publish(path, owner, None)
        exact_publish(authority_path, dict(schema=SCHEMA+'-owner-authority', owner=binding), None)
    return owner, binding, digest(dict(schema=SCHEMA, owner=binding))


def material_budget(directory, owner, caps):
    existing = owned_tree(directory)
    if any(Path(name).stat().st_nlink != 1 for name in existing):
        raise ValueError('owned source-local output cannot contain borrowed hardlink aliases')
    return ContinuationBudget(directory, started_at=owner['started_unix'],
        max_actions=caps['max_actions'], max_seconds=caps['max_seconds'],
        max_bytes=caps['max_output_bytes'], max_memory_bytes=caps['max_memory_bytes'],
        worker_memory_bytes=caps['worker_memory_bytes'], existing_bytes=sum(existing.values()))


def reserve_terminal_and_ledger(budget):
    remaining = max(0, budget.max_actions*budget.RECORD.size-budget.ledger.stat().st_size)
    budget.reserve(TERMINAL_BYTES+remaining)


def context_bindings(args, directory):
    """Authenticate owned restart hints without changing the scientific request."""
    result = {}
    for side in ('left', 'right'):
        path = getattr(args, side+'_context')
        expected = getattr(args, side+'_context_sha256')
        if (path is None) != (expected is None):
            raise ValueError('both --'+side+'-context and its SHA256 are required')
        if path is not None:
            if args.mode != 'pair':
                raise ValueError('source-local contexts are meaningful only in pair mode')
            if not args.resume:
                raise ValueError('explicit source-local contexts require the original --resume owner')
            expected_path = Path(directory)/'scientific-work'/'pair'/side/'complete.json'
            if path.resolve() != expected_path.resolve():
                raise ValueError('source-local context must belong to this persistent pair namespace')
            binding = file_binding(path.resolve())
            if binding['sha256'] != expected:
                raise ValueError('explicit '+side+' source-local context digest differs')
            record = json.loads(authenticate(binding).read_text())
            if (record.get('schema') != 'source-local-q-pullback-v1'
                    or record.get('complete') is not True):
                raise ValueError('complete source-local pullback required, not a mapped-Q checkpoint')
            result[side+'_context'] = binding
    return result


@contextmanager
def terminal_signals():
    previous = {}
    def interrupted(signum, frame):
        raise InterruptedError('source-local material interrupted by signal '+str(signum))
    try:
        for number in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
            previous[number] = signal.signal(number, interrupted)
        yield
    finally:
        for number, handler in previous.items():
            signal.signal(number, handler)


def error_record(error):
    return dict(type=type(error).__name__, message=str(error)[:16384])


def finish_receipt(directory, invocation, result, budget):
    """A reserved terminal receipt must not consult an exhausted action/deadline cap."""
    if budget is not None:
        try:
            result['budget'] = budget.snapshot()
        except BaseException as error:
            result['budget_snapshot_error'] = error_record(error)
    try:
        result['measured_output_file_bytes_before_terminal'] = sum(owned_tree(directory).values())
        stat = os.statvfs(directory)
        result['filesystem_free_bytes_before_terminal'] = stat.f_bavail*stat.f_frsize
    except BaseException as error:
        result['storage_measurement_error'] = error_record(error)
    value = encode(result)
    data = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode()
    if len(data) > TERMINAL_BYTES:
        raise RuntimeError('source-local terminal receipt exceeds reserved bytes')
    return exact_publish(Path(directory)/'invocations'/invocation/'receipt.json', value, None)

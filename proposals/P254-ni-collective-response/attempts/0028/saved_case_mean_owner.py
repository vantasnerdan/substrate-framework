"""Original-clock ownership and existing durable budgets for saved-case means."""
from contextlib import contextmanager
from hashlib import sha256
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import signal
from time import time

from mixed_material_driver_support import (
    ContinuationBudget, canonical, exact_publish, owned_tree)
from substrate_framework.streamed_material_state import authenticate, file_binding

SCHEMA = 'P254-saved-case-mean-owner-v1'
TERMINAL_BYTES = 2*1024**2
CAP_NAMES = ('max_actions', 'max_seconds', 'max_output_bytes',
             'max_memory_bytes', 'worker_memory_bytes')
_fork_terminal_scope = False


def _reset_fork_signals():
    if _fork_terminal_scope:
        for number in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGALRM):
            signal.signal(number, signal.SIG_DFL)
        signal.setitimer(signal.ITIMER_REAL, 0)


os.register_at_fork(after_in_child=_reset_fork_signals)


def reap_owned_workers():
    """Do not publish a terminal or release ownership with live slot children."""
    children = mp.active_children()
    for child in children:
        if child.is_alive():
            child.kill()
    for child in children:
        child.join()
    return [dict(pid=str(child.pid), exitcode=child.exitcode) for child in children]


def isolated_output(args, repository):
    output = args.output.absolute()
    if any(part.is_symlink() for part in (output, *output.parents)):
        raise ValueError('owned output cannot traverse a symlink')
    output = output.resolve()
    repository = repository.resolve()
    if output.is_relative_to(repository) or repository.is_relative_to(output):
        raise ValueError('saved mean output must be outside and disjoint from the repository')
    borrowed = [args.initialization, args.qualification, args.oracle, args.column_oracle]
    if args.reuse_initial_contractions:
        borrowed.append(args.reuse_initial_contractions)
    if any(path.resolve().is_relative_to(output) for path in borrowed):
        raise ValueError('owned output cannot contain original physical inputs or proofs')
    return output


def execution_caps(args):
    if (type(args.max_actions) is not int or args.max_actions < 1
            or not math.isfinite(args.max_seconds) or args.max_seconds <= 0
            or args.max_output_bytes < 1
            or not 0 < args.worker_memory_bytes <= args.max_memory_bytes):
        raise ValueError('positive explicit action/wall/output/aggregate/worker allocations required')
    return {name: getattr(args, name) for name in CAP_NAMES}


def original_owner(directory, request, *, resume, started_unix):
    path, authority = directory/'owner.json', directory/'owner-authority.json'
    request = canonical(request)
    if resume:
        binding = file_binding(path)
        owner = json.loads(authenticate(binding).read_text())
        if (json.loads(authority.read_text()) != dict(schema=SCHEMA+'-authority', owner=binding)
                or owner.get('schema') != SCHEMA or owner.get('request') != request
                or type(owner.get('started_unix')) not in (int, float)
                or not math.isfinite(owner['started_unix']) or owner['started_unix'] > time()):
            raise ValueError('same original mean owner, sources, inputs, limits and caps required')
    else:
        owner = dict(schema=SCHEMA, request=request, started_unix=started_unix)
        data = (json.dumps(owner, indent=2, allow_nan=False)+'\n').encode()
        prospective = dict(path=str(path.resolve()), bytes=len(data), sha256=sha256(data).hexdigest())
        authority_data = (json.dumps(dict(schema=SCHEMA+'-authority', owner=prospective),
                                    indent=2, allow_nan=False)+'\n').encode()
        caps = request['caps']
        required = (sum(owned_tree(directory).values())+len(data)+len(authority_data)
                    +TERMINAL_BYTES+caps['max_actions']*ContinuationBudget.RECORD.size)
        if required > caps['max_output_bytes']:
            raise ValueError('output allocation cannot admit original metadata, terminal and action ledger')
        binding = exact_publish(path, owner, None)
        exact_publish(authority, dict(schema=SCHEMA+'-authority', owner=binding), None)
    return owner, binding


def continuation_budget(directory, owner):
    caps = owner['request']['caps']
    budget = ContinuationBudget(directory, started_at=owner['started_unix'],
        max_actions=caps['max_actions'], max_seconds=caps['max_seconds'],
        max_bytes=caps['max_output_bytes'], max_memory_bytes=caps['max_memory_bytes'],
        worker_memory_bytes=caps['worker_memory_bytes'], existing_bytes=sum(owned_tree(directory).values()))
    remaining = max(0, budget.max_actions*budget.RECORD.size-budget.ledger.stat().st_size)
    budget.reserve(TERMINAL_BYTES+remaining)
    budget.saved_case_terminal_reserved = True
    return budget


@contextmanager
def terminal_signals():
    global _fork_terminal_scope
    previous_scope = _fork_terminal_scope
    _fork_terminal_scope = True
    previous = {}
    def interrupted(signum, frame):
        signal.setitimer(signal.ITIMER_REAL, 0)
        raise InterruptedError('saved-case mean interrupted by signal '+str(signum))
    try:
        for number in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
            previous[number] = signal.signal(number, interrupted)
        yield
    finally:
        _fork_terminal_scope = previous_scope
        for number, handler in previous.items():
            signal.signal(number, handler)


def finish_receipt(directory, invocation, result, budget):
    if budget is None or not getattr(budget, 'saved_case_terminal_reserved', False):
        raise RuntimeError('terminal output was not admitted; report on stderr only')
    result['budget'] = budget.snapshot()
    result['finished_unix'] = time()
    data = (json.dumps(result, indent=2, allow_nan=False)+'\n').encode()
    if len(data) > TERMINAL_BYTES:
        raise RuntimeError('saved-case mean terminal report exceeds reserved bytes')
    return exact_publish(directory/'invocations'/invocation/'receipt.json', result, None)

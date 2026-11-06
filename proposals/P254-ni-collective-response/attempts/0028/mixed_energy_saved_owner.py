"""Independent saved-energy ownership, exclusion and observed cost accounting.

Producer files are borrowed read-only. Only this consumer's owner, durable action
ledger, isolated work and immutable invocation receipts belong to its budget.
"""
from contextlib import contextmanager
import fcntl
import json
import math
import os
from pathlib import Path
import signal
from time import time
from uuid import uuid4

from mixed_material_driver_support import (
    ContinuationBudget, canonical, exact_publish, owned_tree)
from substrate_framework.partitioned_joint_slots import digest
from substrate_framework.streamed_material_state import authenticate, file_binding

SCHEMA = 'P254-saved-native-mixed-energy-consumer-v1'
TERMINAL_BYTES = 2*1024**2
CAP_NAMES = ('max_actions', 'max_seconds', 'max_output_bytes',
             'max_memory_bytes', 'worker_memory_bytes')


def execution_caps(args):
    """Retain the existing material executor's resource ceilings."""
    if (len(set(args.cpus)) != len(args.cpus)
            or not set(args.cpus) <= os.sched_getaffinity(0)):
        raise ValueError('unique explicitly available CPUs6-10 required')
    if (args.max_actions < 1 or not math.isfinite(args.max_seconds)
            or args.max_seconds <= 0 or args.max_output_bytes < 1
            or not 0 < args.worker_memory_bytes <= args.max_memory_bytes):
        raise ValueError('positive explicit action/wall/live-output/aggregate/worker caps required')
    if (args.max_memory_bytes > 18*1024**3
            or args.worker_memory_bytes > 4*1024**3
            or args.max_output_bytes > 150*1024**3):
        raise ValueError('original aggregate18GiB/worker4GiB/live150GiB ceilings cannot be enlarged')
    return {name: getattr(args, name) for name in CAP_NAMES}


def isolated_output(args, repository):
    """Reject overlap before creating anything, including symlink traversal."""
    output = args.output.absolute()
    if any(part.is_symlink() for part in (output, *output.parents)):
        raise ValueError('owned output cannot traverse a symlink')
    output = output.resolve()
    producer, state = args.producer.resolve(), args.state.resolve()
    if (output == repository or output.is_relative_to(repository)
            or repository.is_relative_to(output)
            or output.is_relative_to(producer) or producer.is_relative_to(output)
            or state.is_relative_to(output)):
        raise ValueError('consumer output must be outside the repository and disjoint from producer/input')
    return output


def consumer_owner(directory, request, *, resume, started_unix):
    """Persist the first start before admitting work; a replay cannot refresh it."""
    path = Path(directory)/'owner.json'
    authority_path = Path(directory)/'owner-authority.json'
    request = canonical(request)
    if resume:
        binding = file_binding(path)
        owner = json.loads(authenticate(binding).read_text())
        authority = json.loads(authority_path.read_text())
        if authority != dict(schema=SCHEMA+'-owner-authority', owner=binding):
            raise ValueError('consumer owner bytes changed; original identity/clock cannot be replaced')
        if (owner.get('schema') != SCHEMA+'-owner'
                or owner.get('request') != request
                or type(owner.get('started_unix')) not in (int, float)
                or not math.isfinite(owner['started_unix'])
                or owner['started_unix'] > time()
                or not isinstance(owner.get('consumer_id'), str)):
            raise ValueError('resume requires identical consumer request/caps/sources and original start')
    else:
        owner = dict(schema=SCHEMA+'-owner', consumer_id=uuid4().hex,
                     request=request, started_unix=started_unix)
        binding = exact_publish(path, owner, None)
        exact_publish(authority_path, dict(schema=SCHEMA+'-owner-authority', owner=binding), None)
    return owner, binding, digest(dict(schema=SCHEMA, owner=binding))


def consumer_budget(directory, owner, caps):
    existing = owned_tree(directory)
    for name in existing:
        if Path(name).stat().st_nlink != 1:
            raise ValueError('owned consumer output cannot contain borrowed hardlink aliases')
    budget = ContinuationBudget(directory, started_at=owner['started_unix'],
        max_actions=caps['max_actions'], max_seconds=caps['max_seconds'],
        max_bytes=caps['max_output_bytes'], max_memory_bytes=caps['max_memory_bytes'],
        worker_memory_bytes=caps['worker_memory_bytes'], existing_bytes=sum(existing.values()))
    return budget


def reserve_terminal_and_ledger(budget):
    remaining = max(0, budget.max_actions*budget.RECORD.size-budget.ledger.stat().st_size)
    budget.reserve(TERMINAL_BYTES+remaining)


@contextmanager
def producer_exclusion(producer_directory, state):
    """Hold the actual existing producer lock before positive-time admission.

    The initial immutable physical boundary is a deliberately read-only smoke:
    it neither takes away the live producer's lock nor borrows its writable work.
    """
    if state['kind'] == 'physical-local-pulse-raw-electronic-zero' and state['time'] == 0:
        yield dict(policy='isolated read-only initial physical-boundary engineering readout',
                   producer_may_be_running=True, positive_time_claim=False)
        return
    if state['kind'] != 'mixed-rk4-physical-endpoint' or state['time'] <= 0:
        raise ValueError('unrecognized physical-boundary exclusion policy')
    path = Path(producer_directory)/'owner.lock'
    if path.is_symlink() or not path.is_file():
        raise ValueError('existing authentic producer owner.lock required for positive-time readout')
    with path.open('rb') as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError('positive-time saved readout refuses a competing producer/writer') from error
        yield dict(policy='existing producer owner.lock exclusive nonblocking; held through consumer exit',
                   lock_path=str(path), producer_may_be_running=False, positive_time_claim=True)


def producer_cost_snapshot(admission):
    """Observe the original producer's clock/admissions without constructing its budget."""
    directory = Path(admission['producer_directory'])
    owner = admission['owner']
    caps = owner['request']['caps']
    observed = time()
    result = dict(owner=admission['owner_binding'], caps=caps,
        original_started_unix=owner['started_unix'],
        original_deadline_unix=owner['started_unix']+caps['max_seconds'],
        observed_at_unix=observed,
        elapsed_since_original_start_seconds=max(0., observed-owner['started_unix']),
        attribution='historical producer only; never charged to or restarted by this consumer',
        action_ledger=None, immutable_invocation_receipts=[])
    ledger, lock_path = directory/'action-admissions.bin', directory/'action-admissions.lock'
    if ledger.exists():
        if ledger.is_symlink() or lock_path.is_symlink() or not lock_path.is_file():
            raise ValueError('producer action-ledger ownership is not authentic')
        with lock_path.open('rb') as lock:
            fcntl.flock(lock.fileno(), fcntl.LOCK_SH)
            binding = file_binding(ledger)
            size = binding['bytes']
            if size % ContinuationBudget.RECORD.size:
                raise ValueError('producer durable action ledger is incomplete')
            actions, columns = 0, 0
            if size:
                with authenticate(binding).open('rb') as stream:
                    stream.seek(-ContinuationBudget.RECORD.size, os.SEEK_END)
                    actions, columns = ContinuationBudget.RECORD.unpack(stream.read(ContinuationBudget.RECORD.size))
            result['action_ledger'] = dict(binding=binding, admitted_actions=actions,
                admitted_columns=columns, accounting='original durable pre-action admissions; not completed work')
    for path in sorted((directory/'invocations').glob('*/receipt.json')):
        binding = file_binding(path)
        receipt = json.loads(authenticate(binding).read_text())
        if receipt.get('owner') == admission['owner_binding']:
            result['immutable_invocation_receipts'].append(dict(binding=binding,
                status=receipt.get('status'), observed_budget=receipt.get('budget'),
                invocation_wall_seconds=receipt.get('invocation_wall_seconds'),
                invocation_parent_cpu_seconds=receipt.get('invocation_parent_cpu_seconds')))
    result['missing_measurements'] = (
        'A running producer may have no terminal receipt. Ledger admissions are observed, '
        'not inferred completed actions, worker CPU, RSS, or device I/O.')
    return result


def execution_snapshot(directory):
    """Path identities suffice: QWork complete and execution receipts are immutable."""
    directory = Path(directory)
    return {str(path) for pattern in ('execution-*.json', 'complete.json')
            for path in directory.rglob(pattern)}


def execution_costs(directory, before):
    """Separate newly measured worker costs from adopted older consumer receipts."""
    directory = Path(directory)
    current, adopted, journals = [], [], []
    for path in sorted(directory.rglob('execution-*.json')):
        binding = file_binding(path)
        journals.append(dict(binding=binding, current_invocation=str(path) not in before))
    for path in sorted(directory.rglob('complete.json')):
        binding = file_binding(path)
        value = json.loads(authenticate(binding).read_text())
        if not {'cpu_seconds', 'hamiltonian_actions', 'process_peak_rss_bytes'} <= value.keys():
            continue
        record = dict(binding=binding, index=value.get('index'), cpu_seconds=value['cpu_seconds'],
            hamiltonian_actions=value['hamiltonian_actions'],
            process_peak_rss_bytes=value['process_peak_rss_bytes'])
        (current if str(path) not in before else adopted).append(record)
    return dict(current_invocation_completed_worker_receipts=current,
        adopted_prior_consumer_completed_worker_receipts=adopted, execution_journals=journals,
        current_completed_worker_cpu_seconds=sum(row['cpu_seconds'] for row in current),
        current_completed_worker_hamiltonian_actions=sum(row['hamiltonian_actions'] for row in current),
        scope='Only new immutable completion receipts are current completed worker cost; '
              'failed/killed/incomplete work stays charged in the durable budget and journals, '
              'but has no fabricated CPU measurement. Producer work is never included.')


def energy_attempt_evidence(directory):
    """Bind the real API's request/results/failure attempts even when it raises."""
    directory = Path(directory)
    paths = [path for path in (directory/'request.json', directory/'result.json') if path.is_file()]
    paths.extend(sorted(directory.glob('attempt-*.json')))
    result = []
    for path in paths:
        binding = file_binding(path)
        value = json.loads(authenticate(binding).read_text())
        result.append(dict(binding=binding, schema=value.get('schema'),
                           status=value.get('status'), error=value.get('error')))
    return result


@contextmanager
def terminal_signals():
    """Convert catchable termination into an honest immutable exit receipt."""
    previous = {}
    def interrupted(signum, frame):
        raise InterruptedError('consumer interrupted by signal '+str(signum))
    for number in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        previous[number] = signal.signal(number, interrupted)
    try:
        yield
    finally:
        for number, handler in previous.items():
            signal.signal(number, handler)


def terminal_receipt(directory, invocation, result):
    """Terminal metadata uses the precharged reserve, not an exhausted budget."""
    data = (json.dumps(result, indent=2, allow_nan=False)+'\n').encode()
    if len(data) > TERMINAL_BYTES:
        raise RuntimeError('consumer terminal receipt exceeds reserved bytes')
    return exact_publish(Path(directory)/'invocations'/invocation/'receipt.json', result, None)

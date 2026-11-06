"""Fork-shared action/storage/deadline accounting for the physical joint executor.

Caps apply to this explicitly allocated invocation; existing output consumes the
storage cap. Attempt journals retain earlier invocation costs, including failures.
The parent enforces wall/RSS caps even while workers are inside native libraries.
"""
import multiprocessing as mp
from pathlib import Path
from time import monotonic

from .streamed_joint_state import JointBudget, JointBudgetExceeded


def resident_bytes(pid):
    try:
        for line in Path('/proc/%d/status' % pid).read_text().splitlines():
            if line.startswith('VmRSS:'):
                return int(line.split()[1])*1024
    except FileNotFoundError:
        return 0
    return 0


class PartitionBudget:
    array = JointBudget.array
    json = JointBudget.json

    def __init__(self, *, max_actions, max_seconds, max_bytes, max_memory_bytes,
                 worker_memory_bytes, existing_bytes=0):
        validated = JointBudget(max_actions=max_actions, max_seconds=max_seconds, max_bytes=max_bytes)
        if (type(max_memory_bytes) is not int or type(worker_memory_bytes) is not int
                or not 0 < worker_memory_bytes <= max_memory_bytes
                or not 0 <= existing_bytes <= max_bytes):
            raise ValueError('positive aggregate/worker RSS and sufficient existing-output allocation required')
        self.max_actions, self.max_seconds, self.max_bytes = max_actions, max_seconds, max_bytes
        self.max_memory_bytes, self.worker_memory_bytes = max_memory_bytes, worker_memory_bytes
        self.started = validated.started
        ctx = mp.get_context('fork')
        self._counts = ctx.Array('q', [0, 0, existing_bytes], lock=True)
        self.peak_aggregate_rss = 0
        self.local_actions = self.local_columns = 0

    @property
    def actions(self):
        return self._counts[0]

    @property
    def columns(self):
        return self._counts[1]

    @property
    def bytes_reserved(self):
        return self._counts[2]

    def check(self):
        if monotonic()-self.started >= self.max_seconds:
            raise JointBudgetExceeded('invocation wall-time allocation exhausted')

    def action(self, vectors):
        self.check()
        with self._counts.get_lock():
            if self._counts[0] >= self.max_actions:
                raise JointBudgetExceeded('shared Hamiltonian action allocation exhausted')
            self._counts[0] += 1
            self._counts[1] += 1 if vectors.ndim == 1 else vectors.shape[1]
        self.local_actions += 1
        self.local_columns += 1 if vectors.ndim == 1 else vectors.shape[1]

    def reserve(self, size):
        self.check()
        with self._counts.get_lock():
            if self._counts[2]+size > self.max_bytes:
                raise JointBudgetExceeded('shared output allocation exhausted')
            self._counts[2] += size

    def check_processes(self, parent, children):
        self.check()
        sizes = [resident_bytes(pid) for pid in children]
        total = resident_bytes(parent)+sum(sizes)
        self.peak_aggregate_rss = max(self.peak_aggregate_rss, total)
        if total > self.max_memory_bytes or any(n > self.worker_memory_bytes for n in sizes):
            raise JointBudgetExceeded('RSS allocation exhausted (conservative sum includes shared pages)')

    def snapshot(self):
        # A watchdog may kill a child while it owns the semaphore. Terminal
        # reads must not acquire that abandoned lock; counters are machine-word
        # integers and reservations precede every actual action/write.
        actions, columns, size = self._counts.get_obj()[:]
        return dict(hamiltonian_actions=actions, action_columns=columns,
                    output_bytes_reserved=size, elapsed_seconds=monotonic()-self.started,
                    max_actions=self.max_actions, max_seconds=self.max_seconds, max_bytes=self.max_bytes,
                    max_memory_bytes=self.max_memory_bytes, worker_memory_bytes=self.worker_memory_bytes,
                    peak_aggregate_rss_bytes=self.peak_aggregate_rss,
                    accounting='invocation-wide shared counters; existing output charged; attempts retained',
                    memory_accounting='parent plus child RSS; shared pages conservatively counted repeatedly',
                    time_cap='parent watchdog at 0.1 second intervals; native calls may be terminated')

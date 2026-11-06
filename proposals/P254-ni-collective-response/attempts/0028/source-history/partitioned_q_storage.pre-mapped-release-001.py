"""Q-only live reservations, durable barriers and payload retirement.

Retirement is parent-only, after workers have joined. Receipts and source-bound
checkpoints are never treated as mutable scratch. Existing scientific artifacts
outside the explicitly owned work requests are not eligible for retirement.
"""
import json
import multiprocessing as mp
import os
from pathlib import Path

from .partitioned_joint_budget import PartitionBudget
from .partitioned_joint_slots import publish
from .streamed_joint_state import JointBudgetExceeded
from .streamed_material_state import authenticate, file_binding


def sync_file(path):
    with Path(path).open('rb') as stream:
        os.fsync(stream.fileno())


def sync_directory(path):
    fd = os.open(Path(path), os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def durable_publish(path, value, budget):
    binding = publish(path, value, budget)
    sync_file(path)
    # Newly created nested work directories must themselves survive a crash.
    for parent in Path(path).resolve().parents:
        sync_directory(parent)
    return binding


class QBudget(PartitionBudget):
    """The byte cap is live storage, not cumulative write traffic.

    Reservations include conservative NPY headers and pending terminal receipts.
    Only measured unlinked bytes are credited back, so header slack remains
    charged. Action/column costs and cumulative reservations never decrease.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        existing = kwargs.get('existing_bytes', 0)
        self._storage = mp.get_context('fork').Array('q', [existing, 0, 0, 0], lock=True)

    def reserve(self, size):
        self.check()
        with self._counts.get_lock(), self._storage.get_lock():
            if self._counts[2]+size > self.max_bytes:
                raise JointBudgetExceeded('q live output allocation exhausted')
            self._counts[2] += size
            self._storage[0] = max(self._storage[0], self._counts[2])
            self._storage[1] += size

    def release(self, size):
        with self._counts.get_lock(), self._storage.get_lock():
            if not 0 <= size <= self._counts[2]:
                raise ValueError('invalid q retirement accounting')
            self._counts[2] -= size
            self._storage[2] += size

    def record_write(self, size):
        """Charge planned mapped overwrites without reserving another live file."""
        with self._storage.get_lock():
            self._storage[3] += size

    def snapshot(self):
        result = super().snapshot()
        peak, allocated, retired, overwritten = self._storage.get_obj()[:]
        result.update(live_bytes_reserved=result['output_bytes_reserved'],
            peak_live_bytes_reserved=peak, cumulative_write_bytes_reserved=allocated+overwritten,
            cumulative_allocation_bytes_reserved=allocated, mapped_overwrite_bytes_reserved=overwritten,
            retired_payload_bytes=retired,
            accounting='q live cap; cumulative allocation/mapped-write/action costs never credited back; measured unlinks only; header/reservation slack retained; logical writes, not physical device I/O')
        return result


def retire_payloads(receipt_path, barrier, paths, owner, budget):
    """Write-ahead retirement; replay safely after any interrupted unlink.

    The barrier is authenticated before every replay. Paths are restricted to
    this owned request; compact attempts/completion/retirement receipts survive.
    A missing payload is accepted ONLY after the durable intent exists.
    """
    receipt_path, owner = Path(receipt_path), Path(owner).resolve()
    authenticate(barrier)
    if receipt_path.exists():
        record = json.loads(receipt_path.read_text())
        if record['barrier'] != barrier or record['owner'] != str(owner):
            raise ValueError('retirement barrier/owner differs')
    else:
        entries = []
        for path in paths:
            path = Path(path)
            if path.is_symlink():
                raise ValueError('retirement cannot follow a symlink')
            path = path.resolve()
            if not path.is_relative_to(owner) or path == owner:
                raise ValueError('retirement outside owned q work')
            entries.append(file_binding(path))
        record = dict(schema='q-retirement-v1', barrier=barrier, owner=str(owner), payloads=entries)
        durable_publish(receipt_path, record, budget)
    for binding in record['payloads']:
        path = Path(binding['path'])
        if path.is_symlink() or not path.resolve().is_relative_to(owner):
            raise ValueError('retirement payload escaped owner')
        if path.exists():
            authenticate(binding)
            size = path.stat().st_size
            path.unlink()
            sync_directory(path.parent)
            budget.release(size)
    return file_binding(receipt_path)

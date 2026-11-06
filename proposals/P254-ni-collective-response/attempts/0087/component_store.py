"""A76 additive storage with durable per-cell ownership and authenticated resume.

The atomic checkpoint's packed bitmap is authoritative. A crash between bitmap
fsync and checkpoint rename can leave only an uncommitted tile ahead of it; that
tile is rolled back, never adopted without its checkpoint. Numeric payload and
per-cell digests are fsynced before any ownership publication.
"""
from pathlib import Path
import base64
import hashlib
import json
import fcntl
import os
import time

import numpy as np
from history_banks import ORDERS, authenticate, binding, save
from substrate_framework.streamed_material_state import authenticate as authenticate_source


def durable_save(path, record):
    save(path, record)
    fd = os.open(Path(path).parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def authenticate_tree(value):
    """Use the owning loader's exact integer/string byte-count descriptor parser."""
    if isinstance(value, dict):
        if {'path', 'bytes', 'sha256'} <= value.keys():
            authenticate_source(value)
        else:
            for child in value.values():
                authenticate_tree(child)
    elif isinstance(value, (tuple, list)):
        for child in value:
            authenticate_tree(child)


def _cell_digest(arrays, a, b, slot):
    digest = hashlib.sha256()
    for order in ORDERS:
        digest.update(np.ascontiguousarray(arrays[order][a, b, slot]).tobytes())
    return np.frombuffer(digest.digest(), dtype=np.uint8)


class ComponentStore:
    """One writer; output cells are assigned, never accumulated or refilled."""
    def __init__(self, output, record, *, resume=False, inputs=None):
        self.output = Path(output).resolve()
        self.record = record
        self.arrays, self.fds = {}, {}
        self.seen = self.digests = None
        self.poisoned = False
        self.closed = False
        self.shape = (31, 31, 128, 12, 12)
        self.coverage_shape = self.shape[:3]
        try:
            if not resume:
                self.output.mkdir(parents=True, exist_ok=False)
            self.fds['writer-lock'] = os.open(self.output / '.writer.lock', os.O_CREAT | os.O_RDWR, 0o600)
            fcntl.flock(self.fds['writer-lock'], fcntl.LOCK_EX | fcntl.LOCK_NB)
            if resume:
                latest = json.loads((self.output / 'checkpoint.json').read_bytes())
                if (latest['checkpoint_bitmap_base64'] != record['checkpoint_bitmap_base64']
                        or latest['input_bindings'] != record['input_bindings']):
                    raise ValueError('checkpoint advanced while preparing resume; invoke explicit resume again')
            mode = 'r+' if resume else 'w+'
            for order in ORDERS:
                path = self.output / ('sigma-' + order + '.npy')
                self.arrays[order] = np.lib.format.open_memmap(
                    path, mode=mode, dtype=np.complex128, shape=self.shape)
                self.fds[order] = os.open(path, os.O_RDWR)
                if not resume:
                    self.arrays[order][:] = 0
                if self.arrays[order].shape != self.shape or self.arrays[order].dtype != np.complex128:
                    raise ValueError('saved sigma shape/dtype differs')
            for name, shape, dtype in (
                    ('slot-coverage', self.coverage_shape, bool),
                    ('cell-digests', (*self.coverage_shape, 32), np.uint8)):
                path = self.output / (name + '.npy')
                arr = np.lib.format.open_memmap(path, mode=mode, dtype=dtype, shape=shape)
                self.fds[name] = os.open(path, os.O_RDWR)
                if arr.shape != shape or arr.dtype != np.dtype(dtype):
                    raise ValueError('saved ownership/digest axes differ')
                if not resume:
                    arr[:] = 0
                if name == 'slot-coverage':
                    self.seen = arr
                else:
                    self.digests = arr
            if resume:
                packed = base64.b64decode(record['checkpoint_bitmap_base64'], validate=True)
                expected_bytes = (int(np.prod(self.coverage_shape)) + 7) // 8
                if len(packed) != expected_bytes:
                    raise ValueError('checkpoint bitmap byte count differs')
                committed = np.unpackbits(np.frombuffer(packed, dtype=np.uint8))[:int(np.prod(self.coverage_shape))]
                committed = committed.reshape(self.coverage_shape).astype(bool)
                if np.any(committed & ~self.seen):
                    raise ValueError('checkpoint claims cells absent from durable bitmap')
                for a, b, slot in np.argwhere(committed):
                    if not np.array_equal(_cell_digest(self.arrays, a, b, slot), self.digests[a, b, slot]):
                        raise ValueError('changed committed scientific cell: %s' % ((int(a), int(b), int(slot)),))
                self.recovered_uncommitted_cells = int(np.count_nonzero(self.seen & ~committed))
                self.seen[:] = committed
                self._flush('slot-coverage', self.seen)
            else:
                self.recovered_uncommitted_cells = 0
                for order in ORDERS:
                    self._flush(order, self.arrays[order])
                self._flush('cell-digests', self.digests)
                self._flush('slot-coverage', self.seen)
                if inputs is None:
                    raise ValueError('fresh component requires its actual input bindings')
                durable_save(self.output / 'inputs.json', inputs)
                self.record['input_bindings'] = binding(self.output / 'inputs.json')
                self.checkpoint()
        except BaseException:
            self.close()
            raise

    def _flush(self, name, arr):
        arr.flush()
        os.fsync(self.fds[name])

    def checkpoint(self):
        self.record['checkpoint_bitmap_base64'] = base64.b64encode(
            np.packbits(np.asarray(self.seen).reshape(-1)).tobytes()).decode('ascii')
        self._coverage()
        durable_save(self.output / 'checkpoint.json', self.record)

    def _coverage(self):
        complete = [[a, b] for a in range(31) for b in range(31) if self.seen[a, b].all()]
        count = int(self.seen.sum())
        q, = self.record['owned_q']
        done = count == 31 * 31 * 128
        self.record['coverage'].update(completed_pairs=complete, completed_slot_cells=count)
        self.record.update(status=('complete' if done else 'partial') + '-owned-component-q%d' % q,
                           completed_q=[q] if done else [], full_time_square=done)

    def commit(self, pairs, external, values):
        if self.poisoned:
            raise RuntimeError('failed commit requires explicit resume')
        pairs, external = tuple(pairs), tuple(int(s) for s in external)
        if (not pairs or not external or len(set(pairs)) != len(pairs)
                or len(set(external)) != len(external)):
            raise ValueError('nonempty unique cell ownership required')
        for a, b in pairs:
            if self.seen[a, b, list(external)].any():
                raise ValueError('committed scientific cell must not be reevaluated')
        expected = (len(external), len(pairs), 12, 12)
        if set(values) != set(ORDERS) or any(
                value.shape != expected or value.dtype != np.complex128 or not np.isfinite(value).all()
                for value in values.values()):
            raise ValueError('complete finite complex128 output tile required')
        for j, (a, b) in enumerate(pairs):
            retarded = values['greater'][:, j] - values['lesser'][:, j] if a >= b else np.zeros_like(values['retarded'][:, j])
            if not np.array_equal(values['retarded'][:, j], retarded):
                raise ValueError('owning retarded convention differs')
        tick = time.monotonic()
        try:
            for j, (a, b) in enumerate(pairs):
                for order in ORDERS:
                    self.arrays[order][a, b, list(external)] = values[order][:, j]
                for slot in external:
                    self.digests[a, b, slot] = _cell_digest(self.arrays, a, b, slot)
            for order in ORDERS:
                self._flush(order, self.arrays[order])
            self._flush('cell-digests', self.digests)
            for a, b in pairs:
                self.seen[a, b, list(external)] = True
            self._flush('slot-coverage', self.seen)
            telemetry = self.record['runs'][-1]
            cells = len(pairs) * len(external)
            telemetry['commits'] += 1
            telemetry['committed_slot_cells'] += cells
            telemetry['commit_flush_seconds'] += time.monotonic() - tick
            telemetry['last_commit'] = dict(pairs=[list(p) for p in pairs], external=list(external),
                                          cells=cells, monotonic_seconds=time.monotonic())
            telemetry['last_commit']['run_elapsed_seconds'] = time.monotonic() - telemetry['started_monotonic']
            telemetry['checkpoint_elapsed_seconds'] = (
                telemetry.get('static_rebuild_seconds', 0.)
                + telemetry['last_commit']['run_elapsed_seconds'])
            telemetry['checkpoint_process_cpu_seconds'] = time.process_time() - telemetry['started_process_cpu']
            telemetry['elapsed_accounting'] = 'measured through last durable commit; uncommitted tail before abrupt kill is not inferred'
            self.checkpoint()
        except BaseException:
            # Never republish an in-memory bitmap after a failed publication.
            self.poisoned = True
            raise

    def publish(self):
        if self.poisoned:
            return None
        self.checkpoint()
        self.record.pop('saved_catalog', None)
        descriptors = {}
        for order in ORDERS:
            self._flush(order, self.arrays[order])
            descriptors[order] = dict(binding(self.output / ('sigma-' + order + '.npy')),
                                      shape=list(self.shape), dtype='complex128')
        self._flush('slot-coverage', self.seen)
        self._flush('cell-digests', self.digests)
        self.record['sigma'] = descriptors
        self.record['coverage']['array'] = dict(binding(self.output / 'slot-coverage.npy'),
                                                shape=list(self.coverage_shape), dtype='bool')
        self.record['cell_digests'] = dict(binding(self.output / 'cell-digests.npy'),
                                         shape=[*self.coverage_shape, 32], dtype='uint8')
        durable_save(self.output / 'catalog.json', self.record)
        self.record['saved_catalog'] = binding(self.output / 'catalog.json')
        self.checkpoint()
        return self.output / 'catalog.json'

    def close(self):
        if self.closed:
            return
        for fd in self.fds.values():
            os.close(fd)
        self.fds.clear()
        for arr in (*self.arrays.values(), self.seen, self.digests):
            if arr is not None:
                arr._mmap.close()
        self.arrays.clear()
        self.closed = True

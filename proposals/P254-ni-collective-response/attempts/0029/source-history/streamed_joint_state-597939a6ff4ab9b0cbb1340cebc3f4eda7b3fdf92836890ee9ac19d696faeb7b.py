"""Authenticated new-model wrappers and immutable completed joint boundaries.

Original electronic frame/D bindings are reused, never relabeled as a scalar-
equivalent legacy restart. Trial stages are not adoptable states. Byte budgets
include incomplete output, and completed generations are never deleted.
"""
import json
from pathlib import Path
from time import monotonic

import numpy as np

from .streamed_material_state import DiskMean, StreamSlot, authenticate, file_binding, save_array


class JointBudgetExceeded(RuntimeError):
    """A cumulative allocation was exhausted; only prior boundaries are states."""


class JointBudget:
    def __init__(self, *, max_actions, max_seconds, max_bytes):
        if (type(max_actions) is not int or max_actions < 1 or
                type(max_bytes) is not int or max_bytes < 1 or
                not np.isfinite(max_seconds) or max_seconds <= 0):
            raise ValueError('positive cumulative action/time/storage allocation required')
        self.max_actions, self.max_seconds, self.max_bytes = max_actions, max_seconds, max_bytes
        self.actions = self.columns = self.bytes_reserved = 0
        self.started = monotonic()

    def check(self):
        if monotonic()-self.started >= self.max_seconds:
            raise JointBudgetExceeded('cumulative wall-time cap exhausted (cooperative action/slot boundaries)')

    def action(self, vectors):
        self.check()
        if self.actions >= self.max_actions:
            raise JointBudgetExceeded('cumulative Hamiltonian action cap exhausted')
        self.actions += 1
        self.columns += 1 if vectors.ndim == 1 else vectors.shape[1]

    def reserve(self, size):
        self.check()
        if self.bytes_reserved+size > self.max_bytes:
            raise JointBudgetExceeded('cumulative output byte cap exhausted')
        self.bytes_reserved += size

    def snapshot(self):
        return dict(hamiltonian_actions=self.actions, action_columns=self.columns,
                    output_bytes_reserved=self.bytes_reserved, elapsed_seconds=monotonic()-self.started,
                    max_actions=self.max_actions, max_seconds=self.max_seconds, max_bytes=self.max_bytes,
                    time_cap='cooperative, cannot interrupt a running BLAS/FFT/action')

    def array(self, path, value):
        value = np.asarray(value)
        # Numeric arrays only: the NPY v1/v2 header fits this conservative bound.
        if value.dtype.hasobject or value.ndim > 4:
            raise ValueError('numeric arrays with at most four axes required')
        self.reserve(value.nbytes+512)
        return save_array(path, value)

    def json(self, path, value):
        encoded = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode()
        self.reserve(len(encoded))
        path = Path(path)
        temporary = path.with_name(path.name+'.incomplete')
        with temporary.open('xb') as stream:
            stream.write(encoded)
        if path.exists():
            raise FileExistsError(path)
        temporary.rename(path)
        return file_binding(path)


class DiskJointMean:
    """Lazy full raw-column state plus homogeneous canonical ions/global field."""
    arrays = DiskMean.arrays

    def __init__(self, boundary, *, expected_identity=None, verify=True):
        self.path = Path(boundary).resolve()
        payload = json.loads(self.path.read_text())
        if payload['schema'] != 'streamed-joint-mean-v1' or payload['complete'] is not True:
            raise ValueError('not a completed joint boundary')
        self.identity = payload['model_identity']
        if not self.identity or (expected_identity is not None and self.identity != expected_identity):
            raise ValueError('joint model identity differs')
        self.time = float(payload['time_ry_inverse'])
        for name, key in [('u', 'u_bohr'), ('p', 'p_canonical'),
                          ('coordinate', 'field_coordinate'), ('momentum', 'field_momentum')]:
            value = np.asarray(payload[key], dtype=float)
            if value.ndim != 1 or not np.isfinite(value).all():
                raise ValueError('finite canonical coordinates required')
            value.setflags(write=False)
            setattr(self, name, value)
        self.slots = tuple(StreamSlot(**dict(s, shape=tuple(s['shape']))) for s in payload['slots'])
        self.origin, self.evidence = payload['original_electronic_boundary'], payload['evidence']
        self.history = payload['history']
        if (not np.isfinite(self.time) or self.p.shape != self.u.shape or self.u.size % 3
                or self.momentum.shape != self.coordinate.shape or not 1 <= self.coordinate.size <= 3
                or not self.slots or len({s.identity for s in self.slots}) != len(self.slots)):
            raise ValueError('incomplete joint coordinates or coverage')
        for s in self.slots:
            if (len(s.shape) != 2 or not 0 < s.shape[1] <= s.shape[0] or s.component < 0
                    or not np.isfinite(s.weight) or s.weight <= 0
                    or not all((s.identity, s.basis_identity, s.gauge_identity, s.source_identity))):
                raise ValueError('invalid original slot metadata')
        if verify:
            authenticate(self.origin)
            seen = set()
            bindings = [b for s in self.slots for b in (s.frame, s.coefficient)]
            for record in self.history:
                bindings.extend(record['arrays'].values())
                bindings.append(record['prior_boundary_binding'])
            for b in bindings:
                key = (b['path'], b['sha256'])
                if key not in seen:
                    authenticate(b)
                    seen.add(key)

    @classmethod
    def wrap_original(cls, original_binding, directory, *, identity, coordinate, momentum,
                      expected_slots, evidence, budget):
        original = DiskMean(authenticate(original_binding))
        if original.time != 0 or tuple(s.identity for s in original.slots) != tuple(expected_slots):
            raise ValueError('only authenticated complete original t0 frames may initialize new scalar')
        if identity == original.identity or not evidence.get('new_scalar') or not evidence.get('field_preparation'):
            raise ValueError('explicit DIFFERENT scalar identity and field preparation required')
        directory = Path(directory).resolve()
        directory.mkdir(parents=True, exist_ok=False)
        slots = [dict(s.metadata(), frame=s.frame, coefficient=s.coefficient) for s in original.slots]
        cls._publish(directory, slots, identity, 0., original.u, original.p, coordinate, momentum,
                     original_binding, [], dict(evidence, legacy_scalar_equivalence_claimed=False), budget)
        return cls(directory/'boundary.json', verify=False)

    @staticmethod
    def _publish(directory, slots, identity, time, u, p, coordinate, momentum, origin, history, evidence, budget):
        budget.json(directory/'boundary.json', dict(schema='streamed-joint-mean-v1', complete=True,
            model_identity=identity, time_ry_inverse=float(time), u_bohr=np.asarray(u).tolist(),
            p_canonical=np.asarray(p).tolist(), field_coordinate=np.asarray(coordinate).tolist(),
            field_momentum=np.asarray(momentum).tolist(), original_electronic_boundary=origin,
            slots=slots, history=history, evidence=evidence))

    def successor(self, directory, frames, *, time, u, p, coordinate, momentum, stage_record, evidence, budget):
        directory = Path(directory).resolve()
        directory.mkdir(parents=True, exist_ok=False)
        slots, iterator = [], iter(frames)
        for index, slot in enumerate(self.slots):
            budget.check()
            psi = np.asarray(next(iterator))
            if psi.shape != slot.shape or psi.dtype != np.complex128 or not np.isfinite(psi).all():
                raise ValueError('successor changed original full rows/columns')
            slots.append(dict(slot.metadata(), coefficient=slot.coefficient,
                frame=budget.array(directory/('frame-%04d.npy' % index), psi)))
        if next(iterator, None) is not None:
            raise ValueError('extra electronic sector')
        # The generator must finish every reduction before publication.
        evidence['budget_at_publication'] = budget.snapshot()
        self._publish(directory, slots, self.identity, time, u, p, coordinate, momentum,
                      self.origin, self.history+[stage_record], evidence, budget)
        return DiskJointMean(directory/'boundary.json', expected_identity=self.identity, verify=False)

"""Full raw-column mean checkpoints; no saved-band projection or ket/bra copy.

The JSON boundary is published last. Incomplete directories are not states.
NPY frames are mmap-readable; hashes bind source/basis/gauge and actual bytes.
No automatic checkpoint deletion, compression, or restart interpolation occurs.
"""
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


def file_binding(path):
    path = Path(path).resolve()
    digest = sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            digest.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def authenticate(binding):
    actual = file_binding(binding['path'])
    if actual['sha256'] != binding['sha256'] or actual['bytes'] != int(binding['bytes']):
        raise ValueError('checkpoint/source bytes differ: ' + binding['path'])
    return Path(binding['path'])


def save_array(path, value):
    with Path(path).open('xb') as stream:
        np.save(stream, value, allow_pickle=False)
    return file_binding(path)


@dataclass(frozen=True)
class StreamSlot:
    identity: str
    component: int
    weight: float
    shape: tuple
    basis_identity: str
    gauge_identity: str
    source_identity: str
    frame: dict
    coefficient: dict

    def metadata(self):
        return {name: getattr(self, name) for name in (
            'identity', 'component', 'weight', 'shape', 'basis_identity',
            'gauge_identity', 'source_identity')}


class DiskMean:
    """One immutable time boundary, with lazy per-slot full frames and D."""
    def __init__(self, boundary, *, expected_identity=None, verify=True):
        self.path = Path(boundary).resolve()
        payload = json.loads(self.path.read_text())
        if payload['schema'] != 'streamed-material-mean-v1':
            raise ValueError('unknown streamed boundary')
        self.identity = payload['model_identity']
        if not self.identity or (expected_identity is not None and self.identity != expected_identity):
            raise ValueError('different model/coverage identity')
        self.time = float(payload['time_ry_inverse'])
        self.u = np.asarray(payload['u_bohr'], dtype=float)
        self.p = np.asarray(payload['p_canonical'], dtype=float)
        self.slots = tuple(StreamSlot(**dict(s, shape=tuple(s['shape']))) for s in payload['slots'])
        self.evidence = payload['evidence']
        if (not self.slots or len({s.identity for s in self.slots}) != len(self.slots)
                or self.u.ndim != 1 or self.p.shape != self.u.shape
                or not np.isfinite(self.time) or not np.isfinite(self.u).all()
                or not np.isfinite(self.p).all()):
            raise ValueError('invalid complete state boundary')
        for s in self.slots:
            if (len(s.shape) != 2 or not 0 < s.shape[1] <= s.shape[0]
                    or s.component < 0 or not np.isfinite(s.weight) or s.weight <= 0
                    or not all((s.identity, s.source_identity, s.basis_identity, s.gauge_identity))):
                raise ValueError('invalid explicit full-sector identity')
        if verify:
            seen = set()
            for s in self.slots:
                for b in (s.frame, s.coefficient):
                    if (b['path'], b['sha256']) not in seen:
                        authenticate(b)
                        seen.add((b['path'], b['sha256']))

    def arrays(self, slot):
        psi = np.load(slot.frame['path'], mmap_mode='r', allow_pickle=False)
        d = np.load(slot.coefficient['path'], mmap_mode='r', allow_pickle=False)
        if psi.shape != slot.shape or d.shape != (slot.shape[1], slot.shape[1]):
            raise ValueError('checkpoint frame/D dimensions differ')
        if psi.dtype != np.complex128 or d.dtype != np.complex128:
            raise ValueError('full complex128 frame and original complex D required')
        return psi, d

    @classmethod
    def create(cls, directory, sectors, *, identity, time, u, p, evidence=None):
        """Consume (metadata, full Psi, original D) once, one sector at a time.

        Metadata must enumerate every explicit slot; caller supplies and binds
        its expected coverage. D files deduplicate only byte-identical matrices.
        """
        directory = Path(directory).resolve()
        directory.mkdir(parents=True, exist_ok=False)
        slots, coefficients = [], {}
        for index, (metadata, psi, d) in enumerate(sectors):
            psi, d = np.asarray(psi), np.asarray(d)
            if (psi.ndim != 2 or d.shape != (psi.shape[1], psi.shape[1])
                    or psi.dtype != np.complex128 or d.dtype != np.complex128
                    or not np.isfinite(psi).all() or not np.isfinite(d).all()):
                raise ValueError('finite full complex raw state required')
            key = sha256(np.ascontiguousarray(d).view(np.uint8)).hexdigest()
            if key not in coefficients:
                coefficients[key] = save_array(directory / ('D-' + key + '.npy'), d)
            record = dict(metadata, shape=list(psi.shape),
                          frame=save_array(directory / ('frame-%04d.npy' % index), psi),
                          coefficient=coefficients[key])
            slots.append(record)
        cls._publish(directory, slots, identity, time, u, p, evidence or {})
        return cls(directory / 'boundary.json', verify=False)

    @staticmethod
    def _publish(directory, slots, identity, time, u, p, evidence):
        payload = dict(schema='streamed-material-mean-v1', model_identity=identity,
                       time_ry_inverse=float(time), u_bohr=np.asarray(u).tolist(),
                       p_canonical=np.asarray(p).tolist(), slots=slots, evidence=evidence)
        with (directory / 'boundary.json').open('x') as stream:
            json.dump(payload, stream, indent=2, allow_nan=False)
            stream.write('\n')

    def successor(self, directory, frames, *, time, u, p, evidence):
        directory = Path(directory).resolve()
        directory.mkdir(parents=True, exist_ok=False)
        slots = []
        iterator = iter(frames)
        for index, slot in enumerate(self.slots):
            psi = np.asarray(next(iterator))
            if psi.shape != slot.shape or psi.dtype != np.complex128 or not np.isfinite(psi).all():
                raise ValueError('successor changed full-frame shape or contains nonfinite values')
            slots.append(dict(slot.metadata(), coefficient=slot.coefficient,
                              frame=save_array(directory / ('frame-%04d.npy' % index), psi)))
        if next(iterator, None) is not None:
            raise ValueError('extra successor sector')
        self._publish(directory, slots, self.identity, time, u, p, evidence)
        return DiskMean(directory / 'boundary.json', expected_identity=self.identity, verify=False)

    def requirements(self, steps, midpoint_iterations):
        frames = sum(16 * s.shape[0] * s.shape[1] for s in self.slots)
        largest = max(16 * s.shape[0] * s.shape[1] for s in self.slots)
        fixed = sum(int(b['bytes']) for b in {s.coefficient['path']: s.coefficient for s in self.slots}.values())
        return dict(sectors=len(self.slots), all_frame_bytes=frames,
                    largest_frame_bytes=largest, fixed_coefficient_file_bytes=fixed,
                    retained_frame_payload_bytes=(steps + 1) * frames,
                    maximum_exponential_sector_calls=steps * (midpoint_iterations + 1) * len(self.slots),
                    memory_scaling='O(max_s(npw_s*ncol_s+ncol_s²)+density_grid+operator_workspace); SciPy exponential workspace is additional, not a hard RSS cap',
                    storage_note='NPY headers, boundary JSON and diagnostic density grids are additional; all completed generations retained')

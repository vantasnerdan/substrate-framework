"""Specific E-I Sigma prefix and immutable Pi-q output; no joined-bank duplication."""
import json
import os
from pathlib import Path
import numpy as np
from reduced_projection_cache import ProjectionCache, atomic_json
from substrate_framework.streamed_material_state import authenticate, file_binding

CHANNELS = ('greater', 'lesser', 'retarded')


class EIOutput:
    def __init__(self, directory, specification, sigma_shape, pi_shape, nq):
        self.directory = Path(directory).resolve()
        self.directory.mkdir(parents=True, exist_ok=True)
        self.specification = specification
        self.sigma_shape, self.pi_shape = tuple(sigma_shape), tuple(pi_shape)
        self.nq, self.recoveries = int(nq), []
        self.spec_path = self.directory / 'specification.json'
        if self.spec_path.exists():
            if json.loads(self.spec_path.read_bytes()) != specification:
                raise ValueError('original source/state/time/implementation differs')
        else:
            atomic_json(self.spec_path, specification)
        self.pi = ProjectionCache(self.directory / 'pi-q', specification)

    @staticmethod
    def key(q):
        return 'q/%02d' % q

    def _validate(self, arrays, shape):
        if set(arrays) != set(CHANNELS) or any(
                x.shape != shape or x.dtype != np.dtype('complex128')
                or not np.isfinite(x).all() for x in arrays.values()):
            raise ValueError('literal ordered/retarded numeric axes/type/data differ')

    def pi_values(self, q):
        arrays = self.pi.get(self.key(q))
        self._validate(arrays, self.pi_shape)
        return arrays

    def restore(self):
        candidates = []
        for slot in range(2):
            sidecar = self.directory / ('sigma-%d.json' % slot)
            if not sidecar.exists():
                continue
            try:
                record = json.loads(sidecar.read_bytes())
                if (record['schema'] != 'p254-a48-sigma-prefix-v1'
                        or record['specification'] != file_binding(self.spec_path)
                        or not 1 <= record['completed_q'] <= self.nq):
                    raise ValueError('Sigma prefix identity/count differs')
                authenticate(record['archive'])
                candidates.append(record)
            except (OSError, ValueError, KeyError) as error:
                self.recoveries.append(dict(sidecar=str(sidecar), rejected=str(error)))
        for record in sorted(candidates, key=lambda x: x['completed_q'], reverse=True):
            try:
                with np.load(authenticate(record['archive']), allow_pickle=False) as saved:
                    arrays = {key: saved[key] for key in saved.files}
                self._validate(arrays, self.sigma_shape)
                for q in range(record['completed_q']):
                    self.pi_values(q)
                return arrays, record['completed_q'], record
            except (OSError, ValueError, KeyError) as error:
                self.recoveries.append(dict(archive=record['archive'], rejected=str(error)))
        if candidates or self.recoveries:
            raise ValueError('no authenticated complete Sigma/Pi prefix survives')
        return {key: np.zeros(self.sigma_shape, complex) for key in CHANNELS}, 0, None

    def commit(self, sigma, completed_q, pi_values, metadata):
        self._validate(sigma, self.sigma_shape)
        self._validate(pi_values, self.pi_shape)
        if not 1 <= completed_q <= self.nq:
            raise ValueError('prefix outside declared complete q range')
        q = completed_q - 1
        if self.pi.has(self.key(q)):
            existing = self.pi_values(q)
            if any(not np.array_equal(existing[k], pi_values[k]) for k in CHANNELS):
                raise ValueError('immutable original Pi-q values changed')
        else:
            self.pi.put(self.key(q), pi_values, metadata=metadata)
        slot = completed_q % 2
        archive = self.directory / ('sigma-%d.npz' % slot)
        stage = self.directory / 'sigma-stage.npz'
        with stage.open('wb') as stream:
            np.savez(stream, **sigma)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(stage, archive)
        record = dict(schema='p254-a48-sigma-prefix-v1', completed_q=completed_q,
                      specification=file_binding(self.spec_path), archive=file_binding(archive))
        atomic_json(self.directory / ('sigma-%d.json' % slot), record)
        return record

    def publish(self, record, result):
        if record['completed_q'] != self.nq or set(self.pi.manifest['items']) != {
                self.key(q) for q in range(self.nq)}:
            raise ValueError('whole original q coverage is incomplete')
        self.pi.manifest['status'] = 'complete'
        self.pi.publish()
        result.update(schema='p254-a48-full-external-EI-v1', status='complete',
                      specification=self.specification, sigma_prefix=record,
                      pi_manifest=file_binding(self.pi.path), recoveries=self.recoveries)
        atomic_json(self.directory / 'result.json', result)
        return result

    def consume(self):
        path = self.directory / 'result.json'
        if not path.exists():
            return None
        result = json.loads(path.read_bytes())
        if (result['schema'] != 'p254-a48-full-external-EI-v1'
                or result['status'] != 'complete'
                or result['specification'] != self.specification
                or result['sigma_prefix']['completed_q'] != self.nq
                or result['pi_manifest'] != file_binding(self.pi.path)
                or self.pi.manifest['status'] != 'complete'
                or set(self.pi.manifest['items']) != {self.key(q) for q in range(self.nq)}):
            raise ValueError('complete original Sigma/Pi identity/coverage differs')
        record = result['sigma_prefix']
        with np.load(authenticate(record['archive']), allow_pickle=False) as saved:
            sigma = {key: saved[key] for key in saved.files}
        self._validate(sigma, self.sigma_shape)
        sigma_norms = {key: float(np.linalg.norm(value)) for key, value in sigma.items()}
        pi_norms = [{key: float(np.linalg.norm(value)) for key, value in self.pi_values(q).items()}
                    for q in range(self.nq)]
        if sigma_norms != result['sigma_norms'] or pi_norms != result['pi_norms']:
            raise ValueError('actual saved Sigma/Pi decoded values differ')
        return dict(record=result, sigma=sigma, sigma_norms=sigma_norms, pi_norms=pi_norms)

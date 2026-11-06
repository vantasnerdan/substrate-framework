"""Atomic immutable coefficient items, partial manifests, and one-q lazy loader."""
import json
import os
from pathlib import Path
import numpy as np
from substrate_framework.streamed_material_state import file_binding, authenticate
from reduced_projection_algebra import q_feature_geometry
from reduced_projection_storage import ImmutableItems, MappedQBanks

SCHEMA = 'p254-physical-activeband-cache-v1'


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.partial')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
        stream.flush(); os.fsync(stream.fileno())
    os.replace(temporary, path)


class ProjectionCache:
    def __init__(self, directory, specification):
        self.directory = Path(directory).resolve()
        self.directory.mkdir(parents=True, exist_ok=True)
        self.path = self.directory/'manifest.json'
        if self.path.exists():
            self.manifest = json.loads(self.path.read_text())
            if self.manifest['schema'] != SCHEMA or self.manifest['specification'] != specification:
                raise ValueError('cache source/basis/configuration binding differs')
        else:
            self.manifest = dict(schema=SCHEMA, status='partial', specification=specification,
                                 items={}, timings={}, failures=[])
            self.publish()
        self.reader = ImmutableItems(self.manifest['items'])

    def publish(self):
        atomic_json(self.path, self.manifest)

    def has(self, key):
        record = self.manifest['items'].get(key)
        if record is None:
            return False
        self.reader.authenticate(record['archive'])
        return True

    def put(self, key, arrays, *, metadata=None, timings=None):
        if self.has(key):
            raise ValueError('completed item is immutable: '+key)
        path = self.directory/(key.replace('/', '-')+'.npz')
        temporary = path.with_suffix('.partial.npz')
        # A stale unreceipted file is not a completed item and is replaced only
        # after the new archive is durably closed. No partial archive is loaded.
        with temporary.open('wb') as stream:
            np.savez(stream, **arrays)
            stream.flush(); os.fsync(stream.fileno())
        os.replace(temporary, path)
        self.manifest['items'][key] = dict(archive=file_binding(path), metadata=metadata or {})
        if timings is not None:
            self.manifest['timings'][key] = timings
        self.publish()

    def get(self, key):
        record = self.manifest['items'][key]
        with np.load(self.reader.authenticate(record['archive']), allow_pickle=False) as data:
            return {name: data[name] for name in data.files}


def load_projected_coefficients(manifest_path, bands=8, features=7, *, allow_partial=False):
    """Load only fixed arrays now; supply physical q banks lazily without copies.

    Default requires the full all-q manifest. Explicit allow_partial accepts only
    a Gamma-complete model; attempting an unexported q raises with its exact key.
    Returned numeric object carries no private material-loader dependency.
    """
    from reduced_coefficients import ReducedCoefficients, QCoefficients
    from substrate_framework.material_clamped_array import ClampedArray
    path = Path(manifest_path).resolve()
    manifest = json.loads(path.read_text())
    spec, items = manifest['specification'], manifest['items']
    allowed = ('complete', 'gamma-complete') if allow_partial else ('complete',)
    if manifest['schema'] != SCHEMA or manifest['status'] not in allowed:
        raise ValueError('complete all128-slot/all64-q projection manifest required; '
                         'explicit allow_partial admits a Gamma-complete manifest')
    if bands not in (8, 12) or bands > spec['bands'] or features not in (7, 27) or features > spec['features']:
        raise ValueError('requested nested basis/features not present')
    for binding in spec['implementation']+spec['physical_implementation']:
        authenticate(binding)
    authenticate(spec['original_boundary'])
    if 'passive/%02d' % bands not in items:
        raise ValueError('explicit passive initial split is not complete for requested band size')
    reader = ImmutableItems(items)
    get = reader.get
    fixed = [get('slot/%03d' % s, ('h0', 'rho0', 'field', 'field_reverse')) for s in range(128)]
    h0 = np.asarray([x['h0'][:bands, :bands] for x in fixed])
    rho0 = np.asarray([x['rho0'][:bands, :bands] for x in fixed])
    field = np.asarray([x['field'][:, :bands, :bands] for x in fixed])
    field_reverse = np.asarray([x['field_reverse'][:, :bands, :bands] for x in fixed])
    q_provider = MappedQBanks(path, manifest, bands, features, QCoefficients, q_feature_geometry)
    mechanical_arrays = get('mechanical')
    mechanical = ClampedArray(**mechanical_arrays)
    indexing = spec['indexing']
    return ReducedCoefficients(mesh=tuple(spec['mesh']), k_indices=np.asarray(indexing['k_indices']),
        q_indices=np.asarray(indexing['q_indices']), slot_spin=np.asarray(indexing['slot_spin']),
        slot_k=np.asarray(indexing['slot_k']), target=np.asarray(indexing['target']),
        minus=np.asarray(indexing['minus']), h0=h0, rho0=rho0,
        field=field, field_reverse=field_reverse, mechanical=mechanical, q_provider=q_provider,
        field_polarizations=mechanical_arrays['polarizations'])

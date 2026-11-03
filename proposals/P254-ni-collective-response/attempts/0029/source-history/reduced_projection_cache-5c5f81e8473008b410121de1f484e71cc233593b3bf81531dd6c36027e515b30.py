"""Atomic immutable coefficient items, partial manifests, and one-q lazy loader."""
import json
import os
from pathlib import Path
import numpy as np
from substrate_framework.streamed_material_state import file_binding, authenticate
from reduced_projection_algebra import q_feature_geometry

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

    def publish(self):
        atomic_json(self.path, self.manifest)

    def has(self, key):
        record = self.manifest['items'].get(key)
        if record is None:
            return False
        authenticate(record['archive'])
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
        with np.load(authenticate(record['archive']), allow_pickle=False) as data:
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
    def get(key, names=None):
        if key not in items:
            raise ValueError('physical projection item is not exported: '+key)
        with np.load(authenticate(items[key]['archive']), allow_pickle=False) as archive:
            return {name: archive[name] for name in (archive.files if names is None else names)}
    fixed = [get('slot/%03d' % s, ('h0', 'rho0', 'field', 'field_reverse')) for s in range(128)]
    h0 = np.asarray([x['h0'][:bands, :bands] for x in fixed])
    rho0 = np.asarray([x['rho0'][:bands, :bands] for x in fixed])
    field = np.asarray([x['field'][:, :bands, :bands] for x in fixed])
    field_reverse = np.asarray([x['field_reverse'][:, :bands, :bands] for x in fixed])
    def q_provider(i):
        if not 0 <= i < 64:
            raise IndexError(i)
        if 'kernel/%02d' % i not in items or any('edge/%02d/%03d' % (i, s) not in items for s in range(128)):
            raise ValueError('physical q bank is incomplete: q='+str(i))
        gs, partner = q_feature_geometry(i, features)
        full_g, _ = q_feature_geometry(i, spec['features'])
        lookup = {tuple(g): j for j, g in enumerate(full_g)}
        within = np.asarray([lookup[tuple(g)] for g in gs])
        selection = np.r_[within, len(full_g)+within]
        # Allocate precisely one q bank, never an entire all-q tensor.
        d = spec['ndof']; r = 2*len(gs)
        ionic = np.empty((128, d, bands, bands), complex)
        reverse = np.empty_like(ionic)
        density = np.empty((128, r, bands, bands), complex)
        potential = np.empty((128, bands, bands, r), complex)
        for s in range(128):
            data = get('edge/%02d/%03d' % (i, s))
            ionic[s] = data['ionic'][:, :bands, :bands]
            reverse[s] = data['ionic_reverse'][:, :bands, :bands]
            density[s] = data['density'][selection, :bands, :bands]
            potential[s] = data['potential'][:bands, :bands, :][:, :, selection]
        kernel = get('kernel/%02d' % i)['kernel'][np.ix_(selection, selection)]
        return QCoefficients(ionic=ionic, ionic_reverse=reverse, density=density,
                             potential=potential, kernel=kernel, feature_partner=partner,
                             feature_g=gs, feature_spin=np.repeat([1, 2], len(gs)))
    mechanical_arrays = get('mechanical')
    mechanical = ClampedArray(**mechanical_arrays)
    indexing = spec['indexing']
    return ReducedCoefficients(mesh=tuple(spec['mesh']), k_indices=np.asarray(indexing['k_indices']),
        q_indices=np.asarray(indexing['q_indices']), slot_spin=np.asarray(indexing['slot_spin']),
        slot_k=np.asarray(indexing['slot_k']), target=np.asarray(indexing['target']),
        minus=np.asarray(indexing['minus']), h0=h0, rho0=rho0,
        field=field, field_reverse=field_reverse, mechanical=mechanical, q_provider=q_provider,
        field_polarizations=mechanical_arrays['polarizations'])

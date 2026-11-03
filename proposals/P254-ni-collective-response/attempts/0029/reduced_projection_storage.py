"""Authenticated immutable items and persistent bounded-memory q assembly."""
from collections import OrderedDict
from hashlib import sha256
import json
import os
from pathlib import Path
import numpy as np
from substrate_framework.streamed_material_state import authenticate, file_binding


class ImmutableItems:
    """Authenticate an immutable binding once during this loader lifetime."""
    def __init__(self, items):
        self.items = items
        self.checked = set()

    def authenticate(self, binding):
        identity = (binding['path'], int(binding['bytes']), binding['sha256'])
        if identity not in self.checked:
            authenticate(binding)
            self.checked.add(identity)
        return Path(binding['path'])

    def get(self, key, names=None):
        if key not in self.items:
            raise ValueError('physical projection item is not exported: '+key)
        with np.load(self.authenticate(self.items[key]['archive']), allow_pickle=False) as archive:
            return {name: archive[name] for name in (archive.files if names is None else names)}


def _publish(path, payload):
    temporary = path.with_suffix('.partial.json')
    with temporary.open('w') as stream:
        json.dump(payload, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())
    os.replace(temporary, path)


class MappedQBanks:
    """Keep Gamma and one most-recent non-Gamma record; payloads are read-only mmap.

    All-q storage is on disk, not an avoidable in-memory copy. Caller-held q
    records may outlive eviction; their mappings remain valid independently.
    """
    def __init__(self, manifest_path, manifest, bands, family, q_type, geometry):
        self.manifest = manifest; self.items = ImmutableItems(manifest['items'])
        self.bands = bands; self.family = family; self.q_type = q_type; self.geometry = geometry
        root_binding = json.dumps(manifest['specification'], sort_keys=True, separators=(',', ':')).encode()
        identity = sha256(root_binding).hexdigest()[:20]
        self.root = Path(manifest_path).resolve().parent/('assembled-'+identity)/('b%02d-f%02d' % (bands, family))
        self.root.mkdir(parents=True, exist_ok=True)
        self.gamma = None; self.recent = OrderedDict()

    def __call__(self, iq):
        if not 0 <= iq < 64:
            raise IndexError(iq)
        if iq == 0 and self.gamma is not None:
            return self.gamma
        if iq in self.recent:
            self.recent.move_to_end(iq)
            return self.recent[iq]
        keys = ['edge/%02d/%03d' % (iq, s) for s in range(128)]+['kernel/%02d' % iq]
        if any(k not in self.manifest['items'] for k in keys):
            raise ValueError('physical q bank is incomplete: q='+str(iq))
        gs, partner = self.geometry(iq, self.family)
        full_g, _ = self.geometry(iq, self.manifest['specification']['features'])
        lookup = {tuple(g): j for j, g in enumerate(full_g)}
        within = np.asarray([lookup[tuple(g)] for g in gs])
        selection = np.r_[within, len(full_g)+within]
        source = [self.manifest['items'][key]['archive'] for key in keys]
        directory = self.root/('q%02d' % iq); directory.mkdir(exist_ok=True)
        receipt_path = directory/'completed.json'
        if receipt_path.exists():
            receipt = json.loads(receipt_path.read_text())
            if receipt['source_archives'] != source:
                raise ValueError('assembled bank source binding differs')
        else:
            receipt = self._assemble(directory, iq, selection, len(gs), source)
            _publish(receipt_path, receipt)
        arrays = {name: np.load(self.items.authenticate(binding), mmap_mode='r', allow_pickle=False)
                  for name, binding in receipt['arrays'].items()}
        value = self.q_type(**arrays, feature_partner=partner, feature_g=gs,
                            feature_spin=np.repeat([1, 2], len(gs)))
        if iq == 0:
            self.gamma = value
        else:
            self.recent[iq] = value
            while len(self.recent) > 1:
                self.recent.popitem(last=False)
        return value

    def _assemble(self, directory, iq, selection, ng, source):
        m, d, r = self.bands, self.manifest['specification']['ndof'], 2*ng
        shapes = dict(ionic=(128, d, m, m), ionic_reverse=(128, d, m, m),
                      density=(128, r, m, m), potential=(128, m, m, r), kernel=(r, r))
        maps = {name: np.lib.format.open_memmap(directory/(name+'.partial.npy'), mode='w+',
                                              dtype=np.complex128, shape=shape)
                for name, shape in shapes.items()}
        for s in range(128):
            data = self.items.get('edge/%02d/%03d' % (iq, s), ('ionic', 'ionic_reverse', 'density', 'potential'))
            maps['ionic'][s] = data['ionic'][:, :m, :m]
            maps['ionic_reverse'][s] = data['ionic_reverse'][:, :m, :m]
            maps['density'][s] = data['density'][selection, :m, :m]
            maps['potential'][s] = data['potential'][:m, :m, :][:, :, selection]
        kernel = self.items.get('kernel/%02d' % iq, ('kernel',))['kernel']
        maps['kernel'][:] = kernel[np.ix_(selection, selection)]
        bindings = {}
        for name in shapes:
            maps[name].flush()
            del maps[name]
            temporary = directory/(name+'.partial.npy')
            with temporary.open('rb') as stream:
                os.fsync(stream.fileno())
            path = directory/(name+'.npy')
            os.replace(temporary, path)
            bindings[name] = file_binding(path)
        return dict(schema='p254-assembled-immutable-q-v1', source_archives=source, arrays=bindings)

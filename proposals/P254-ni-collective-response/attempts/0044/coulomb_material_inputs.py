"""Read only original density members into a bounded RAM bank; never copy banks to disk."""
from dataclasses import dataclass
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0029'), str(HERE.parent/'0041')]
from substrate_framework.streamed_material_state import authenticate, file_binding
from reduced_projection_storage import ImmutableItems
from reduced_projection_algebra import q_feature_geometry
from material_electron_ion_inputs import FrozenReference, MANIFEST_SHA, MANIFEST_BYTES, DEFAULT_MANIFEST


@dataclass
class CoulombInputs:
    h0: np.ndarray
    rho0: np.ndarray
    target: np.ndarray
    minus: np.ndarray
    slot_spin: np.ndarray
    slot_k: np.ndarray
    k_indices: np.ndarray
    q_indices: np.ndarray
    density: tuple
    feature_g: tuple
    feature_partner: tuple
    reciprocal: np.ndarray
    omega: float
    tpiba: float
    bindings: dict
    n: int = 64
    ns: int = 128
    nq: int = 64
    m: int = 12


def _density_bank(reader, spec, cap_bytes):
    geometries = [q_feature_geometry(q, 27) for q in range(64)]
    required = sum(128*len(gs)*12*12*16 for gs, _ in geometries)
    if required > cap_bytes:
        raise MemoryError(f'compact all-q density bank needs {required} bytes; cap={cap_bytes}')
    banks = []
    for q, (gs, partner) in enumerate(geometries):
        if not np.array_equal(gs, spec['feature_g'][q]) or not np.array_equal(partner, spec['feature_partner'][q]):
            raise ValueError('producer q closure differs')
        bank = np.empty((128, len(gs), 12, 12), complex)
        for source in range(128):
            key = f'edge/{q:02d}/{source:03d}'
            meta = reader.items[key]['metadata']
            raw = reader.get(key, ('density', 'feature_g', 'feature_partner'))
            if (meta['source'] != source or meta['q_index'] != q
                    or meta['target'] != spec['indexing']['target'][q][source]
                    or not np.array_equal(raw['feature_g'], gs)
                    or not np.array_equal(raw['feature_partner'], partner)):
                raise ValueError('directed density identity differs: '+key)
            d = raw['density']
            spin = source//64
            if d.shape != (2*len(gs), 12, 12) or np.any(d[(1-spin)*len(gs):(2-spin)*len(gs)] != 0):
                raise ValueError('explicit spin support differs: '+key)
            bank[source] = d[spin*len(gs):(spin+1)*len(gs)]
        bank.flags.writeable = False
        banks.append(bank)
    return tuple(banks), tuple(g[0] for g in geometries), tuple(g[1] for g in geometries)


def load_coulomb_inputs(path=DEFAULT_MANIFEST, *, density_cap_bytes=900*1024**2):
    binding = file_binding(path)
    if (binding['sha256'], binding['bytes']) != (MANIFEST_SHA, MANIFEST_BYTES):
        raise ValueError('not the original full material manifest')
    manifest = json.loads(Path(path).read_bytes())
    spec = manifest['specification']
    if (manifest['status'] != 'complete' or manifest['failures'] or len(manifest['items']) != 8548
            or (spec['case'], spec['bands'], spec['features'], spec['mesh']) != ('loaded-ni31h6',12,27,[4,4,4])
            or spec['original_source_weight_in_rho'] is not False):
        raise ValueError('original full material scope differs')
    authenticate(spec['original_boundary'])
    for source in spec['implementation']+spec['physical_implementation']:
        authenticate(source)
    reader = ImmutableItems(manifest['items'])
    fixed = [reader.get(f'slot/{s:03d}', ('h0','rho0')) for s in range(128)]
    receipt_binding = spec['context']['context']['background_receipt']
    receipt = json.loads(Path(authenticate(receipt_binding)).read_bytes())
    with np.load(authenticate(receipt['archive']), allow_pickle=False) as archive:
        reciprocal = archive['reciprocal_bohr_inv']
        cell = archive['cell_bohr']
    density, gs, partners = _density_bank(reader, spec, density_cap_bytes)
    idx = {k: np.asarray(spec['indexing'][k], dtype=int) for k in
           ('target','minus','slot_spin','slot_k','k_indices','q_indices')}
    expected = np.indices((4,4,4)).reshape(3,-1).T
    if not np.array_equal(idx['k_indices'], expected) or not np.array_equal(idx['q_indices'], expected):
        raise ValueError('original complete momentum ordering differs')
    for q, v in enumerate(expected):
        tk = np.ravel_multi_index(tuple(((expected[idx['slot_k']]+v)%4).T), (4,4,4))
        if not np.array_equal(idx['target'][q], np.arange(128)//64*64+tk):
            raise ValueError('source/target/spin map differs')
        qm = int(np.ravel_multi_index(tuple(-v%4), (4,4,4)))
        if idx['minus'][q] != qm or not np.array_equal(partners[qm][partners[q]], np.arange(len(partners[q]))):
            raise ValueError('minus/partner map differs')
    return CoulombInputs(h0=np.asarray([x['h0'] for x in fixed]),
        rho0=np.asarray([x['rho0'] for x in fixed]), **idx,
        density=density, feature_g=gs, feature_partner=partners,
        reciprocal=reciprocal, omega=abs(float(np.linalg.det(cell))),
        tpiba=2*np.pi/receipt['metadata']['raw_reference_metadata']['alat_bohr'],
        bindings=dict(manifest=binding, original_boundary=spec['original_boundary'],
            background_receipt=receipt_binding, background_archive=receipt['archive'],
            authenticated_archives=[dict(path=p, bytes=n, sha256=h) for p,n,h in sorted(reader.checked)],
            density_resident_bytes=sum(x.nbytes for x in density),
            producer_implementation=spec['implementation']+spec['physical_implementation']))

"""Actual Ni31H6 sources; no fixed-fixture inputs or physical ensemble.

Uses settled QE 7.5 saved record/XML and NC readers. Geometry, full37 atom
order, both spin-local Gamma identities and raw MV metadata are retained.
Source authority is implementation-registration -> source-access-inventory;
read_bound_bytes verifies the same immutable bytes subsequently decoded.
"""
from pathlib import Path
import json
import numpy as np
from substrate_framework.qe_saved_native import read_bound_bytes, load_wavefunctions, load_density
from substrate_framework.qe_saved_reference import load_saved_reference
from substrate_framework.qe_nc_upf import read_nc_upf
from substrate_framework.qe_nc_radial import build_radial_tables
from substrate_framework.qe_material_ionic_action import GammaIonicAction

HERE = Path(__file__).resolve().parent
DIRECTIONS = [{'atom_one_based': 32, 'symbol': 'H', 'Cartesian_axis': 'x'},
              {'atom_one_based': 1, 'symbol': 'Ni', 'Cartesian_axis': 'x'}]


def bound_json(record):
    raw, receipt = read_bound_bytes(record)
    return json.loads(raw), dict(receipt)


def load_material(inventory):
    records = {Path(x['path']).name: x for x in inventory['material_inputs']}
    chain, _ = bound_json(inventory['material_frozen_chain'])
    if chain['case'] != 'loaded-ni31h6' or chain['files']['copied_scf_xml'] != records['data-file-schema.xml']:
        raise ValueError('original frozen material source chain differs')
    if chain['files']['snapshot'] != inventory['material_historical_snapshot_verified']:
        raise ValueError('original material snapshot binding differs')
    bound_json(inventory['material_historical_snapshot_verified'])
    for symbol in ('Ni', 'H'):
        original = chain['files']['pseudo_'+symbol]
        actual = records[symbol+'_ONCV_PBE-1.2.upf']
        if (original['bytes'], original['sha256']) != (actual['bytes'], actual['sha256']):
            raise ValueError('UPF differs from original material chain')
    ref = load_saved_reference(records['data-file-schema.xml'])
    if (ref.nat != 37 or ref.atom_names != ('Ni',)*31 + ('H',)*6 or
            ref.atom_indices != tuple(range(1, 38)) or ref.nbnd_by_spin != (346, 346) or
            ref.nks != 1 or not ref.lsda or ref.noncolin or ref.spinorbit or ref.gamma_only):
        raise ValueError('not the declared full37atom collinear two-spin material')
    if ref.fft_grid != ref.fft_smooth or ref.ngm != ref.ngms:
        raise ValueError('declared NC material requires equal saved dense/smooth support')
    waves = tuple(load_wavefunctions(records[name]) for name in ('wfcup1.dat', 'wfcdw1.dat'))
    density = load_density(records['charge-density.dat'])
    for spin, w in enumerate(waves, 1):
        point, energies, occupations = ref.select(1, spin)
        if (w.ik != 1 or w.ispin != spin or w.gamma_only or w.npol != 1 or
                w.nbnd != 346 or w.igwx != 18613 or point.npw != w.igwx or
                np.any(w.k_cart_bohr_inv != 0) or np.any(point.k_cart_bohr_inv != 0)):
            raise ValueError('exact full-complex Gamma/spin/PW identity mismatch')
        if not np.allclose(w.reciprocal_bohr_inv, ref.reciprocal_bohr_inv, rtol=2e-13, atol=2e-14):
            raise ValueError('saved physical reciprocal columns disagree')
    if not np.array_equal(waves[0].miller, waves[1].miller):
        raise ValueError('spin PW rows differ; cannot reuse projector factors')
    if (density.gamma_only or density.ngm != ref.ngm or density.nspin != 2 or
            not np.allclose(density.reciprocal_bohr_inv, ref.reciprocal_bohr_inv, rtol=2e-13, atol=2e-14)):
        raise ValueError('full local Fourier support identity mismatch')
    qmax = max(np.linalg.norm((ref.reciprocal_bohr_inv @ m).T, axis=1).max()
               for m in (waves[0].miller, density.miller))
    omega = abs(float(np.linalg.det(ref.cell_bohr)))
    tables = {}
    for species in ref.species:
        record = records[species.pseudo_file]
        upf = read_nc_upf(record['path'], expected_sha256=record['sha256'], expected_bytes=record['bytes'])
        tables[species.name] = build_radial_tables(upf, omega, float(qmax))
    action = GammaIonicAction(ref.fft_grid, waves[0].miller, density.miller,
                              ref.reciprocal_bohr_inv, ref.positions_bohr, ref.atom_names, tables)
    directions = []
    for spec in DIRECTIONS:
        u = np.zeros((37, 3), dtype=np.complex128)
        u[spec['atom_one_based']-1, 0] = 1
        directions.append(u)
    return ref, waves, density, action, directions


def coordinate_arrays(ref, waves, density):
    """Original rows/atoms, physical units, actual energies and unmodified MV."""
    return {'pw_miller': waves[0].miller, 'local_miller': density.miller,
            'reciprocal_bohr_inv': ref.reciprocal_bohr_inv,
            'positions_bohr': ref.positions_bohr, 'cell_bohr': ref.cell_bohr,
            'atom_one_based': np.asarray(ref.atom_indices, dtype=np.int64),
            'atom_symbols': np.asarray(ref.atom_names),
            'fft_shape': np.asarray(ref.fft_grid, dtype=np.int64),
            'energies_ry': np.stack([ref.select(1, s)[1] for s in (1, 2)]),
            'raw_mv_occupations_not_probabilities': np.stack([ref.select(1, s)[2] for s in (1, 2)]),
            'bands_one_based': np.arange(1, 347, dtype=np.int64),
            'spins_one_based': np.asarray([1, 2], dtype=np.int64)}

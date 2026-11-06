"""Bound original0024 intake and numeric-only factor checkpoint reconstruction.

No ionic action/radial setup or material state is constructed. The saved
functional consumer requires a caller-owned receipt (normally from a bound
terminal); it verifies the immutable archive through unchanged0024 IO before
restoring FactorFunctional directly, without FFTs, contractions or eigensolves.
"""
import hashlib
from pathlib import Path
from types import MappingProxyType

import numpy as np
from material_run_contract import bound_json, canonical, recover_block
from substrate_framework.qe_saved_reference import load_saved_reference
from substrate_framework.qe_saved_native import load_wavefunctions
from substrate_framework.physical_factor_preparation import FactorFunctional
from substrate_framework.material_observable_factors import NumberProbeGrid, CompleteGridNormalization

ORIGINAL_NAME = 'terminal-cab1323dcfbe4b38a50e1eadcfa231b9.json'
ORIGINAL_BYTES = '867367'
ORIGINAL_SHA256 = '729f17467f67b330e5b68838e66886e74c6a908237cbe7a70fe02f7892de543d'
DIRECTIONS = [{'atom_one_based': 32, 'symbol': 'H', 'Cartesian_axis': 'x'},
              {'atom_one_based': 1, 'symbol': 'Ni', 'Cartesian_axis': 'x'}]
OPERATORS = ('force_H32_x', 'force_Ni1_x', 'number_cos_b1', 'number_sin_b1')
SCOPE = 'selected k=q=Gamma; both collinear spins separately; NOT complete BvK'


def original_arrays(receipt, identity):
    """Verify against the receipt embedded in the bound ORIGINAL terminal."""
    archive = Path(receipt['archive']['path'])
    result = recover_block(archive.parent, archive.stem, identity, receipt['metadata'])
    if result != receipt:
        raise ValueError('original terminal/committed receipt mismatch: '+archive.stem)
    return archive


def load_original(inputs):
    if set(inputs) != {'xml', 'wavefunctions', 'ionic_terminal'} or set(inputs['wavefunctions']) != {'1', '2'}:
        raise ValueError('inputs require xml, wavefunctions {1,2}, ionic_terminal only')
    record = inputs['ionic_terminal']
    if (Path(record['path']).name != ORIGINAL_NAME or record['bytes'] != ORIGINAL_BYTES
            or record['sha256'] != ORIGINAL_SHA256):
        raise ValueError('only immutable original0024 terminal is admitted')
    terminal = bound_json(record)
    recipe, identity = terminal['recipe'], terminal['recipe_sha256']
    if (terminal['status'] != 'complete' or terminal['exit_code'] != 0
            or hashlib.sha256(canonical(recipe)).hexdigest() != identity
            or recipe['directions'] != DIRECTIONS or recipe['bands'] != 346
            or recipe['npw'] != 18613 or recipe['nat'] != 37 or recipe['spins'] != [1, 2]):
        raise ValueError('original ionic recipe/terminal identity mismatch')
    for item in (inputs['xml'], *inputs['wavefunctions'].values()):
        if item not in recipe['inputs']:
            raise ValueError('material input is not an original0024 input')
    ref = load_saved_reference(inputs['xml'])
    if (ref.nat != 37 or ref.atom_names != ('Ni',)*31+('H',)*6
            or ref.atom_indices != tuple(range(1, 38)) or ref.nbnd_by_spin != (346, 346)
            or not ref.lsda or ref.noncolin or ref.spinorbit or ref.gamma_only
            or ref.fft_grid != ref.fft_smooth or ref.ngm != ref.ngms):
        raise ValueError('not the original full37atom collinear saved material')
    archive = original_arrays(terminal['coordinates'], identity)
    with np.load(archive, allow_pickle=False) as saved:
        coordinates = {name: saved[name] for name in saved.files}
    directions = np.zeros((2, 37, 3), dtype=np.complex128)
    directions[0, 31, 0] = directions[1, 0, 0] = 1
    expected = {'positions_bohr': ref.positions_bohr, 'cell_bohr': ref.cell_bohr,
                'reciprocal_bohr_inv': ref.reciprocal_bohr_inv,
                'atom_one_based': ref.atom_indices, 'atom_symbols': ref.atom_names,
                'fft_shape': ref.fft_grid, 'bands_one_based': np.arange(1, 347),
                'spins_one_based': [1, 2], 'directions_dimensionless': directions,
                'energies_ry': np.stack([ref.select(1, s)[1] for s in (1, 2)]),
                'raw_mv_occupations_not_probabilities': np.stack([ref.select(1, s)[2] for s in (1, 2)])}
    for name, value in expected.items():
        if not np.array_equal(coordinates[name], value):
            raise ValueError('original coordinate identity differs: '+name)
    if terminal['coordinates']['metadata']['direction_order'] != DIRECTIONS:
        raise ValueError('original coordinate direction order differs')
    return ref, coordinates, terminal


def spin_receipts(terminal, spin):
    """Require every original block exactly once, in its original band recipe."""
    blocks = terminal['blocks']
    if len(blocks) != 4*((346+terminal['recipe']['block_bands']-1)//terminal['recipe']['block_bands']):
        raise ValueError('original block count differs')
    selected = []
    for direction in DIRECTIONS:
        for start in range(1, 347, terminal['recipe']['block_bands']):
            bands = list(range(start, min(start+terminal['recipe']['block_bands'], 347)))
            matches = [r for r in blocks if r['metadata']['spin_one_based'] == spin
                       and r['metadata']['direction'] == direction
                       and r['metadata']['source_bands_one_based'] == bands]
            if len(matches) != 1:
                raise ValueError('missing/duplicate original full-PW band block')
            receipt = matches[0]
            metadata = receipt['metadata']
            if (metadata['k_spin_local_one_based'] != 1 or metadata['units'] != 'Ry/bohr'
                    or metadata['target_saved_bands_one_based'] != list(range(1, 347))
                    or metadata['full_pw_rows'] != 'coordinates.npz:pw_miller original column order'
                    or metadata['atom_order'] != 'coordinates.npz:atom_one_based/atom_symbols'):
                raise ValueError('original block coordinate semantics differ')
            selected.append(receipt)
    return selected


def load_spin(inputs, ref, coordinates, terminal, spin, receipts, budget, *, assemble):
    wave = load_wavefunctions(inputs['wavefunctions'][str(spin)])
    point, _, _ = ref.select(1, spin)
    if (wave.ik != 1 or wave.ispin != spin or wave.gamma_only or wave.npol != 1
            or wave.nbnd != 346 or wave.igwx != 18613 or point.npw != wave.igwx
            or np.any(wave.k_cart_bohr_inv != 0) or np.any(point.k_cart_bohr_inv != 0)
            or wave.coefficients.shape != (18613, 346)
            or not np.array_equal(wave.miller, coordinates['pw_miller'])
            or not np.allclose(wave.reciprocal_bohr_inv, ref.reciprocal_bohr_inv, rtol=2e-13, atol=2e-14)):
        raise ValueError('actual Gamma band/PW/spin identity mismatch')
    actions = {atom: np.empty_like(wave.coefficients) for atom in (32, 1)} if assemble else None
    for receipt in receipts:
        budget()
        archive = original_arrays(receipt, terminal['recipe_sha256'])
        metadata = receipt['metadata']
        bands = metadata['source_bands_one_based']
        with np.load(archive, allow_pickle=False) as saved:
            raw = saved['gpsi_ry_per_bohr']
            if (not np.array_equal(saved['source_bands_one_based'], bands)
                    or raw.shape != (18613, len(bands)) or raw.dtype != np.dtype('complex128')
                    or not np.isfinite(raw).all()
                    or saved['saved_band_overlap_ry_per_bohr'].shape != (346, len(bands))):
                raise ValueError('original gPsi array identities differ')
            if assemble:
                actions[metadata['direction']['atom_one_based']][:, bands[0]-1:bands[-1]] = raw
    budget()
    return wave, actions


def factor_provenance(inputs, settings, coordinates, terminal, spin, receipts):
    if (set(settings) != {'complete_grid_shape', 'normalization_authority'}
            or settings['complete_grid_shape'] != [4, 4, 4]):
        raise ValueError('settings require complete_grid_shape=[4,4,4] and normalization_authority')
    normalization = CompleteGridNormalization(tuple(settings['complete_grid_shape']), settings['normalization_authority'])
    basis = hashlib.sha256(canonical({name: terminal['coordinates']['arrays'][name]
        for name in ('pw_miller', 'cell_bohr', 'reciprocal_bohr_inv')})).hexdigest()
    grid = NumberProbeGrid(tuple(int(n) for n in coordinates['fft_shape']), coordinates['pw_miller'],
                           coordinates['cell_bohr'], coordinates['reciprocal_bohr_inv'], basis)
    provenance = {'inputs': inputs, 'spin_one_based': spin, 'k_spin_local_one_based': 1,
                  'source_bands_one_based': list(range(1, 347)), 'full_pw_dimension': 18613,
                  'original_recipe_sha256': terminal['recipe_sha256'],
                  'original_coordinates': terminal['coordinates'], 'original_blocks': receipts,
                  'number_grid': grid.metadata(), 'normalization': normalization.metadata(),
                  'scope': SCOPE, 'state': None, 'D': None, 'density_coefficients': None,
                  'temperature': None, 'ionic': None, 'cross': None, 'total': None}
    return provenance, grid


def load_saved_functional(target, key, identity, expected_receipt):
    """Root readback API: trusted terminal receipt in; restored functional out.

    No preparation, population mapping, rank recheck or action recomputation.
    Exact array integrity is established by recover_block; structural checks
    reject incompatible serialized layouts before the immutable object exists.
    """
    target = Path(target)
    receipt = recover_block(target, key, identity, expected_receipt['metadata'])
    if receipt is None or receipt != expected_receipt:
        raise ValueError('saved functional does not match caller-owned receipt')
    metadata = receipt['metadata']['functional']
    layout = metadata['array_layout']
    names = tuple(layout['overlaps'])
    if (set(names) != set(OPERATORS) or metadata['factor_count'] != 346
            or metadata['full_pw_dimension'] != 18613
            or metadata['state'] is not None or metadata['density_coefficients'] is not None):
        raise ValueError('not the full state-independent material functional')
    products = {(item['left'], item['right']): item['array'] for item in layout['ordered_products']}
    if len(layout['ordered_products']) != 16 or set(products) != {(a, b) for a in names for b in names}:
        raise ValueError('directed T_AB coverage differs')
    keys = [layout['gram'], *layout['overlaps'].values(), *products.values()]
    if len(set(keys)) != 21 or set(keys) != set(receipt['arrays']):
        raise ValueError('Gram/M/T array layout differs')
    with np.load(target/(key+'.npz'), allow_pickle=False) as saved:
        arrays = {name: saved[name] for name in keys}
    for value in arrays.values():
        if value.shape != (346, 346) or value.dtype != np.dtype('complex128') or not np.isfinite(value).all():
            raise ValueError('invalid numeric functional array')
        value.flags.writeable = False
    functional = FactorFunctional(arrays[layout['gram']],
        MappingProxyType({name: arrays[layout['overlaps'][name]] for name in names}),
        MappingProxyType({pair: arrays[key] for pair, key in products.items()}),
        MappingProxyType({name: MappingProxyType(metadata['operators'][name]) for name in names}),
        metadata['factor_identity'], metadata['pw_basis_identity'], metadata['full_pw_dimension'])
    if functional.archive_payload()[1] != metadata:
        raise ValueError('restored functional metadata differs from serialized definition')
    return functional

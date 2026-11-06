"""Authenticated complete neutral ledger and separately selected Gamma raw D.

Readback never refills populations, regenerates factors/H0, or loads a WFC.
It authenticates completed Bloch coverage, original factor provenance, current
source bindings and numeric archives; manager post-exit success is separate.
"""
from dataclasses import dataclass
import hashlib
from pathlib import Path
from types import MappingProxyType
import numpy as np

from material_run_contract import binding, bound_json, canonical, recover_block, verify_consumers
from material_reference_output import load_arrays
from .material_bloch_output import load_saved_bloch_coverage
from .material_preparation_inputs import (KEYS, validate_settings, preparation_policy,
                                         load_factor_terminal)
from .material_preparation_gamma import gamma_metadata

RECEIPT_FIELDS = ('recipe_sha256', 'archive', 'arrays', 'metadata', 'readback')


def population_arrays(specification):
    sectors, levels = specification.sectors, specification.filling.levels
    return dict(populations=np.array(specification.filling.populations).reshape(128, 346),
        energies_ry=np.array([l.energy for l in levels]).reshape(128, 346),
        measures=np.array([levels[start].measure for start, _ in specification.band_slices]),
        spin=np.array([s.spin for s in sectors], dtype=np.int64),
        mesh_index=np.array([s.mesh_index for s in sectors], dtype=np.int64),
        source_ik=np.array([s.source_ik for s in sectors], dtype=np.int64),
        source_sha256=np.array([s.source_binding.sha256 for s in sectors]),
        frame_identity=np.array([s.frame_identity for s in sectors]),
        basis_identity=np.array([s.pw_basis_identity for s in sectors]))


def validate_population_arrays(arrays, metadata, coverage):
    shapes = {'populations': (128, 346), 'energies_ry': (128, 346), 'measures': (128,),
              'spin': (128,), 'mesh_index': (128, 3), 'source_ik': (128,),
              'source_sha256': (128,), 'frame_identity': (128,), 'basis_identity': (128,)}
    _shapes(arrays, shapes)
    for i, record in enumerate(coverage.slot_records):
        expected = {'spin': record.ispin, 'mesh_index': record.mesh_index, 'source_ik': record.source_ik,
                    'source_sha256': record.source_binding['sha256'], 'frame_identity': record.gauge_identity,
                    'basis_identity': record.basis_identity,
                    'energies_ry': coverage.grid.reference.select(record.source_ik, record.ispin)[1]}
        if any(not np.array_equal(arrays[key][i], value) for key, value in expected.items()):
            raise ValueError('saved population source/slot/frame/energy identity differs')
    p = arrays['populations']
    if (np.any(p < 0) or np.any(p > 1) or not np.array_equal(arrays['measures'], np.full(128, 1/64))
            or abs(np.sum(p * arrays['measures'][:, None]) - metadata['neutral_target']) > 1e-10
            or metadata['neutral_target'] != coverage.grid.reference.nelec):
        raise ValueError('saved full-grid mean-number/population measure differs')
    order = np.argsort(arrays['energies_ry'].ravel(), kind='stable')
    energy, population = arrays['energies_ry'].ravel()[order], p.ravel()[order]
    if (np.any(np.diff(population) > 0)
            or np.any((np.diff(energy) == 0) & (np.diff(population) != 0))):
        raise ValueError('saved populations are not scalar exact-energy zero-T blocks')


def _shapes(arrays, shapes):
    if set(arrays) != set(shapes):
        raise ValueError('preparation archive array keys differ')
    for key, shape in shapes.items():
        value = arrays[key]
        if value.shape != shape or (value.dtype.kind in 'fc' and not np.isfinite(value).all()):
            raise ValueError('preparation archive shape/finite-value mismatch: ' + key)


def validate_selected_arrays(arrays, population_arrays, coverage):
    shapes = {'density_coefficients': (346, 346), 'gram': (346, 346), 'populations': (346,),
              'energies_ry': (346,), 'density_eigenvalues': (346,), 'selected_number_pair': (2,),
              'means': (2,), 'ordered_connected': (2,), 'pauli_roundoff_allowance': (1,),
              'gram_min_eigenvalue': (1,)}
    _shapes(arrays, shapes)
    index = next(i for i, r in enumerate(coverage.slot_records) if r.ispin == 1 and r.mesh_index == (0, 0, 0))
    for key in ('populations', 'energies_ry'):
        if not np.array_equal(arrays[key], population_arrays[key][index]):
            raise ValueError('selected sector populations differ from complete ledger')
    allowance = arrays['pauli_roundoff_allowance'][0]
    spectrum = arrays['density_eigenvalues']
    if (not 0 < allowance < 1e-8 or arrays['gram_min_eigenvalue'][0] <= 0
            or spectrum.min() < -allowance or spectrum.max() > 1 + allowance
            or abs(arrays['selected_number_pair'][0] - arrays['selected_number_pair'][1]) > 1e-8
            or abs(arrays['selected_number_pair'][1] - arrays['populations'].sum()) > 1e-10):
        raise ValueError('saved selected number/Pauli diagnostics differ')
    if any(arrays[key].dtype != np.dtype('complex128') for key in
           ('density_coefficients', 'gram', 'selected_number_pair', 'means', 'ordered_connected')):
        raise ValueError('raw complex128 coefficient/ordered data required')


def load_saved_preparation_block(target, key, identity, expected_receipt):
    target = Path(target).resolve()
    expected = {field: expected_receipt[field] for field in RECEIPT_FIELDS}
    if key not in KEYS:
        raise ValueError('unknown completed preparation boundary')
    actual = recover_block(target, key, identity, expected['metadata'])
    if actual != expected or actual['archive']['path'] != str(target / (key + '.npz')):
        raise ValueError('preparation block differs from authenticated terminal receipt')
    arrays = load_arrays(actual)
    for value in arrays.values():
        value.flags.writeable = False
    return MappingProxyType(arrays)


@dataclass(frozen=True)
class SavedMaterialPreparation:
    population_ledger: object
    selected_gamma: object
    receipts: object
    coverage_identity: str
    recipe: object


def load_saved_material_preparation(terminal_binding):
    """Restore authenticated numeric results, not a refilled or evolved state."""
    terminal = bound_json(terminal_binding)
    recipe, identity = terminal['recipe'], terminal['recipe_sha256']
    if (terminal['status'] != 'complete' or terminal['exit_code'] != 0
            or recipe['kind'] != 'material-neutral-preparation'
            or hashlib.sha256(canonical(recipe)).hexdigest() != identity):
        raise ValueError('complete authenticated material preparation terminal required')
    validate_settings(recipe['inputs'], recipe['settings'])
    final = terminal['final_implementation_sources']
    verify_consumers(recipe, final)
    sources = {r['path']: r for r in recipe['implementation_sources']}
    sources.update({r['path']: r for r in final})
    for source in sources.values():
        if binding(source['path']) != source:
            raise ValueError('preparation consumer source changed: ' + source['path'])
    if sources.get(str(Path(__file__).resolve())) != binding(__file__):
        raise ValueError('current preparation reader is not source-authenticated')
    outputs = terminal['outputs']
    if len(outputs) != 2 or {entry['key'] for entry in outputs} != set(KEYS):
        raise ValueError('both completed neutral ledger and actual selected Gamma result required')
    entries = {entry['key']: entry for entry in outputs}
    coverage = load_saved_bloch_coverage(recipe['inputs']['bloch_terminal'])
    factor_terminal = load_factor_terminal(recipe['inputs'], coverage)
    _, _, expected_gamma = gamma_metadata(recipe['inputs'], recipe['settings'], coverage, factor_terminal)
    if entries[KEYS[1]]['metadata'] != expected_gamma:
        raise ValueError('saved Gamma original-frame/source policy differs')
    metadata = entries[KEYS[0]]['metadata']
    from dataclasses import asdict
    if (metadata['schema'] != 'material-neutral-populations-v1'
            or metadata['inputs'] != recipe['inputs'] or metadata['settings'] != recipe['settings']
            or metadata['coverage_identity'] != coverage.coverage_identity
            or metadata['policy'] != asdict(preparation_policy(recipe['settings'], coverage.coverage_identity))):
        raise ValueError('saved neutral policy/source coverage differs')
    restored = {key: load_saved_preparation_block(Path(entry['archive']['path']).parent, key, identity, entry)
                for key, entry in entries.items()}
    validate_population_arrays(restored[KEYS[0]], metadata, coverage)
    validate_selected_arrays(restored[KEYS[1]], restored[KEYS[0]], coverage)
    return SavedMaterialPreparation(restored[KEYS[0]], restored[KEYS[1]], MappingProxyType(entries),
                                    coverage.coverage_identity, MappingProxyType(recipe))

"""Join one ORIGINAL Gamma raw frame to its admitted A25 functional honestly.

Bloch and A25 basis labels are retained separately. Their equivalence is tied
to original XML/WFC bytes, identity transport, complete ordered row coordinates
and all original columns, not inferred from equal coefficient arrays alone.
"""
from dataclasses import replace
from pathlib import Path
import numpy as np

from material_factor_output import load_saved_functional, original_arrays
from material_reference_output import load_arrays
from material_run_contract import bound_json, canonical
from material_pilot_io import array_record
from .material_preparation import SourceBinding, instantiate_material_sector
from .qe_saved_native import load_wavefunctions

PROBES = ('force_H32_x', 'number_cos_b1')


def gamma_metadata(inputs, settings, coverage, terminal):
    matches = [r for r in coverage.slot_records if r.ispin == 1 and r.mesh_index == (0, 0, 0)]
    if len(matches) != 1:
        raise ValueError('one original spin1 Gamma slot required')
    record = matches[0]
    entry = next(o for o in terminal['outputs'] if o['key'] == 'spin1-functional')
    receipt = entry['receipt']
    if bound_json(entry['receipt_binding']) != receipt or entry['metadata'] != receipt['metadata']:
        raise ValueError('accepted factor terminal/receipt mismatch')
    provenance, metadata = receipt['metadata']['provenance'], receipt['metadata']['functional']
    original = bound_json(terminal['recipe']['inputs']['ionic_terminal'])
    if (provenance['inputs'] != terminal['recipe']['inputs']
            or provenance['original_coordinates'] != original['coordinates']
            or provenance['original_recipe_sha256'] != original['recipe_sha256']
            or provenance['spin_one_based'] != 1 or provenance['k_spin_local_one_based'] != 1
            or provenance['source_bands_one_based'] != list(range(1, 347))
            or provenance['full_pw_dimension'] != 18613):
        raise ValueError('accepted full-row/full-column source provenance differs')
    if (not record.original_representative or record.source_ik != 1
            or record.raw_band_order != 'original-all-saved-columns'
            or record.row_order != 'inherited-source-order'
            or tuple(record.reciprocal_shift) != (0, 0, 0)
            or tuple(record.target_k_bohr_inv) != (0., 0., 0.)
            or tuple(record.source_k_bohr_inv) != (0., 0., 0.)
            or record.source_binding != terminal['recipe']['inputs']['wavefunctions']['1']
            or record.source_frame_identity != metadata['factor_identity']
            or record.source_frame_identity != record.source_binding['sha256']
            or record.gamma_basis_join_identity != metadata['pw_basis_identity']
            or (record.n_pw, record.n_bands) != (18613, 346)):
        raise ValueError('no authenticated original Gamma/A25 raw-frame basis join')
    coordinates = provenance['original_coordinates']
    original_arrays(coordinates, provenance['original_recipe_sha256'])
    coordinate_arrays = load_arrays(coordinates)
    source_arrays = coverage.source_receipts[1, 1]['arrays']
    reference_reciprocal = source_arrays['reference_reciprocal_bohr_inv']
    if array_record(coordinate_arrays['reciprocal_bohr_inv']) != reference_reciprocal:
        raise ValueError('historical coordinate reciprocal record differs from admitted XML reference')
    # Native int32 rows and canonical int64 transport have distinct byte records.
    # Identity transport preserves their exact integer values and ordering.
    target_miller = coordinate_arrays['pw_miller'].astype(np.int64, copy=False)
    if (array_record(coordinate_arrays['pw_miller']) != record.source_miller_record
            or array_record(target_miller) != record.target_miller_record
            or not np.array_equal(coordinate_arrays['cell_bohr'], coverage.grid.reference.cell_bohr)
            or not np.array_equal(coordinate_arrays['reciprocal_bohr_inv'], coverage.grid.reference.reciprocal_bohr_inv)):
        raise ValueError('original coordinate order/cell/reciprocal identity does not join')
    join = dict(bloch_basis_identity=record.basis_identity, bloch_gauge_identity=record.gauge_identity,
        accepted_basis_identity=metadata['pw_basis_identity'], original_frame_identity=metadata['factor_identity'],
        original_coordinates=coordinates['archive'], original_source=record.source_binding,
        source_coefficients_record=record.source_coefficients_record,
        source_miller_record=record.source_miller_record, original_representative=True,
        target_miller_record=record.target_miller_record,
        miller_scope='distinct Native int32 and canonical int64 records; exact ordered integer identity',
        binary_reciprocal_record=source_arrays['reciprocal_bohr_inv'],
        reference_reciprocal_record=reference_reciprocal,
        reciprocal_scope='historical basis uses XML reference; binary B separately admitted at atol=2e-13, not asserted equal',
        justification='same bound original XML/WFC; identity operation/no fold; identical complete ordered row coordinates and original all-column frame')
    result = dict(schema='material-selected-gamma-preparation-v1', inputs=inputs, settings=settings,
        coverage_identity=coverage.coverage_identity, spin_one_based=1, mesh_index=[0, 0, 0],
        source_ik=1, source_bands_one_based=list(range(1, 347)), full_pw_dimension=18613,
        factor_receipt_binding=entry['receipt_binding'], frame_join=join, probe_names=list(PROBES),
        ordered_pairs=[[PROBES[0], PROBES[1]], [PROBES[1], PROBES[0]]],
        probe_units=[metadata['operators'][name]['units'] for name in PROBES],
        scope='ONE selected spin1 Gamma sector: raw unweighted means/connected products, NOT neutral total',
        number_scope='Tr(S D) for this sector versus its assigned full-grid populations; no extra spin factor',
        complement='exactly empty initially, full-PW ordered products retain its contributions',
        stationarity='not imposed or established; no H0 propagation or zero-drift inference',
        ionic=None, cross=None, held_state=None, total=None, physical_temperature=None)
    return record, receipt, result


def instantiate_gamma(specification, record, receipt, metadata):
    archive = Path(receipt['archive']['path'])
    functional = load_saved_functional(archive.parent, archive.stem, receipt['recipe_sha256'], receipt)
    wave = load_wavefunctions(record.source_binding)
    if (wave.ispin != 1 or wave.ik != 1 or wave.coefficients.shape != (18613, 346)
            or array_record(wave.coefficients) != record.source_coefficients_record
            or array_record(wave.miller) != record.source_miller_record
            or not np.array_equal(wave.k_cart_bohr_inv, [0., 0., 0.])):
        raise ValueError('actual all-row/all-column original WFC differs from admitted coverage')
    # The explicit equivalence above authorizes use of the ORIGINAL labels.
    # It does not replace the distinct Bloch basis identity in the saved ledger.
    sectors = tuple(replace(s, frame_identity=functional.factor_identity,
                           pw_basis_identity=functional.pw_basis_identity)
                    if s.key == (1, (0, 0, 0)) else s for s in specification.sectors)
    joined = replace(specification, sectors=sectors)
    preparation = instantiate_material_sector(joined, (1, (0, 0, 0)), functional,
        wave.coefficients, source_binding=SourceBinding.from_mapping(record.source_binding),
        reference_binding=specification.reference_binding,
        coverage_authority=specification.policy.coverage_authority,
        frame_authority=canonical(metadata['frame_join']).decode())
    index = next(i for i, s in enumerate(specification.sectors) if s.key == (1, (0, 0, 0)))
    start, stop = specification.band_slices[index]
    populations = np.array(specification.filling.populations[start:stop])
    number = np.trace(functional.gram @ preparation.density_coefficients)
    if abs(number - populations.sum()) > 1e-8:
        raise ValueError('actual Gram-metric selected number differs from assigned populations')
    return dict(density_coefficients=preparation.density_coefficients, gram=functional.gram,
        populations=populations, energies_ry=np.array([l.energy for l in specification.filling.levels[start:stop]]),
        density_eigenvalues=np.array(preparation.admissibility['density_eigenvalues']),
        selected_number_pair=np.array([number, populations.sum()], dtype=complex),
        means=np.array([preparation.mean(name) for name in PROBES]),
        ordered_connected=np.array([preparation.ordered(*PROBES), preparation.ordered(*PROBES[::-1])])), preparation.admissibility

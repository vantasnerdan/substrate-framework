"""Authenticated original-source checkpoints for the full conditional reference."""
from pathlib import Path

import numpy as np

from material_run_contract import bound_json, recover_block
from material_reference_output import load_arrays
from substrate_framework.material_reference_state import restore_source, source_population_indices
from substrate_framework.material_reference_transport import SourceObservables


SOURCE_FIELDS = ('source', 'preparation_terminal', 'coverage_identity', 'spin', 'source_ik',
                 'state_authority', 'original_frame_identity', 'original_basis_identity', 'all_raw_columns')


def checkpoint(binding):
    expected = bound_json(binding)
    archive = Path(expected['archive']['path'])
    actual = recover_block(archive.parent, archive.stem, expected['recipe_sha256'], expected['metadata'])
    if actual != expected:
        raise ValueError('completed source checkpoint differs from its bound receipt')
    return expected, load_arrays(actual)


def source_observables(manifest, saved, grid, diagnostics):
    """Stream each original source once; no WFC, Gram or density recomputation."""
    seen = set()
    for record in manifest['sources']:
        identity = record['spin'], record['source_ik']
        if identity in seen:
            raise ValueError('duplicate original source')
        seen.add(identity)
        populations, pairs = source_population_indices(grid, saved.population_ledger, *identity)
        originals = [(i, slot) for i, slot in pairs if slot.operation.index == 1
                     and not np.any(slot.reciprocal_shift)]
        if len(originals) != 1:
            raise ValueError('unique original source representative required')
        index, representative = originals[0]
        density_receipt, moments = checkpoint(record['density_gradient'])
        metadata = density_receipt['metadata']
        if (metadata['spin'], metadata['source_ik']) != identity or not metadata['all_raw_columns']:
            raise ValueError('source moments have different spin/source/column coverage')
        if (metadata['preparation_terminal'] != manifest['preparation_terminal']
                or metadata['coverage_identity'] != saved.coverage_identity
                or metadata['source'] != dict(representative.source_binding)
                or metadata['mesh_index'] != list(representative.mesh_index)
                or metadata['population_index'] != index or metadata['measure'] != representative.weight):
            raise ValueError('moments differ from the accepted source/slot/preparation')
        if record['reuse_selected_gamma']:
            if record['state'] is not None or identity != (1, 1) or representative.mesh_index != (0, 0, 0):
                raise ValueError('A27 cache belongs only to original spin1 Gamma')
            arrays = saved.selected_gamma
        else:
            state_receipt, arrays = checkpoint(record['state'])
            if any(state_receipt['metadata'].get(key) != metadata[key] for key in SOURCE_FIELDS):
                raise ValueError('state and moments refer to different raw preparations')
        if not np.array_equal(arrays['populations'], populations):
            raise ValueError('source populations differ from the complete ledger')
        state = restore_source(arrays, source_identity=representative.source_binding['sha256'])
        curvature_receipt, curvature = checkpoint(record['curvature'])
        if record['curvature_kind'] == 'expansion':
            cm = curvature_receipt['metadata']['source_metadata']
            if (not np.array_equal(curvature['source_gram'], state.gram)
                    or not np.array_equal(curvature['source_density_coefficients'], state.density_coefficients)):
                raise ValueError('selected expansion contains a different raw state')
        elif record['curvature_kind'] == 'source':
            cm = curvature_receipt['metadata']
        else:
            raise ValueError('unknown original curvature checkpoint kind')
        if any(cm.get(key) != metadata[key] for key in SOURCE_FIELDS):
            raise ValueError('curvature and density have different source/preparation authority')
        if not np.array_equal(curvature['conditional_density_bohr3'], moments['density_bohr3']):
            raise ValueError('curvature did not consume the exact same conditional density')
        h_receipt, h = checkpoint(record['hamiltonian'])
        if any(h_receipt['metadata'].get(key) != metadata[key] for key in SOURCE_FIELDS):
            raise ValueError('H0 moments belong to a different source preparation')
        diagnostics.append(dict(spin=identity[0], source_ik=identity[1],
            source_identity=state.source_identity, slots=len(pairs),
            total_slot_measure=sum(slot.weight for _, slot in pairs),
            population_number=float(np.sum(populations)),
            energy=h['one_body_energy_ry'].item(),
            h0_commutator_norm_squared=h['density_commutator_norm_squared_ry2'].item(),
            saved_label_residual_min=float(np.min(h['saved_label_residual_norm_ry'])),
            saved_label_residual_max=float(np.max(h['saved_label_residual_norm_ry']))))
        yield SourceObservables(representative, state, moments['density_bohr3'],
            moments['positive_gradient_ry_bohr'], curvature['conditional_nonlocal_ry_bohr2'],
            curvature['original_mv_nonlocal_ry_bohr2'], curvature['mv_density_bohr3'])
        del arrays, moments, curvature, h, state
    expected = {(slot.ispin, slot.source_ik) for slot in grid.slots}
    if seen != expected:
        raise ValueError('original-source coverage is incomplete')

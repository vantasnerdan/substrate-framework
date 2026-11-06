"""Actual full-column source preparation and source-slot reference moments."""
import time
from pathlib import Path
import numpy as np

from material_run_contract import bound_json, implementation, verify_consumers, recover_block
from material_reference_output import load_arrays
from substrate_framework.material_preparation_output import load_saved_material_preparation
from substrate_framework.material_bloch_admission import restore_accepted, wave_bindings
from substrate_framework.qe_saved_reference import load_saved_reference
from substrate_framework.qe_saved_native import load_wavefunctions
from substrate_framework.qe_bloch_grid import build_bloch_grid
from substrate_framework.qe_bloch_frame import transport_bloch_frame
from substrate_framework.qe_bloch_dyads import BlochDyadCoupling
from substrate_framework.qe_material_hamiltonian import MaterialHamiltonian
from substrate_framework.material_reference_state import (
    prepare_source, restore_source, source_population_indices,
    sector_density_moments, sector_hamiltonian_moments)
from substrate_framework.material_frozen_state_curvature import conditional_minus_mv


def produce(run, recipe, paths):
    settings = recipe['settings']
    bound_json(settings['construction_registration'])
    for record in settings['oracles']:
        if bound_json(record)['status'] != 'passed':
            raise ValueError('new state/dyad equation block has not passed')
    if settings['slot_mode'] not in ('original-representative', 'all-source-slots'):
        raise ValueError('explicit complete source-slot selection required')
    verify_consumers(run.grant, implementation(paths))
    saved = load_saved_material_preparation(recipe['inputs']['preparation_terminal'])
    bloch = bound_json(saved.recipe['inputs']['bloch_terminal'])['recipe']
    original = bloch['inputs']
    reference = load_saved_reference(original['xml'])
    grid = build_bloch_grid(reference, wave_bindings(original), symmetry_policy=bloch['settings']['symmetry_policy'])
    populations, pairs = source_population_indices(grid, saved.population_ledger,
                                                   settings['spin'], settings['source_ik'])
    representative = next((index, slot) for index, slot in pairs
                          if slot.operation.index == 1 and not np.any(slot.reciprocal_shift))
    selected = (representative,) if settings['slot_mode'] == 'original-representative' else pairs
    source = load_wavefunctions(representative[1].source_binding)
    frame = transport_bloch_frame(reference, representative[1], source)
    authority = 'A27 neutral source recipe '+saved.recipe['settings']['preparation_recipe']+'; '+saved.coverage_identity
    metadata = dict(schema='material-common-reference-source-v1', source=dict(source.source),
        preparation_terminal=recipe['inputs']['preparation_terminal'],
        coverage_identity=saved.coverage_identity, spin=settings['spin'], source_ik=settings['source_ik'],
        state_authority=authority, original_frame_identity=frame.gauge_identity,
        original_basis_identity=frame.basis.identity, all_raw_columns=True,
        scope='conditional source state; not H0 equilibrium or a full weighted material total')
    previous = recipe['inputs'].get('prepared_source_receipt')
    if previous is None:
        state = prepare_source(frame.coefficients, populations, source_identity=source.source['sha256'],
                               basis_identity=frame.basis.identity, authority=authority)
    else:
        expected = bound_json(previous)
        archive = Path(expected['archive']['path'])
        actual = recover_block(archive.parent, archive.stem, expected['recipe_sha256'], metadata)
        if actual != expected:
            raise ValueError('prior completed source state differs from its bound receipt')
        restored = load_arrays(actual)
        if not np.array_equal(restored['populations'], populations):
            raise ValueError('completed source state belongs to different populations')
        state = restore_source(restored, source_identity=source.source['sha256'])
        run.summary['source_state_reused_from'] = previous
    if previous is None:
        source_arrays = dict(gram=state.gram, density_coefficients=state.density_coefficients,
            populations=state.populations,
            density_eigenvalues=np.asarray(state.admissibility['density_eigenvalues']))
        run.commit('source-state', source_arrays, metadata)
    run.event('background_restore_start')
    background, background_receipt = restore_accepted(original)
    run.summary['background_receipt'] = background_receipt
    run.event('background_restore_complete')
    for index, slot in selected:
        began = time.monotonic()
        frame = transport_bloch_frame(reference, slot, source)
        run.event('slot_coupling_start', mesh_index=list(slot.mesh_index))
        coupling = BlochDyadCoupling(background.inputs, frame.basis, block_size=settings['block_size'])
        hamiltonian = MaterialHamiltonian(background, coupling.factors, slot.ispin,
                                          block_size=settings['block_size'])
        slot_metadata = dict(metadata, boundary='full-column-source-slot-moments',
            population_index=index, mesh_index=list(slot.mesh_index), basis_identity=frame.basis.identity,
            gauge_identity=frame.gauge_identity, measure=slot.weight,
            units=dict(density='bohr^-3', positive_gradient='Ry/bohr', H0='Ry'),
            held_density=None, held_gradient=None, total_force=None)
        key = 'slot-'+str(slot.mesh_slot)
        previous_moments = recipe['inputs'].get('density_gradient_receipts', {}).get(key)
        if previous_moments is not None:
            expected = bound_json(previous_moments)
            if any(expected['metadata'].get(name) != value for name, value in slot_metadata.items()):
                raise ValueError('completed density/gradient belongs to a different source slot')
            archive = Path(expected['archive']['path'])
            actual = recover_block(archive.parent, archive.stem,
                                   expected['recipe_sha256'], expected['metadata'])
            if actual != expected:
                raise ValueError('completed density/gradient differs from its bound receipt')
            run.summary.setdefault('density_gradient_reused_from', {})[key] = previous_moments
            source_density = load_arrays(actual)['density_bohr3']
            run.event('slot_density_gradient_reused', mesh_index=list(slot.mesh_index))
        else:
            run.event('slot_density_gradient_start', mesh_index=list(slot.mesh_index))
            direction = np.zeros(3*reference.nat, complex)
            direction[3*31], direction[1] = 1., 1j
            arrays = sector_density_moments(frame, state, coupling, check_direction=direction)
            arrays['check_direction'] = direction
            arrays['mesh_index'] = np.asarray(slot.mesh_index, dtype=np.int64)
            arrays['measure'] = np.asarray(slot.weight)
            for name in ('number_pair', 'gradient_action_pair'):
                left, right = arrays[name]
                if abs(left-right) > 1e-8*(1+abs(right)):
                    raise ValueError('actual full-column trace pairing failed: '+name)
            run.commit(key+'-density-gradient', arrays, slot_metadata)
            source_density = arrays['density_bohr3']
            del arrays
        run.event('slot_hamiltonian_start', mesh_index=list(slot.mesh_index))
        arrays = sector_hamiltonian_moments(frame, state, hamiltonian)
        run.commit(key+'-hamiltonian', arrays, dict(slot_metadata,
            stationarity='actual full-row commutator norm retained; not imposed'))
        run.event('slot_frozen_curvature_start', mesh_index=list(slot.mesh_index))
        arrays = conditional_minus_mv(frame.coefficients, state, coupling,
            frame.raw_occupations, conditional_density=source_density)
        run.commit(key+'-frozen-curvature', arrays, dict(slot_metadata,
            boundary='full-column-source-frozen-electronic-curvature',
            units=dict(slot_metadata['units'], curvature='Ry/bohr^2'),
            curvature_scope='conditional and original numerical-MV nonlocal blocks separate; orbital local difference is diagnostic'))
        run.summary.setdefault('slot_seconds', {})[key] = time.monotonic()-began
        del coupling, hamiltonian, arrays, source_density
        run.check()
    run.summary.update(status='complete', selected_slots=len(selected),
        full_material_total=False, final_implementation_sources=implementation(paths))
    verify_consumers(run.grant, run.summary['final_implementation_sources'])

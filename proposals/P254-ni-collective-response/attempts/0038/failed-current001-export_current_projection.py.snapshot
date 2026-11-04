"""Export three raw1038 physical-current matrices for the independent A38 consumer.

No temporal integration occurs here. The only physical actions are three charged
J(original346) re-evaluations and three new J(extra692) evaluations. Full-PW
current action arrays are transient; only the small matrix bank is persisted.
"""
import argparse
import json
import os
from pathlib import Path
import resource
from time import perf_counter
import traceback

from current_projection_inputs import (CONTRACT, CURRENT_NAMES, HERE, ROOT,
    evaluate_currents, field_map, initial_comparison, load_handoff,
    load_projection_bank, restore_inputs)
from checkpoint_directories import make_directory
from expanded_bath_inputs import SOURCE_COMPATIBILITY
from expanded_influence_readout import save_json, save_npz
from export_expanded_bath import WorkJournal
from substrate_framework.streamed_material_state import file_binding


SCOPE = (
    'Original source0/Gamma/first Ni1-x output; unchanged A37 original346 plus692 '
    'empty directions and nonstationary original S/D/R. Raw kinetic+physical-UPF '
    'J=-L/2 in the same saved finite-bath coordinates, with qualified Jdagger=J '
    'aliases, all25 initial complex F/J covariances and an explicit initial '
    'finite-bath current projection deficit. No new H/F actions, basis selection, '
    'trajectory, temporal integration or spectral reconstruction. This producer '
    'does not establish joint temporal influence, late-time accuracy, physical '
    'loss, convergence, corrected A32 totals, full moving/backreaction/RPA, '
    'all128-source/64q closure or whole physical E/C/N completion.'
)


def implementation_bindings():
    """Bind the implementation actually used, separately from cached producers."""
    own = [HERE/name for name in ('export_current_projection.py', 'current_projection_inputs.py')]
    reused = {
        '0037': ('expanded_bath_inputs.py', 'export_expanded_bath.py', 'expanded_influence_readout.py'),
        '0036': ('canonical_current_noise.py',),
        '0031': ('cached_bath_inputs.py', 'export_frozen_bath_comparator.py'),
        '0030': ('packet_hamiltonian_actions.py', 'original_noise_source.py',
                 'export_output_directed_actions.py', 'canonical_force_filter.py'),
        '0034': ('checkpoint_directories.py',),
        '0028': ('clamped_scientific_readout.py',)}
    dependencies = [HERE.parent/attempt/name for attempt, names in reused.items() for name in names]
    dependencies += [ROOT/'src'/'substrate_framework'/name for name in (
        'physical_factor_preparation.py', 'streamed_material_state.py',
        'material_physical_operators.py', 'material_frozen_action.py')]
    return [file_binding(path) for path in own+dependencies]


def work_record(inputs, calls, phases, started):
    if (len(calls) != 6 or any(call['status'] != 'complete' for call in calls)
            or [call['columns'] for call in calls] != [346, 692]*3
            or inputs.extra_calls):
        raise ValueError('exactly six completed current evaluations and zero H/F callbacks required')
    return dict(phases, total_seconds=perf_counter()-started,
        peak_process_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*1024,
        cpu_affinity=sorted(os.sched_getaffinity(0)),
        thread_environment={name: os.environ.get(name) for name in (
            'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')},
        original_context_restoration=inputs.context.operator_work,
        unused_ionic_direction_factory_setup=inputs.direction_preparations,
        unused_ionic_direction_factory_setup_seconds=sum(
            row['elapsed_seconds'] for row in inputs.direction_preparations),
        current_actions=calls, current_action_evaluations=6,
        original_current_action_reevaluations=3, original_current_columns_reevaluated=3*346,
        original_current_reevaluation_reason='A35 did not retain action arrays; necessary new matrix elements, intentionally charged, not cache reuse',
        new_extra_current_action_evaluations=3, new_extra_current_action_columns=3*692,
        total_current_action_columns=sum(row['columns'] for row in calls),
        current_action_seconds=sum(row['elapsed_seconds'] for row in calls),
        qualified_adjoint_reuses_by_hermiticity=6, independent_current_adjoint_actions=0,
        native_current_internal_column_blocks=sum(row['internal_column_blocks'] for row in calls),
        counter_scope='six logical current evaluations; 64-column native blocking, not six separate adjoint evaluations',
        reused_cached_original_H_F_action_legs=6, original_H_action_replays=0,
        original_force_action_replays=0, extra_H_actions=0, extra_force_actions=0,
        material_preparation_replays=0, population_refills=0, original_gram_recomputations=0,
        basis_selection_replays=0, hamiltonian_eigensolver_replays=0,
        bath_spectral_reconstructions=0, trajectory_actions=0, temporal_integrations=0,
        initial_car_admissibility_evaluations=2, initial_car_small_gram_eigh_calls=2,
        initial_car_small_density_eigvalsh_calls=2,
        admissibility_scope='Unchanged A31 CAR consumer validates supplied S/D; no H/basis/spectral eigensolver replay',
        current_action_arrays_written=False, original_banks_copied=False)


def complete_manifest(contract, bindings, data, inputs, bank_binding, specifications,
                      initial, field, work, journal_binding):
    original = data['basis_complete']
    context = inputs.context
    polarizations, matrix, amplitude = field
    result = dict(schema='P254-A38-current-projection-v1', status='complete',
        contract=dict(CONTRACT), **{key: bindings[key] for key in (
            'basis_complete', 'basis', 'a37_filtered_complete', 'a37_filtered_operators',
            'a37_temporal_result', 'current_reference')},
        current_bank=bank_binding, basis_identity=original['basis_identity'],
        original_state_authority=context.state_authority,
        original_factor_identity=context.factor_identity,
        original_pw_basis_identity=context.pw_basis_identity,
        original_source=inputs.original_source, restoration_source=inputs.restoration_source,
        source_compatibility=dict(SOURCE_COMPATIBILITY), original_frame=context.frame,
        original_coordinates=context.source['coordinates'],
        cached_original_actions=inputs.cached_action_records,
        cached_basis_implementation=original['implementation'],
        cached_basis_operator_specifications=original['operator_specifications'],
        source_record=0, quadrature_weight=context.quadrature_weight, ncells=inputs.ncells,
        polarizations=polarizations, current_force_matrix=matrix, source_field_amplitude=amplitude,
        a36_field_reference=bindings['a36_field_reference'],
        initial_covariance=initial, current_specifications=specifications,
        current_adjoint_aliases={name: name for name in CURRENT_NAMES},
        current_array_contract=dict(keys=list(CURRENT_NAMES), dtype='complex128', shape=[1038,1038],
            units='electronic_current_hartree', basis='B=[Psi S^-1/2,E], original346 then extra692',
            adjoints='same actual arrays under checked real direction/Hermitian original UPF D, not manufactured matrix.conj().T',
            raw_matrix_policy='No Hermitization, clipping, rank/quality cutoff, averaging or defect removal'),
        normalization=dict(current='raw J_Ha=-L_Ry/2 including kinetic and physical UPF',
            ionic='reuse paid integrated A37 F/sqrt(N) operators without further scaling',
            field='source amplitude sqrt(N*w0)*M; N64*w0=1 so amplitude is M',
            current_bank_quadrature_applied=False, current_bank_canonical_field_map_applied=False,
            additional_spin_factor=False, additional_mass_or_oscillator_factor=False,
            independent_marginal_noise_addition=False),
        work=work, implementation=implementation_bindings(), work_journal=journal_binding,
        recovery=dict(complete_current_bank_and_manifest_reusable_without_current_reevaluation=True,
            unfinished_current_arrays_resumable=False,
            maximum_lost_current_evaluations_before_complete_marker=6,
            scope='Only complete authenticated current bank plus manifest is reusable. '
                'Partial arrays are not checkpoints; original A37 banks stay referenced, never copied.'),
        scope=SCOPE)
    missing = set(contract['handoff']['required_manifest_fields'])-set(result)
    if missing:
        raise ValueError('missing fixed current-projection metadata: '+repr(sorted(missing)))
    return result


def run(args):
    directory = make_directory(Path(args.output).resolve(), exist_ok=False)
    started = perf_counter()
    calls, phases = [], {}
    inputs = None
    stage = 'authenticate_fixed_handoff'
    try:
        with WorkJournal(directory/'physical-work.jsonl') as progress:
            progress(dict(stage=stage, status='started'))
            t = perf_counter()
            contract, bindings, data = load_handoff(args)
            phases['authentication_seconds'] = perf_counter()-t
            stage = 'restore_original_context_and_cached_H_F_legs'
            progress(dict(stage=stage, status='started'))
            t = perf_counter()
            inputs = restore_inputs(data['basis_complete'], data['current_reference'], progress)
            phases['restoration_and_cached_loading_seconds'] = perf_counter()-t
            bank = load_projection_bank(bindings['basis'], data['basis_complete'], inputs)
            field = field_map(data['a36_field_reference'], inputs)
            progress(dict(stage=stage, status='complete', cached_action_legs=6, H_F_actions=0))
            stage = 'evaluate_and_project_six_current_actions'
            t = perf_counter()
            matrices, currents, specifications = evaluate_currents(inputs, bank,
                data['current_reference'], progress, calls)
            phases['current_factory_actions_and_projection_seconds'] = perf_counter()-t
            stage = 'save_current_matrix_bank'
            progress(dict(stage=stage, status='started'))
            t = perf_counter()
            bank_binding = save_npz(directory/'current-operators.npz', matrices)
            phases['current_bank_saving_seconds'] = perf_counter()-t
            progress(dict(stage=stage, status='complete', current_bank=bank_binding))
            stage = 'initial_full_pw_and_finite_bath_25_pair_CAR'
            progress(dict(stage=stage, status='started'))
            t = perf_counter()
            initial = initial_comparison(inputs, bank, matrices, currents,
                data['basis_complete']['basis_identity'], data['current_reference'])
            phases['initial_25_pair_covariance_seconds'] = perf_counter()-t
            del currents
            progress(dict(stage=stage, status='complete', ordered_pairs=25,
                          elapsed_seconds=phases['initial_25_pair_covariance_seconds']))
            work = work_record(inputs, calls, phases, started)
        stage = 'publish_complete_current_projection'
        manifest = complete_manifest(contract, bindings, data, inputs, bank_binding,
            specifications, initial, field, work, file_binding(directory/'physical-work.jsonl'))
        binding = save_json(directory/'current-projection-complete.json', manifest)
    except BaseException as error:
        failure = dict(schema='P254-A38-current-projection-failure-v1', status='failed',
            stage=stage, error=repr(error), traceback=traceback.format_exc(),
            elapsed_seconds=perf_counter()-started, contract=dict(CONTRACT),
            current_actions=calls, completed_phases=phases, unfinished_current_arrays_resumable=False,
            original_context_restoration=inputs.context.operator_work if inputs is not None else None,
            scope='Failure retained in unique output directory. No complete current projection '
                'or resumable current-action arrays claimed; no automatic retry or old cost reset.')
        try:
            save_json(directory/'current-projection-failed.json', failure)
        except BaseException as persistence_error:
            error.add_note('Failure persistence also failed: '+repr(persistence_error))
        raise
    print(json.dumps(dict(current_projection_complete=binding, current_bank=bank_binding,
                          total_seconds=perf_counter()-started), allow_nan=False), flush=True)
    return binding


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('basis-complete', 'a37-filtered-complete', 'a37-temporal-result',
                 'current-reference', 'output'):
        parser.add_argument('--'+name, required=True)
    run(parser.parse_args())


if __name__ == '__main__':
    main()

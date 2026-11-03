"""Run the existing original-column mean engine after saved-case restoration.

No numerical equation or historical source file is replaced here. The original
model authority and current execution/source identity remain distinct.
"""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from joint_column_mean import ColumnLimits, ColumnPartitionedMean
from gramian_column_mean import bounded_source_columns
from matched_joint_initial_reuse import import_initial_contractions
from run_joint_column_mean_oracle import sources as column_sources
from run_joint_partitioned_material import (
    implementation_bindings, persistent_json, diagnostics, source_identity)
from saved_case_mean_authority import restore_completed_mean, restore_saved_case
from substrate_framework.partitioned_joint_mean import save_stage
from substrate_framework.streamed_material_state import authenticate, file_binding

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def numerical_authority(midpoint_binding, column_binding):
    midpoint = json.loads(authenticate(midpoint_binding).read_text())
    if (midpoint.get('schema') != 'P254-joint-partitioned-oracle-v1'
            or midpoint.get('status') != 'passed'):
        raise ValueError('executed partitioned physical midpoint proof required')
    for binding in midpoint['source_bindings']:
        if Path(binding['path']).parent == ROOT/'src/substrate_framework':
            authenticate(binding)
    column = json.loads(authenticate(column_binding).read_text())
    if (column.get('schema') != 'P254-joint-column-controller-oracle-v1'
            or column.get('status') != 'passed'
            or column['source_bindings'] != column_sources()):
        raise ValueError('current original-column propagation proof required')


def seed_sources():
    return implementation_bindings()+[file_binding(HERE/name) for name in (
        'saved_case_mean_authority.py', 'saved_case_mean_execution.py',
        'saved_case_mean_owner.py', 'run_saved_case_mean.py',
        'gramian_column_mean.py', 'gramian_krylov_bounds.py')]


def resource_plan(context, state, args):
    sets = 3 if args.mode == 'refine' else args.steps if args.mode == 'evolve' else 0
    return dict(workers=len(args.cpus), simultaneous_slots=len(args.cpus),
        original_frame_bytes=sum(int(slot.frame['bytes']) for slot in state.slots),
        max_single_frame_bytes=max(int(slot.frame['bytes']) for slot in state.slots),
        retained_new_frame_sets=sets,
        complex_contraction_payload_bytes_per_complete_barrier=
            len(state.slots)*16*(context.reference.nbar.shape[0]+len(state.u)+5),
        maximum_reduction_barriers_before_cache_reuse=1+sets*(args.max_iterations+2),
        storage_estimate_note='Existing mean-engine payload formula; excludes headers, JSON, stages, reservations and failed attempts. Hard allocation remains authoritative.',
        midpoint_storage='per-slot contractions, not trial full-frame sets')


def refinement(first, second, budget):
    error = 0.
    for left, right in zip(first.slots, second.slots, strict=True):
        budget.check()
        error = max(error, float(np.linalg.norm(first.arrays(left)[0]-second.arrays(right)[0])))
    return dict(max_frame_frobenius=error,
        coordinate_max=max(float(np.max(np.abs(getattr(first,name)-getattr(second,name))))
                           for name in ('u', 'p', 'coordinate', 'momentum')),
        accuracy_certified=False, interval_not_an_accuracy_license=True)


def execute(args, directory, inputs, sources, budget, result):
    numerical_authority(inputs['oracle'], inputs['column_oracle'])
    context, state, provenance = restore_saved_case(args.case, inputs['initialization'],
        inputs['qualification'], reality_tolerance=args.reality_tolerance)
    if inputs.get('mean_continuation') is not None:
        state, provenance['continuation'] = restore_completed_mean(
            args.case, state, inputs['initialization'], inputs['mean_continuation'],
            branch_index=args.start_from_mean_branch)
    limits = ColumnLimits(args.column_dimension, args.column_workspace_bytes, args.column_error_target)
    source = sha256((source_identity(sources)+json.dumps(limits.__dict__, sort_keys=True)).encode()).hexdigest()
    descriptor = dict(identity=context.identity, specification=context.specification,
        implementation=sources, source_identity=source, initialization=inputs['initialization'],
        column_limits=limits.__dict__, block_size=args.block_size,
        reality_tolerance=args.reality_tolerance, current_execution_qualification=inputs['qualification'])
    persistent_json(directory/'model.json', descriptor, budget)
    result.update(case=args.case, oracle=inputs['oracle'], column_oracle=inputs['column_oracle'],
        initialization=inputs['initialization'], initial_boundary=file_binding(state.path),
        qualified_saved_authority=provenance, source_identity=source,
        coefficient_defect_upper_bound=None, positive_time_symmetry=False,
        resource_plan=resource_plan(context, state, args))
    serial = context.make_mean(budget, block_size=args.block_size)
    model = ColumnPartitionedMean.from_serial(serial, directory=directory/'work',
        cpus=args.cpus, source_identity=source, limits=limits)
    model.exponential = bounded_source_columns
    if inputs['initial_reuse'] is not None:
        result['initial_contraction_reuse'] = import_initial_contractions(
            model, state, inputs['initial_reuse'], block_size=args.block_size)
    initial = model.reduce(state)
    result.update(initial_reduction=save_stage(initial,
        directory/('initial-reduction-'+file_binding(state.path)['sha256']), budget,
        metadata=dict(kind='full-initial-reduction', boundary=str(state.path), time_ry_inverse=state.time)),
        initial_energy_total_ry=initial.energy, initial_number_per_cell=initial.number)
    if args.mode == 'reduce':
        return
    finals, result['branches'] = [], []
    for steps in ((1, 2) if args.mode == 'refine' else (args.steps,)):
        final, _ = model.evolve(state, directory/('interval-%s-steps-%d' % (args.interval, steps)),
            physical_interval=args.interval, steps=steps, max_steps=args.max_steps,
            max_iterations=args.max_iterations, tolerance=args.tolerance,
            refinement_label='steps-'+str(steps))
        stage = model.reduce(final)
        result['branches'].append(dict(steps=steps, boundary=file_binding(final.path),
            diagnostics=diagnostics(state, final, initial, stage, budget)))
        finals.append(final)
    if len(finals) == 2:
        result['refinement'] = refinement(*finals, budget)

#!/usr/bin/env python3
"""One original full-column Native source: direct split and H0 action feasibility.

This does not propagate the coupled material, extend its mean, or establish
spatial/coherence coverage. The fixed original-start allocation is never retried
or enlarged by this driver. All existing physical inputs remain read-only.
"""
import argparse
import json
import os
from pathlib import Path
import sys
from time import monotonic, time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'

import numpy as np
from joint_column_mean import ColumnLimits
from matched_joint_case import load_case
from moving_exponential_h0 import PhysicalH0Exponential
from moving_exponential_material import DirectMeanResidual
from run_joint_q_material import bind_saved_model
from run_mixed_q_material import watchdog
from saved_case_mean_owner import terminal_signals
from substrate_framework.material_control_reference import bound_json
from substrate_framework.material_joint_response import RealBlochEdge
from substrate_framework.material_physical_operators import (
    physical_hamiltonian, physical_coupling, physical_finiteq,
    PhysicalFieldVertices, PhysicalFrozenAction)
from substrate_framework.partitioned_joint_budget import PartitionBudget
from substrate_framework.streamed_joint_reduction import JointStage as MeanStage
from substrate_framework.streamed_joint_state import DiskJointMean
from substrate_framework.streamed_material_state import authenticate, file_binding


class ColumnBudget(PartitionBudget):
    def __init__(self, *, max_columns, **caps):
        super().__init__(**caps)
        self.max_columns = max_columns

    def action(self, vectors):
        columns = 1 if vectors.ndim == 1 else vectors.shape[1]
        if self.columns+columns > self.max_columns:
            raise RuntimeError('original cumulative H0 action-column allocation exhausted')
        super().action(vectors)


def reuse_completed_split(plan, result):
    prior = bound_json(plan['cumulative']['prior'])
    if prior.get('H0_propagations') != []:
        raise ValueError('this continuation must not replay a completed H0 propagation')
    for key in ('model_identity', 'source_index', 'original_frame', 'original_D', 'saved_mean_stage'):
        if prior[key] != result[key]:
            raise ValueError('prior successful split refers to different '+key)
    old_plan = bound_json(prior['registration'])
    changed = {str(HERE/'moving_exponential_h0.py'), str(Path(__file__).resolve())}
    for binding in old_plan['source_bindings']:
        if binding['path'] not in changed:
            authenticate(binding)
    comparisons = prior['split_comparisons']
    if any(comparisons[key]['maximum_difference'] > comparisons[key]['comparison_threshold']
           for key in ('False', 'True')):
        raise ValueError('prior split comparison was not successful')
    result.update(split_comparisons=comparisons,
                  split_reuse=plan['cumulative']['prior'])


def run(plan, directory, budget, result):
    initialization = bound_json(plan['initialization'])
    saved = initialization['saved']
    qualification = bound_json(plan['physical_qualification'])
    context = bind_saved_model(load_case(plan['case']), saved,
                               qualification.get('source_compatibility'))
    initial = DiskJointMean(authenticate(saved['joint_boundary']),
                           expected_identity=context.identity, verify=False)
    mean = bound_json(plan['mean_receipt'])
    final = DiskJointMean(authenticate(mean['branches'][0]['boundary']),
                         expected_identity=context.identity, verify=False)
    index = plan['source_index']
    slot = initial.slots[index]
    original = context.grid.load_frame(context.ordered[index])
    basis = original.basis
    if basis.identity != slot.basis_identity or original.gauge_identity != slot.gauge_identity:
        raise ValueError('original source basis/gauge differs')
    del original
    psi = np.load(authenticate(slot.frame), mmap_mode='r', allow_pickle=False)
    authenticate(slot.coefficient)
    result.update(model_identity=context.identity, source_index=index,
        full_frame_shape=list(psi.shape), original_frame=slot.frame,
        original_D=slot.coefficient, all_original_columns=True,
        saved_mean_stage=final.history[-1])
    h0 = physical_hamiltonian(context.background, basis, basis.slot.ispin, plan['block_size'])
    action = PhysicalH0Exponential(h0, budget=budget,
        limits=ColumnLimits(**plan['column_limits']), identity=context.identity+':'+basis.identity,
        receipt_directory=directory/'exponentials', fft_block_size=plan['block_size'])
    if 'cumulative' in plan:
        reuse_completed_split(plan, result)
    else:
        verify_split(context, basis, stage=MeanStage.restore(final.history[-1]),
                     h0=h0, psi=psi, action=action, budget=budget, result=result,
                     block_size=plan['block_size'])
    result['H0_propagations'] = []
    for interval in plan['intervals']:
        before = budget.snapshot()
        started = monotonic()
        evolved = action.propagate(psi, interval)
        result['H0_propagations'].append(dict(interval_ry_inverse=interval,
            full_columns=True, elapsed_seconds=monotonic()-started,
            actions_before=before, actions_after=budget.snapshot()))
        del evolved


def verify_split(context, basis, *, stage, h0, psi, action, budget, result, block_size):
    gamma = physical_coupling(context.background.inputs, basis, h0.factors, block_size)
    edge = RealBlochEdge.create(context.background.inputs.reference, basis, basis, np.zeros(3))
    rectangular = physical_finiteq(context.background.inputs, edge, h0.factors, h0.factors, block_size)
    field = PhysicalFieldVertices(h0).prepare(stage.shift) if np.any(stage.shift) else None
    residual = DirectMeanResidual(rectangular, stage.potential[:, basis.slot.ispin-1], stage.u, field)
    full = PhysicalFrozenAction(h0, gamma, potential=stage.potential[:, basis.slot.ispin-1],
                                direction=stage.u, shift=stage.shift)
    result['split_comparisons'] = {}
    for adjoint in (False, True):
        budget.action(psi)
        full_value = full.adjoint_apply(psi) if adjoint else full.apply(psi)
        budget.action(psi)
        bare = action.frozen.adjoint_apply(psi) if adjoint else action.frozen.apply(psi)
        part = residual.adjoint_apply(psi) if adjoint else residual.apply(psi)
        difference = float(np.max(np.abs(full_value-(bare+part))))
        threshold = 1e-11*(1+float(np.max(np.abs(full_value))))
        result['split_comparisons'][str(adjoint)] = dict(maximum_difference=difference,
            comparison_threshold=threshold, direct_residual_norm=float(np.linalg.norm(part)),
            scope='operator decomposition/roundoff comparison, not material accuracy')
        if difference > threshold:
            raise AssertionError('direct physical mean residual differs from settled full H')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads(args.registration.read_bytes())
    for binding in plan['source_bindings']:
        authenticate(binding)
    directory = Path(plan['directory'])
    directory.mkdir(parents=True, exist_ok=False)
    caps = dict(plan['budget'])
    cumulative = plan.get('cumulative')
    prior = bound_json(cumulative['prior']) if cumulative else None
    if cumulative:
        caps['max_seconds'] = cumulative['original_started_unix']+cumulative['max_seconds']-time()
        caps['max_actions'] = cumulative['max_actions']-prior['budget']['hamiltonian_actions']
        caps['max_bytes'] -= int(cumulative['prior']['bytes'])
        max_columns = cumulative['max_columns']-prior['budget']['action_columns']
    else:
        max_columns = plan['max_columns']
    budget = ColumnBudget(max_columns=max_columns, **caps)
    budget.reserve(2*1024**2)
    result = dict(schema='P254-moving-exponential-Native-smoke-v1', status='partial',
        registration=file_binding(args.registration), pid=str(os.getpid()),
        affinity=sorted(os.sched_getaffinity(0)),
        scope='One original physical source/operator feasibility; not a coupled trajectory, mean extension or coherence window')
    print('P254_MOVING_EXP_NATIVE_START '+json.dumps(result), flush=True)
    try:
        if result['affinity'] != [plan['cpu']]:
            raise ValueError('Native probe does not own its declared singleton CPU')
        with terminal_signals(), watchdog(budget):
            run(plan, directory, budget, result)
        result['status'] = 'complete'
    except BaseException as error:
        result['error'] = dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc())
    finally:
        result['budget'] = budget.snapshot()
        if cumulative:
            result['cumulative_budget'] = dict(
                original_started_unix=cumulative['original_started_unix'],
                elapsed_seconds=time()-cumulative['original_started_unix'],
                hamiltonian_actions=prior['budget']['hamiltonian_actions']+budget.actions,
                action_columns=prior['budget']['action_columns']+budget.columns,
                max_seconds=cumulative['max_seconds'], max_actions=cumulative['max_actions'],
                max_columns=cumulative['max_columns'], prior=cumulative['prior'])
        result['source_bindings'] = [file_binding(HERE/name) for name in (
            'moving_exponential_h0.py', 'moving_exponential_material.py', Path(__file__).name)]
        with (directory/'result.json').open('x') as stream:
            json.dump(result, stream, indent=2, allow_nan=False)
            stream.write('\n')
        print(('[dft-done]' if result['status'] == 'complete' else '[dft-fail]')
              +' moving-exponential-native '+json.dumps(result), flush=True)
    return 0 if result['status'] == 'complete' else 2


if __name__ == '__main__':
    raise SystemExit(main())

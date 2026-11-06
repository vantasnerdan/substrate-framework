#!/usr/bin/env python3
"""Finish the wall-limited storage oracle without replaying admitted physics.

The original failed receipt remains unchanged. Its passed checks are admitted
only against identical source bytes; the remaining non-Gamma scheduled pair,
worker recovery, forward/dual stage recovery and live accounting run here.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from run_joint_partitioned_q_oracle import (
    CHECKS, DenseTangent, DiskFiniteHistory, FiniteSystem, FockOracle,
    PartitionedQResponse, PartitionedTwoSidedCovariance, JointProductPreparation,
    QBudget, MappedResponse, bindings, check, encode, file_binding, mapped,
    pack_legs, recovery, stage_recovery)
from substrate_framework.streamed_material_state import authenticate


def remaining_pair(system, fock, directory, budget, cpus):
    response = system.response(False)
    dense = DenseTangent(system, response)
    executor = PartitionedQResponse(response, DiskFiniteHistory(system, response, budget),
        kernel_adjoint=lambda n, **kw: system.k*n, grid_measure=1., budget=budget,
        identity='nonGamma-storage-completion', cpus=cpus,
        directory=directory/'executor', work_directory=directory/'work')
    preparation = JointProductPreparation(response, system.stage(0),
        ionic_ordered=system.wi, ionic_minus_ordered=system.wi,
        field_ordered=np.zeros((0, 0)),
        electronic_identity='raw-singular-D', ionic_identity='initial-ionic-W',
        field_identity='initial-field-W', product_identity='original-product')
    left = dense.probe(system.random((dense.size,)), 'completion-A')
    right = dense.probe(system.random((dense.size,)), 'completion-B')
    protected = []
    retire = executor.retire_frames_except

    def checked_retirement(needed, barrier):
        keep = []
        for candidate in executor.frames.directory.glob('*/generation-*'):
            if (candidate/'retire.json').exists() or not (candidate/'request.json').exists():
                continue
            request = json.loads((candidate/'request.json').read_text())
            if request['input']['time'] in needed:
                for receipt in candidate.glob('item-*/complete.json'):
                    keep.extend(json.loads(receipt.read_text())['arrays'].values())
        retire(needed, barrier)
        for binding in keep:
            executor.frames.verify(binding)
        protected.append(len(keep))

    executor.retire_frames_except = checked_retirement
    t, s = .027, .019
    actual = PartitionedTwoSidedCovariance(preparation, executor).pair(
        left, right, np.linspace(0, t, 5), np.linspace(0, s, 5), directory/'pair')
    pa, pb = dense.pullback(t, left), dense.pullback(s, right)
    for ordering in ('AB', 'BA'):
        expected = sum(fock.parts(response, pa, pb, ordering).values())
        check('nonGamma_scheduled_frame_retirement_'+ordering,
              actual[ordering]['total'], expected, tolerance=2e-8)
    if not sum(protected) or list(executor.frames.directory.rglob('*.npy')):
        raise AssertionError('future frame protection or last-use retirement not exercised')
    CHECKS.append(dict(name='nonGamma_both_leg_frame_last_use_schedule', passed=True,
        protected_future_bindings=sum(protected), final_live_frame_files=0))
    checkpoint = actual['left_checkpoints'][-1]
    original, _ = MappedResponse.restore(response, checkpoint, budget, executor.identity)
    fresh = PartitionedQResponse(response, DiskFiniteHistory(system, response, budget),
        kernel_adjoint=lambda n, **kw: system.k*n, grid_measure=1., budget=budget,
        identity=executor.identity, cpus=cpus, directory=directory/'fresh',
        work_directory=directory/'work')
    before = budget.actions
    resumed, _ = PartitionedTwoSidedCovariance(preparation, fresh).pullback(
        left, np.linspace(0, t, 5), directory/'resumed', checkpoint=checkpoint)
    check('nonGamma_disk_dual_restart', pack_legs(resumed), pack_legs(original), tolerance=0.)
    check('nonGamma_restart_reuses_completed_H_work', budget.actions, before, tolerance=0.)
    writes = budget.snapshot()['cumulative_write_bytes_reserved']
    automatic, _ = PartitionedTwoSidedCovariance(preparation, fresh).pullback(
        left, np.linspace(0, t, 5), directory/'automatic')
    check('nonGamma_automatic_interval_restart', pack_legs(automatic), pack_legs(original), tolerance=0.)
    check('nonGamma_automatic_restart_reuses_completed_H_work', budget.actions, before, tolerance=0.)
    check('nonGamma_automatic_restart_without_seed_copy',
          budget.snapshot()['cumulative_write_bytes_reserved'], writes, tolerance=0.)
    size = len(pack_legs(response.zero_tangent()))
    x = mapped(response, directory/'mutable-x', budget, system.random((size,)))
    ax = mapped(response, directory/'Ax', budget)
    recovered = mapped(response, directory/'recovered-Ax', budget)
    executor.rhs(.0173, x, ax)
    fresh.rhs(.0173, x, recovered)
    check('nonGamma_fresh_disk_stage_recovery', pack_legs(recovered), pack_legs(ax), tolerance=0.)
    x.ket[0][0, 0] += .1j
    fresh.rhs(.0173, x, recovered)
    check('nonGamma_mutable_scratch_not_completed_input',
          pack_legs(recovered), pack_legs(ax), mutation=True)
    return response


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous', type=Path, required=True)
    parser.add_argument('--previous-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--cpus', nargs='+', type=int, choices=range(6, 11), required=True)
    parser.add_argument('--max-seconds', type=float, required=True)
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpus[0]})
    prior_binding = file_binding(args.previous)
    if prior_binding['sha256'] != args.previous_sha256:
        raise ValueError('explicit prior proof digest differs')
    prior = json.loads(args.previous.read_text())
    sources = bindings()
    if (prior['source_bindings'] != sources or prior['status'] != 'failed'
            or prior.get('error') != 'JointBudgetExceeded: invocation wall-time allocation exhausted'
            or not prior['checks'] or not all(row['passed'] for row in prior['checks'])
            or prior['checks'][-1]['name'] != 'nonGamma_N_weight_sensitive'):
        raise ValueError('exact completed prefix and unchanged implementation required')
    for binding in sources:
        authenticate(binding)
    CHECKS.extend(prior['checks'])
    args.work.mkdir(parents=True, exist_ok=False)
    budget = QBudget(max_actions=1000000, max_seconds=args.max_seconds,
        max_bytes=1000000000, max_memory_bytes=8589934592,
        worker_memory_bytes=1610612736)
    result = dict(schema='P254-joint-partitioned-q-oracle-v1', status='failed',
        source_bindings=sources+[file_binding(__file__)], checks=CHECKS,
        cpus=args.cpus, retained_prefix=prior_binding,
        retained_check_count=len(prior['checks']),
        continuation_scope='No replay of passed original dense/Fock/forward/dual checks')
    try:
        system = FiniteSystem()
        fock = FockOracle(system)
        response = remaining_pair(system, fock, args.work/'nonGamma-scheduled', budget, args.cpus)
        recovery(system, response, args.work/'recovery', budget, args.cpus)
        stage_recovery(system, response, args.work/'stage-recovery', budget, args.cpus)
        counters = budget.snapshot()
        measured = sum(path.stat().st_size for path in args.work.rglob('*') if path.is_file())
        if (counters['retired_payload_bytes'] <= 0
                or counters['cumulative_write_bytes_reserved'] <= counters['live_bytes_reserved']
                or counters['peak_live_bytes_reserved'] < counters['live_bytes_reserved']
                or measured > counters['live_bytes_reserved']):
            raise AssertionError('live/peak/cumulative storage accounting differs from actual ownership')
        CHECKS.append(dict(name='bounded_live_storage_accounting', passed=True,
                          measured_file_bytes=measured, counters=counters))
        result['status'] = 'passed'
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result['budget'] = budget.snapshot()
    with args.output.open('x') as stream:
        json.dump(encode(result), stream, indent=2, allow_nan=False)
        stream.write('\n')
    if result['status'] != 'passed':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()

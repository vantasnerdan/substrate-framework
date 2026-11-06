#!/usr/bin/env python3
"""Compare existing column Arnoldi with an actual saved physical midpoint slot.

A bounded zero-field first-stage pilot, not a new integrator or trajectory.
The unchanged floating analytic envelopes are not interval certificates.
"""
import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
from time import perf_counter
from types import SimpleNamespace

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np
import run_joint_partitioned_material as execution
from run_column_krylov_material import operator_norm, ambient_growth
from run_column_krylov_oracle import sources as column_sources
from substrate_framework.column_frame_krylov import ColumnFrameKrylov
from substrate_framework.material_physical_operators import PhysicalFrozenAction
from substrate_framework.partitioned_joint_slots import stage_identity, digest, array_identity
from substrate_framework.streamed_joint_reduction import JointStage
from substrate_framework.streamed_joint_state import DiskJointMean, JointBudget
from substrate_framework.streamed_material_state import authenticate, file_binding


def bound_json(binding):
    return json.loads(Path(authenticate(binding)).read_text())


def restore(args, budget):
    registration = json.loads(args.registration.read_text())
    if registration['scope'] != 'exact-saved-physical-first-midpoint-slot':
        raise ValueError('registered physical first-stage comparison required')
    for source in registration['sources']:
        authenticate(source)
    oracle = bound_json(registration['column_oracle'])
    if oracle['status'] != 'passed' or oracle['source_bindings'] != column_sources():
        raise ValueError('unchanged checked column mathematics required')
    request = bound_json(registration['request'])
    prior = bound_json(registration['initial_result'])
    comparison = bound_json(registration['comparison'])
    model_record = bound_json(registration['model'])
    if (prior['status'] != 'complete' or not comparison['complete']
            or comparison['request'] != digest(request)
            or comparison['index'] != registration['slot_index']
            or execution.implementation_bindings() != model_record['implementation']):
        raise ValueError('actual complete source-bound physical comparison differs')
    restored, _, ions, joint, coefficients, field, _, identity = execution.settled.load_settled()
    options = SimpleNamespace(block_size=model_record['block_size'],
                              reality_tolerance=model_record['reality_tolerance'])
    model = execution.settled.build_model(restored, ions, joint, coefficients, field, identity, options, budget)
    state = DiskJointMean(authenticate(prior['initial_boundary']), expected_identity=identity, verify=False)
    stage = JointStage.restore(prior['initial_reduction'])
    slot = state.slots[registration['slot_index']]
    if (request['boundary'] != file_binding(state.path) or request['model_identity'] != identity
            or stage_identity(stage) != request['frozen_stage']
            or json.loads(json.dumps(asdict(slot))) != request['slots'][registration['slot_index']]
            or comparison['slot'] != slot.identity
            or request['reduction_coordinates'] != [array_identity(x) for x in
                model.midpoint_coordinates(state, stage, 2*request['interval'])]):
        raise ValueError('operator/state/stage differs from exact stalled midpoint request')
    authenticate(slot.frame)
    authenticate(slot.coefficient)
    return registration, request, comparison, model, state, stage, slot


def probe(args, budget, result):
    registration, request, comparison, model, state, stage, slot = restore(args, budget)
    psi, d = state.arrays(slot)
    initial_gram = psi.conj().T @ psi
    result.update(registration=file_binding(args.registration), original_slot=asdict(slot),
                  interval=request['interval'], scipy_elapsed_seconds=comparison['elapsed_seconds'])
    with model.operators(slot) as (coupling, h0, vertex):
        frozen = PhysicalFrozenAction(h0, coupling, potential=stage.potential[:, slot.component],
                                      direction=stage.u, shift=stage.shift)
        if frozen.field is not None:
            raise ValueError('this registered first-stage norm pilot requires exactly zero field shift')
        norm, norm_record = operator_norm(frozen.frozen)
        growth, growth_record = ambient_growth(frozen.frozen)
        started = perf_counter()
        column = ColumnFrameKrylov(model.counted(frozen.apply), psi,
            max_dimension=args.max_dimension, action_cap=args.max_dimension*psi.shape[1],
            max_memory_bytes=args.workspace_bytes, ambient_growth_bound=growth,
            operator_norm_bound=norm, operator_identity=digest(dict(request=request, slot=slot.identity)))
        result.update(norm_bound=norm_record, growth_bound=growth_record, memory_plan=column.memory_plan)
        result['stages'] = []
        for index in range(args.max_dimension):
            evolved, envelope = column.propagate(request['interval'])
            row = dict(generation=index, elapsed_seconds=perf_counter()-started,
                       calls=column.calls, envelope=envelope,
                       dimensions=column.sizes.tolist())
            budget.json(args.output/('stage-%02d.json' % index), row)
            result['stages'].append(file_binding(args.output/('stage-%02d.json' % index)))
            if envelope['absolute_error_envelope_estimate'] <= args.target:
                break
            expansion = column.expand()
            if expansion['accepted_columns'] == 0:
                raise RuntimeError('bounded original-column expansion did not reach requested estimate')
        else:
            raise RuntimeError('bounded original-column estimate remains above target')
        result['propagation_seconds'] = perf_counter()-started
        result['estimated_error'] = envelope['absolute_error_envelope_estimate']
        result['frame'] = budget.array(args.output/'frame.npy', evolved)
        result['gram_change_frobenius'] = float(np.linalg.norm(evolved.conj().T @ evolved-initial_gram))
        n, g = coupling.prepared_moments(evolved, d)
        ell = vertex.moments(evolved, d)['field_gradient']
        hpsi = model.counted(PhysicalFrozenAction(h0).apply)(evolved)
        values = dict(density=n, gradient=g, ell=ell,
            one_body=np.asarray(np.einsum('ij,ji->', d, evolved.conj().T @ hpsi)),
            number=np.asarray(model.measure*np.sum(n)))
        result['observable_comparison'] = {}
        for name, value in values.items():
            expected = np.load(authenticate(comparison['arrays'][name]), allow_pickle=False)
            if value.shape != expected.shape:
                raise ValueError('saved physical observable shape differs: '+name)
            result['observable_comparison'][name] = dict(
                max_absolute_difference=float(np.max(np.abs(value-expected))),
                frobenius_difference=float(np.linalg.norm(value-expected)),
                reference_norm=float(np.linalg.norm(expected)),
                actual=budget.array(args.output/(name+'.npy'), value),
                scipy=comparison['arrays'][name])
        result['comparison_scope'] = 'actual full raw physical slot; saved SciPy observables, not a retained SciPy frame'
        result['status'] = 'complete'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpu', type=int, choices=range(6, 11), required=True)
    parser.add_argument('--max-seconds', type=float, required=True)
    parser.add_argument('--max-dimension', type=int, required=True)
    parser.add_argument('--workspace-bytes', type=int, required=True)
    parser.add_argument('--target', type=float, required=True)
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    args.output.mkdir(parents=True, exist_ok=False)
    budget = JointBudget(max_actions=32, max_seconds=args.max_seconds, max_bytes=2*1024**3)
    budget.reserve(1024**2)
    result = dict(schema='P254-physical-column-probe-v1', status='partial', cpu=args.cpu,
                  target=args.target, material_trajectory=False, qualification='floating envelope, not interval certified')
    try:
        probe(args, budget, result)
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result['budget'] = budget.snapshot()
    (args.output/'result.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    if result['status'] != 'complete':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()

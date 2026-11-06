#!/usr/bin/env python3
"""Recover actual midpoint residuals from settled contractions, without H0 work.

Use the original stage reducer and residual function. Every reconstructed frozen
stage must match the next actual worker request exactly before its residual is
reported as an executed-stage result. No live mean files are written.
"""
import argparse
import json
import os
from pathlib import Path

from matched_joint_case import load_case
from substrate_framework.partitioned_joint_budget import PartitionBudget
from substrate_framework.partitioned_joint_slots import SlotExecutor, array_identity, stage_identity
from substrate_framework.streamed_joint_state import DiskJointMean
from substrate_framework.streamed_material_state import authenticate, file_binding


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    registration = json.loads(args.registration.read_bytes())
    args.output.mkdir(parents=True, exist_ok=False)
    budget = PartitionBudget(max_actions=1, max_seconds=300, max_bytes=8 << 20,
                             max_memory_bytes=4 << 30, worker_memory_bytes=4 << 30)
    result = dict(schema='P254-saved-midpoint-residual-intake-v1', status='incomplete',
                  registration=file_binding(args.registration), reader=file_binding(__file__),
                  stages=[], residuals=[], scope='Actual homogeneous midpoint iteration; '
                  'not a positive-time boundary, convergence forecast or local-q response.')
    try:
        model = json.loads(authenticate(registration['model']).read_bytes())
        for record in model['implementation']:
            authenticate(record)
        step = json.loads(authenticate(registration['step']).read_bytes())
        requests = [json.loads(authenticate(record).read_bytes())
                    for record in registration['requests']]
        context = load_case(registration['case'], reality_tolerance=model['reality_tolerance'])
        if context.identity != model['identity'] or context.specification != model['specification']:
            raise ValueError('current matched case differs from the executed model')
        state = DiskJointMean(authenticate(requests[0]['boundary']),
                              expected_identity=context.identity, verify=False)
        serial = context.make_mean(budget, block_size=model['block_size'])
        executor = SlotExecutor(serial, args.output/'reader-work',
                                cpus=[registration['cpu']], source_identity=model['source_identity'])
        previous = None
        for index, request in enumerate(requests[:-1]):
            budget.check_processes(os.getpid(), [])
            if (request['source_identity'] != model['source_identity']
                    or request['boundary'] != requests[0]['boundary']
                    or request['operation'] != 'reduce'
                    or request['interval'] != (0. if index == 0 else step['dt']/2)):
                raise ValueError('saved request is not this original midpoint sequence')
            coordinates = ((state.u, state.p, state.coordinate, state.momentum)
                           if previous is None else serial.midpoint_coordinates(state, previous, step['dt']))
            if [array_identity(value) for value in coordinates] != request['reduction_coordinates']:
                raise ValueError('reconstructed coordinates differ from the actual request')
            directory = Path(registration['requests'][index]['path']).parent
            records = executor.collect(directory, directory.name, state)
            stage = executor.stage(directory, directory.name, state, coordinates, records)
            reconstructed = stage_identity(stage)
            expected = requests[index+1]['frozen_stage']
            differences = [key for key in reconstructed if reconstructed[key] != expected[key]]
            evidence = dict(request=registration['requests'][index],
                            next_request=registration['requests'][index+1],
                            exact_original_stage_match=not differences, differing_fields=differences,
                            reconstructed_stage=reconstructed,
                            contraction_receipts=[file_binding(directory/f'slot-{i:04d}'/'complete.json')
                                                  for i in range(len(records))])
            result['stages'].append(evidence)
            if differences:
                raise ValueError('reconstructed stage is not the exact executed stage: '+repr(differences))
            if previous is not None:
                values = serial.residual(previous, stage)
                result['residuals'].append(dict(iteration=index, values=values,
                    maximum=max(values.values()), requested_tolerance=step['tolerance'],
                    accepted=max(values.values()) <= step['tolerance']))
            previous = stage
        if budget.actions != 0:
            raise ValueError('saved-contraction intake unexpectedly performed H0 actions')
        result['status'] = 'complete'
    except Exception as error:
        result['error'] = repr(error)
    finally:
        result['budget'] = budget.snapshot()
        budget.json(args.output/'result.json', result)
        print(json.dumps(dict(status=result['status'], output=file_binding(args.output/'result.json'),
                              residuals=result['residuals'], budget=result['budget'])), flush=True)
        print('[dft-done] A28 saved midpoint intake' if result['status'] == 'complete' else
              '[dft-fail] A28 saved midpoint intake', flush=True)
    return 0 if result['status'] == 'complete' else 1


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        raise SystemExit(main())

#!/usr/bin/env python3
"""Initialize/read back any registered real matched case, without an H0 run.

See matched-joint-case-registration.json for exact commands and limitations.
Arrays must live outside Git. Initialize stores all128 full raw frames and D,
or source-verifies/reuses --original-boundary without copying those arrays.
It writes the local/q/mean ionic covariance and global field covariance.
Read-back independently reloads the bound original source states and compares
every saved frame/D entry; it is not a metadata-only or toy-material proof.
"""
import argparse
import json
import os
from pathlib import Path

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
import matched_reference_io as io
from matched_joint_case import load_case
from matched_joint_sources import registration
from matched_joint_state import initialize, read_back
from substrate_framework.material_control_reference import bound_json
from substrate_framework.streamed_joint_state import JointBudget


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    parser.add_argument('--mode', choices=('initialize', 'read-back'), required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--oracle', type=Path, required=True)
    parser.add_argument('--initialization', type=Path)
    parser.add_argument('--initialization-sha256')
    parser.add_argument('--original-boundary', type=Path)
    parser.add_argument('--original-boundary-sha256')
    parser.add_argument('--cpu', type=int, choices=range(6, 11), required=True)
    parser.add_argument('--max-seconds', type=float, required=True)
    parser.add_argument('--max-output-bytes', type=int, required=True)
    args = parser.parse_args()
    if args.mode == 'read-back' and (args.initialization is None or not args.initialization_sha256):
        parser.error('read-back requires exact initialization path and digest')
    if bool(args.original_boundary) != bool(args.original_boundary_sha256):
        parser.error('existing original boundary requires its exact digest')
    if args.original_boundary and args.mode != 'initialize':
        parser.error('existing original boundary is an initialize input')
    os.sched_setaffinity(0, {args.cpu})
    # Check outside-Git path without creating it: exclusive directory is below.
    output = args.output.resolve()
    if output == io.ROOT or io.ROOT in output.parents:
        parser.error('original numerical state storage must remain outside Git')
    io.runtime()
    budget = JointBudget(max_actions=1, max_seconds=args.max_seconds, max_bytes=args.max_output_bytes)
    budget.reserve(65536)
    output.mkdir(parents=True, exist_ok=False)
    result = dict(schema='P254-matched-joint-case-v1', status='partial', mode=args.mode,
                  case=args.case, cpu=args.cpu, blas_threads=1, material_trajectory=False)
    try:
        reg, _, sources = registration()
        oracle_binding = io.binding(args.oracle)
        oracle = bound_json(oracle_binding)
        # This oracle exercises scalar/preparation joins, not this CLI or state
        # serialization. Their current bytes and actual save/read proof are
        # recorded separately; do not replay unchanged mathematics to bless I/O.
        consumer_paths = {str(Path(__file__).resolve()), str(io.HERE/'matched_joint_state.py')}
        proved = {b['path']: b['sha256'] for b in oracle.get('implementation_sources', [])}
        if (oracle.get('schema') != 'P254-matched-joint-oracle-v1' or oracle.get('status') != 'passed'
                or oracle['registration'] != reg
                or any(proved.get(b['path']) != b['sha256'] for b in sources
                       if b['path'] not in consumer_paths)):
            raise ValueError('executed independent oracle bound to current numerical sources required')
        context = load_case(args.case)
        # Ordinary constructors, no physical operator evaluation or trajectory.
        context.make_mean(budget)
        budget.json(output/'model.json', dict(identity=context.identity, specification=context.specification))
        result.update(oracle=oracle_binding, implementation_sources=sources,
            approximation=context.aggregate['approximation'], coefficient_defect_upper_bound=None,
            nonstationary_initial_state=True, positive_time_symmetry=False,
            conditional_preparation='case-specific sudden-release product; not coupled equilibrium')
        if args.mode == 'initialize':
            original_binding = io.binding(args.original_boundary) if args.original_boundary else None
            if original_binding and original_binding['sha256'] != args.original_boundary_sha256:
                raise ValueError('explicit original boundary digest differs')
            saved = initialize(context, output, budget, original_binding=original_binding)
            result['saved'] = saved
            # Disk constructors authenticate full arrays during wrapping; the
            # separate read-back command additionally regenerates original frames.
            result['proof'] = dict(complete_original_slots=128, disk_boundary_read=True,
                                  exhaustive_source_comparison_command='--mode read-back')
        else:
            record = io.binding(args.initialization)
            if record['sha256'] != args.initialization_sha256:
                raise ValueError('explicit initialization digest differs')
            initialization = bound_json(record)
            if initialization['status'] != 'complete' or initialization['mode'] != 'initialize':
                raise ValueError('completed actual initialization required')
            result['initialization'] = record
            result['proof'] = read_back(context, initialization['saved'], budget)
        result['status'] = 'complete'
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result['budget'] = budget.snapshot()
    raw = (json.dumps(result, indent=2, allow_nan=False)+'\n').encode()
    if len(raw) > 65536:
        raise RuntimeError('reserved terminal capacity exceeded')
    with (output/'result.json').open('xb') as stream:
        stream.write(raw)
    if result['status'] != 'complete':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()

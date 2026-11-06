#!/usr/bin/env python3
"""Case-selected physical joint mean with bounded original-column propagation.

Use an actual run_matched_joint_case.py initialization, its exact digest, the
retained partitioned-midpoint oracle and the new column-controller oracle.
Resume reduce/evolve/refine with the SAME output and unchanged solver limits.
Completed per-slot work is retained; incomplete Arnoldi bases are not resumable.
The existing scalar/midpoint equations and full raw frame/D interface are used.
No legacy loaded-scalar fallback or positive-time symmetry is permitted.
"""
import argparse
import fcntl
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
from uuid import uuid4

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np
from matched_joint_case import load_case
from matched_joint_sources import registration as matched_registration
from joint_column_mean import ColumnLimits, ColumnPartitionedMean
from run_joint_column_mean_oracle import sources as column_sources
from matched_joint_initial_reuse import import_initial_contractions
from substrate_framework.partitioned_joint_budget import PartitionBudget
from substrate_framework.partitioned_joint_mean import save_stage
from substrate_framework.partitioned_joint_slots import publish
from substrate_framework.streamed_joint_state import DiskJointMean
from substrate_framework.streamed_material_state import authenticate, file_binding


def implementation_bindings():
    paths = [Path(__file__), HERE/'joint_column_mean.py', HERE/'matched_joint_initial_reuse.py']
    paths += [ROOT/'src/substrate_framework'/name for name in (
        'partitioned_joint_budget.py', 'partitioned_joint_slots.py', 'partitioned_joint_mean.py',
        'streamed_joint_mean.py', 'streamed_joint_reduction.py', 'streamed_joint_state.py',
        'streamed_material_state.py', 'streamed_material_evolution.py')]
    bindings = {b['path']: file_binding(b['path']) for b in matched_registration()[2]}
    bindings.update({b['path']: b for b in column_sources()})
    for path in paths:
        b = file_binding(path)
        bindings[b['path']] = b
    return [bindings[path] for path in sorted(bindings)]


def source_identity(bindings):
    return sha256(json.dumps({b['path']: b['sha256'] for b in bindings}, sort_keys=True).encode()).hexdigest()


def persistent_json(path, value, budget):
    if path.exists():
        if json.loads(path.read_text()) != value:
            raise ValueError('existing execution identity differs: '+str(path))
    else:
        publish(path, value, budget)


def diagnostics(initial, final, initial_stage, final_stage, budget):
    gram_error = frame_drift = 0.
    for old, new in zip(initial.slots, final.slots):
        budget.check()
        a, da = initial.arrays(old)
        b, db = final.arrays(new)
        if old.coefficient != new.coefficient or old.metadata() != new.metadata():
            raise ValueError('original D/basis/weight changed during mean evolution')
        gram_error = max(gram_error, float(np.linalg.norm(a.conj().T @ a-b.conj().T @ b)))
        frame_drift = max(frame_drift, float(np.linalg.norm(a-b)))
    return dict(energy_total_ry=final_stage.energy, energy_drift_ry=final_stage.energy-initial_stage.energy,
        number_per_cell=final_stage.number, number_drift=final_stage.number-initial_stage.number,
        max_original_column_gram_error=gram_error, max_raw_frame_frobenius_drift=frame_drift,
        density_drift_max=float(np.max(np.abs(final_stage.density-initial_stage.density))),
        canonical_coordinate_drift={k:float(np.max(np.abs(getattr(final,k)-getattr(initial,k))))
                                    for k in ('u','p','coordinate','momentum')})


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--mode', choices=('initialize','reduce','evolve','refine'), required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--oracle', type=Path, required=True)
    p.add_argument('--checkpoint', type=Path)
    p.add_argument('--case', required=True)
    p.add_argument('--initialization', type=Path, required=True)
    p.add_argument('--initialization-sha256', required=True)
    p.add_argument('--column-oracle', type=Path, required=True)
    p.add_argument('--column-dimension', type=int, required=True)
    p.add_argument('--column-workspace-bytes', type=int, required=True)
    p.add_argument('--column-error-target', type=float, required=True)
    p.add_argument('--reuse-initial-contractions', type=Path)
    p.add_argument('--reuse-initial-sha256')
    p.add_argument('--checkpoint-sha256')
    p.add_argument('--cpus', type=int, nargs='+', choices=range(6,11), required=True)
    p.add_argument('--max-actions', type=int, required=True)
    p.add_argument('--max-seconds', type=float, required=True)
    p.add_argument('--max-output-bytes', type=int, required=True)
    p.add_argument('--max-memory-bytes', type=int, required=True)
    p.add_argument('--worker-memory-bytes', type=int, required=True)
    p.add_argument('--interval', type=float, default=.001)
    p.add_argument('--steps', type=int, default=1)
    p.add_argument('--max-steps', type=int, default=2)
    p.add_argument('--max-iterations', type=int, default=30)
    p.add_argument('--tolerance', type=float, default=1e-11)
    p.add_argument('--reality-tolerance', type=float, default=1e-10)
    p.add_argument('--block-size', type=int, default=8)
    args = p.parse_args()
    if len(set(args.cpus)) != len(args.cpus):
        p.error('CPUs must be unique explicit allocations')
    os.sched_setaffinity(0, {args.cpus[0]})
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    lock = (args.output/'executor.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    existing = sum(path.stat().st_size for path in args.output.rglob('*') if path.is_file())
    budget = PartitionBudget(max_actions=args.max_actions, max_seconds=args.max_seconds,
        max_bytes=args.max_output_bytes, max_memory_bytes=args.max_memory_bytes,
        worker_memory_bytes=args.worker_memory_bytes, existing_bytes=existing)
    budget.reserve(65536)
    result = dict(schema='P254-joint-partitioned-material-v1', status='partial', mode=args.mode,
        cpus=args.cpus, blas_threads=1, requested_interval=args.interval if args.mode in ('evolve','refine') else 0.,
        scope='all-slot homogeneous physical mean; not whole local-q/quantum completion')
    try:
        bindings = implementation_bindings()
        oracle = json.loads(args.oracle.read_text())
        if oracle.get('schema') != 'P254-joint-partitioned-oracle-v1' or oracle.get('status') != 'passed':
            raise ValueError('executed partitioned saved-output oracle required')
        # The unchanged numerical midpoint/partitioning proof is not a proof of
        # a former CLI or scalar constructor. Current case I/O is exercised below.
        for binding in oracle['source_bindings']:
            if Path(binding['path']).parent == ROOT/'src/substrate_framework':
                authenticate(binding)
        column_oracle = json.loads(args.column_oracle.read_text())
        if (column_oracle.get('schema') != 'P254-joint-column-controller-oracle-v1'
                or column_oracle.get('status') != 'passed'
                or column_oracle['source_bindings'] != column_sources()):
            raise ValueError('current changed column-controller/field-bound proof required')
        initialization_binding = file_binding(args.initialization)
        if initialization_binding['sha256'] != args.initialization_sha256:
            raise ValueError('explicit actual case initialization differs')
        initialization = json.loads(args.initialization.read_text())
        context = load_case(args.case, reality_tolerance=args.reality_tolerance)
        saved = initialization.get('saved', {})
        if (initialization.get('status') != 'complete' or initialization.get('mode') != 'initialize'
                or initialization.get('case') != args.case or saved.get('model_identity') != context.identity):
            raise ValueError('complete initialization of this exact current case required')
        initial = DiskJointMean(authenticate(saved['joint_boundary']), expected_identity=context.identity)
        if initial.time != 0 or initial.origin != saved['original_boundary']:
            raise ValueError('original case state/preparation differs')
        state = initial
        if args.checkpoint:
            binding = file_binding(args.checkpoint)
            if binding['sha256'] != args.checkpoint_sha256:
                raise ValueError('checkpoint requires exact explicit digest')
            state = DiskJointMean(args.checkpoint, expected_identity=context.identity)
            if (state.origin != initial.origin or (state.time > 0 and
                    state.history[0]['prior_boundary_binding'] != file_binding(initial.path))):
                raise ValueError('checkpoint does not continue this actual case initialization')
        limits = ColumnLimits(args.column_dimension, args.column_workspace_bytes, args.column_error_target)
        source = sha256((source_identity(bindings)+json.dumps(limits.__dict__, sort_keys=True)).encode()).hexdigest()
        budget.check_processes(os.getpid(), [])
        persistent_json(args.output/'model.json', dict(identity=context.identity,
            specification=context.specification, implementation=bindings, source_identity=source,
            initialization=initialization_binding, column_limits=limits.__dict__,
            block_size=args.block_size, reality_tolerance=args.reality_tolerance), budget)
        result.update(case=args.case, oracle=file_binding(args.oracle), column_oracle=file_binding(args.column_oracle),
            initialization=initialization_binding, initial_boundary=file_binding(state.path),
            coefficient_defect_upper_bound=None, positive_time_symmetry=False)
        frame_sets = 3 if args.mode == 'refine' else args.steps if args.mode == 'evolve' else 0
        nr = context.reference.nbar.shape[0]
        contraction_payload = len(state.slots)*16*(nr+len(state.u)+5)
        result['resource_plan'] = dict(workers=len(args.cpus), simultaneous_slots=len(args.cpus),
            original_frame_bytes=sum(int(s.frame['bytes']) for s in state.slots),
            max_single_frame_bytes=max(int(s.frame['bytes']) for s in state.slots),
            retained_new_frame_sets=frame_sets,
            complex_contraction_payload_bytes_per_complete_barrier=contraction_payload,
            maximum_reduction_barriers_before_cache_reuse=1+frame_sets*(args.max_iterations+2),
            storage_estimate_note='payload excludes NPY/JSON/stage arrays and failed-attempt storage; hard reservation remains authoritative',
            midpoint_storage='per-slot n/G/ell/H0-energy receipts, not trial full frames',
            output_accounting='all incomplete and complete files retained; hard cumulative invocation reservation',
            memory='mmap original frames/D; one operator/slot per worker; aggregate RSS watchdog')
        if args.mode != 'initialize':
            serial = context.make_mean(budget, block_size=args.block_size)
            model = ColumnPartitionedMean.from_serial(serial, directory=args.output/'work',
                cpus=args.cpus, source_identity=source, limits=limits)
            if args.reuse_initial_contractions:
                reused = file_binding(args.reuse_initial_contractions)
                if reused['sha256'] != args.reuse_initial_sha256:
                    raise ValueError('explicit initial contraction reuse proof differs')
                result['initial_contraction_reuse'] = import_initial_contractions(
                    model, state, reused, block_size=args.block_size)
            initial = model.reduce(state)
            initial_record = save_stage(initial, args.output/('initial-reduction-'+file_binding(state.path)['sha256']), budget,
                metadata=dict(kind='full-initial-reduction', boundary=str(state.path), time_ry_inverse=state.time))
            result.update(initial_reduction=initial_record, initial_energy_total_ry=initial.energy,
                          initial_number_per_cell=initial.number)
            if args.mode in ('evolve','refine'):
                branches = (1,2) if args.mode == 'refine' else (args.steps,)
                results = []
                result['branches'] = []
                for steps in branches:
                    final, history = model.evolve(state, args.output/('interval-%s-steps-%d' % (args.interval,steps)),
                        physical_interval=args.interval, steps=steps, max_steps=args.max_steps,
                        max_iterations=args.max_iterations, tolerance=args.tolerance, refinement_label='steps-'+str(steps))
                    final_stage = model.reduce(final)
                    result['branches'].append(dict(steps=steps, boundary=file_binding(final.path),
                        diagnostics=diagnostics(state,final,initial,final_stage,budget)))
                    results.append(final)
                if len(results) == 2:
                    a,b = results
                    error = 0.
                    for s,t in zip(a.slots,b.slots):
                        budget.check()
                        error = max(error,float(np.linalg.norm(a.arrays(s)[0]-b.arrays(t)[0])))
                    result['refinement'] = dict(max_frame_frobenius=error,
                        coordinate_max=max(float(np.max(np.abs(getattr(a,k)-getattr(b,k)))) for k in ('u','p','coordinate','momentum')),
                        accuracy_certified=False, interval_not_an_accuracy_license=True)
        result['status'] = 'complete'
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result['completed_boundaries'] = [str(path) for path in args.output.glob('**/boundary.json')]
    result['budget'] = budget.snapshot()
    encoded = (json.dumps(result,indent=2,allow_nan=False)+'\n').encode()
    if len(encoded) > 65536:
        raise RuntimeError('terminal report exceeds reserved storage')
    output = args.output/('result-'+uuid4().hex+'.json')
    with output.open('xb') as stream:
        stream.write(encoded)
    print(str(output), flush=True)
    if result['status'] != 'complete':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()

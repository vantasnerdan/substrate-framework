#!/usr/bin/env python3
"""Bounded resumable all-slot physical joint mean; no MPI, CPU allocation explicit.

Use --mode initialize/reduce/evolve/refine with the SAME --output to resume.
Each invocation has a new explicitly granted action/time allocation; all existing
output (including abandoned attempts) consumes --max-output-bytes. --cpus names
only currently allocated cores in 6..10, never an implicit request for all five.
Example after the source-bound oracle passes:
  python run_joint_partitioned_material.py --mode reduce --output /var/tmp/joint-partitioned \
    --oracle /var/tmp/joint-partitioned-oracle.json --cpus 6 7 \
    --max-actions 1000000 --max-seconds 7200 --max-output-bytes 4000000000 \
    --max-memory-bytes 28000000000 --worker-memory-bytes 7000000000
Then --mode refine --interval .001 --max-steps 2 with a new sufficient storage
allocation (three complete new frame sets alone require about 39.4 GB).
Completed per-slot contractions are inspectable below work/<request>/slot-*/;
only complete.json receipts are reusable. History is the serial joint interface.
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
import run_joint_streamed_material as settled
from substrate_framework.partitioned_joint_budget import PartitionBudget
from substrate_framework.partitioned_joint_mean import PartitionedJointMean, save_stage
from substrate_framework.partitioned_joint_slots import publish
from substrate_framework.streamed_joint_state import DiskJointMean
from substrate_framework.streamed_material_state import file_binding


def implementation_bindings():
    paths = [Path(__file__), HERE/'run_joint_partitioned_oracle.py', HERE/'joint-partitioned-registration.json']
    paths += [ROOT/'src/substrate_framework'/name for name in (
        'partitioned_joint_budget.py', 'partitioned_joint_slots.py', 'partitioned_joint_mean.py',
        'prepared_reference.py', 'physical_field_environment.py', 'qe_bloch_frame.py',
        'qe_real_bloch_grid.py', 'qe_bloch_dyads.py', 'qe_bloch_finiteq_ionic.py',
        'qe_material_grid.py', 'qe_material_hamiltonian.py', 'qe_material_field_derivatives.py',
        'qe_nc_projectors.py', 'material_mixed_field_ionic.py', 'material_frozen_action.py')]
    bindings = {b['path']: b for b in settled.implementation_bindings()}
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
        expected = {b['path']: b['sha256'] for b in oracle['source_bindings']}
        if any(expected.get(b['path']) != b['sha256'] for b in bindings):
            raise ValueError('partitioned oracle does not bind the current complete implementation')
        result['oracle'] = file_binding(args.oracle)
        restored, original, ions, joint, coefficients, field, spec, identity = settled.load_settled()
        budget.check_processes(os.getpid(), [])
        source = source_identity(bindings)
        model_record = dict(identity=identity, specification=spec, implementation=bindings,
                            source_identity=source, block_size=args.block_size,
                            reality_tolerance=args.reality_tolerance)
        persistent_json(args.output/'model.json', model_record, budget)
        if args.checkpoint:
            binding = file_binding(args.checkpoint)
            if binding['sha256'] != args.checkpoint_sha256:
                raise ValueError('checkpoint requires exact explicit digest')
            state = DiskJointMean(args.checkpoint, expected_identity=identity)
        elif list(args.output.glob('initial*/boundary.json')):
            state = DiskJointMean(sorted(args.output.glob('initial*/boundary.json'))[0], expected_identity=identity)
        else:
            # Incomplete wrapper directories contain no electronic arrays. Keep
            # them as evidence and create a new wrapper generation on retry.
            initial_dir = args.output/'initial'
            if initial_dir.exists():
                initial_dir = args.output/('initial-'+uuid4().hex)
            state = DiskJointMean.wrap_original(original, initial_dir, identity=identity,
                coordinate=joint['field_preparation']['coordinate_mean'],
                momentum=joint['field_preparation']['momentum_mean'],
                expected_slots=tuple(settled.slot_identity(s) for s in restored[6]),
                evidence=dict(new_scalar=spec, field_preparation=spec['joint']), budget=budget)
        if state.time == 0 and (not np.array_equal(state.u, ions['u_bohr']) or not np.array_equal(state.p, ions['p_canonical'])):
            raise ValueError('initial ionic preparation changed')
        result['initial_boundary'] = file_binding(state.path)
        frame_sets = 3 if args.mode == 'refine' else args.steps if args.mode == 'evolve' else 0
        nr = restored[0]['saved_density_bohr3'].shape[0]
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
            serial = settled.build_model(restored, ions, joint, coefficients, field, identity, args, budget)
            model = PartitionedJointMean.from_serial(serial, directory=args.output/'work',
                cpus=args.cpus, source_identity=source)
            initial = model.reduce(state)
            initial_record = save_stage(initial, args.output/('initial-reduction-'+file_binding(state.path)['sha256']), budget,
                metadata=dict(kind='full-initial-reduction', boundary=str(state.path), time_ry_inverse=state.time))
            result.update(initial_reduction=initial_record, initial_energy_total_ry=initial.energy,
                          initial_number_per_cell=initial.number)
            if state.time == 0:
                result['initial_force_aggregate_difference_max'] = float(np.max(np.abs(initial.force-np.asarray(joint['initial_join']['force_per_cell_ry_bohr']).reshape(-1))))
                result['initial_field_rhs_aggregate_difference_max'] = float(np.max(np.abs(initial.field_force-joint['initial_join']['field_momentum_rhs'])))
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

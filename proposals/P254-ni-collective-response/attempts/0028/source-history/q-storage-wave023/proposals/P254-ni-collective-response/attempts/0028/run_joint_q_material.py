#!/usr/bin/env python3
"""Selected actual moving-baseline q response; never a whole-local reconstruction.

Supply a completed physical joint mean boundary AND its explicit digest, an
executed source-bound partitioned q oracle, --cpus 6 7 8 9 10 and explicit
action/wall/disk/aggregate-RSS/per-worker-RSS caps. --mode all computes forward
and BOTH genuine initial-source pullbacks. --resume keeps the same scientific
request/output root and authenticates completed edge/frame work; caps and CPU
allocation may change. Mutable scratch is never a restart input. Invocation
receipts are result-<id>.json; checkpoints remain immutable under invocations/.
Disk retains per-time baseline frames and per-stage edge work, including failed
attempts, under --max-output-bytes; plan mode reports payload estimates. No cache
eviction or all-q resident allocation. Global K/K-adjoint runs only after the
complete edge barrier, in a watchdog-supervised worker; local feedback builds
only the finite-q bridge. First targets: q=0, +.25 and -.25 separately, H37x.
No material accuracy license is invented: every result requires refinement.
"""
import argparse
from functools import partial
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
from uuid import uuid4

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from run_joint_streamed_material import load_settled, slot_identity
from run_clamped_array_material import authenticated_arrays
from material_run_contract import bound_json
from substrate_framework.material_bloch_admission import restore_accepted
from substrate_framework.streamed_material_state import file_binding
from substrate_framework.partitioned_joint_budget import PartitionBudget
from substrate_framework.partitioned_joint_slots import publish
from substrate_framework.material_joint_response import MaterialJointHistory
from substrate_framework.joint_kernel_adjoint import physical_kernel_adjoint
from substrate_framework.streamed_joint_response import MappedResponse
from substrate_framework.partitioned_q_response import PartitionedQResponse, partition_plan
from substrate_framework.streamed_joint_covariance import StreamedTwoSidedCovariance
from substrate_framework.joint_ordered_sources import JointProductPreparation
from substrate_framework.joint_response_observables import total_force_probe


MODULES = ('partitioned_q_response.py', 'partitioned_q_work.py', 'partitioned_joint_budget.py',
           'partitioned_joint_slots.py', 'streamed_joint_response.py', 'streamed_joint_covariance.py', 'material_joint_response.py',
           'joint_kernel_adjoint.py', 'joint_bloch_response.py', 'joint_ordered_sources.py',
           'joint_response_observables.py', 'material_physical_operators.py',
           'qe_bloch_finiteq.py', 'qe_bloch_finiteq_ionic.py', 'qe_real_bloch_grid.py',
           'spin_pbe_kernel.py', 'streamed_joint_state.py', 'streamed_joint_reduction.py',
           'streamed_joint_mean.py', 'qe_material_grid.py', 'qe_material_reference.py',
           'qe_bloch_edge.py', 'qe_bloch_frame.py', 'qe_bloch_dyads.py',
           'qe_material_hamiltonian.py', 'material_frozen_action.py',
           'streamed_material_evolution.py', 'streamed_material_state.py',
           'material_joint_field_scalar.py')


def bindings():
    return [file_binding(p) for p in [HERE/'run_joint_q_material.py', HERE/'run_joint_q_oracle.py',
        HERE/'run_joint_partitioned_q_oracle.py', HERE/'joint-partitioned-q-registration.json',
        HERE/'joint-q-registration.json', HERE/'run_joint_response_oracle.py',
        HERE/'run_joint_streamed_material.py']+[ROOT/'src/substrate_framework'/name for name in MODULES]]


def encode(value):
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    if isinstance(value, np.ndarray):
        return encode(value.tolist())
    if isinstance(value, (complex, np.complexfloating)):
        return dict(real=float(value.real), imag=float(value.imag))
    if isinstance(value, np.generic):
        return value.item()
    return value


def selected_readout(executor, state, probe, time):
    if hasattr(executor.edge_stage, 'clear_interval'):
        executor.edge_stage.clear_interval()
    if hasattr(executor, 'prepare_time'):
        executor.prepare_time(time)
    value = sum(np.vdot(getattr(probe, key), getattr(state, key)) for key in ('u','p','coordinate','momentum'))
    for i, edge in enumerate(executor.r.edges):
        executor.budget.check()
        ps = executor.edge_stage.frame(edge.source, time)
        pt = executor.edge_stage.frame(edge.target, time)
        value += executor.r.sectors[edge.source].weight*(np.vdot(probe.actions[i](ps), state.ket[i])
            +np.einsum('ij,ji->', state.bra[i], probe.adjoints[i](pt)))
    return complex(value)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--history', type=Path, required=True)
    p.add_argument('--history-sha256', required=True)
    p.add_argument('--oracle', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--mode', choices=('plan','forward','pair','all'), required=True)
    p.add_argument('--q', nargs=3, type=float, required=True)
    p.add_argument('--cpus', nargs='+', type=int, choices=range(6,11), required=True)
    p.add_argument('--max-memory-bytes', type=int, required=True)
    p.add_argument('--worker-memory-bytes', type=int, required=True)
    p.add_argument('--resume', action='store_true', help='Reuse authenticated immutable edge/frame work in --output')
    p.add_argument('--max-actions', type=int, required=True)
    p.add_argument('--max-seconds', type=float, required=True)
    p.add_argument('--max-output-bytes', type=int, required=True)
    p.add_argument('--time', type=float, required=True)
    p.add_argument('--second-time', type=float, required=True)
    p.add_argument('--subdivisions', type=int, required=True)
    p.add_argument('--atom', type=int, default=37)
    p.add_argument('--axis', choices=('x','y','z'), default='x')
    p.add_argument('--cell', nargs=3, type=int, default=(0,0,0))
    p.add_argument('--amplitude', type=float, default=1.)
    p.add_argument('--dynamics', choices=('screened-joint','bare-electronic'), default='screened-joint')
    p.add_argument('--block-size', type=int, default=8)
    p.add_argument('--checkpoint', type=Path)
    p.add_argument('--checkpoint-sha256')
    for side in ('left', 'right'):
        p.add_argument('--'+side+'-checkpoint', type=Path)
        p.add_argument('--'+side+'-checkpoint-sha256')
    args = p.parse_args()
    if len(set(args.cpus)) != len(args.cpus):
        raise ValueError('allocated CPUs must be unique')
    os.sched_setaffinity(0, {args.cpus[0]})
    args.output.mkdir(parents=True, exist_ok=args.resume)
    existing = sum(p.stat().st_size for p in args.output.rglob('*') if p.is_file())
    budget = PartitionBudget(max_actions=args.max_actions, max_seconds=args.max_seconds,
        max_bytes=args.max_output_bytes, max_memory_bytes=args.max_memory_bytes,
        worker_memory_bytes=args.worker_memory_bytes, existing_bytes=existing)
    invocation = uuid4().hex
    run_directory = args.output/'invocations'/invocation
    run_directory.mkdir(parents=True)
    budget.reserve(131072)
    result = dict(schema='P254-joint-q-material-v2', status='partial', q=args.q, cpus=args.cpus,
        blas_threads=1, dynamics=args.dynamics, whole_spatial_reconstruction=False,
        required_q_for_whole_local=64, numerical_accuracy_certified=False)
    try:
        if args.checkpoint and args.mode != 'forward':
            raise ValueError('forward restart uses --mode forward; covariance always pulls both probes to original t0')
        source_bindings = bindings()
        oracle = bound_json(file_binding(args.oracle))
        expected = {x['path']: x['sha256'] for x in oracle['source_bindings']}
        if oracle['schema'] != 'P254-joint-partitioned-q-oracle-v1' or oracle['status'] != 'passed':
            raise ValueError('executed source-bound pinned q forward/dual/Fock/recovery oracle required')
        if any(expected.get(x['path']) != x['sha256'] for x in source_bindings):
            raise ValueError('oracle does not bind current q implementation')
        hb = file_binding(args.history)
        if hb['sha256'] != args.history_sha256:
            raise ValueError('explicit completed mean history digest differs')
        restored, original_binding, ions, joint, coeff, field, spec, model_identity = load_settled()
        arrays, receipt, manifest, saved, original, grid, ordered, pairs = restored
        background, _ = restore_accepted(original)
        if not np.array_equal(arrays['saved_density_bohr3'], background.rho_updown):
            raise ValueError('physical H0/K saved common reference differs')
        history = MaterialJointHistory(hb, model_identity=model_identity, background=background,
            grid=grid, ordered=ordered, slot_identity=slot_identity, budget=budget, block_size=args.block_size)
        if (not np.array_equal(history.initial.u, ions['u_bohr'])
                or not np.array_equal(history.initial.p, ions['p_canonical'])
                or not np.array_equal(history.initial.coordinate, joint['field_preparation']['coordinate_mean'])
                or not np.array_equal(history.initial.momentum, joint['field_preparation']['momentum_mean'])):
            raise ValueError('initial canonical means differ from original product preparation')
        q = np.asarray(args.q)
        mesh = np.array([grid.reference.monkhorst_pack['nk%d' % i] for i in (1,2,3)])
        qi = tuple(np.rint(q*mesh).astype(int) % mesh)
        ndof = 3*grid.reference.nat
        r = history.bind_response(q=q, coefficients=coeff,
            common_curvature=arrays['common_q_grid_ry_bohr2'][qi].reshape(ndof,ndof),
            masses=joint['original_atoms']['masses_ry'], charges=joint['original_atoms']['valence_charges_e'],
            polarizations=field.polarizations, field_mass=field.mass, ncells=64,
            density_shape=arrays['saved_density_bohr3'].shape, scalar_identity=model_identity,
            coverage_identity=saved.coverage_identity)
        times, second_times = history.time_grid(args.time,args.subdivisions), history.time_grid(args.second_time,args.subdivisions)
        if not 1 <= args.atom <= grid.reference.nat or not np.isfinite(args.amplitude):
            raise ValueError('original atom and finite coherent amplitude required')
        if args.atom == 37 and grid.reference.atom_names[36] != 'H':
            raise ValueError('original H37 authority differs')
        direction = np.zeros(ndof)
        direction[3*(args.atom-1)+'xyz'.index(args.axis)] = 1
        probe = total_force_probe(r, direction, identity='original-%s%d%s-total-force' % (grid.reference.atom_names[args.atom-1],args.atom,args.axis))
        result.update(source_bindings=source_bindings, baseline=hb, model_identity=model_identity,
            resource_plan=partition_plan(r,times,second_times,len(args.cpus)), time_grid=times, second_time_grid=second_times,
            observable=probe.identity, all_explicit_edges=len(r.edges), fixed_original_D=True,
            mean_history_method='authenticated actual midpoint stages; full-PW exponential frame reconstruction')
        if args.mode != 'plan':
            identity = sha256(json.dumps(dict(model=model_identity,q=args.q,dynamics=args.dynamics,
                sources=source_bindings),sort_keys=True).encode()).hexdigest()
            request = dict(identity=identity, baseline=hb, q=args.q, mode=args.mode,
                times=times.tolist(), second_times=second_times.tolist(), atom=args.atom, axis=args.axis,
                amplitude=args.amplitude, cell=list(args.cell), subdivisions=args.subdivisions,
                block_size=args.block_size, checkpoint=str(args.checkpoint), checkpoint_sha256=args.checkpoint_sha256,
                left_checkpoint=str(args.left_checkpoint), left_checkpoint_sha256=args.left_checkpoint_sha256,
                right_checkpoint=str(args.right_checkpoint), right_checkpoint_sha256=args.right_checkpoint_sha256)
            request_path = args.output/'request.json'
            if request_path.exists():
                if json.loads(request_path.read_text()) != request:
                    raise ValueError('resume requires identical scientific request and source bytes')
            else:
                publish(request_path,request,budget)
            measure = abs(np.linalg.det(grid.reference.cell_bohr))/int(np.prod(grid.reference.fft_smooth))
            executor = PartitionedQResponse(r, history, kernel_adjoint=partial(physical_kernel_adjoint,background),
                grid_measure=measure, budget=budget, directory=run_directory/'executor', identity=identity,
                dynamics=args.dynamics, cpus=args.cpus, work_directory=args.output/'scientific-work')
            executor.provenance = dict(baseline=hb, history=history.records)
            if args.mode in ('forward','all'):
                if args.checkpoint:
                    cb = file_binding(args.checkpoint)
                    if cb['sha256'] != args.checkpoint_sha256:
                        raise ValueError('explicit q checkpoint digest differs')
                    initial, record = MappedResponse.restore(r,cb,budget,identity)
                    previous = record['provenance']['history']
                    if history.records[:len(previous)] != previous or record['kind'] != 'forward':
                        raise ValueError('checkpoint is not a forward state on this exact extended history')
                    start = record['time']
                    times = np.r_[start, times[times > start]]
                    if len(times) < 2:
                        raise ValueError('extension must request a later covered time')
                else:
                    initial = MappedResponse(r,run_directory/'coherent-initial',budget,create=True)
                    initial.p[:] = args.amplitude*np.exp(-2j*np.pi*q @ np.asarray(args.cell))*direction/r.sqrt_n
                    initial.publish(time=0.,identity=identity,budget=budget,kind='coherent-initial',provenance=executor.provenance)
                final, checkpoints = executor.propagate(initial,times,run_directory/'forward')
                result['coherent'] = dict(readout=selected_readout(executor,final,probe,args.time),
                    checkpoints=checkpoints, amplitude=args.amplitude, cell=args.cell,
                    convention='local canonical momentum kick q component; 1/sqrt(N), not whole local response')
            if args.mode in ('pair','all'):
                ionic = bound_json(spec['ionic'])
                ia = authenticated_arrays(ionic['receipt'])
                prep = JointProductPreparation(r,history.initial_stage(),ionic_ordered=ia['local_ordered'],
                    ionic_minus_ordered=ia['local_ordered'],
                    field_ordered=joint['field_preparation']['ordered_covariance'] if r.gamma else np.zeros((0,0)),
                    electronic_identity=original_binding['sha256'],ionic_identity=spec['ionic']['sha256'],
                    field_identity=spec['joint']['sha256'],product_identity='original-product-'+model_identity)
                resumes = {}
                for side in ('left', 'right'):
                    path = getattr(args, side+'_checkpoint')
                    if path:
                        binding = file_binding(path)
                        if binding['sha256'] != getattr(args, side+'_checkpoint_sha256'):
                            raise ValueError('explicit '+side+' dual checkpoint digest differs')
                        resumes[side+'_checkpoint'] = binding
                result['quantum'] = StreamedTwoSidedCovariance(prep,executor).pair(
                    probe,probe,times,second_times,run_directory/'quantum',**resumes)
        result['status'] = 'complete'
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result['budget'] = budget.snapshot()
    encoded = (json.dumps(encode(result),indent=2,allow_nan=False)+'\n').encode()
    if len(encoded) > 131072:
        raise RuntimeError('reserved terminal receipt exceeded')
    with (args.output/('result-'+invocation+'.json')).open('xb') as stream:
        stream.write(encoded)
    if result['status'] != 'complete':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()

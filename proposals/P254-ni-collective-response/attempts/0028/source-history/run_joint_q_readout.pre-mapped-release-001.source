#!/usr/bin/env python3
"""Read one authenticated committed material endpoint; never propagate a trajectory.

All scientific selections and resource caps are explicit. The historical producer
identity is pinned independently of current consumer sources. Only owned scalar
work/metadata are charged; borrowed checkpoint, mean, and preparation files stay
read-only. --resume retains the original deadline, admitted action costs, and
actual completed QWork worker receipts. No synthetic case/model/frame fallback.
"""
import argparse
from functools import partial
import json
import os
from pathlib import Path
import resource
import sys
from time import monotonic, process_time, time
from uuid import uuid4

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE), str(HERE.parent/'0025'), str(HERE.parent/'0024')]
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from matched_joint_case import load_case
from matched_joint_state import preparation_arrays
from run_joint_q_material import (
    bindings, admit_oracle, bind_saved_model, selected_initialization, encode,
    producer_identity)
from mixed_material_driver_support import ContinuationBudget, owned_directory
from partitioned_material_readout import partitioned_selected_readout
from substrate_framework.material_control_reference import bound_json, bound_arrays
from substrate_framework.streamed_material_state import authenticate, file_binding
from substrate_framework.streamed_joint_response import MappedResponse
from substrate_framework.partitioned_q_response import PartitionedQResponse
from substrate_framework.partitioned_q_work import QWork
from substrate_framework.partitioned_q_storage import durable_publish, sync_directory
from substrate_framework.joint_kernel_adjoint import physical_kernel_adjoint
from substrate_framework.joint_response_observables import total_force_probe
from substrate_framework.joint_field_probes import OBSERVABLES, global_field_probe, validate_selection


SCHEMA = 'P254-committed-joint-q-readout-v1'
TERMINAL_BYTES = 2*1024**2
READOUT_SOURCES = ('partitioned_material_readout.py', 'run_joint_q_readout.py',
                   'mixed_material_driver_support.py')


def _bound(path, expected, label):
    binding = file_binding(Path(path).resolve())
    if binding['sha256'] != expected:
        raise ValueError('explicit '+label+' SHA256 differs: '+binding['path'])
    return binding, bound_json(binding)


def _owned_bytes(directory):
    total = 0
    for path in directory.rglob('*'):
        if path.is_symlink():
            raise ValueError('owned readout output cannot contain symlink aliases: '+str(path))
        if path.is_file():
            if path.stat().st_nlink != 1:
                raise ValueError('owned readout output cannot contain borrowed hardlink aliases: '+str(path))
            total += path.stat().st_size
    return total


def _admit_readout(oracle, sources):
    if (oracle.get('schema') != 'P254-joint-q-readout-oracle-v1'
            or oracle.get('status') != 'passed' or not oracle.get('checks')
            or any(check.get('passed') is not True for check in oracle['checks'])):
        raise ValueError('executed independent full-density complex readout/recovery oracle required')
    expected = {item['path']: item['sha256'] for item in oracle['source_bindings']}
    required = [source for source in sources if Path(source['path']).name in
        READOUT_SOURCES+('partitioned_q_work.py', 'partitioned_q_storage.py',
                         'partitioned_q_response.py', 'streamed_joint_response.py')]
    if any(expected.get(source['path']) != source['sha256'] for source in required):
        raise ValueError('readout oracle does not bind the current consumer/scalar/recovery implementation')


def _caps(args):
    if len(set(args.cpus)) != len(args.cpus):
        raise ValueError('allocated CPUs6-10 must be unique')
    limits = ((args.max_actions, 1024, 'Hamiltonian actions'),
              (args.max_seconds, 7200., 'seconds'),
              (args.max_memory_bytes, 16*1024**3, 'aggregate Q RSS bytes'),
              (args.worker_memory_bytes, 4*1024**3, 'worker RSS bytes'),
              (args.max_output_bytes, 256*1024**2, 'owned output bytes'),
              (args.max_open_files, 8192, 'open files'))
    for value, ceiling, label in limits:
        if not np.isfinite(value) or not 0 < value <= ceiling:
            raise ValueError('explicit '+label+' cap must be positive and <= '+str(ceiling))
    if args.max_output_bytes < TERMINAL_BYTES:
        raise ValueError('owned output cap must also reserve the complete terminal provenance report')
    if args.worker_memory_bytes > args.max_memory_bytes:
        raise ValueError('worker memory cap cannot exceed aggregate Q memory cap')
    if not set(args.cpus).issubset(os.sched_getaffinity(0)):
        raise ValueError('requested CPUs are outside the current allocated affinity')
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    if hard != resource.RLIM_INFINITY and args.max_open_files > hard:
        raise ValueError('requested file cap exceeds the inherited hard file limit')
    resource.setrlimit(resource.RLIMIT_NOFILE, (args.max_open_files, hard))
    os.sched_setaffinity(0, {args.cpus[0]})
    return dict(cpus=args.cpus, blas_threads=1, max_actions=args.max_actions,
        max_seconds=args.max_seconds, aggregate_q_rss_bytes=args.max_memory_bytes,
        per_worker_rss_bytes=args.worker_memory_bytes, owned_output_bytes=args.max_output_bytes,
        open_files_soft=resource.getrlimit(resource.RLIMIT_NOFILE)[0],
        open_files_hard=resource.getrlimit(resource.RLIMIT_NOFILE)[1],
        manager_affinity=sorted(os.sched_getaffinity(0)),
        manager_container_limit='external launcher authority; this CLI enforces aggregate parent+workers Q RSS')


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', required=True)
    for name in ('initialization', 'history', 'oracle', 'readout-oracle', 'checkpoint', 'producer-request'):
        parser.add_argument('--'+name, type=Path, required=True)
        parser.add_argument('--'+name+'-sha256', required=True)
    parser.add_argument('--producer-identity', required=True,
                        help='PINNED historical checkpoint producer, not a current source hash')
    parser.add_argument('--q', nargs=3, type=float, required=True)
    parser.add_argument('--time', type=float, required=True)
    parser.add_argument('--atom', type=int, required=True, help='original one-based atom order')
    parser.add_argument('--axis', choices=('x', 'y', 'z'), required=True)
    parser.add_argument('--observable', choices=OBSERVABLES, required=True)
    parser.add_argument('--field-axis', choices=('x', 'y', 'z'))
    parser.add_argument('--dynamics', choices=('screened-joint', 'bare-electronic'), required=True)
    parser.add_argument('--block-size', type=int, required=True)
    parser.add_argument('--cpus', nargs='+', type=int, choices=range(6, 11), required=True)
    parser.add_argument('--max-memory-bytes', type=int, required=True)
    parser.add_argument('--worker-memory-bytes', type=int, required=True)
    parser.add_argument('--max-actions', type=int, required=True)
    parser.add_argument('--max-seconds', type=float, required=True)
    parser.add_argument('--max-output-bytes', type=int, required=True)
    parser.add_argument('--max-open-files', type=int, required=True)
    parser.add_argument('--output', type=Path, required=True, help='exclusively owned readout directory')
    parser.add_argument('--report', type=Path,
                        help='new report path inside --output; default result-INVOCATION.json')
    parser.add_argument('--resume', action='store_true')
    return parser


def main(arguments=None):
    args = _parser().parse_args(arguments)
    validate_selection(args.observable, args.field_axis, args.q)
    if (not np.isfinite(args.q).all() or not np.isfinite(args.time) or args.time < 0
            or args.block_size < 1 or len(args.producer_identity) != 64
            or any(c not in '0123456789abcdef' for c in args.producer_identity)):
        raise ValueError('finite explicit q/time, positive block size and pinned SHA256 producer identity required')
    caps = _caps(args)
    args.output = args.output.absolute()
    if args.output.is_symlink():
        raise ValueError('owned output cannot be a symlink')
    args.output = args.output.resolve()
    invocation = uuid4().hex
    report = args.report.absolute().resolve() if args.report else args.output/('result-'+invocation+'.json')
    if not report.is_relative_to(args.output) or report == args.output or report.exists():
        raise ValueError('report must be a new exclusively owned path inside --output')
    protected = [getattr(args, key) for key in
                 ('initialization', 'history', 'oracle', 'readout_oracle', 'checkpoint', 'producer_request')]
    if (any(Path(path).resolve().is_relative_to(args.output) for path in protected)
            or args.output.is_relative_to(args.checkpoint.resolve().parent)):
        raise ValueError('owned output must be disjoint from borrowed scientific inputs/checkpoint directory')
    started, cpu_started = monotonic(), process_time()
    with owned_directory(args.output, args.resume):
        owner_path = args.output/'owner.json'
        prior_owner = json.loads(owner_path.read_text()) if args.resume else None
        started_at = prior_owner['started_unix'] if prior_owner is not None else time()
        budget = ContinuationBudget(args.output, started_at=started_at,
            max_actions=args.max_actions, max_seconds=args.max_seconds,
            max_bytes=args.max_output_bytes, max_memory_bytes=args.max_memory_bytes,
            worker_memory_bytes=args.worker_memory_bytes, existing_bytes=_owned_bytes(args.output))
        budget.reserve(TERMINAL_BYTES)
        result = dict(schema=SCHEMA, status='failed', complete=False, invocation=invocation,
            execution_caps=caps, historical_producer_identity=args.producer_identity,
            requested_inputs=dict(case=args.case, q=args.q, time=args.time, atom=args.atom,
                axis=args.axis, observable=args.observable, field_axis=args.field_axis,
                dynamics=args.dynamics, block_size=args.block_size,
                files={name: dict(path=str(getattr(args, name).resolve()),
                                 sha256=getattr(args, name+'_sha256'))
                       for name in ('initialization', 'history', 'oracle', 'readout_oracle',
                                    'checkpoint', 'producer_request')}),
            scope=dict(selected_q_only=True, whole_spatial_reconstruction=False,
                numerical_accuracy_certified=False, propagation_rerun=False,
                not_claimed=['full64q closure', 'noise/covariance', 'energy or energy transfer',
                             'ranking', 'finite-amplitude/size/refinement convergence']))
        work = None
        before = budget.snapshot()
        try:
            operator_sources = bindings()
            source_bindings = list({item['path']: item for item in operator_sources+
                [file_binding(HERE/name) for name in READOUT_SOURCES]}.values())
            for source in source_bindings:
                authenticate(source)
            owner = dict(request=result['requested_inputs'], sources=source_bindings,
                         producer_identity=args.producer_identity, started_unix=started_at,
                         execution_caps=caps)
            if prior_owner is not None:
                if prior_owner != owner:
                    raise ValueError('resume requires the exact original checkpoint/source/operator request')
            else:
                durable_publish(owner_path, owner, budget)
            result['source_bindings'] = source_bindings
            ob, oracle = _bound(args.oracle, args.oracle_sha256, 'physical operator oracle')
            rb, readout_oracle = _bound(args.readout_oracle, args.readout_oracle_sha256, 'readout oracle')
            admit_oracle(oracle, operator_sources, args.observable)
            _admit_readout(readout_oracle, source_bindings)
            hb, _ = _bound(args.history, args.history_sha256, 'completed mean history')
            cb, checkpoint = _bound(args.checkpoint, args.checkpoint_sha256, 'committed checkpoint')
            pb, producer_request = _bound(args.producer_request, args.producer_request_sha256,
                                           'historical producer request')
            initialization_binding, saved = selected_initialization(args)
            if (checkpoint.get('identity') != args.producer_identity
                    or checkpoint.get('kind') != 'forward' or checkpoint.get('complete') is not True
                    or checkpoint.get('time') != args.time or checkpoint.get('q') != args.q):
                raise ValueError('pinned checkpoint must be a complete physical forward endpoint at exact q/time')
            if (producer_request.get('identity') != args.producer_identity
                    or producer_request.get('case') != args.case
                    or producer_request.get('initialization') != initialization_binding
                    or producer_request.get('baseline') != hb or producer_request.get('q') != args.q
                    or producer_request.get('mode') not in ('forward', 'all')
                    or args.time not in producer_request.get('times', [])
                    or producer_request.get('block_size') != args.block_size):
                raise ValueError('explicit historical producer request differs from selected checkpoint/case/history')
            context = bind_saved_model(load_case(args.case), saved, oracle.get('source_compatibility'))
            budget.check_processes(os.getpid(), ())
            if (context.identity != saved['model_identity']
                    or context.specification != saved['specification']):
                raise ValueError('saved initialization differs from the exact selected matched scalar/preparation')
            history = context.make_history(hb, budget, block_size=args.block_size)
            budget.check_processes(os.getpid(), ())
            if (file_binding(history.initial.path) != saved['joint_boundary']
                    or history.initial.origin != saved['original_boundary']):
                raise ValueError('completed history must start at the selected actual original initialization')
            if (not np.array_equal(history.initial.u, np.zeros_like(context.reference.masses))
                    or not np.array_equal(history.initial.p, np.zeros_like(context.reference.masses))
                    or not np.array_equal(history.initial.coordinate, context.field_preparation['coordinate_mean'])
                    or not np.array_equal(history.initial.momentum, context.field_preparation['momentum_mean'])):
                raise ValueError('initial canonical means differ from the original matched product preparation')
            prepared = bound_arrays(saved['preparation'])
            expected = preparation_arrays(context)
            if set(prepared) != set(expected) or any(not np.array_equal(prepared[k], v) for k, v in expected.items()):
                raise ValueError('saved matched ionic/field preparation differs')
            if checkpoint.get('provenance') != dict(baseline=hb, history=history.records):
                raise ValueError('committed checkpoint is not bound to this exact authenticated completed history')
            # Production readout borrows actual saved endpoint frames. It must
            # not spend this owned scalar budget reconstructing missing means.
            saved_frames = history._saved_frames(args.time)
            if saved_frames is None:
                raise ValueError('actual saved mean frames at the exact checkpoint time required; no propagation/reconstruction permitted')
            grid = context.grid
            if not 1 <= args.atom <= grid.reference.nat:
                raise ValueError('selected atom is outside the original one-based material atom order')
            r = context.bind_response(history, np.asarray(args.q))
            budget.check_processes(os.getpid(), ())
            direction = np.zeros(len(context.reference.masses))
            direction[3*(args.atom-1)+'xyz'.index(args.axis)] = 1.
            observable_metadata = None
            if args.observable == 'total-force':
                probe = total_force_probe(r, direction, identity='original-%s%d%s-total-force' % (
                    grid.reference.atom_names[args.atom-1], args.atom, args.axis))
            else:
                probe, observable_metadata = global_field_probe(r, context.field,
                    observable=args.observable, axis=args.field_axis)
            consumer_identity = producer_identity(context.identity, args.q, args.dynamics,
                                                   source_bindings, observable_metadata)
            operator_binding = dict(case=args.case, initialization=initialization_binding,
                model_identity=context.identity, specification=context.specification, history=hb,
                q=args.q, time=args.time, atom=args.atom, axis=args.axis,
                observable=args.observable, field_axis=args.field_axis, probe_identity=probe.identity,
                observable_metadata=observable_metadata, block_size=args.block_size,
                dynamics=args.dynamics, sources=source_bindings)
            request = dict(schema=SCHEMA, checkpoint=cb, historical_producer_identity=args.producer_identity,
                historical_producer_request=pb, physical_operator_oracle=ob, readout_oracle=rb,
                consumer_identity=consumer_identity, operator=operator_binding)
            run_directory = args.output/'invocations'/invocation
            run_directory.mkdir(parents=True)
            measure = abs(np.linalg.det(grid.reference.cell_bohr))/int(np.prod(grid.reference.fft_smooth))
            executor = PartitionedQResponse(r, history,
                kernel_adjoint=partial(physical_kernel_adjoint, context.background), grid_measure=measure,
                budget=budget, directory=run_directory/'executor', identity=consumer_identity,
                dynamics=args.dynamics, cpus=args.cpus, work_directory=args.output/'scientific-work')
            executor.provenance = dict(baseline=hb, history=history.records)
            executor.readout_source_bindings = source_bindings
            executor.readout_operator_binding = operator_binding
            state, state_record = MappedResponse.restore(r, cb, budget, args.producer_identity)
            work = QWork(args.output/'scalar-work', budget, args.cpus, consumer_identity)
            journals_before = {str(path) for path in work.directory.glob('*/execution-*.json')}
            scalar_before = budget.snapshot()
            value = partitioned_selected_readout(executor, state, probe, args.time, work, state_binding=cb)
            key, directory, count = work.last
            receipts = work.collect(directory, key, count)
            journals = sorted(directory.glob('execution-*.json'))
            current_journals = [path for path in journals if str(path) not in journals_before]
            current_items = set()
            for path in current_journals:
                execution = json.loads(path.read_text())
                if execution.get('status') == 'allocated':
                    current_items.update(execution['items'])
            result.update(status='complete', complete=True, request=request,
                source_bindings=source_bindings, physical_operator_oracle=ob, readout_oracle=rb,
                case=args.case, initialization=initialization_binding, baseline=hb,
                committed_checkpoint=cb, checkpoint_provenance=state_record['provenance'],
                historical_producer_request=dict(binding=pb, metadata=producer_request),
                consumer_identity=consumer_identity, model_identity=context.identity,
                scalar=context.specification['scalar'], matched_aggregate=context.aggregate_binding,
                observable=probe.identity, observable_selection=args.observable, units=probe.units,
                observable_metadata=observable_metadata, atom=args.atom,
                species=grid.reference.atom_names[args.atom-1], axis=args.axis,
                q=args.q, time_ry_inverse=args.time, dynamics=args.dynamics,
                readout=value, saved_mean_frame_boundary=file_binding(saved_frames.path),
                all_explicit_edges=len(r.edges), fixed_original_D=True,
                original_source_weights=[float(s.weight) for s in r.sectors],
                original_D=[slot.coefficient for slot in history.initial.slots],
                contraction='original canonical vdot sum once, then weighted independent ket/bra edge terms in original index order; no extra D/conjugation/real/cell factor',
                scalar_work=dict(key=key, directory=str(directory), count=count,
                    request=file_binding(directory/'request.json'), receipts=receipts,
                    execution_journals=[file_binding(path) for path in journals],
                    current_invocation_execution_journals=[file_binding(path) for path in current_journals],
                    current_invocation_allocated_indices=sorted(current_items),
                    receipt_bindings=[file_binding(directory/('item-%04d' % i)/'complete.json') for i in range(count)]),
                measurements=dict(budget_before=before, budget_before_scalar=scalar_before,
                    current_completed_worker_cpu_seconds=sum(x['cpu_seconds'] for x in receipts if x['index'] in current_items),
                    historical_propagation_costs='not re-executed, reset or included in this owned scalar budget; immutable original producer journals remain historical',
                    receipt_worker_cpu_seconds=sum(x['cpu_seconds'] for x in receipts),
                    receipt_worker_hamiltonian_actions=sum(x['hamiltonian_actions'] for x in receipts),
                    receipt_worker_peak_rss_bytes=max(x['process_peak_rss_bytes'] for x in receipts),
                    receipt_scope='actual immutable original worker costs; replay adopts completed items, never relabels their cost as current execution'),
                budget_scope='durable original-start wall deadline and pre-action counters; storage snapshots cover this invocation plus remeasured existing owned bytes; prior execution/allocation snapshots remain immutable journals, never attributed to a replay',
                storage_scope='exclusively owned scalar outputs/request/cost journals/reports; borrowed checkpoint/mean/preparation bytes not recharged or retired')
        except Exception as error:
            result['error'] = type(error).__name__+': '+str(error)
            result['recovery'] = (
                'immutable completed scalar QWork items retained; rerun identical --output with --resume and a new --report; original wall deadline and pre-action admissions remain charged'
                if owner_path.exists() else
                'input admission failed before a scientific owner was committed; no scalar work was launched; use a new output directory with corrected explicit inputs')
            if work is not None and work.last is not None:
                key, directory, count = work.last
                result['incomplete_work'] = dict(key=key, directory=str(directory), count=count)
                result['incomplete_work']['execution_journals'] = [
                    file_binding(path) for path in sorted(directory.glob('execution-*.json'))]
                try:
                    completed = [work.receipt(directory, key, i) for i in range(count)]
                    result['incomplete_work']['completed_receipts'] = [item for item in completed if item is not None]
                except Exception as receipt_error:
                    result['incomplete_work']['receipt_authentication_error'] = type(receipt_error).__name__+': '+str(receipt_error)
        result['budget'] = budget.snapshot()
        result['invocation_measurement'] = dict(elapsed_seconds=monotonic()-started,
            manager_cpu_seconds=process_time()-cpu_started,
            manager_peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
            owned_file_bytes_before_report=_owned_bytes(args.output))
        data = (json.dumps(encode(result), indent=2, allow_nan=False)+'\n').encode()
        if len(data) > TERMINAL_BYTES:
            raise RuntimeError('complete terminal provenance exceeds reserved report bytes')
        report.parent.mkdir(parents=True, exist_ok=True)
        with report.open('xb') as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        sync_directory(report.parent)
        print(json.dumps(dict(status=result['status'], report=file_binding(report)), allow_nan=False))
        if result['status'] != 'complete':
            raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()

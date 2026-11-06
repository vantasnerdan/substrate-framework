#!/usr/bin/env python3
"""Root-only pinned q oracle: original dense/Fock physics plus real disk recovery.

No material trajectory is inferred. Run with explicit allocated CPUs6-10; all
caps are invocation-wide. Existing q tolerances remain unchanged. Small positive
times reduce fork count, not the scientific constructions or two-time sources.
"""
import argparse
from contextlib import contextmanager
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from run_joint_response_oracle import (FiniteSystem, DenseTangent, FockOracle, check,
    reject, CHECKS, symbolic, pack_legs)
from run_joint_q_oracle import FiniteHistory, mapped, kernel_oracle
from run_joint_q_material import bindings, encode
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.partitioned_joint_slots import array_identity
from substrate_framework.partitioned_q_work import QWork, arrays
from substrate_framework.partitioned_q_response import PartitionedQResponse, PartitionedTwoSidedCovariance
from substrate_framework.streamed_joint_response import MappedResponse, StreamedJointResponse
from substrate_framework.streamed_material_state import file_binding
import substrate_framework.partitioned_q_response as q_response_module
from substrate_framework.streamed_joint_covariance import initial_contraction
from substrate_framework.joint_ordered_sources import JointProductPreparation


class DiskFiniteHistory(FiniteHistory):
    """Actual finite moving frames stored through production immutable workers."""
    def __init__(self, *args):
        super().__init__(*args)
        self.frame_work = None
        self.records = {}

    def stage_identity(self, time):
        frames, *_, hs = self.system.mean_values(self.system.mean.sol(time))
        return dict(time=float(time), frames=[array_identity(p) for p in frames],
            h=[array_identity(h) for h in hs], wrong_adjoint=self.wrong_adjoint,
            coupling=[dict(g=array_identity(e.coupling.g),v=array_identity(e.coupling.v))
                      for e in self.response.edges])

    def prepare_time(self, time):
        identity = dict(time=float(time), frames=self.stage_identity(time)['frames'])
        self.records[float(time)] = self.frame_work.execute('finite-actual-baseline', identity,
            len(self.response.sectors), lambda index: dict(frame=self.system.stage(time).frames[index]))

    def frame(self, index, time):
        record = self.records[float(time)][index]
        self.frame_work.verify(record['arrays']['frame'])
        return arrays(record)['frame']

    @contextmanager
    def __call__(self, index, time):
        with super().__call__(index, time) as data:
            e = self.response.edges[index]
            yield self.frame(e.source,time), self.frame(e.target,time), *data[2:]

    @contextmanager
    def feedback(self, index, time):
        e = self.response.edges[index]
        yield self.frame(e.source,time), self.frame(e.target,time), e.coupling.bridge


def exercise(system, fock, gamma, directory, budget, cpus):
    r = system.response(gamma)
    tag = 'Gamma' if gamma else 'nonGamma'
    dense = DenseTangent(system,r)
    history = DiskFiniteHistory(system,r,budget)
    kwargs = dict(kernel_adjoint=lambda n,**kw:system.k*n, grid_measure=1., budget=budget,
                  identity=tag+'-partitioned-original', cpus=cpus, work_directory=directory/'work')
    executor = PartitionedQResponse(r,history,directory=directory/'executor',**kwargs)
    length = len(pack_legs(r.zero_tangent()))
    x = mapped(r,directory/'x',budget,.07*system.random((length,)))
    y = mapped(r,directory/'y',budget,system.random((length,)))
    ax, ay = mapped(r,directory/'Ax',budget), mapped(r,directory/'AstarY',budget)
    time = .0173
    executor.rhs(time,x,ax)
    executor.rhs(time,y,ay,dual=True)
    check(tag+'_actual_discrete_generator_adjoint',y.dot(ax),ay.dot(x),tolerance=2e-11)
    # Independent full-density derivative includes moving frames' derivatives,
    # not simply the leg generator mapped at a frozen time.
    stage = system.stage(time)
    physical = dense.physical(stage,ax)
    blocks = [l @ (-1j*stage.h_actions[e.source](stage.frames[e.source])).conj().T
              +(-1j*stage.h_actions[e.target](stage.frames[e.target])) @ b
              for e,l,b in zip(r.edges,x.ket,x.bra)]
    physical += dense.pack(blocks,np.zeros(r.ndof),np.zeros(r.ndof),np.zeros(r.nfield),np.zeros(r.nfield))
    check(tag+'_original_full_density_RHS',physical,dense.rhs(time,dense.physical(stage,x)),tolerance=2e-11)
    history.wrong_adjoint = True
    executor.rhs(time,y,ay,dual=True)
    check(tag+'_wrong_H_adjoint_sensitive',y.dot(ax),ay.dot(x),mutation=True)
    history.wrong_adjoint = False
    executor.dynamics = 'bare-electronic'
    executor.rhs(time,x,ay)
    check(tag+'_feedback_sensitive',pack_legs(ay),pack_legs(ax),mutation=True)
    executor.dynamics = 'screened-joint'
    if gamma:
        old = r.sqrt_n
        r.sqrt_n = 1.
        executor.rhs(time,x,ay)
        check(tag+'_sqrtN_sensitive',pack_legs(ay),pack_legs(ax),mutation=True)
        r.sqrt_n = old
    else:
        old = r.edges[0].coupling.v.copy()
        r.edges[0].coupling.v[:] = old.conj()
        executor.rhs(time,x,ay)
        check(tag+'_directed_q_not_minus_q_sensitive',pack_legs(ay),pack_legs(ax),mutation=True)
        r.edges[0].coupling.v[:] = old
    errors, finals = [], []
    t, s = .027, .019
    expected = dense.transition(t) @ dense.physical(system.stage(0),x)
    for steps in (4,8):
        final, checkpoints = executor.propagate(x,np.linspace(0,t,steps+1),directory/('forward-%d'%steps))
        finals.append(final)
        errors.append(float(np.linalg.norm(dense.physical(system.stage(t),final)-expected)))
        restored,_ = MappedResponse.restore(r,checkpoints[-1],budget,executor.identity)
        check(tag+'_immutable_checkpoint_%d'%steps,pack_legs(restored),pack_legs(final),tolerance=0.)
    check(tag+'_full_density_exact_consumer',dense.physical(system.stage(t),finals[-1]),expected,tolerance=2e-8)
    if not errors[1] < errors[0]/8:
        raise AssertionError('RK4 refinement failed: '+repr(errors))
    CHECKS.append(dict(name=tag+'_RK4_refinement',passed=True,errors=errors))
    dual,_ = executor.propagate(y,np.linspace(0,t,5),directory/'dual',dual=True)
    check(tag+'_RK4_actual_reverse_graph',y.dot(finals[0]),dual.dot(x),tolerance=2e-11)
    preparation = JointProductPreparation(r,system.stage(0),ionic_ordered=system.wi,
        ionic_minus_ordered=system.wi,field_ordered=system.wf if gamma else np.zeros((0,0)),
        electronic_identity='raw-singular-D',ionic_identity='initial-ionic-W',
        field_identity='initial-field-W',product_identity='original-product')
    left, right = dense.probe(system.random((dense.size,)),tag+'-A'), dense.probe(system.random((dense.size,)),tag+'-B')
    covariance = PartitionedTwoSidedCovariance(preparation,executor)
    a, ar = covariance.pullback(left,np.linspace(0,t,5),directory/'left')
    b, _ = covariance.pullback(right,np.linspace(0,s,5),directory/'right')
    result = initial_contraction(preparation,a,b)
    pa,pb = dense.pullback(t,left),dense.pullback(s,right)
    for ordering in ('AB','BA'):
        parts = fock.parts(r,pa,pb,ordering)
        for origin in ('electronic','ionic','field'):
            check(tag+'_Fock_'+ordering+'_'+origin,result[ordering][origin],parts[origin],tolerance=2e-8)
        ca = dense.transition(t).conj().T @ dense.coefficient(left)
        cb = dense.transition(s).conj().T @ dense.coefficient(right)
        check(tag+'_dense_two_time_'+ordering,result[ordering]['total'],np.vdot(ca,dense.covariance(ordering) @ cb),tolerance=2e-8)
    b0 = mapped(r,directory/'b0',budget)
    executor.terminal_probe(0.,right,b0)
    check(tag+'_genuine_second_pullback',initial_contraction(preparation,a,b0)['AB']['total'],result['AB']['total'],mutation=True)
    check(tag+'_N_weight_sensitive',result['AB']['electronic']/r.ncells,result['AB']['electronic'],mutation=True)
    scheduled = PartitionedQResponse(r,DiskFiniteHistory(system,r,budget),directory=directory/'scheduled',
        **dict(kwargs,work_directory=directory/'scheduled-work'))
    retire_frames = scheduled.retire_frames_except
    protected = []
    def checked_frame_retirement(needed, barrier):
        keep = []
        for candidate in scheduled.frames.directory.glob('*/generation-*'):
            if not (candidate/'retire.json').exists() and (candidate/'request.json').exists():
                request = json.loads((candidate/'request.json').read_text())
                if request['input']['time'] in needed:
                    for receipt in candidate.glob('item-*/complete.json'):
                        keep.extend(json.loads(receipt.read_text())['arrays'].values())
        retire_frames(needed,barrier)
        for binding in keep:
            scheduled.frames.verify(binding)
        protected.append(len(keep))
    scheduled.retire_frames_except = checked_frame_retirement
    scheduled_result = PartitionedTwoSidedCovariance(preparation,scheduled).pair(
        left,right,np.linspace(0,t,5),np.linspace(0,s,5),directory/'scheduled-pair')
    for ordering in ('AB','BA'):
        check(tag+'_scheduled_frame_retirement_'+ordering,
              scheduled_result[ordering]['total'],result[ordering]['total'],tolerance=0.)
    if not sum(protected) or list(scheduled.frames.directory.rglob('*.npy')):
        raise AssertionError('future frame protection or last-use retirement was not exercised')
    CHECKS.append(dict(name=tag+'_both_leg_frame_last_use_schedule',passed=True,
                      protected_future_bindings=sum(protected),final_live_frame_files=0))
    # Fresh executor/fresh frame adapter on the same immutable work directory.
    fresh = PartitionedQResponse(r,DiskFiniteHistory(system,r,budget),directory=directory/'fresh',**kwargs)
    before = budget.actions
    resumed,_ = PartitionedTwoSidedCovariance(preparation,fresh).pullback(left,np.linspace(0,t,5),directory/'resumed',checkpoint=ar[-1])
    check(tag+'_disk_dual_restart',pack_legs(resumed),pack_legs(a),tolerance=0.)
    check(tag+'_restart_reuses_completed_H_work',budget.actions,before,tolerance=0.)
    writes_before = budget.snapshot()['cumulative_write_bytes_reserved']
    automatic,_ = PartitionedTwoSidedCovariance(preparation,fresh).pullback(left,np.linspace(0,t,5),directory/'automatic')
    check(tag+'_automatic_interval_restart',pack_legs(automatic),pack_legs(a),tolerance=0.)
    check(tag+'_automatic_restart_reuses_completed_H_work',budget.actions,before,tolerance=0.)
    check(tag+'_automatic_restart_without_seed_copy',
          budget.snapshot()['cumulative_write_bytes_reserved'],writes_before,tolerance=0.)
    fresh.rhs(time,x,ay)
    check(tag+'_fresh_disk_stage_recovery',pack_legs(ay),pack_legs(ax),tolerance=0.)
    x.ket[0][0,0] += .1j
    fresh.rhs(time,x,ay)
    check(tag+'_mutable_scratch_not_completed_input',pack_legs(ay),pack_legs(ax),mutation=True)
    return r, executor, x, ax


def recovery(system, r, directory, budget, cpus):
    history = DiskFiniteHistory(system,r,budget)
    calls = mp.get_context('fork').Value('q',0)
    old_kernel = r.kernel
    def counted_kernel(n, **kw):
        with calls.get_lock():
            calls.value += 1
        return old_kernel(n, **kw)
    r.kernel = counted_kernel
    kwargs = dict(kernel_adjoint=lambda n,**kw:system.k*n,grid_measure=1.,budget=budget,
                  identity='barrier-recovery',work_directory=directory/'work')
    executor = PartitionedQResponse(r,history,directory=directory/'first',cpus=cpus[:1],**kwargs)
    x = mapped(r,directory/'x',budget,system.random((len(pack_legs(r.zero_tangent())),)))
    out = mapped(r,directory/'out',budget)
    actual = executor.independent_edge
    def interrupted(index,time,state,dual):
        if index == 1:
            raise RuntimeError('oracle deliberate interruption after completed edge zero')
        return actual(index,time,state,dual)
    executor.independent_edge = interrupted
    try:
        executor.rhs(.011,x,out)
    except RuntimeError as error:
        if 'q worker failed' not in str(error):
            raise
        CHECKS.append(dict(name='real_worker_interruption',passed=True,error=str(error)))
    else:
        raise AssertionError('deliberate worker interruption did not stop the stage')
    key,path,count = executor.edge_work.last
    first = executor.edge_work.receipt(path,key,0)
    if first is None:
        raise AssertionError('completed first edge was not retained')
    reject('incomplete_barrier_rejects_feedback',lambda:executor.edge_work.collect(path,key,count))
    check('no_partial_global_kernel',calls.value,0,tolerance=0.)
    fresh = PartitionedQResponse(r,DiskFiniteHistory(system,r,budget),directory=directory/'second',cpus=cpus,**kwargs)
    fresh.rhs(.011,x,out)
    check('one_complete_global_feedback',calls.value,1,tolerance=0.)
    if json.loads((path/'item-0000'/'complete.json').read_text()) != first:
        raise AssertionError('restart replaced the compact completed scientific edge receipt')
    reject('retired_work_not_checkpoint',lambda:fresh.edge_work.receipt(path,key,0))
    if not (path/'retire.json').exists() or Path(first['arrays']['ket']['path']).exists():
        raise AssertionError('obsolete RHS partial was not retired after durable full RHS')
    check('recovered_actual_equations',pack_legs(out),pack_legs(r.rhs(system.stage(.011),x)),tolerance=2e-11)
    CHECKS.append(dict(name='immutable_completed_edge_reused',passed=True,receipt=first,obsolete_payloads_retired=True))
    # All five CPUs perform real original Hamiltonian actions, not echo tasks.
    worker = QWork(directory/'cpu-proof',budget,cpus,'original-finite-H-actions')
    hs = system.mean_values(system.mean.sol(.013))[-1]
    vectors = [system.random((h.shape[1],128)) for h in hs]
    def action(index):
        h,v = hs[index%len(hs)],vectors[index%len(hs)]
        budget.action(v)
        return dict(value=h @ v)
    records = worker.execute('independent-original-H',dict(h=[array_identity(h) for h in hs],
        v=[array_identity(v) for v in vectors]),2*len(cpus),action)
    for index,record in enumerate(records):
        check('CPU_action_%d'%index,arrays(record)['value'],hs[index%len(hs)] @ vectors[index%len(hs)],tolerance=0.)
        if record['affinity'] != [record['cpu']] or record['pid'] == os.getpid() or record['cpu_seconds'] <= 0:
            raise AssertionError('actual pinned child CPU execution missing')
    if set(record['cpu'] for record in records) != set(cpus):
        raise AssertionError('not all allocated workers performed real CPU work')
    CHECKS.append(dict(name='all_allocated_CPUs_actual_work',passed=True,workers=records))
    # Deliberately corrupt only this disposable oracle work product.
    rhs_record = json.loads(Path(fresh.last_rhs['path']).read_text())
    binding = rhs_record['records'][0]['arrays']['ket']
    with Path(binding['path']).open('ab') as stream:
        stream.write(b'corruption')
    reject('modified_completed_array_rejected',lambda:fresh.rhs(.011,x,out))
    r.kernel = old_kernel


def stage_recovery(system, r, directory, budget, cpus):
    """The existing ordered RK4 graph, including an interrupted mixed edge prefix."""
    for dual in (False, True):
        name = 'dual' if dual else 'forward'
        root = directory/name
        kwargs = dict(kernel_adjoint=lambda n,**kw:system.k*n, grid_measure=1., budget=budget,
            identity='bounded-stage-'+name, cpus=cpus, work_directory=root/'work')
        executor = PartitionedQResponse(r,DiskFiniteHistory(system,r,budget),directory=root/'first',**kwargs)
        length = len(pack_legs(r.zero_tangent()))
        x = mapped(r,root/'x',budget,system.random((length,)))
        expected = mapped(r,root/'expected',budget)
        # Execute the preexisting graph/order with exactly the same partitioned
        # RHS, isolating the storage refactor from alternative RHS arithmetic.
        reference = StreamedJointResponse(r,executor.edge_stage,
            kernel_adjoint=kwargs['kernel_adjoint'],grid_measure=1.,budget=budget,
            directory=root/'reference',identity=executor.identity)
        reference.rhs = executor.rhs
        reference.step(x,expected,0.,.009,dual=dual)
        original_retire = q_response_module.retire_payloads
        interrupted = []
        def interrupt_after_edge(receipt, barrier, paths, owner, allocation):
            result = original_retire(receipt,barrier,paths,owner,allocation)
            if (Path(receipt).name == 'retire-edge-0000.json'
                    and Path(owner).name == 'accumulator' and Path(owner).parent.name == 'combine-0'
                    and not interrupted):
                interrupted.append(barrier)
                raise RuntimeError('oracle interruption after durable edge replacement')
            return result
        q_response_module.retire_payloads = interrupt_after_edge
        try:
            executor.propagate(x,[0.,.009],root/'first-result',dual=dual)
        except RuntimeError as error:
            if 'after durable edge replacement' not in str(error):
                raise
        else:
            raise AssertionError('edge-prefix interruption was not exercised')
        finally:
            q_response_module.retire_payloads = original_retire
        prefix = interrupted[0]
        prefix_before = file_binding(prefix['path'])
        # Corrupting uncommitted scratch cannot alter the recovered answer.
        prefix_directory = Path(prefix['path']).parent
        scratch = prefix_directory/'accumulator'/'leg-0002.npy'
        with scratch.open('xb') as stream:
            stream.write(b'not an NPY checkpoint')
        budget.reserve(scratch.stat().st_size)
        fresh = PartitionedQResponse(r,DiskFiniteHistory(system,r,budget),directory=root/'fresh',**kwargs)
        calls = []
        actual_rhs = fresh.rhs
        def counted_rhs(time, *args, **options):
            calls.append(time)
            return actual_rhs(time,*args,**options)
        fresh.rhs = counted_rhs
        final, checkpoints = fresh.propagate(x,[0.,.009],root/'recovered',dual=dual)
        check(name+'_stagewise_exact_existing_RK4_order',pack_legs(final),pack_legs(expected),tolerance=0.)
        if file_binding(prefix['path']) != prefix_before or len(calls) != 2:
            raise AssertionError('completed RHS stages/edge prefix were not reused')
        CHECKS.append(dict(name=name+'_completed_prefix_and_RHS_reused',passed=True,
                          durable_prefix=prefix_before,rhs_calls_after_restart=calls))
        check(name+'_only_live_final_checkpoint_returned',len(checkpoints),1,tolerance=0.)
        live_legs = list((root/'work'/'propagation').rglob('*.npy'))
        check(name+'_superseded_stage_payloads_retired',len(live_legs),2*len(r.edges),tolerance=0.)
        frame_files = list((root/'work'/'frames').rglob('*.npy'))
        if not frame_files:
            raise AssertionError('baseline frames needed by the second leg were silently retired')
        CHECKS.append(dict(name=name+'_frame_cache_retained_across_legs',passed=True,files=len(frame_files)))
        before = budget.actions
        again,_ = fresh.propagate(x,[0.,.009],root/'again',dual=dual)
        check(name+'_complete_propagation_reused',pack_legs(again),pack_legs(final),tolerance=0.)
        check(name+'_complete_recovery_no_H_recomputation',budget.actions,before,tolerance=0.)
        record = json.loads(Path(checkpoints[0]['path']).read_text())
        with Path(record['legs'][0]['path']).open('ab') as stream:
            stream.write(b'corruption')
        reject(name+'_corrupt_stage_checkpoint_rejected',
               lambda:MappedResponse.restore(r,checkpoints[0],budget,fresh.identity))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--cpus',nargs='+',type=int,choices=range(6,11),required=True)
    parser.add_argument('--max-seconds',type=float,default=1800.)
    parser.add_argument('--max-output-bytes',type=int,default=1000000000)
    parser.add_argument('--max-memory-bytes',type=int,default=12000000000)
    parser.add_argument('--worker-memory-bytes',type=int,default=2000000000)
    args = parser.parse_args()
    os.sched_setaffinity(0,{args.cpus[0]})
    args.work.mkdir(parents=True,exist_ok=False)
    budget = QBudget(max_actions=1000000,max_seconds=args.max_seconds,
        max_bytes=args.max_output_bytes,max_memory_bytes=args.max_memory_bytes,
        worker_memory_bytes=args.worker_memory_bytes)
    registration = json.loads((HERE/'joint-partitioned-q-registration.json').read_text())
    if registration['schema'] != 'P254-joint-partitioned-q-registration-v1':
        raise ValueError('prospective partitioned q registration required')
    result = dict(schema='P254-joint-partitioned-q-oracle-v1',status='failed',
                  source_bindings=bindings(),checks=CHECKS,cpus=args.cpus)
    try:
        symbolic()
        kernel_oracle()
        system = FiniteSystem()
        fock = FockOracle(system)
        exercise(system,fock,True,args.work/'Gamma',budget,args.cpus)
        r,*_ = exercise(system,fock,False,args.work/'nonGamma',budget,args.cpus)
        recovery(system,r,args.work/'recovery',budget,args.cpus)
        stage_recovery(system,r,args.work/'stage-recovery',budget,args.cpus)
        counters = budget.snapshot()
        measured = sum(path.stat().st_size for path in args.work.rglob('*') if path.is_file())
        if (counters['retired_payload_bytes'] <= 0 or counters['cumulative_write_bytes_reserved'] <= counters['live_bytes_reserved']
                or counters['peak_live_bytes_reserved'] < counters['live_bytes_reserved']
                or measured > counters['live_bytes_reserved']):
            raise AssertionError('live/peak/cumulative storage accounting differs from actual ownership')
        CHECKS.append(dict(name='bounded_live_storage_accounting',passed=True,
                          measured_file_bytes=measured,counters=counters))
        result['status'] = 'passed'
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result['budget'] = budget.snapshot()
    with args.output.open('x') as stream:
        json.dump(encode(result),stream,indent=2,allow_nan=False)
        stream.write('\n')
    if result['status'] != 'passed':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()

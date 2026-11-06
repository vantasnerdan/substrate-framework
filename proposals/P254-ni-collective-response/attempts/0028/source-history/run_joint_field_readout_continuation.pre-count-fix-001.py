#!/usr/bin/env python3
"""Finish only the two uncompleted field pairs from the authentic 129-check prefix.

The existing 56 producer sources, original result and completed work receipts are
not edited. Production resume authentication owns unfinished per-pair work;
new invocation paths avoid collisions with the timed-out executor's scratch.
"""
import argparse
import json
import os
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from run_joint_field_readout_oracle import (
    FiniteSystem, DenseTangent, FockOracle, CHECKS, check, reject,
    DiskFiniteHistory, JointProductPreparation, UniformMaxwellField,
    PartitionedQResponse, PartitionedTwoSidedCovariance, QBudget, MappedResponse,
    global_field_probe, bindings, encode, producer_identity, admit_oracle,
    file_binding, bound_json, retain_original)

REGISTRATION = HERE/'joint-field-readout-continuation-registration-001.json'
PREFIX_SHA = '255a17b14c4acd2569b1f4a260d22c6d617b141228ffd3a92dcd8abf78fc7ef2'


def authenticate_prefix(registration):
    if (registration['schema'] != 'P254-joint-field-readout-continuation-registration-v1'
            or registration['prefix']['sha256'] != PREFIX_SHA):
        raise ValueError('exact actual failed-prefix preregistration required')
    prefix = bound_json(registration['prefix'])
    current = bindings()
    if (prefix['schema'] != 'P254-joint-partitioned-q-oracle-v1'
            or prefix['extension'] != 'physical-global-field-readouts'
            or prefix['status'] != 'failed'
            or prefix['error'] != 'RuntimeError: q worker failed; immutable completed work retained'
            or prefix['new_check_count'] != 129 or len(prefix['checks']) != 129
            or any(c.get('passed') is not True for c in prefix['checks'])
            or [c['name'] for c in prefix['checks']] != registration['prefix_check_names']
            or prefix['source_bindings'] != registration['producer_sources']
            or current != prefix['source_bindings'] or len(current) != 56):
        raise ValueError('actual prefix checks or exact producer source bindings differ')
    original_registration = json.loads((HERE/'joint-field-readout-registration-001.json').read_text())
    original = retain_original(original_registration)
    if (prefix['retained_check_count'] != 100 or prefix['retained_checks'] != original['checks']
            or prefix['retained_proof'] != original_registration['retained_proof']):
        raise ValueError('original settled 100-check retention differs')
    if prefix['checks'][-1]['name'] != 'field-coordinate__field-coordinate_foreign_observable_dual_checkpoint_rejected':
        raise ValueError('continuation boundary differs')
    return prefix, current


def complete_remaining(original_work, invocation, budget, cpus, producer_sources):
    system = FiniteSystem()
    response = system.response(True)
    dense = DenseTangent(system,response)
    # Rebuilding reference data is necessary, not replaying the completed pair.
    # FockOracle performs one constructor trace assertion: disclose it as setup,
    # never count it again among the 129 retained or 39 remaining checks.
    fock = FockOracle(system)
    setup_checks = list(CHECKS)
    CHECKS.clear()
    field = UniformMaxwellField(8*np.pi*system.mf,system.e,'finite physical Maxwell volume')
    probes, metadata, independent = {}, {}, {}
    zero_blocks = [np.zeros(shape,complex) for shape in dense.shapes]
    for observable in ('field-coordinate','field-momentum'):
        probes[observable],metadata[observable] = global_field_probe(response,field,observable=observable,axis='x')
        coefficient = dense.pack(zero_blocks,np.zeros(3),np.zeros(3),
            np.array([.8 if observable == 'field-coordinate' else 0.]),
            np.array([.8 if observable == 'field-momentum' else 0.]))
        independent[observable] = dense.probe(coefficient,'independent-'+observable)
    _,y_metadata = global_field_probe(response,field,observable='field-momentum',axis='y')
    identity = lambda meta:producer_identity('finite-field-model',[0.,0.,0.],
        'screened-joint',producer_sources,meta)
    preparation = JointProductPreparation(response,system.stage(0),ionic_ordered=system.wi,
        ionic_minus_ordered=system.wi,field_ordered=system.wf,electronic_identity='original-finite-D',
        ionic_identity='original-finite-ion-W',field_identity='original-finite-field-W',
        product_identity='original-finite-product')
    t,s = .09,.06
    pair_receipts = []
    for left_name,right_name in (('field-momentum','field-momentum'),('field-coordinate','field-momentum')):
        tag = left_name+'__'+right_name
        work = original_work/tag/'work'
        executor = PartitionedQResponse(response,DiskFiniteHistory(system,response,budget),
            kernel_adjoint=lambda n,**kw:system.k*n,grid_measure=1.,budget=budget,
            directory=invocation/tag/'executor',work_directory=work,cpus=cpus,identity=identity(metadata[left_name]))
        before = budget.snapshot()
        pair = PartitionedTwoSidedCovariance(preparation,executor).pair(probes[left_name],probes[right_name],
            np.linspace(0,t,17),np.linspace(0,s,17),invocation/(tag+'-pair'))
        a,b = dense.pullback(t,independent[left_name]),dense.pullback(s,independent[right_name])
        ca,cb = dense.coefficient(a),dense.coefficient(b)
        for ordering in ('AB','BA'):
            parts = fock.parts(response,a,b,ordering)
            for origin in ('electronic','ionic','field'):
                check(tag+'_'+ordering+'_Fock_'+origin,pair[ordering][origin],parts[origin])
                check(tag+'_'+ordering+'_actual_initial_origin_'+origin,pair[ordering][origin],0.,
                      mutation=True,tolerance=1e-14)
            expected = np.vdot(ca,dense.covariance(ordering) @ cb)
            check(tag+'_'+ordering+'_independent_dense',pair[ordering]['total'],expected)
            check(tag+'_'+ordering+'_missing_field_origin_sensitive',
                  pair[ordering]['total']-pair[ordering]['field'],expected,mutation=True)
            check(tag+'_'+ordering+'_wrong_N_normalization_sensitive',
                  response.ncells*pair[ordering]['total'],expected,mutation=True)
        reject(tag+'_foreign_observable_dual_checkpoint_rejected',lambda:MappedResponse.restore(
            response,pair['left_checkpoints'][-1],budget,identity(y_metadata)))
        pair_receipts.append(dict(pair=tag,immutable_work=str(work),result=pair,
            invocation_budget_before=before,invocation_budget_after=budget.snapshot()))
    check('coordinate_momentum_covariance_swap_sensitive',
          np.vdot(dense.coefficient(independent['field-coordinate']),
                  dense.covariance('AB') @ dense.coefficient(independent['field-coordinate'])),
          np.vdot(dense.coefficient(independent['field-momentum']),
                  dense.covariance('AB') @ dense.coefficient(independent['field-momentum'])),mutation=True)
    return setup_checks,pair_receipts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--invocation',type=Path,required=True)
    parser.add_argument('--cpus',nargs='+',type=int,choices=range(6,11),required=True)
    parser.add_argument('--max-seconds',type=float,default=1800.)
    parser.add_argument('--max-output-bytes',type=int,default=1000000000)
    parser.add_argument('--max-memory-bytes',type=int,default=4294967296)
    parser.add_argument('--worker-memory-bytes',type=int,default=1610612736)
    args = parser.parse_args()
    if len(set(args.cpus)) != len(args.cpus):
        raise ValueError('unique allocated CPUs required')
    registration = json.loads(REGISTRATION.read_text())
    prefix,producer_sources = authenticate_prefix(registration)
    original_work = Path(registration['original_work']).resolve()
    if not original_work.is_dir() or args.invocation.resolve().is_relative_to(original_work):
        raise ValueError('actual original work and separate fresh invocation root required')
    if args.output.exists():
        raise ValueError('never overwrite an oracle result')
    os.sched_setaffinity(0,{args.cpus[0]})
    args.invocation.mkdir(parents=True,exist_ok=False)
    live_work_bytes = sum(path.stat().st_size for path in original_work.rglob('*') if path.is_file())
    existing_bytes = live_work_bytes+int(registration['prefix']['bytes'])
    budget = QBudget(max_actions=1000000,max_seconds=args.max_seconds,max_bytes=args.max_output_bytes,
        max_memory_bytes=args.max_memory_bytes,worker_memory_bytes=args.worker_memory_bytes,existing_bytes=existing_bytes)
    continuation_sources = [file_binding(Path(__file__)),file_binding(REGISTRATION)]
    result = dict(schema=prefix['schema'],status='failed',extension=prefix['extension'],material_result=False,
        source_bindings=producer_sources+continuation_sources,producer_source_bindings=producer_sources,
        continuation_source_bindings=continuation_sources,cpus=args.cpus,
        retained_proof=prefix['retained_proof'],retained_checks=prefix['retained_checks'],retained_check_count=100,
        retained_prefix=registration['prefix'],retained_prefix_checks=prefix['checks'],retained_prefix_check_count=129,
        continuation_checks=CHECKS,original_invocation_budget=prefix['budget'],
        original_external_cost=registration['original_external_cost'],
        existing_live_bytes_charged=dict(work=live_work_bytes,prefix_receipt=int(registration['prefix']['bytes']),total=existing_bytes),
        scope='Only two unfinished pairs and final swap; exact retained original100 and field129; no material result')
    try:
        setup,pairs = complete_remaining(original_work,args.invocation,budget,args.cpus,producer_sources)
        result.update(reference_setup_checks=setup,completed_pairs=pairs)
        if len(CHECKS) != registration['remaining_check_count'] or any(c.get('passed') is not True for c in CHECKS):
            raise ValueError('remaining finite check coverage differs')
        if set(c['name'] for c in CHECKS) & set(registration['prefix_check_names']):
            raise ValueError('a retained field check was replayed as new evidence')
        authenticate_prefix(registration)
        if [file_binding(Path(__file__)),file_binding(REGISTRATION)] != continuation_sources:
            raise ValueError('continuation source moved during execution')
        result['checks'] = prefix['checks']+CHECKS
        result['new_check_count'] = len(result['checks'])
        completed = dict(result,status='passed')
        for observable in ('total-force','field-coordinate','field-momentum'):
            admit_oracle(completed,producer_sources,observable)
        result['status'] = 'passed'
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result['checks'] = prefix['checks']+CHECKS
    result['new_check_count'] = len(result['checks'])
    result['continuation_check_count'] = len(CHECKS)
    result['budget'] = budget.snapshot()
    additive = ('hamiltonian_actions','action_columns','elapsed_seconds','cumulative_write_bytes_reserved',
                'cumulative_allocation_bytes_reserved','mapped_overwrite_bytes_reserved','retired_payload_bytes')
    result['cumulative_cost'] = {key:prefix['budget'][key]+result['budget'][key] for key in additive}
    result['cumulative_cost'].update(live_bytes_reserved=result['budget']['live_bytes_reserved'],
        peak_live_bytes_reserved=max(prefix['budget']['peak_live_bytes_reserved'],result['budget']['peak_live_bytes_reserved']),
        note='Original and continuation invocation costs added, not reset; initial existing bytes are live occupancy, not new writes. Budget elapsed excludes unmeasured process overhead; original external wall/CPU recorded separately.')
    with args.output.open('x') as stream:
        json.dump(encode(result),stream,indent=2,allow_nan=False)
        stream.write('\n')
    if result['status'] != 'passed':
        raise SystemExit(2)


if __name__ == '__main__':
    from threadpoolctl import threadpool_limits
    with threadpool_limits(limits=1):
        main()

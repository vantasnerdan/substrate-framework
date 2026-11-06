#!/usr/bin/env python3
"""Changed-path physical field oracle; retain authenticated proof003, never replay it.

Run only by the integration owner. Finite moving frames are genuine oracle
trajectories, not material results or evidence of field energy transfer.
"""
import argparse
from copy import deepcopy
from hashlib import sha256
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
from run_joint_response_oracle import FiniteSystem, DenseTangent, FockOracle, CHECKS, check, reject, pack_legs
from run_joint_partitioned_q_oracle import DiskFiniteHistory
from run_joint_q_oracle import mapped
from run_joint_q_material import bindings, encode, producer_identity, admit_oracle
from partitioned_material_readout import partitioned_selected_readout
from substrate_framework.joint_field_probes import global_field_probe, validate_selection
from substrate_framework.joint_ordered_sources import JointProductPreparation
from substrate_framework.joint_response_observables import local_momentum_kick
from substrate_framework.physical_field_environment import UniformMaxwellField
from substrate_framework.partitioned_q_response import PartitionedQResponse, PartitionedTwoSidedCovariance
from substrate_framework.partitioned_q_storage import QBudget
from substrate_framework.streamed_joint_response import MappedResponse
from substrate_framework.streamed_material_state import file_binding
from substrate_framework.material_control_reference import bound_json


RETAINED_SHA = 'd1c495cba5ae8124f83658a00023a2b1305218f8ab270f565dc1b0a1851e6d9c'


def exact_binding(actual, expected):
    if actual['sha256'] != expected['sha256'] or int(actual['bytes']) != int(expected['bytes']):
        raise ValueError('retained numerical source differs: '+expected['path'])


def retain_original(registration):
    binding = registration['retained_proof']
    if binding['sha256'] != RETAINED_SHA:
        raise ValueError('only the actual settled proof003 may be retained')
    prior = bound_json(binding)
    if (prior.get('schema') != 'P254-joint-partitioned-q-oracle-v1'
            or prior.get('status') != 'passed' or len(prior['checks']) != 100
            or any(c.get('passed') is not True for c in prior['checks'])):
        raise ValueError('complete actual 100-check proof003 required')
    original = bound_json(registration['original_sources'])
    if original['proof'] != binding or original['source_bindings'] != [
            dict(x, bytes=int(x['bytes'])) for x in prior['source_bindings']]:
        raise ValueError('pre-edit source inventory does not exactly bind proof003')
    producer = str(HERE/'run_joint_q_material.py')
    snapshot = registration['producer_snapshot']
    if (original['producer_snapshot']['snapshot'] != snapshot['path']
            or original['producer_snapshot']['original']['path'] != producer):
        raise ValueError('original producer snapshot authority differs')
    for source in prior['source_bindings']:
        path = snapshot['path'] if source['path'] == producer else source['path']
        actual = file_binding(path)
        exact_binding(actual, source)
        if source['path'] == producer:
            exact_binding(actual, snapshot)
    return prior


def new_executor(system, response, directory, budget, cpus, identity):
    return PartitionedQResponse(response, DiskFiniteHistory(system,response,budget),
        kernel_adjoint=lambda n,**kw:system.k*n, grid_measure=1., budget=budget,
        directory=directory/'executor', work_directory=directory/'work', cpus=cpus, identity=identity)


def exercise(directory, budget, cpus):
    system = FiniteSystem()
    response = system.response(True)
    dense = DenseTangent(system,response)
    fock = FockOracle(system)
    field = UniformMaxwellField(8*np.pi*system.mf, system.e, 'finite physical Maxwell volume')
    probes, metadata, independent = {}, {}, {}
    zero_blocks = [np.zeros(shape,complex) for shape in dense.shapes]
    for observable in ('field-coordinate','field-momentum'):
        probe, meta = global_field_probe(response,field,observable=observable,axis='x')
        probes[observable], metadata[observable] = probe, meta
        # Independent Cartesian x projection: original finite E=(.8,.6,0).
        coefficient = dense.pack(zero_blocks,np.zeros(3),np.zeros(3),
            np.array([.8 if observable == 'field-coordinate' else 0.]),
            np.array([.8 if observable == 'field-momentum' else 0.]))
        independent[observable] = dense.probe(coefficient,'independent-'+observable)
        check(observable+'_terminal_coefficient',dense.coefficient(probe),coefficient,tolerance=0.)
        check(observable+'_no_sqrtN_normalization',dense.coefficient(probe)*response.sqrt_n,
              coefficient,mutation=True)
        for index, edge in enumerate(response.edges):
            source, target = response.sectors[edge.source], response.sectors[edge.target]
            check(observable+'_forward_zero_full_PW_'+str(index),
                  probe.actions[index](np.ones((source.shape[0],2))),np.zeros((target.shape[0],2)),tolerance=0.)
            check(observable+'_independent_adjoint_zero_full_PW_'+str(index),
                  probe.adjoints[index](np.ones((target.shape[0],3))),np.zeros((source.shape[0],3)),tolerance=0.)
    check('physical_Maxwell_mass',field.mass,system.mf,tolerance=2e-14)
    check('canonical_coordinate_momentum_not_interchangeable',dense.coefficient(probes['field-coordinate']),
          dense.coefficient(probes['field-momentum']),mutation=True)
    for observable, axis, q in (('field-coordinate','x',[0,0,.25]),('field-momentum',None,[0,0,0]),
                                ('total-force','x',[0,0,0]),('energy','x',[0,0,0]),
                                ('field-coordinate','w',[0,0,0])):
        reject('unsupported_selection_'+str((observable,axis,q)),
               lambda o=observable,a=axis,q=q:validate_selection(o,a,q))
    reject('absent_physical_polarization_rejected',lambda:global_field_probe(
        response,field,observable='field-coordinate',axis='z'))
    reject('nonGamma_response_rejected',lambda:global_field_probe(
        system.response(False),field,observable='field-coordinate',axis='x'))
    reject('wrong_physical_mass_rejected',lambda:global_field_probe(response,
        UniformMaxwellField(field.volume_bohr3*response.ncells,system.e,'wrong volume'),
        observable='field-coordinate',axis='x'))
    reject('wrong_physical_polarizations_rejected',lambda:global_field_probe(response,
        UniformMaxwellField(field.volume_bohr3,np.array([[1.,0.,0.]]),'wrong polarization'),
        observable='field-coordinate',axis='x'))
    sources = bindings()
    identity = lambda meta:producer_identity('finite-field-model',[0.,0.,0.],'screened-joint',sources,meta)
    y_probe, y_metadata = global_field_probe(response,field,observable='field-momentum',axis='y')
    identities = [identity(None),identity(metadata['field-coordinate']),identity(metadata['field-momentum']),identity(y_metadata)]
    mutated = deepcopy(metadata['field-coordinate'])
    mutated['normalization']['mass'] *= response.ncells
    identities.append(identity(mutated))
    if len(set(identities)) != len(identities):
        raise AssertionError('observable/axis/normalization identity collision')
    CHECKS.append(dict(name='observable_axis_normalization_request_identity_separation',passed=True))
    preparation = JointProductPreparation(response,system.stage(0),ionic_ordered=system.wi,
        ionic_minus_ordered=system.wi,field_ordered=system.wf,electronic_identity='original-finite-D',
        ionic_identity='original-finite-ion-W',field_identity='original-finite-field-W',
        product_identity='original-finite-product')
    initial_executor = new_executor(system,response,directory/'canonical',budget,cpus,identity(metadata['field-coordinate']))
    covariance = PartitionedTwoSidedCovariance(preparation,initial_executor)
    for axis, right, factor in (('x',probes['field-momentum'],.64),('y',y_probe,.48)):
        pair = covariance.pair(probes['field-coordinate'],right,[0.],[0.],directory/('canonical-'+axis))
        check('Gamma_canonical_commutator_'+axis,pair['commutator'],1j*factor,tolerance=2e-12)
        for ordering in ('AB','BA'):
            for origin in ('electronic','ionic'):
                check('initial_field_only_'+axis+'_'+ordering+'_'+origin,pair[ordering][origin],0,tolerance=0.)
    t,s = .09,.06
    # Ionic seed stays separate from the measured Cartesian field axis.
    kick = local_momentum_kick(response,[0,0,0],[0,1,0],.7)
    seed_physical = dense.physical(system.stage(0),kick)
    expected_state = dense.transition(t) @ seed_physical
    errors = []
    for steps in (8,16):
        executor = new_executor(system,response,directory/('coherent-'+str(steps)),budget,cpus,'finite-coherent-'+str(steps))
        executor.frame_tail_times = {t}
        initial = mapped(response,directory/('kick-'+str(steps)),budget,pack_legs(kick))
        final,checkpoints = executor.propagate(initial,np.linspace(0,t,steps+1),directory/('forward-'+str(steps)))
        errors.append(float(np.linalg.norm(dense.physical(system.stage(t),final)-expected_state)))
        for observable,probe in probes.items():
            executor.readout_source_bindings = sources
            executor.readout_operator_binding = dict(model='finite-field-model',
                observable=observable, coefficient=encode(dense.coefficient(independent[observable])))
            actual = partitioned_selected_readout(
                executor,final,probe,t,executor.edge_work,state_binding=checkpoints[-1])
            expected = np.vdot(dense.coefficient(independent[observable]),expected_state)
            check(observable+'_actual_coherent_dense_'+str(steps),actual,expected)
            check(observable+'_coherent_not_zero_'+str(steps),actual,0.,mutation=True,tolerance=1e-12)
    if not errors[1] < errors[0]/8:
        raise AssertionError('new field coherent RK4 refinement failed: '+repr(errors))
    CHECKS.append(dict(name='new_field_coherent_RK4_refinement',passed=True,errors=errors))
    for left_name,right_name in (('field-coordinate','field-coordinate'),('field-momentum','field-momentum'),
                                 ('field-coordinate','field-momentum')):
        tag = left_name+'__'+right_name
        executor = new_executor(system,response,directory/tag,budget,cpus,identity(metadata[left_name]))
        pair = PartitionedTwoSidedCovariance(preparation,executor).pair(probes[left_name],probes[right_name],
            np.linspace(0,t,17),np.linspace(0,s,17),directory/(tag+'-pair'))
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
    check('coordinate_momentum_covariance_swap_sensitive',
          np.vdot(dense.coefficient(independent['field-coordinate']),
                  dense.covariance('AB') @ dense.coefficient(independent['field-coordinate'])),
          np.vdot(dense.coefficient(independent['field-momentum']),
                  dense.covariance('AB') @ dense.coefficient(independent['field-momentum'])),mutation=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--cpus',nargs='+',type=int,choices=range(6,11),required=True)
    parser.add_argument('--max-seconds',type=float,default=1800.)
    parser.add_argument('--max-output-bytes',type=int,default=1000000000)
    parser.add_argument('--max-memory-bytes',type=int,default=8589934592)
    parser.add_argument('--worker-memory-bytes',type=int,default=1610612736)
    args = parser.parse_args()
    if len(set(args.cpus)) != len(args.cpus):
        raise ValueError('unique allocated CPUs required')
    os.sched_setaffinity(0,{args.cpus[0]})
    args.work.mkdir(parents=True,exist_ok=False)
    budget = QBudget(max_actions=1000000,max_seconds=args.max_seconds,max_bytes=args.max_output_bytes,
        max_memory_bytes=args.max_memory_bytes,worker_memory_bytes=args.worker_memory_bytes)
    registration = json.loads((HERE/'joint-field-readout-registration-001.json').read_text())
    result = dict(schema='P254-joint-partitioned-q-oracle-v1',status='failed',source_bindings=bindings(),
        checks=CHECKS,cpus=args.cpus,extension='physical-global-field-readouts',material_result=False)
    try:
        if registration['schema'] != 'P254-joint-field-readout-registration-v1':
            raise ValueError('focused preregistration required')
        prior = retain_original(registration)
        result.update(retained_proof=registration['retained_proof'],retained_checks=prior['checks'],
            retained_check_count=len(prior['checks']),
            retention_scope='actual original proof003; every unchanged numerical source exact; changed CLI preserved separately')
        CHECKS.append(dict(name='exact_unchanged_source_retention',passed=True,count=len(prior['source_bindings'])-1))
        wrong = dict(prior['source_bindings'][1],sha256='0'*64)
        reject('retained_source_mutation_rejected',lambda:exact_binding(file_binding(wrong['path']),wrong))
        wrong_registration = deepcopy(registration)
        wrong_registration['retained_proof']['sha256'] = '0'*64
        reject('wrong_original_proof_rejected',lambda:retain_original(wrong_registration))
        candidate = dict(schema=result['schema'],status='passed',source_bindings=result['source_bindings'])
        # This only tests the admission function; it is NOT an emitted proof.
        admit_oracle(candidate,result['source_bindings'])
        reject('unextended_q_proof_cannot_admit_field',lambda:admit_oracle(
            candidate,result['source_bindings'],'field-coordinate'))
        for index,source in enumerate(result['source_bindings']):
            mutated = deepcopy(result['source_bindings'])
            mutated[index]['sha256'] = '0'*64
            reject('admission_source_mutation_'+Path(source['path']).name,
                   lambda mutated=mutated:admit_oracle(candidate,mutated))
        reject('failed_proof_rejected',lambda:admit_oracle(dict(candidate,status='failed'),result['source_bindings']))
        reject('old_proof_cannot_admit_changed_CLI',lambda:admit_oracle(prior,result['source_bindings']))
        exercise(args.work,budget,args.cpus)
        # Do not publish an admissible proof if sources moved during execution.
        if bindings() != result['source_bindings']:
            raise ValueError('sources changed during finite oracle execution')
        completed = dict(result,status='passed',new_check_count=len(CHECKS))
        for observable in ('total-force','field-coordinate','field-momentum'):
            admit_oracle(completed,result['source_bindings'],observable)
        result['status'] = 'passed'
    except Exception as error:
        result['error'] = type(error).__name__+': '+str(error)
    result['new_check_count'] = len(CHECKS)
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

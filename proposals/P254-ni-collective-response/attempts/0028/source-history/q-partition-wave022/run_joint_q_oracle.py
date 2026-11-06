#!/usr/bin/env python3
"""Root-only independent numerical oracle for the real mapped q execution path.

The reused fixture owns an independent full-density DOP853 moving mean/tangent
and finite Fock comparator. Production is exercised through actual mapped files,
streamed per-edge actions and the implemented reverse RK4 graph, never a mock
pullback callback. Dense matrices occur only in this finite oracle.
"""
import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace as NS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE)]
for variable in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np
from run_joint_response_oracle import (FiniteSystem, DenseTangent, FockOracle, check,
    CHECKS, symbolic, pack_legs, unpack_legs, MatrixAction)
from run_joint_q_material import bindings, encode
from substrate_framework.streamed_joint_response import StreamedJointResponse, MappedResponse
from substrate_framework.streamed_joint_covariance import StreamedTwoSidedCovariance, initial_contraction
from substrate_framework.streamed_joint_state import JointBudget
from substrate_framework.streamed_material_state import file_binding
from substrate_framework.joint_ordered_sources import JointProductPreparation
from substrate_framework.joint_response_observables import readout
from substrate_framework.joint_kernel_adjoint import physical_kernel_adjoint
from substrate_framework.spin_pbe_kernel import apply_spin_pbe_arrays
from substrate_framework.qe_material_grid import PhysicalRhoFFT


def mapped(r, path, budget, vector=None):
    out = MappedResponse(r,path,budget,create=True)
    if vector is not None:
        source = unpack_legs(r,vector)
        for a,b in zip(out.values(), [*source.ket,*source.bra,source.u,source.p,source.coordinate,source.momentum]):
            a[:] = b
    return out


class FiniteHistory:
    def __init__(self, system, response, budget):
        self.system, self.response, self.budget = system, response, budget
        self.wrong_adjoint = False

    def counted(self, action):
        def apply(vectors):
            self.budget.action(vectors)
            return action(vectors)
        return apply

    @contextmanager
    def __call__(self, index, time):
        r, s = self.response, self.system
        e, stage = r.edges[index], s.stage(time)
        hs, ht = stage.h_actions[e.source], stage.h_actions[e.target]
        hsa, hta = stage.h_adjoints[e.source], stage.h_adjoints[e.target]
        if self.wrong_adjoint:
            hsa = lambda x: stage.h_actions[e.source](x.conj()).conj()
            hta = lambda x: stage.h_actions[e.target](x.conj()).conj()
        yield (stage.frames[e.source],stage.frames[e.target],self.counted(hs),self.counted(hsa),
               self.counted(ht),self.counted(hta),e.coupling,
               r.field_vertices[e.source].axes if r.gamma else ())


def kernel_oracle():
    rng = np.random.default_rng(2252832)
    miller = np.array(list(np.ndindex((3,3,3))),int).T-1
    fft = PhysicalRhoFFT.create((3,3,3),miller,np.diag([.8,1.1,1.3]),13.)
    # Nonsymmetric arbitrary coefficients deliberately forbid assuming K=K†.
    coefficients = {name:rng.normal(size=(27,2,2)) for name in ('dmuxc','dvxc_rr','dvxc_sr','dvxc_ss','dvxc_s')}
    grho = rng.normal(size=(3,27,2))
    background = NS(fft=fft,coefficients=coefficients,grho=grho,inputs=NS(reference=NS(alat_bohr=5.)))
    q = np.array([.07,-.13,.21])
    def original(x):
        return apply_spin_pbe_arrays(fft.at_q(q),x.reshape(27,2),grho=grho,coefficients=coefficients,tpiba=2*np.pi/5.)
    matrix = np.column_stack([original(x).ravel() for x in np.eye(54,dtype=complex)])
    y = rng.normal(size=(27,2))+1j*rng.normal(size=(27,2))
    expected = (matrix.conj().T @ y.ravel()).reshape(27,2)
    actual = physical_kernel_adjoint(background,y,q_bohr_inv=q)
    check('actual_Native_kernel_adjoint_dense',actual,expected,tolerance=2e-11)
    check('kernel_selfadjoint_assumption_sensitive',original(y),expected,mutation=True)
    check('kernel_q_sign_sensitive',physical_kernel_adjoint(background,y,q_bohr_inv=-q),expected,mutation=True)


def exercise(system, fock, gamma, directory, budget):
    r = system.response(gamma)
    tag = 'Gamma' if gamma else 'nonGamma'
    dense = DenseTangent(system,r)
    history = FiniteHistory(system,r,budget)
    executor = StreamedJointResponse(r,history,kernel_adjoint=lambda n,**kw:system.k*n,
        grid_measure=1.,budget=budget,directory=directory/'executor',identity=tag+'-finite-real-streamed')
    preparation = JointProductPreparation(r,system.stage(0),ionic_ordered=system.wi,ionic_minus_ordered=system.wi,
        field_ordered=system.wf if gamma else np.zeros((0,0)),electronic_identity='original-raw-singular-D',
        ionic_identity='original-oscillator-W',field_identity='original-global-field-W',product_identity='original-product')
    left = dense.probe(system.random((dense.size,)),tag+'-A')
    right = dense.probe(system.random((dense.size,)),tag+'-B')
    length = len(pack_legs(r.zero_tangent()))
    x = mapped(r,directory/'random-tangent',budget,.07*system.random((length,)))
    y = mapped(r,directory/'random-cotangent',budget,system.random((length,)))
    ax, ay = mapped(r,directory/'Ax',budget), mapped(r,directory/'AstarY',budget)
    executor.rhs(.173,x,ax)
    executor.rhs(.173,y,ay,dual=True)
    check(tag+'_instantaneous_actual_dual',y.dot(ax),ay.dot(x),tolerance=2e-11)
    check(tag+'_settled_forward_equation',pack_legs(ax),pack_legs(r.rhs(system.stage(.173),x)),tolerance=2e-11)
    # A wrong H-transpose and omitted feedback must actually separate.
    history.wrong_adjoint = True
    executor.rhs(.173,y,ay,dual=True)
    check(tag+'_adjoint_mutation_sensitive',y.dot(ax),ay.dot(x),mutation=True)
    history.wrong_adjoint = False
    executor.dynamics = 'bare-electronic'
    executor.rhs(.173,x,ay)
    check(tag+'_feedback_mutation_sensitive',pack_legs(ay),pack_legs(ax),mutation=True)
    executor.dynamics = 'screened-joint'
    if gamma:
        original = r.sqrt_n
        r.sqrt_n = 1.
        executor.rhs(.173,x,ay)
        check(tag+'_sqrtN_mutation_sensitive',pack_legs(ay),pack_legs(ax),mutation=True)
        r.sqrt_n = original
    else:
        original = r.edges[0].coupling.v.copy()
        r.edges[0].coupling.v[:] = original.conj()
        executor.rhs(.173,x,ay)
        check(tag+'_directed_q_support_mutation_sensitive',pack_legs(ay),pack_legs(ax),mutation=True)
        r.edges[0].coupling.v[:] = original
    # A full physical-density DOP853 fundamental solution is independent of
    # production leg propagation and all its discretization choices.
    errors = []
    final_time = .27
    expected = dense.transition(final_time) @ dense.physical(system.stage(0),x)
    finals = []
    for steps in (4,8,32):
        final, records = executor.propagate(x,np.linspace(0,final_time,steps+1),directory/('forward-%d' % steps))
        actual = dense.physical(system.stage(final_time),final)
        errors.append(float(np.linalg.norm(actual-expected)))
        finals.append(final)
        restored, record = MappedResponse.restore(r,records[-1],budget,executor.identity)
        check(tag+'_checkpoint_restore_%d' % steps,pack_legs(restored),pack_legs(final),tolerance=0.)
    check(tag+'_full_density_DOP853',dense.physical(system.stage(final_time),finals[-1]),expected,tolerance=2e-8)
    if not errors[1] < errors[0]/8 or not errors[2] < errors[1]/50:
        raise AssertionError(tag+' forward RK4 refinement '+repr(errors))
    CHECKS.append(dict(name=tag+'_forward_refinement',passed=True,errors=errors))
    # Exact adjoint of the numerical map must hold for arbitrary independent
    # legs, not only a source/probe subset or a numerically fitted operator.
    dual, records = executor.propagate(y,np.linspace(0,final_time,9),directory/'random-dual',dual=True)
    check(tag+'_exact_discrete_duality',y.dot(finals[1]),dual.dot(x),tolerance=2e-11)
    covariance = StreamedTwoSidedCovariance(preparation,executor)
    a, ar = covariance.pullback(left,np.linspace(0,.27,33),directory/'pullback-A')
    b, _ = covariance.pullback(right,np.linspace(0,.19,33),directory/'pullback-B')
    equal_b, _ = covariance.pullback(right,np.linspace(0,.27,33),directory/'pullback-B-equal')
    resumed, _ = covariance.pullback(left,np.linspace(0,.27,33),directory/'recovered-A',checkpoint=ar[15])
    check(tag+'_partial_dual_recovery',pack_legs(resumed),pack_legs(a),tolerance=0.)
    from run_joint_response_oracle import reject
    reject(tag+'_wrong_checkpoint_identity',lambda: MappedResponse.restore(r,ar[-1],budget,'wrong-identity'))
    reject(tag+'_wrong_dual_probe',lambda: covariance.pullback(right,np.linspace(0,.27,33),
           directory/'rejected-wrong-probe',checkpoint=ar[15]))
    a0 = mapped(r,directory/'A-initial',budget)
    b0 = mapped(r,directory/'B-initial',budget)
    executor.terminal_probe(0.,left,a0)
    executor.terminal_probe(0.,right,b0)
    initial = initial_contraction(preparation,a0,b0)
    fine = initial_contraction(preparation,a,b)
    equal = initial_contraction(preparation,a,equal_b)
    check(tag+'_readout_duality_physical_probe',readout(r,system.stage(.27),finals[-1],left)['total'],a.dot(x),tolerance=2e-11)
    for label,result,t,s in (('initial',initial,0.,0.),('unequal',fine,.27,.19),('equal',equal,.27,.27)):
        pa, pb = dense.pullback(t,left), dense.pullback(s,right)
        for ordering in ('AB','BA'):
            fock_parts = fock.parts(r,pa,pb,ordering)
            for origin in ('electronic','ionic','field'):
                check(tag+'_'+label+'_'+ordering+'_'+origin,result[ordering][origin],fock_parts[origin],tolerance=2e-8)
            ca = dense.transition(t).conj().T @ dense.coefficient(left)
            cb = dense.transition(s).conj().T @ dense.coefficient(right)
            check(tag+'_'+label+'_'+ordering+'_dense_two_sided',result[ordering]['total'],
                np.vdot(ca,dense.covariance(ordering) @ cb),tolerance=2e-8)
    check(tag+'_nonzero_second_time_sensitive',initial_contraction(preparation,a,b0)['AB']['total'],fine['AB']['total'],mutation=True)
    check(tag+'_weights_and_N_sensitive',fine['AB']['electronic']/r.ncells,fine['AB']['electronic'],mutation=True)
    defects = max(max(d['left'],d['right']) for d in fine['raw_frame_compatibility'])
    check(tag+'_raw_nonorthogonal_compatibility',defects,0.,tolerance=2e-8)
    coarse_a, _ = covariance.pullback(left,np.linspace(0,.27,5),directory/'coarse-A')
    coarse_b, _ = covariance.pullback(right,np.linspace(0,.19,5),directory/'coarse-B')
    coarse = initial_contraction(preparation,coarse_a,coarse_b)
    coarse_defect = max(max(d['left'],d['right']) for d in coarse['raw_frame_compatibility'])
    if not defects < coarse_defect/100:
        raise AssertionError('raw frame compatibility did not refine')
    CHECKS.append(dict(name=tag+'_compatibility_refinement',passed=True,coarse=coarse_defect,fine=defects))
    # Compare the directly contracted initial cotangent to the settled ordered
    # source readout without constructing a full operator from cotangent columns.
    from substrate_framework.joint_response_observables import source_readout
    for ordering in ('AB','BA'):
        check(tag+'_initial_source_contract_'+ordering,initial[ordering]['total'],
              source_readout(r,system.stage(0),preparation.source(right,ordering=ordering),left)['total'],tolerance=2e-11)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--work',type=Path,required=True)
    args = parser.parse_args()
    args.work.mkdir(parents=True,exist_ok=False)
    registration = json.loads((HERE/'joint-q-registration.json').read_text())
    if registration['schema'] != 'P254-joint-q-registration-v1':
        raise ValueError('prospective registration missing')
    result = dict(schema='P254-joint-q-oracle-v1',status='failed',source_bindings=bindings(),checks=CHECKS)
    budget = JointBudget(max_actions=10000000,max_seconds=7200.,max_bytes=1000000000)
    try:
        symbolic()
        kernel_oracle()
        system = FiniteSystem()
        fock = FockOracle(system)
        exercise(system,fock,True,args.work/'Gamma',budget)
        exercise(system,fock,False,args.work/'nonGamma',budget)
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

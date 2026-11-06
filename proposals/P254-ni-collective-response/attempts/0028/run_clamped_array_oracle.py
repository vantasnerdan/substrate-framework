#!/usr/bin/env python3
"""Independent SymPy scalar and SciPy dense-real-space oracle; Root executes."""
from pathlib import Path
import argparse
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(HERE))

import numpy as np
import sympy as sp
from scipy.linalg import expm
from substrate_framework.material_clamped_array import ClampedArray
from run_clamped_array_material import source_bindings, BoundedOutputs


def symbolic_checks():
    n,m,mf,z,k,c,t,f,ell = sp.symbols('N M Mf Z K C T F ell',positive=True)
    q,p,a,pi = sp.symbols('Q P a Pi',real=True)
    # Independently substitute the homogeneous real-cell scalar, then Hessian.
    energy = n*((p/sp.sqrt(n)-z*a)**2/(2*m)+k*(q/sp.sqrt(n))**2/2
                +t*q*a/sp.sqrt(n)-f*q/sp.sqrt(n)+ell*a+c*a*a/2)+pi*pi/(2*mf)
    expected = sp.Matrix([[k,0,sp.sqrt(n)*t,0],[0,1/m,-sp.sqrt(n)*z/m,0],
                          [sp.sqrt(n)*t,-sp.sqrt(n)*z/m,n*(c+z*z/m),0],[0,0,0,1/mf]])
    hessian = sp.hessian(energy,(q,p,a,pi))
    if sp.simplify(hessian-expected) != sp.zeros(4):
        raise AssertionError('symbolic Gamma Hessian differs')
    gradient = sp.Matrix([sp.diff(energy,x).subs({q:0,p:0,a:0,pi:0}) for x in (q,p,a,pi)])
    if sp.simplify(gradient-sp.Matrix([-sp.sqrt(n)*f,0,n*ell,0])) != sp.zeros(4,1):
        raise AssertionError('symbolic affine force normalization differs')
    root = (-1+sp.sqrt(3)*sp.I)/2
    fourier = sp.Matrix(3,3,lambda i,j:root**(-i*j)/sp.sqrt(3))
    pairing = sp.Matrix([[1,0,0],[0,0,1],[0,1,0]])
    if sp.simplify(fourier*fourier.T-pairing) != sp.zeros(3):
        raise AssertionError('ordered Fourier covariance must pair q with minus q')
    if sp.simplify(fourier*fourier.conjugate().T-sp.eye(3)) != sp.zeros(3):
        raise AssertionError('orthonormal canonical Fourier scaling differs')
    return dict(gamma_hessian=str(expected),affine_gradient=str(gradient),
                ordered_q_pairing=str(pairing),p_minus_Za_cross='-sqrt(N)*Z/M')


def fixture():
    # Three distinct periodic cells expose q/-q, oriented cross-cell transport.
    lattice = np.zeros((3,1,1,3,3))
    lattice[0,0,0] = np.diag([2.4,2.9,3.3])
    lattice[1,0,0] = [[-.21,.07,.03],[-.04,-.17,.05],[.02,-.06,-.13]]
    lattice[2,0,0] = lattice[1,0,0].T
    mass = np.array([1.7,2.1,2.6])
    omega = np.sqrt(np.diag(lattice[0,0,0])/mass)
    w = np.block([[np.diag(1/(2*mass*omega)),.5j*np.eye(3)],
                  [-.5j*np.eye(3),np.diag(mass*omega/2)]])
    field_mass = 2.3
    field_width = np.array([1.9,2.2])
    wf = np.block([[np.diag(1/(2*field_width)),.5j*np.eye(2)],
                   [-.5j*np.eye(2),np.diag(field_width/2)]])
    return dict(force_lattice=lattice,masses=mass,charges=np.array([.8]),
        delta_curvature=np.array([[.19,.02,-.01],[.02,.16,.03],[-.01,.03,.12]]),
        contact=np.diag([.6,.7,.8]),mixed=np.array([[.11,.03,-.04],[-.02,.09,.02],[.06,-.03,.08]]),
        force=np.array([.17,-.12,.08]),field_gradient=np.array([.07,-.05,.03]),
        polarizations=np.array([[1.,0.,0.],[0.,.6,.8]]),field_mass=field_mass,
        ionic_ordered=w,field_ordered=wf,trap_omega=omega)


def dense_reference(data, mutation=None):
    """Assemble the real-cell scalar directly, without any production q helper."""
    n,d,f = 3,3,2
    size = 2*n*d+2*f
    h,j,w = (np.zeros((size,size),complex) for _ in range(3))
    b = np.zeros(size,complex)
    a,pi = slice(18,20),slice(20,22)
    e = data['polarizations']
    zmap = data['charges'][0]*e.T
    mass = data['masses']
    for r in range(n):
        u,p = slice(6*r,6*r+3),slice(6*r+3,6*r+6)
        j[u,p],j[p,u] = np.eye(d),-np.eye(d)
        w[6*r:6*r+6,6*r:6*r+6] = data['ionic_ordered']
        h[p,p] = np.diag(1/mass)
        b[u] = -data['force'] if mutation != 'omit_force' else 0
        for s in range(n):
            h[u,6*s:6*s+3] = data['force_lattice'][(r-s)%n,0,0]
            if mutation == 'omit_cross_cell' and r != s:
                h[u,6*s:6*s+3] = 0
        if mutation != 'omit_deltaC':
            h[u,u] += data['delta_curvature']
        cross = data['mixed']@e.T
        if mutation == 'omit_mixed':
            cross = np.zeros_like(cross)
        if mutation == 'wrong_sqrtN':
            cross = np.sqrt(n)*cross
        h[u,a],h[a,u] = cross,cross.T
        crossp = -zmap/mass[:,None]
        if mutation == 'wrong_pZa_sign':
            crossp = -crossp
        h[p,a],h[a,p] = crossp,crossp.T
        h[a,a] += zmap.T@(zmap/mass[:,None])
        if mutation != 'omit_contact':
            h[a,a] += e@data['contact']@e.T
        b[a] += e@data['field_gradient']
    h[pi,pi] = np.eye(f)/data['field_mass']
    j[a,pi],j[pi,a] = np.eye(f),-np.eye(f)
    w[18:,18:] = data['field_ordered']
    return h,j,w,b


def dense_evolve(h,j,b,time):
    aug = np.zeros((len(h)+1,len(h)+1),complex)
    aug[:-1,:-1],aug[:-1,-1] = j@h,j@b
    result = expm(time*aug)
    return result[:-1,:-1],result[:-1,-1]


def flatten_mean(mean):
    return np.r_[mean['ions'].reshape(-1),mean['field']]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True,help='new directory outside repository')
    args = parser.parse_args()
    symbolic = symbolic_checks()
    data = fixture()
    model = ClampedArray(**data)
    h,j,w,b = dense_reference(data)
    errors,sensitivity = {},{}

    def check(name,actual,expected,tolerance=3e-11):
        error = float(np.linalg.norm(actual-expected)/max(1.,np.linalg.norm(expected)))
        errors[name] = error
        if not np.isfinite(error) or error > tolerance:
            raise AssertionError(f'{name}: normalized error {error} > {tolerance}')

    def sensitive(name,wrong,expected,minimum=1e-6):
        error = float(np.linalg.norm(wrong-expected)/max(1.,np.linalg.norm(expected)))
        sensitivity[name] = error
        if not np.isfinite(error) or error < minimum:
            raise AssertionError(f'insensitive negative control {name}: {error}')

    time,other_time = .71,.23
    ev,other,zero = model.evolve(time),model.evolve(other_time),model.evolve(0.)
    s,drift = dense_evolve(h,j,b,time)
    ss,_ = dense_evolve(h,j,b,other_time)
    local,common = model.kicks((1,0,0),1)
    amplitude = np.sqrt(2*data['masses'][1]*data['trap_omega'][1])
    initial_local,initial_common = np.zeros(len(h)),np.zeros(len(h))
    initial_local[6+3+1] = amplitude
    initial_common[[4,10,16]] = amplitude/np.sqrt(3)
    check('baseline_affine_force_and_current',flatten_mean(model.mean(ev)),drift)
    check('local_one_quantum_response',flatten_mean(model.mean(ev,local)),s@initial_local+drift)
    check('equal_energy_common_response',flatten_mean(model.mean(ev,common)),s@initial_common+drift)
    check('one_quantum_energy',initial_local@h@initial_local/2,data['trap_omega'][1])
    check('equal_total_energy',initial_common@h@initial_common/2,initial_local@h@initial_local/2)
    check('canonical_initial_local',model.ordered_pair(zero,zero,(1,0,0),(1,0,0),range(6)),w[6:12,6:12])
    check('initial_independent_cells',model.ordered_pair(zero,zero,(1,0,0),(0,0,0),range(6)),np.zeros((6,6)))
    check('per_cell_mean_ccr',model.w/3-(model.w/3).T,1j*j[:6,:6]/3)
    for r in range(3):
        for v in range(3):
            check(f'two_time_ordered_{r}_{v}',model.ordered_pair(ev,other,(r,0,0),(v,0,0),range(6)),
                  (s@w@ss.T)[6*r:6*r+6,6*v:6*v+6])
            check(f'equal_time_ordered_{r}_{v}',model.ordered_pair(ev,ev,(r,0,0),(v,0,0),range(6)),
                  (s@w@s.T)[6*r:6*r+6,6*v:6*v+6])
    check('global_field_two_time',model.field_ordered_pair(ev,other),(s@w@ss.T)[18:,18:])
    check('ion_global_field_two_time',model.ion_field_ordered_pair(ev,other,range(6)),(s@w@ss.T)[:6,18:])
    for name,kick,x0 in (('baseline',None,np.zeros(len(h))),('local',local,initial_local),('common',common,initial_common)):
        mean = model.mean(ev,kick)
        diag = model.diagnostics(ev,mean['q'])
        x = s@x0+drift
        energy = x@h@x/2+b@x+np.einsum('ij,ij->',h,s@w@s.T)/2
        initial_energy = x0@h@x0/2+b@x0+np.einsum('ij,ij->',h,w)/2
        check(name+'_dense_energy',diag['energy_ry'],energy)
        check(name+'_energy_conservation',diag['energy_ry'],initial_energy)
        check(name+'_ordered_ccr',diag['ordered_ccr_relative_defect'],0.)
        check(name+'_paired_symplectic',diag['paired_symplectic_relative_defect'],0.)
    for mutation in ('omit_force','omit_cross_cell','omit_deltaC','omit_mixed','wrong_sqrtN','wrong_pZa_sign','omit_contact'):
        hh,jj,ww,bb = dense_reference(data,mutation)
        bad_s,bad_drift = dense_evolve(hh,jj,bb,time)
        sensitive(mutation,np.r_[bad_s@initial_local+bad_drift,(bad_s@ww@bad_s.T).reshape(-1)],
                  np.r_[s@initial_local+drift,(s@w@s.T).reshape(-1)])
    # Wrong q pairing, dagger, and 1/N-vs-orthonormal scaling are separate traps.
    wrong_pair,wrong_dagger = np.zeros((6,6),complex),np.zeros((6,6),complex)
    for k in range(3):
        phase = np.exp(2j*np.pi*k/3)/3
        wrong_pair += phase*(ev.transitions[k][:6]@model.w0[k]@other.transitions[k][:6].T)
        wrong_dagger += phase*(ev.transitions[k][:6]@model.w0[k]@other.transitions[model.minus[k]][:6].conj().T)
    expected_pair = (s@w@ss.T)[6:12,:6]
    sensitive('wrong_q_pairing',wrong_pair,expected_pair)
    sensitive('wrong_covariance_dagger',wrong_dagger,expected_pair)
    sensitive('wrong_mean_CCR_scaling',model.w,model.w/3)
    # No hidden reality extraction or symmetric Hessian repair, even on defective input.
    defective = dict(data)
    defective['delta_curvature'] = data['delta_curvature'].astype(complex).copy()
    defective['delta_curvature'][0,1] += .013j
    raw = ClampedArray(**defective)
    check('raw_complex_coefficient_retained',raw.h[0][0,1]-model.h[0][0,1],.013j)
    if raw.coefficient_diagnostics()['hermiticity_defect'] <= 0:
        raise AssertionError('raw Hermiticity defect hidden')
    out = BoundedOutputs(args.output,4*1024*1024)
    report = dict(schema='P254-clamped-array-oracle-v1',status='passed',source_bindings=source_bindings(),
        symbolic=symbolic,normalized_errors=errors,sensitivity=sensitivity,
        coefficient_diagnostics=model.coefficient_diagnostics(),
        covariance_normalization='norm(error)/max(1,norm(reference)), ordered covariance with q/-q transpose',
        independent_reference='SymPy differentiated real-cell scalar; SciPy expm of separately assembled full real-space Hamiltonian',
        actual_material_execution=False,whole_joint_dynamics_complete=False)
    record = out.json('result.json',report)
    print(json.dumps(dict(status='passed',result=record)))


if __name__ == '__main__':
    main()

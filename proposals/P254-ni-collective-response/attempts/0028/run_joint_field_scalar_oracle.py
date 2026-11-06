#!/usr/bin/env python3
"""Independent SymPy Hamiltonian/unit and dense complex-frame oracle.

Synthetic finite examples only, executed by Root; never material science output.
Production scalar/overlap routines are compared with explicit dense matrices,
independent radial quadrature values, and differentiated scalar expressions.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'),str(HERE.parent/'0025'),str(HERE.parent/'0024')]

import numpy as np
import sympy as sp
from scipy.special import spherical_jn
import run_joint_field_material as consumer
from substrate_framework import material_mixed_field_ionic as mixed
from substrate_framework import material_joint_field_scalar as joint
from substrate_framework.qe_material_field_derivatives import UPFRadial
from substrate_framework.qe_material_hamiltonian import MaterialHamiltonian
from substrate_framework.prepared_reference import PreparedReference
from substrate_framework.prepared_factor_dynamics import PreparedFactorDynamics
from substrate_framework.physical_field_environment import UniformMaxwellField


CHECKS = {}
MUTATIONS = {}


def close(name, actual, expected, tolerance=2e-10):
    actual,expected = np.asarray(actual),np.asarray(expected)
    error = float(np.max(np.abs(actual-expected)))
    scale = max(1.,float(np.max(np.abs(expected))))
    CHECKS[name] = dict(error=error,relative_scale=scale,tolerance=tolerance)
    if not np.isfinite(error) or error > tolerance*scale:
        raise AssertionError(name+': '+repr(CHECKS[name]))


def sensitive(name, wrong, right, minimum=1e-7):
    error = float(np.max(np.abs(np.asarray(wrong)-np.asarray(right))))
    MUTATIONS[name] = error
    if not np.isfinite(error) or error <= minimum:
        raise AssertionError('insensitive mutation: '+name)


def symbolic_oracle():
    N,u,p,a,pi,M,m,z,r,r0,K,h,d0,g,dg,F,C,dC,L,Ca,T = sp.symbols(
        'N u p a pi M m z r r0 K h d0 g dg F C dC L Ca T',real=True)
    scalar = N*((h+d0)*r+K*(r-r0)**2/2+u*((g+dg)*r-F)
                +(C+dC)*u*u/2+a*L*r+Ca*a*a/2+u*T*a+(p-z*a)**2/(2*m))+pi*pi/(2*M)
    expected = {'u':(p-z*a)/m,'p':F-(g+dg)*r-(C+dC)*u-T*a,
                'a':pi/M,'pi':-N*(L*r+Ca*a+T*u-z*(p-z*a)/m),
                'rho':h+d0+K*(r-r0)+u*(g+dg)+a*L}
    actual = {'u':sp.diff(scalar,p)/N,'p':-sp.diff(scalar,u)/N,
              'a':sp.diff(scalar,pi),'pi':-sp.diff(scalar,a),'rho':sp.diff(scalar,r)/N}
    for name in expected:
        assert sp.simplify(actual[name]-expected[name]) == 0,name
    conservation = sum(sp.diff(scalar,x)*actual[name] for x,name in [(u,'u'),(p,'p'),(a,'a'),(pi,'pi')])
    assert sp.simplify(conservation) == 0
    assert sp.diff(actual['rho'],a,2) == 0 and sp.diff(actual['rho'],u,a) == 0
    # Starting SI formulas, not an asserted Ry prefactor.
    e,a0,eps,hbar,me,Omega,adot,mi = sp.symbols('e a0 eps hbar me Omega adot mi',positive=True)
    Ry = e**2/(8*sp.pi*eps*a0)
    Esi = -Ry*adot/(e*a0)
    field_energy = eps*Omega*a0**3*Esi**2/(2*Ry)
    assert sp.simplify(field_energy-Omega*adot**2/(16*sp.pi)) == 0
    assert sp.simplify(Esi/(2*Ry/(e*a0))+adot/2) == 0
    # a0=4*pi*eps*hbar^2/(me*e^2) is substituted before mass reduction.
    ion_mass = mi*a0*a0*Ry/hbar**2
    assert sp.simplify(ion_mass.subs(a0,4*sp.pi*eps*hbar**2/(me*e**2))-mi/(2*me)) == 0
    electron_kinetic = (p+a)**2
    assert sp.simplify(sp.diff(electron_kinetic,a)-2*(p+a)) == 0
    electron_current = -sp.diff(electron_kinetic,a)/2
    ion_current = -sp.diff((p-z*a)**2/(2*m),a)/2
    assert sp.simplify(electron_current+p+a) == 0
    assert sp.simplify(ion_current-z*(p-z*a)/(2*m)) == 0
    CHECKS['sympy'] = dict(hamilton_derivatives=True,total_energy_conservation=True,
        omitted_cubic_rho_terms=True,si_maxwell_mass=True,ry_ion_mass=True,
        electron_negative_charge=True,ion_positive_charge=True,time_and_electric_factor_two=True)


def dense_mixed_oracle(rng):
    # An explicit retained radial table, with independent literal Simpson rule.
    radial_grid = np.linspace(.1,1.3,7)
    beta_r = np.column_stack([np.exp(-radial_grid),radial_grid*np.exp(-.7*radial_grid),
                              radial_grid**2*np.exp(-.4*radial_grid)])
    channel_l = np.array([0,1,1,1,2,2,2,2,2])
    channel_radial = channel_l.copy()
    upf = SimpleNamespace(r=radial_grid,rab=np.full(7,.2),beta_mesh=7,beta_r=beta_r,
        nbeta=3,angular_l=np.arange(3),channel_l=channel_l,channel_lm=np.arange(9),
        channel_radial=channel_radial,nh=9)
    tables = SimpleNamespace(upf=upf,omega_bohr3=19.)
    quadrature = np.array([1,4,2,4,2,4,1])*.2/3
    def centered_literal(p):
        x,y,z = p.T
        q = np.linalg.norm(p,axis=1)
        s = np.column_stack([np.ones(len(p))/np.sqrt(4*np.pi),
            np.sqrt(3/(4*np.pi))*z,-np.sqrt(3/(4*np.pi))*x,-np.sqrt(3/(4*np.pi))*y,
            np.sqrt(5/(4*np.pi))*(z*z-(x*x+y*y)/2),
            -np.sqrt(15/(4*np.pi))*x*z,-np.sqrt(15/(4*np.pi))*y*z,
            np.sqrt(15/(16*np.pi))*(x*x-y*y),np.sqrt(15/(4*np.pi))*x*y])
        value = np.empty((len(p),9),complex)
        for channel,ell in enumerate(channel_l):
            radial_value = spherical_jn(int(ell),q[:,None]*radial_grid) @ (
                quadrature*radial_grid*beta_r[:,ell])
            value[:,channel] = (-1j)**int(ell)*s[:,channel]/q**ell*radial_value*4*np.pi/np.sqrt(19.)
        return value
    p = rng.normal(size=(11,3))*.6+np.array([.4,.2,-.1])
    tau = np.array([.31,-.27,.19])
    frame = (rng.normal(size=(11,4))+1j*rng.normal(size=(11,4)))/4
    dc = rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
    dphysical = dc @ dc.conj().T/8
    dordered = dc/3
    active = np.diag(np.array([.9,-.4,-.4,-.4,.2,.2,.2,.2,.2]))
    phase = np.exp(-1j*p @ tau)[:,None]
    jet = mixed.first_cartesian_jets(UPFRadial(tables),p)
    h = 2e-5
    literal = centered_literal(p)
    close('literal radial/angular value',jet[0],literal)
    for j in range(3):
        e = np.eye(3)[j]*h
        derivative = (centered_literal(p-2*e)-8*centered_literal(p-e)
                      +8*centered_literal(p+e)-centered_literal(p+2*e))/(12*h)
        close('independent first UPF jet '+str(j),jet[j+1],derivative,2e-8)
    placed = [v*phase for v in jet]
    jets = mixed.mixed_projector_jets(placed[0],placed[1:],p)
    overlaps = np.asarray([b.conj().T @ frame for b in jets])
    def dense_mixed(dchannel):
        return np.asarray([[jets[7+3*i+j] @ dchannel @ jets[0].conj().T
            +jets[4+i] @ dchannel @ jets[1+j].conj().T
            +jets[1+j] @ dchannel @ jets[4+i].conj().T
            +jets[0] @ dchannel @ jets[7+3*i+j].conj().T
            for j in range(3)] for i in range(3)])
    general = rng.normal(size=(9,9))+1j*rng.normal(size=(9,9))
    general = (general+general.conj().T)/8
    for label,density in [('physical complex D',dphysical),('independent ordered D',dordered)]:
        rho = frame @ density @ frame.conj().T
        for channel_name,dchannel in [('scalar NC',active),('complex channel',general)]:
            operators = dense_mixed(dchannel)
            expected = np.einsum('ijab,ba->ij',operators,rho)
            observed = mixed.mixed_overlap_coefficient(overlaps,dchannel,density)
            close(label+' '+channel_name+' dense mixed trace',observed,expected)
            if label == 'independent ordered D' and channel_name == 'complex channel':
                sensitive('conjugated D',mixed.mixed_overlap_coefficient(overlaps,dchannel,density.conj()),expected)
                sensitive('transpose instead of adjoint overlaps',mixed.mixed_overlap_coefficient(
                    np.asarray([b.T @ frame for b in jets]),dchannel,density),expected)
                bad = operators.copy()
                bad[0,0] -= jets[7] @ dchannel @ jets[0].conj().T
                sensitive('omitted mixed product leg',np.einsum('ijab,ba->ij',bad,rho),expected)
    def beta(u,a,gauge=False):
        wave = p+a if gauge else p
        return centered_literal(p+a)*np.exp(-1j*wave @ (tau+u))[:,None]
    def v(u,a,gauge=False):
        b = beta(u,a,gauge)
        return b @ active @ b.conj().T
    zero = np.zeros(3)
    close('atom common Peierls phase cancels',v(np.array([.12,.09,-.06]),np.array([.17,-.13,.08])),
          v(np.array([.12,.09,-.06]),np.array([.17,-.13,.08]),True))
    dense = dense_mixed(active)
    step = 2e-4
    for i in range(3):
        for j in range(3):
            eu,ea = np.eye(3)[i]*step,np.eye(3)[j]*step
            finite = (v(eu,ea)-v(eu,-ea)-v(-eu,ea)+v(-eu,-ea))/(4*step*step)
            finite_gauge = (v(eu,ea,True)-v(eu,-ea,True)-v(-eu,ea,True)+v(-eu,-ea,True))/(4*step*step)
            close('independent mixed finite derivative '+str((i,j)),dense[i,j],finite,2e-6)
            close('phase dressed mixed finite derivative '+str((i,j)),dense[i,j],finite_gauge,2e-6)
    sensitive('ionic phase sign',-dense,dense)
    # Exercise the whole blocked contraction and physical moving-state vertices.
    # Synthetic object bypasses saved-material provenance only inside this oracle.
    legacy = placed[0]*(1+.04*np.sum(p*p,axis=1))[:,None]
    species = SimpleNamespace(name='fixture',atom_indices=(1,),active_d=active,beta=legacy)
    ham = object.__new__(MaterialHamiltonian)
    object.__setattr__(ham,'factors',SimpleNamespace(species=(species,),
        fft=SimpleNamespace(reciprocal_bohr_inv=np.eye(3),miller=p.T),k_bohr_inv=np.zeros(3)))
    object.__setattr__(ham,'background',SimpleNamespace(inputs=SimpleNamespace(
        reference=SimpleNamespace(positions_bohr=tau[:,None]),tables={'fixture':tables})))
    observed = mixed.full_column_mixed(ham,frame,dphysical,row_block=3)
    rho = frame @ dphysical @ frame.conj().T
    expected = np.einsum('ijab,ba->ij',dense,rho)
    close('all-column blocked mixed contraction',observed['mixed_ionic_field'][0],expected)
    def legacy_v(u):
        b = legacy*np.exp(-1j*p @ u)[:,None]
        return b @ active @ b.conj().T
    delta0 = v(zero,zero)-legacy_v(zero)
    delta_g = np.asarray([(v(h*e,zero)-legacy_v(h*e)-v(-h*e,zero)+legacy_v(-h*e))/(2*h)
                          for e in np.eye(3)])
    ell = np.asarray([np.diag(2*p[:,j])+(v(zero,h*e)-v(zero,-h*e))/(2*h)
                      for j,e in enumerate(np.eye(3))])
    u,a = np.array([.03,-.02,.04]),np.array([.05,-.07,.02])
    result = mixed.UPFPreparedVertices(ham).evaluate(frame,dphysical,u,a)
    close('unanchored zero-field energy',result['delta_energy'],np.trace(delta0 @ rho))
    close('physical delta ionic gradient',result['delta_gradient'],np.einsum('iab,ba->i',delta_g,rho),2e-7)
    close('physical kinetic plus UPF field gradient',result['field_gradient'],np.einsum('iab,ba->i',ell,rho),2e-7)
    close('full moving-state electronic vertex',result['action'],
          (delta0+np.einsum('i,iab->ab',u,delta_g)+np.einsum('i,iab->ab',a,ell)) @ frame,2e-7)
    sensitive('anchored zero-field correction omission',0,np.trace(delta0 @ rho))
    sensitive('kinetic optical omission',result['field_gradient']-np.einsum('iab,ba->i',
              np.asarray([np.diag(2*p[:,j]) for j in range(3)]),rho),result['field_gradient'])


def joint_dense_oracle(rng):
    nrow,ncol,ncells = 5,3,4
    def hermitian():
        x = rng.normal(size=(nrow,nrow))+1j*rng.normal(size=(nrow,nrow))
        return (x+x.conj().T)/7
    h0,delta0 = hermitian(),hermitian()
    gs,dgs,ls = [np.asarray([hermitian() for _ in range(3)]) for _ in range(3)]
    frames = tuple((rng.normal(size=(nrow,ncol))+1j*rng.normal(size=(nrow,ncol)))/5 for _ in range(2))
    densities = []
    for _ in frames:
        x = rng.normal(size=(ncol,ncol))+1j*rng.normal(size=(ncol,ncol))
        densities.append(x @ x.conj().T/8)
    class Coupling:
        npw,nr,ndof = nrow,1,3
        cell_volume_bohr3,grid_measure_bohr3 = 2.3,2.3
        def moments(self,left,right,coefficient=None):
            d = np.eye(left.shape[1]) if coefficient is None else coefficient
            rho = left @ d @ right.conj().T
            return np.array([np.trace(rho)]),np.einsum('iab,ba->i',gs,rho)
        def potential_action(self,value,psi):
            return value[0]*psi
        def ionic_action(self,u,psi):
            return np.einsum('i,iab->ab',u,gs) @ psi
    # density pairing has dv=2.3; choose n=Tr(rho)/dv for conjugate action.
    coupling = Coupling()
    original_moments = coupling.moments
    coupling.moments = lambda l,r,d=None: (original_moments(l,r,d)[0]/2.3,original_moments(l,r,d)[1])
    reference = PreparedReference(nbar=np.array([[.2]]),ionic_force=np.array([.3,-.4,.2]),
        held_density=np.zeros((1,1)),held_gradient=np.zeros(3),curvature=np.diag([1.2,1.5,1.8]),
        masses=np.array([8.,8.,8.]),support_force=np.array([-.1,.05,.07]),kernel=lambda x:.4*x,
        ncells=ncells,expansion_authority='synthetic exact scalar',ionic_authority='synthetic force',
        coverage_authority='two independent slots',support_authority='constant support',
        coordinate_identity='one synthetic atom',state_authority='complex full frames')
    sectors = tuple(SimpleNamespace(initial_frame=f,d=d,weight=1/ncells,coupling=coupling,
        density_component=0,h0=lambda x:h0 @ x) for f,d in zip(frames,densities))
    base = PreparedFactorDynamics(reference,sectors)
    class Vertex:
        def evaluate(self,psi,d,u,a):
            rho = psi @ d @ psi.conj().T
            return dict(action=(delta0+np.einsum('i,iab->ab',u,dgs)+np.einsum('i,iab->ab',a,ls)) @ psi,
                delta_energy=np.trace(delta0 @ rho),delta_gradient=np.einsum('iab,ba->i',dgs,rho),
                field_gradient=np.einsum('iab,ba->i',ls,rho))
    coeff = joint.JointFieldCoefficients(np.diag([.2,-.1,.3]),
        np.array([[1.4,.2,-.1],[.2,1.8,.3],[-.1,.3,1.6]]),rng.normal(size=(3,3))*.3,'synthetic fixed state')
    field = UniformMaxwellField(2.3*ncells,np.eye(3),'synthetic periodic volume')
    z = np.array([2.])
    model = joint.JointPreparedFieldDynamics(base,(Vertex(),Vertex()),coeff,field,z,scalar_authority='oracle scalar')
    state = model.mean(frames,np.array([.1,-.2,.15]),np.array([.3,-.25,.2]),
                       np.array([.07,-.06,.04]),np.array([.09,-.08,.11]))
    def literal_energy(frames,u,p,a,pi):
        rhos = [psi @ d @ psi.conj().T for psi,d in zip(frames,densities)]
        density = sum(np.trace(rho)/ncells for rho in rhos)/2.3
        potential_energy = sum(np.trace((h0+delta0+np.einsum('i,iab->ab',u,gs+dgs)
                                          +np.einsum('i,iab->ab',a,ls)) @ rho)/ncells for rho in rhos)
        scalar = (potential_energy+2.3*.4*(density-.2)**2/2
            -u @ (reference.ionic_force+reference.support_force)
            +u @ (reference.curvature+coeff.delta_curvature) @ u/2
            +a @ coeff.contact @ a/2+u @ coeff.mixed @ a
            +np.sum((p-2*a)**2/reference.masses)/2)
        return ncells*scalar+pi @ pi/(2*field.mass)
    result = model.evaluate(state)
    u,p,a,pi = state.matter.u,state.matter.p,state.coordinate,state.momentum
    close('literal total scalar',result['energy'],literal_energy(frames,u,p,a,pi))
    h = 1e-5
    args = [u,p,a,pi]
    gradients = []
    for index in range(4):
        gradient = []
        for axis in np.eye(3):
            plus,minus = list(args),list(args)
            plus[index] = args[index]+h*axis
            minus[index] = args[index]-h*axis
            gradient.append((literal_energy(frames,*plus)-literal_energy(frames,*minus))/(2*h))
        gradients.append(np.asarray(gradient))
    close('canonical ionic force',result['matter_rhs'].p,-gradients[0]/ncells,2e-8)
    close('positive ion minimal velocity',result['matter_rhs'].u,gradients[1]/ncells,2e-8)
    close('canonical Maxwell momentum force',result['field_rhs'][1],-gradients[2],2e-8)
    close('canonical Maxwell velocity',result['field_rhs'][0],gradients[3],2e-8)
    close('total charge current sign',result['electronic_current_hartree']+result['ionic_current_hartree'],
          -gradients[2]/2,2e-8)
    for index in range(2):
        variation = rng.normal(size=frames[index].shape)+1j*rng.normal(size=frames[index].shape)
        plus,minus = list(frames),list(frames)
        plus[index] = frames[index]+h*variation
        minus[index] = frames[index]-h*variation
        numerical = (literal_energy(plus,*args)-literal_energy(minus,*args))/(2*h)
        expected = 2*np.real(np.trace(densities[index] @ variation.conj().T @ result['h_frames'][index]))
        close('electronic functional derivative '+str(index),numerical,expected,2e-8)
    dt = 1e-5
    def flowed(sign):
        psis = [psi+sign*dt*rhs for psi,rhs in zip(frames,result['matter_rhs'].frames)]
        return literal_energy(psis,u+sign*dt*result['matter_rhs'].u.real,
            p+sign*dt*result['matter_rhs'].p.real,a+sign*dt*result['field_rhs'][0].real,
            pi+sign*dt*result['field_rhs'][1].real)
    close('same-scalar directional energy conservation',(flowed(1)-flowed(-1))/(2*dt),0,2e-7)
    sensitive('wrong electron charge current sign',-result['electronic_current_hartree'],result['electronic_current_hartree'])
    sensitive('wrong ion minimal coupling sign',(p+2*a)/reference.masses,result['matter_rhs'].u)
    sensitive('missing N field normalization',result['field_rhs'][1]/ncells,result['field_rhs'][1])
    sensitive('fixed contact omitted',result['field_rhs'][1]+ncells*coeff.contact @ a,result['field_rhs'][1])
    sensitive('fixed mixed force omitted',result['matter_rhs'].p+coeff.mixed @ a,result['matter_rhs'].p)
    # No field-contact/mixed operator is inserted in rho derivative.
    wrong_h = result['h_frames'][0]+(a @ coeff.contact @ a/2+u @ coeff.mixed @ a)*frames[0]
    sensitive('spurious rho derivative of fixed contact or mixed coefficient',wrong_h,result['h_frames'][0])
    prepared = joint.conditional_field_gaussian(field,coeff.contact,reference.masses,z,ncells,state_authority='oracle product')
    close('Gaussian canonical uncertainty',prepared['covariance_qq'] @ prepared['covariance_pp'],np.eye(3)/4)
    close('Gaussian ordered qp sign',prepared['ordered_covariance'][:3,3:],.5j*np.eye(3))
    literal_clamped = ncells*(coeff.contact+np.diag(4/reference.masses))
    close('canonical ion diamagnetic stiffness retained',prepared['clamped_stiffness_raw'],literal_clamped)
    close('Gaussian same-scalar quadratic energy',prepared['quadratic_zero_point_energy_ry'],
          np.trace(literal_clamped @ prepared['covariance_qq'])/2
          +np.trace(prepared['covariance_pp'])/(2*field.mass))
    sensitive('omitted canonical ion diamagnetism',ncells*coeff.contact,prepared['clamped_stiffness_raw'])
    eig = np.linalg.eigvalsh(prepared['ordered_covariance'])
    if eig.min() < -1e-10:
        raise AssertionError('ordered Gaussian covariance not positive')
    rejected = False
    try:
        joint.conditional_field_gaussian(field,-10*np.eye(3),reference.masses,z,ncells,state_authority='negative fixture')
    except joint.MissingFieldPreparation:
        rejected = True
    assert rejected,'negative clamped stiffness must never be clipped'
    CHECKS['nonpositive preparation rejection'] = True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    registration,_ = consumer.authenticate_sources()
    implementations = consumer.implementation_bindings()
    rng = np.random.default_rng(254225)
    symbolic_oracle()
    dense_mixed_oracle(rng)
    joint_dense_oracle(rng)
    result = dict(schema='P254-joint-field-scalar-oracle-v1',status='passed',
        registration=registration,implementation_sources=implementations,checks=CHECKS,
        mutation_sensitivities=MUTATIONS,material_executed=False,
        scope='independent finite symbolic/dense oracle; no physical material coefficients or propagated dynamics')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(status='passed',output=str(args.output),checks=len(CHECKS),mutations=len(MUTATIONS))))


if __name__ == '__main__':
    main()

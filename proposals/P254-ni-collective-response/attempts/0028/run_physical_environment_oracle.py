#!/usr/bin/env python3
"""Independent exact algebra/Fock-space and physical-action oracles; Root executes.

No material evolution. Synthetic fixtures are explicitly finite mathematical
operators, not fictitious material checkpoints. All mutations reach public
production consumers. A report is emitted only after every assertion passes.
"""
from pathlib import Path
from types import SimpleNamespace
from dataclasses import replace
import argparse
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import sympy as sp
from scipy.linalg import expm
from substrate_framework.physical_field_environment import (
    UniformMaxwellField, UPFUniformCoupling, charged_ion_kinematics, joint_canonical_rhs)
from substrate_framework.nonstationary_force_kernel import (
    NonstationaryForceKernel, HeisenbergForceFactors)
from substrate_framework.physical_factor_preparation import (
    FactorFunctional, FactorPreparation, ActionFactors)
from substrate_framework.qe_material_field_derivatives import MaterialFieldDerivatives, UPFRadial
from run_field_derivative_oracle import fixture, literal_integral


def close(actual, expected, tolerance=2e-11):
    if not np.allclose(actual, expected, rtol=tolerance, atol=tolerance):
        raise AssertionError(('independent oracle mismatch', actual, expected))


def rejected(call):
    try:
        call()
    except ValueError:
        return True
    raise AssertionError('consumer accepted corrupted input')


def canonical_blocks():
    # Derive coefficient from SI, not from production mass formula.
    eps, vol, a0, e, hb, ry, adot = sp.symbols('eps V a0 e hbar Ry adot', positive=True)
    electric_si = -ry/(e*a0)*adot
    energy = eps*vol*a0**3*electric_si**2/(2*ry)
    energy = sp.simplify(energy.subs(ry, e**2/(8*sp.pi*eps*a0)))
    assert energy == vol*adot**2/(16*sp.pi)
    mass = sp.diff(energy, adot, 2)
    a, p, pe, pi, z, mi, k, r, g = sp.symbols('a p pe Pi Z M k R g', real=True)
    he = (pe+a)**2+g*r
    hi = (p-z*a)**2/(2*mi)+k*r*r/2
    hf = pi*pi/(2*mass)
    scalar = he+hi+hf
    av = sp.diff(scalar, pi)
    pv = -sp.diff(scalar, a)
    rv = sp.diff(scalar, p)
    pdot = -sp.diff(scalar, r)
    assert sp.simplify(sp.diff(scalar,a)*av+sp.diff(scalar,pi)*pv
                       +sp.diff(scalar,r)*rv+sp.diff(scalar,p)*pdot) == 0
    assert sp.diff(he,a,2) == 2
    assert sp.diff(he,a)/(-2) == -pe-a
    assert sp.simplify(sp.diff(hi,a)+z*rv) == 0
    assert sp.simplify(sp.diff(hi,a,2)-z*z/mi) == 0
    # Exact constant-offset Ward obstruction: H(a)=G(a)H G†(a)+Delta.
    # Its derivative is i[X,H], NOT i[X,H+Delta] unless [X,Delta]=0.
    x = sp.diag(0,1)
    h = sp.Matrix([[1,2],[2,3]])
    delta = sp.Matrix([[0,1],[1,0]])
    comm = lambda left,right: left*right-right*left
    ward_difference = sp.I*comm(x,delta)
    assert ward_difference != sp.zeros(2)
    assert -comm(x,comm(x,h+delta)) != -comm(x,comm(x,h))
    # Finite matrices cannot satisfy canonical [x,p]=iI: trace obstruction.
    assert sp.trace(sp.I*sp.eye(2)) != 0
    field = UniformMaxwellField.from_cell(np.diag([2.,3.,4.]),
        polarizations=np.eye(3),cell_authority='synthetic exact 24 bohr^3 cell')
    close(field.mass, float(mass.subs(vol,24)))
    momentum = np.array([.2,-.3,.4]); coordinate = np.array([.1,.2,-.1])
    ge = np.array([.6,-.2,.1]); pions = np.array([[.2,.7,-.3],[.4,-.2,.1]])
    masses = np.array([3.,7.]); charges = np.array([1.,2.])
    force = np.array([[.1,.3,-.2],[-.3,.2,.4]])
    out = joint_canonical_rhs(field,coordinate,momentum,electronic_gradient=ge,
        ion_momentum=pions,masses_ry=masses,charges_e=charges,
        conservative_ionic_force=force,scalar_authority='synthetic conservative scalar')
    v = (pions-charges[:,None]*coordinate)/masses[:,None]
    gi = -np.sum(charges[:,None]*v,axis=0)
    close(out['field_coordinate'],momentum/(24/(8*np.pi)))
    close(out['field_momentum'],-ge-gi)
    close(out['ion_position'],v); close(out['ion_momentum'],force)
    close(out['field_power']+(ge+gi)@out['field_coordinate'],0)
    close(field.electric_field_hartree(momentum),-out['field_coordinate']/2)
    ions = charged_ion_kinematics(pions,masses,charges,coordinate)
    close(ions['field_contact'],np.eye(3)*sum(charges*charges/masses))
    mutations = {'Hartree_instead_of_Ry_mass':float(abs(field.mass-24/(4*np.pi))),
                 'reversed_field_current_sign':float(np.linalg.norm(2*out['field_momentum'])),
                 'missing_ionic_minimal_shift':float(np.linalg.norm(v-pions/masses[:,None])),
                 'missing_electric_factor_half':float(np.linalg.norm(field.electric_field_hartree(momentum)+out['field_coordinate']))}
    assert min(mutations.values())>1e-4
    return dict(SI_to_Ry=str(energy),field_mass=str(mass),joint_energy_identity=True,
                undressed_offset_Ward_obstruction=str(ward_difference),mutations=mutations)


def physical_actions():
    # The independent angular integral differentiates exp(-i r p.n), not the
    # production solid-harmonic/Bessel derivative implementation.
    tables=fixture(); radial=UPFRadial(tables)
    p=np.array([[0.,0.,0.],[.2,-.3,.4],[1.1,.7,-.6],[-.31,0.,.5]])
    shift=np.array([.12,-.07,.04]); left=np.array([.7,-.2,.4]); right=np.array([-.3,.8,.1])
    position=np.array([.13,-.21,.09]); phase=np.exp(-1j*(p@position))[:,None]
    b0=literal_integral(tables,p,left,right)[0]*phase
    dim=b0.shape[1]
    active=np.diag(np.linspace(.4,1.2,dim)); active[1,9]=active[9,1]=.23
    old=b0+(.03+.01j)*np.ones_like(b0)
    local=np.diag([.1,.2,-.1,.3])
    h0=np.diag(np.sum(p*p,axis=1))+local+old@active@old.conj().T
    species=SimpleNamespace(name='synthetic',active_d=active,atom_indices=[1],beta=old)
    ham=SimpleNamespace(block_size=2,factors=SimpleNamespace(species=[species]),
        background=SimpleNamespace(inputs=SimpleNamespace(reference=SimpleNamespace(
            positions_bohr=position[:,None]))),apply=lambda v:h0@v)
    def vectors(v):
        v=np.asarray(v); return (v[:,None],True) if v.ndim==1 else (v,False)
    # Assemble a finite mathematical fixture of the existing physical family;
    # no accepted material/background object or checkpoint is claimed here.
    derivatives=object.__new__(MaterialFieldDerivatives)
    derivatives.hamiltonian=ham; derivatives.momenta=p
    derivatives.original=SimpleNamespace(_vectors=vectors)
    derivatives.radial={'synthetic':radial}
    model=UPFUniformCoupling(derivatives)
    b,be,bf,bef=(v*phase for v in literal_integral(tables,p+shift,left,right))
    dense_h=np.diag(np.sum((p+shift)**2,axis=1))+local+b@active@b.conj().T
    dense_l=np.diag(2*((p+shift)@left))+be@active@b.conj().T+b@active@be.conj().T
    dense_c=(2*(left@right)*np.eye(4)+bef@active@b.conj().T
        +be@active@bf.conj().T+bf@active@be.conj().T+b@active@bef.conj().T)
    v=np.array([[1,.2j],[.3+.2j,.6],[.4,-.2j],[.1j,.5]],complex)
    got_l,got_c=model.derivative_actions(v,shift,left,right)
    close(model.action(v,shift),dense_h@v)
    close(got_l,dense_l@v); close(got_c,dense_c@v)
    close(model.current_action(v,shift,left),-dense_l@v/2)
    d=np.array([[.2,.03j],[-.03j,.15]])
    moments=model.moments(v,d,shift,block_size=1,row_block_size=2)
    close(moments['energy_ry'],np.trace(d@v.conj().T@dense_h@v))
    close(left@moments['gradient'],np.trace(d@v.conj().T@dense_l@v))
    close(left@moments['contact']@right,np.trace(d@v.conj().T@dense_c@v))
    # The A25 offset is a different H(0), although current/contact agree.
    anchored=derivatives.field_action(v,shift,anchored=True)
    mutations={'anchored_A25_mislabeled_physical':float(np.linalg.norm(anchored-model.action(v,shift))),
               'omit_contact_cross_legs':float(np.linalg.norm((be@active@bf.conj().T+bf@active@be.conj().T)@v)),
               'reverse_current_sign':float(np.linalg.norm(2*model.current_action(v,shift,left)))}
    assert min(mutations.values())>1e-6
    return dict(physical_UPF_action_current_contact=True,unanchored_zero_field_distinct=True,
                independent_literal_spherical_integral=True,mutations=mutations)


def fock_blocks():
    psi=sp.Matrix([[1,0],[1,1],[0,1]])
    d=sp.Matrix([[2,1],[1,1]])/10
    rho=psi*d*psi.H
    # Independent exterior-power density on the full 8-dimensional Fock space.
    # Gamma(R)_{I,J}=det R[I,J], rho_F=det(1-rho) Gamma(rho/(1-rho)).
    ratio=rho*(sp.eye(3)-rho).inv(); pref=(sp.eye(3)-rho).det()
    occupied=[[i for i in range(3) if state&(1<<i)] for state in range(8)]
    rf=sp.zeros(8)
    for i,ii in enumerate(occupied):
        for j,jj in enumerate(occupied):
            if len(ii)==len(jj):
                rf[i,j]=pref*(ratio.extract(ii,jj).det() if ii else 1)
    assert sp.trace(rf)==1
    annih=[]
    for mode in range(3):
        c=sp.zeros(8)
        for state in range(8):
            if state&(1<<mode):
                c[state^(1<<mode),state]=(-1)**sum(bool(state&(1<<i)) for i in range(mode))
        annih.append(c)
    for i in range(3):
        for j in range(3):
            assert annih[i]*annih[j].H+annih[j].H*annih[i]==(sp.eye(8) if i==j else sp.zeros(8))
    def second(a):
        return sum((a[i,j]*annih[i].H*annih[j] for i in range(3) for j in range(3)),sp.zeros(8))
    h=sp.diag(0,1,2)
    a=sp.Matrix([[1,1,sp.I],[1,0,2],[-sp.I,2,-1]])
    b=sp.Matrix([[0,sp.I,1],[-sp.I,2,1],[1,1,1]])
    assert h*rho-rho*h != sp.zeros(3)
    def heis(o,t):
        u=sp.diag(1,sp.exp(-sp.I*t),sp.exp(-2*sp.I*t))
        return sp.simplify(u.H*o*u)
    exact=[]
    for t,s in [(sp.pi/2,0),(sp.pi,sp.pi/2)]:
        aa,bb=heis(a,t),heis(b,s)
        fa,fb=second(aa),second(bb)
        ma,mb=sp.trace(rf*fa),sp.trace(rf*fb)
        cab=sp.simplify(sp.trace(rf*fa*fb)-ma*mb)
        cba=sp.simplify(sp.trace(rf*fb*fa)-mb*ma)
        assert sp.simplify(cab-sp.trace(rho*aa*(sp.eye(3)-rho)*bb))==0
        assert sp.simplify(cab-cba-sp.trace(rho*(aa*bb-bb*aa)))==0
        exact.append((cab,cba))
    assert exact[0][0]!=exact[1][0]  # same time separation, different absolute times
    npy=lambda x:np.array(x.tolist(),dtype=complex)
    pn,dn,an,bn,hn,rfn=map(npy,(psi,d,a,b,h,rf))
    ops={'A':an,'B':bn,'number':np.eye(3,dtype=complex)}
    probes={name:ActionFactors(op@pn,op@pn,name,name,'Ry/bohr','fixture-basis','exact synthetic Fock oracle')
            for name,op in ops.items()}
    functional=FactorFunctional.from_actions(pn,probes,factor_identity='fixture-frame',pw_basis_identity='fixture-basis')
    prep=FactorPreparation(functional,dn,'explicit nonstationary rational quasifree fixture')
    kernel=NonstationaryForceKernel(prep,row_block_size=1)
    def endpoint(op,t,*,returned=True):
        u=expm(-1j*hn*t)
        act=(u.conj().T@op@u@pn) if returned else op@u@pn
        adj=(u.conj().T@op.conj().T@u@pn) if returned else op.conj().T@u@pn
        return HeisenbergForceFactors(t,act,adj,'fixture force','Ry/bohr',
            'fixture-frame','fixture-basis','same exact diagonal quadratic history')
    first=kernel.pair(pn,endpoint(an,np.pi/2),endpoint(bn,0.))
    close(first.ordered_ab,complex(exact[0][0])); close(first.ordered_ba,complex(exact[0][1]))
    close(first.retarded(perturbation_sign=1),-1j*complex(exact[0][0]-exact[0][1]))
    close(first.retarded(perturbation_sign=-1),1j*complex(exact[0][0]-exact[0][1]))
    early=kernel.pair(pn,endpoint(an,0.),endpoint(bn,1.))
    close(early.retarded(perturbation_sign=1),0.)
    second_pair=kernel.pair(pn,endpoint(an,np.pi),endpoint(bn,np.pi/2))
    close(second_pair.ordered_ab,complex(exact[1][0]))
    # SciPy independent full-Fock expm, arbitrary noncommensurate times.
    hf=npy(second(h)); af=npy(second(a)); bf=npy(second(b))
    t,s=.71,.19
    ut,us=expm(-1j*hf*t),expm(-1j*hf*s)
    fa,fb=ut.conj().T@af@ut,us.conj().T@bf@us
    dense=np.trace(rfn@fa@fb)-np.trace(rfn@fa)*np.trace(rfn@fb)
    pair=kernel.pair(pn,endpoint(an,t),endpoint(bn,s))
    close(pair.ordered_ab,dense)
    # Non-Hermitian probes require independent adjoint actions as well.
    nonhermitian=an+1j*bn
    nh=kernel.pair(pn,endpoint(nonhermitian,t),endpoint(bn,s))
    fnh=ut.conj().T@(af+1j*bf)@ut
    close(nh.ordered_ab,np.trace(rfn@fnh@fb)-np.trace(rfn@fnh)*np.trace(rfn@fb))
    endpoints=[endpoint(op,time) for op,time in [(an,0.),(bn,.23),(an,.71),(bn,1.3)]]
    covariance=np.array([[kernel.pair(pn,l,r).ordered_ab for r in endpoints] for l in endpoints])
    close(covariance,covariance.conj().T)
    assert np.linalg.eigvalsh(covariance).min()>-2e-12
    recurrent=kernel.pair(pn,endpoint(an,t+2*np.pi),endpoint(bn,s+2*np.pi))
    close(recurrent.ordered_ab,pair.ordered_ab)
    initial=kernel.pair(pn,endpoint(an,0.),endpoint(an,0.))
    later=kernel.pair(pn,endpoint(an,np.pi/2),endpoint(an,0.))
    return_pair=kernel.pair(pn,endpoint(an,2*np.pi),endpoint(an,0.))
    close(return_pair.ordered_ab,initial.ordered_ab)
    assert abs(later.ordered_ab-initial.ordered_ab)>1e-3
    # Exercise actual initial-checkpoint consumer with a clearly labeled fixture.
    arrays,metadata=functional.archive_payload(); names=tuple(ops)
    means=np.array([prep.mean(n) for n in names])
    c=np.array([[prep.ordered(l,r) for r in names] for l in names])
    arrays.update(density_coefficients=dn,means=means,ordered_connected=c,
        reverse_ordered_connected=c.T,symmetrized_connected=(c+c.T)/2,commutator=c-c.T,
        ordered_second_moment=c+means[:,None]*means[None,:],
        reverse_ordered_second_moment=(c+means[:,None]*means[None,:]).T)
    metadata.update(schema='material-equal-time-force-noise-v1',all_raw_columns=True,
        full_physical_pw_rows=True,source_slot_measure_applied=False,
        state_authority=prep.state_authority,operator_order=list(names))
    consume=lambda aa:NonstationaryForceKernel.from_initial_checkpoint(aa,metadata,
        checkpoint_authority='SYNTHETIC oracle only; not a material archive')
    restored=consume(arrays)
    close(restored.initial_pair('A','B').ordered_ab,c[0,1])
    mutations={}
    for name,key,replacement in [('transpose_AB','ordered_connected',c.T),
            ('erase_commutator','commutator',np.zeros_like(c)),
            ('symmetrize_quantum_noise','ordered_connected',(c+c.T)/2),
            ('replace_D_by_diagonal','density_coefficients',np.diag(np.diag(dn)))]:
        changed=dict(arrays); changed[key]=replacement
        mutations[name]=rejected(lambda:consume(changed))
    projected=dict(arrays)
    for item in metadata['array_layout']['ordered_products']:
        ml=functional.overlaps[item['left']]; mr=functional.overlaps[item['right']]
        projected[item['array']]=ml@np.linalg.solve(functional.gram,mr)
    mutations['discard_empty_complement']=rejected(lambda:consume(projected))
    bad=kernel.pair(pn,endpoint(an,t,returned=False),endpoint(bn,s,returned=False))
    mutations['omit_return_propagator_exposure']=float(abs(bad.ordered_ab-pair.ordered_ab))
    mutations['stationarity_exposure']=float(abs(first.ordered_ab-second_pair.ordered_ab))
    mutations['reverse_retarded_sign_exposure']=float(abs(2*first.retarded(perturbation_sign=1)))
    assert min(mutations[key] for key in mutations if key.endswith('exposure'))>1e-4
    mutations['truncate_raw_column']=rejected(lambda:kernel.pair(pn,
        replace(endpoint(an,t),action=endpoint(an,t).action[:,:1]),endpoint(bn,s)))
    return dict(exact_Fock_pairs=[[str(v) for v in item] for item in exact],
        CAR=True,nonstationary=True,unequal_time_PSD=True,finite_recurrence=True,
        independent_SciPy_Fock_expm=True,non_Hermitian_adjoint=True,mutations=mutations)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    # Source authenticity precedes every science block; implementation receipts
    # bind this script and both new modules so mutations remain consumer-visible.
    registration=json.loads((HERE/'physical-environment-registration.json').read_text())
    for item in registration['source_bindings']:
        if hashlib.sha256(Path(item['path']).read_bytes()).hexdigest()!=item['sha256']:
            raise ValueError('registered source bytes changed: '+item['path'])
    result=dict(schema='P254-physical-environment-oracle-v1',canonical=canonical_blocks(),
        physical_actions=physical_actions(),ordered_kernel=fock_blocks(),
        scope='independent synthetic checks only; no material two-time evolution')
    sources=[Path(__file__),HERE/'physical-environment-registration.json',
        ROOT/'src/substrate_framework/physical_field_environment.py',
        ROOT/'src/substrate_framework/nonstationary_force_kernel.py',
        HERE/'run_field_derivative_oracle.py']
    result['implementation_sources']=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources]
    text=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if args.output:
        with args.output.open('x') as stream: stream.write(text)
    print(text,end='')


if __name__=='__main__':
    main()

#!/usr/bin/env python3
"""NEW finite mathematical oracle; no Ni/H state or finite-time claim.

Run only by Root after source intake: --result PATH. The first block uses
SymPy on an independently written scalar. Subsequent comparators use original
small density matrices, literal occupation-bit CAR, and full-line canonical
Gaussian moments. Production never receives these dense comparator matrices
except through the explicitly finite action adapter. No old oracle is run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import traceback

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / 'src'))
CHECKS = []
FAILURES = []


def exact(name, difference):
    import sympy as sp
    ok = sp.simplify(difference) == 0
    CHECKS.append(dict(name=name, passed=bool(ok), exact=True))
    if not ok:
        raise AssertionError(name + ': ' + str(difference))


def compare(name, actual, expected):
    import numpy as np
    a, b = np.asarray(actual), np.asarray(expected)
    compatible = a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
    residual = float(np.max(np.abs(a - b), initial=0)) if compatible else None
    ok = bool(compatible and np.allclose(a, b, rtol=1e-10, atol=1e-12))
    CHECKS.append(dict(name=name, passed=ok, max_abs_residual=residual,
                       rtol=1e-10, atol=1e-12))
    if not ok:
        raise AssertionError(name + ': residual=' + str(residual))


def mutation(name, wrong, expected):
    import numpy as np
    a, b = np.asarray(wrong), np.asarray(expected)
    compatible = a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
    residual = float(np.max(np.abs(a - b), initial=0)) if compatible else None
    ok = bool(compatible and not np.allclose(a, b, rtol=1e-10, atol=1e-12))
    CHECKS.append(dict(name='mutation_' + name, passed=bool(ok), detected=bool(ok),
                       max_abs_residual=residual))
    if not ok:
        raise AssertionError('mutation not detected: ' + name)


def symbolic_scalar_first():
    import sympy as s
    n1, n2, b1, b2, h1, h2 = s.symbols('n1 n2 b1 b2 h1 h2', real=True)
    u, v, p, q, m, mass = s.symbols('u v p q m mass', real=True)
    k, l, mix, c, d, off = s.symbols('k l mix c d off', real=True)
    g, j, f, z, support, t = s.symbols('g j f z support t', real=True)
    gh, jh, dv = s.symbols('gh jh dv', real=True)
    x, y = n1 + h1 - b1, n2 + h2 - b2
    # One scalar, written without importing any production derivative.
    energy = (dv * (k*x*x + 2*mix*x*y + l*y*y)/2
              + u*(g+gh-f-support) + v*(j+jh-z-t)
              + (c*u*u+2*off*u*v+d*v*v)/2 + p*p/(2*m)+q*q/(2*mass))
    exact('scalar_mixed_spin_density_1', s.diff(energy, n1) - dv*(k*x+mix*y))
    exact('scalar_mixed_spin_density_2', s.diff(energy, n2) - dv*(mix*x+l*y))
    exact('scalar_nonzero_ionic_baseline', -s.diff(energy, u) - (f+support-g-gh-c*u-off*v))
    exact('scalar_second_ionic_coordinate', -s.diff(energy, v) - (z+t-j-jh-off*u-d*v))
    exact('scalar_distinct_mass_1', s.diff(energy, p)-p/m)
    exact('scalar_distinct_mass_2', s.diff(energy, q)-q/mass)
    exact('scalar_mixed_density_hessian', s.diff(energy, n1, n2)-dv*mix)
    exact('scalar_fixed_curvature', s.diff(energy, u, v)-off)
    a, b = s.symbols('origin_a origin_b', real=True)
    shifted = energy.subs({u:u+a, v:v+b}, simultaneous=True)
    expected = (energy + a*(g+gh-f-support) + b*(j+jh-z-t)
                + u*(c*a+off*b)+v*(off*a+d*b)+(c*a*a+2*off*a*b+d*b*b)/2)
    exact('scalar_coordinate_origin_identity', shifted-expected)


def finite_blocks():
    import numpy as np
    from dataclasses import replace
    from substrate_framework.physical_factor_preparation import ActionFactors, FactorFunctional, FactorPreparation
    from substrate_framework.physical_bilinear_preparation import IonicPreparation
    from substrate_framework.prepared_reference import PreparedReference, PreparedSector
    from substrate_framework.prepared_factor_dynamics import PreparedFactorDynamics
    from substrate_framework.prepared_covariance_source import (
        PreparedProbe, PreparedCovariance, FactorCrossPreparation, readout,
        probe_mean, ordered_readout, total_force_probe)

    rng = np.random.default_rng(22526)
    size, rank, ndof, nr, ncells = 3, 2, 2, 2, 4
    weight, dv = 1/ncells, 1.7

    def hermitian(scale=1):
        z = rng.normal(size=(size, size)) + 1j*rng.normal(size=(size, size))
        return scale*(z+z.conj().T)

    # Actual tiny operators, with nonzero source-to-empty-complement entries.
    hs = (hermitian(.3), hermitian(.2))
    gs = tuple(np.stack([hermitian(.17), hermitian(.13)]) for _ in range(2))
    density_ops = []
    for spin in range(2):
        vec = rng.normal(size=(nr, size)) + 1j*rng.normal(size=(nr, size))
        density_ops.append(np.einsum('ri,rj->rij', vec, vec.conj())/9)
    ks = np.array([[.7,.13,-.11,.08],[.13,.9,.07,-.12],
                   [-.11,.07,.6,.19],[.08,-.12,.19,.8]])
    curvature = np.array([[1.3,.21],[.21,.85]])
    masses = np.array([2.7,5.1])
    qframe = np.eye(size, dtype=complex)[:, :rank]
    transforms = (np.array([[1.2,.17+.09j],[.04j,.8]]),
                  np.array([[.9,.12j],[.19,1.1]]))
    occupations = (np.array([.35,.35,0.]), np.array([.6,.2,0.]))
    rhos = tuple(np.diag(f).astype(complex) for f in occupations)

    class FiniteCoupling:
        npw, nr, ndof = size, 2, 2
        grid_measure_bohr3, cell_volume_bohr3 = dv, dv*nr

        def __init__(self, spin):
            self.spin = spin

        def moments(self, left, right, coefficient=None):
            rho = left @ right.conj().T if coefficient is None else left @ coefficient @ right.conj().T
            return (np.einsum('rij,ji->r', density_ops[self.spin], rho),
                    np.einsum('aij,ji->a', gs[self.spin], rho))

        def potential_action(self, field, vectors):
            return dv*np.einsum('r,rij->ij', field, density_ops[self.spin]) @ vectors

        def ionic_action(self, direction, vectors):
            return np.einsum('a,aij->ij', direction, gs[self.spin]) @ vectors

    def sector(spin, transform):
        psi = qframe @ transform
        inv = np.linalg.inv(transform)
        dcoef = inv @ np.diag(occupations[spin][:rank]) @ inv.conj().T
        label = 'synthetic-finite-spin-' + str(spin)
        identity = ActionFactors(psi, psi, 'I', 'I', '1', label, 'finite mathematical identity')
        functional = FactorFunctional.from_actions(psi, {'I':identity}, factor_identity=label+'-raw', pw_basis_identity=label)
        prep = FactorPreparation(functional, dcoef, 'synthetic independent Bernoulli modes; empty complement, not Ni/H')
        return PreparedSector(prep, psi, FiniteCoupling(spin), lambda x, h=hs[spin]:h@x,
                              spin, weight, label, 'finite Hermitian operators in declared original basis')

    sectors = tuple(sector(i, t) for i,t in enumerate(transforms))
    reference = PreparedReference(
        nbar=np.array([[.31,.23],[.18,.36]]), ionic_force=np.array([.19,-.13]),
        held_density=np.array([[.04,.02],[.03,.06]]),
        held_gradient=np.array([-.03,.05]), curvature=curvature, masses=masses,
        support_force=np.array([.02,.01]), kernel=lambda n:(ks@n.ravel()).reshape(nr,2), ncells=ncells,
        expansion_authority='explicit synthetic joint-second-order scalar, fixed coefficients',
        ionic_authority='synthetic linear ionic scalar; no Native force interpretation',
        coverage_authority='two independent active sectors, explicit held remainder, N=4 convention',
        support_authority='declared synthetic constant conservative force, not a centering counterterm',
        coordinate_identity='two complete synthetic coordinates in Ry convention',
        state_authority='synthetic quasifree independent spins, specified mean number, empty initial complement')
    model = PreparedFactorDynamics(reference, sectors)
    state = model.initial_mean(np.array([.11,-.07]), np.array([.08,-.13]))

    # Direct original-density-matrix defining scalar and its derivatives.
    def original(rho, u, p, h0=hs, ref=reference):
        density = np.column_stack([weight*np.einsum('rij,ji->r', a, x)
                                   for a,x in zip(density_ops,rho)])
        grad = sum(weight*np.einsum('aij,ji->a', g,x) for g,x in zip(gs,rho))
        delta = density + ref.held_density - ref.nbar
        field = (ks@delta.ravel()).reshape(nr,2)
        h = tuple(hbase+dv*np.einsum('r,rij->ij', field[:,i], density_ops[i])
                  + np.einsum('a,aij->ij',u,gs[i]) for i,hbase in enumerate(h0))
        force = ref.ionic_force + ref.support_force-ref.held_gradient-grad-curvature@u
        drho = tuple(-1j*(a@x-x@a) for a,x in zip(h,rho))
        energy = (sum(weight*np.trace(a@x) for a,x in zip(h0,rho))
                  +dv*delta.ravel()@ks@delta.ravel()/2
                  +u@(grad+ref.held_gradient-ref.ionic_force-ref.support_force)
                  +u@curvature@u/2+np.sum(p*p/masses)/2)
        return h, drho, force, energy, density, grad

    ev = model.evaluate_mean(state)
    h, drho, force, energy, density, grad = original(rhos,state.u,state.p)
    compare('mean_nonorthogonal_density', ev.active_density,density)
    compare('mean_all_coordinate_gradient', ev.active_gradient,grad)
    compare('mean_scalar_energy', ev.energy,energy)
    compare('mean_ionic_acceleration', ev.acceleration,force/masses)
    compare('mean_velocity_distinct_masses', ev.rhs.u,state.p/masses)
    for i,s in enumerate(sectors):
        psi, dot = state.frames[i],ev.rhs.frames[i]
        compare('mean_full_space_commutator_'+str(i),dot@s.d@psi.conj().T+psi@s.d@dot.conj().T,drho[i])
        mutation('missing_complement_'+str(i), np.pad(drho[i][:2,:2],((0,1),(0,1))),drho[i])
    mutation('missing_ionic_baseline_and_support', -grad-curvature@state.u,force)
    mutation('missing_density_offset', density,ev.delta_density)
    mutation('force_sign', grad+reference.ionic_force+reference.support_force-reference.held_gradient-curvature@state.u,force)

    # Coordinate origin shift changes scalar constant, never physical drift.
    shift = np.array([.08,-.03])
    shifted_h0 = tuple(a+np.einsum('j,jab->ab',shift,g) for a,g in zip(hs,gs))
    shifted_sectors = tuple(replace(s,h0=lambda x,a=a:a@x) for s,a in zip(sectors,shifted_h0))
    shifted_ref = replace(reference,ionic_force=reference.ionic_force-curvature@shift)
    shifted_model = PreparedFactorDynamics(shifted_ref,shifted_sectors)
    shifted_ev = shifted_model.evaluate_mean(shifted_model.initial_mean(state.u-shift,state.p))
    compare('coordinate_origin_force',shifted_ev.rhs.p,force)
    for i in range(2):
        compare('coordinate_origin_hframe_'+str(i),shifted_ev.h_frames[i],ev.h_frames[i])
    constant = shift@(reference.held_gradient-reference.ionic_force-reference.support_force)+shift@curvature@shift/2
    compare('coordinate_origin_energy_constant',ev.energy-shifted_ev.energy,constant)

    left = tuple(rng.normal(size=(size,rank))+1j*rng.normal(size=(size,rank)) for _ in sectors)
    right = tuple(rng.normal(size=(rank,size))+1j*rng.normal(size=(rank,size)) for _ in sectors)
    tangent = model.tangent(left,right,np.array([.1+.04j,-.03j]),np.array([.07j,-.05+.02j]))

    def reconstruct(st, tan):
        return tuple(l@psi.conj().T+psi@r for psi,l,r in zip(st.frames,tan.ket,tan.bra))

    def reconstruct_dot(st, mean_rhs, tan, rhs):
        return tuple(dl@psi.conj().T+l@dp.conj().T+dp@r+psi@dr
                     for psi,dp,l,r,dl,dr in zip(st.frames,mean_rhs.frames,tan.ket,tan.bra,rhs.ket,rhs.bra))

    def original_tangent(delta, du, dp):
        dn = np.column_stack([weight*np.einsum('rij,ji->r',a,x) for a,x in zip(density_ops,delta)])
        field = (ks@dn.ravel()).reshape(nr,2)
        dg = sum(weight*np.einsum('aij,ji->a',g,x) for g,x in zip(gs,delta))
        drives = tuple(dv*np.einsum('r,rij->ij',field[:,i],density_ops[i])+np.einsum('a,aij->ij',du,gs[i]) for i in range(2))
        out = tuple(-1j*(a@x-x@a+b@rho-rho@b) for a,x,b,rho in zip(h,delta,drives,rhos))
        return out,dp/masses,-curvature@du-dg

    delta = reconstruct(state,tangent)
    expected_delta,expected_u,expected_p = original_tangent(delta,tangent.u,tangent.p)
    rhs = model.tangent_rhs(ev,tangent)
    for i,x in enumerate(reconstruct_dot(state,ev.rhs,tangent,rhs)):
        compare('independent_complex_moving_tangent_'+str(i),x,expected_delta[i])
        mutation('missing_moving_frame_connection_'+str(i),rhs.ket[i]@state.frames[i].conj().T+state.frames[i]@rhs.bra[i],expected_delta[i])
    compare('complex_tangent_ionic_velocity',rhs.u,expected_u)
    compare('complex_tangent_ionic_force',rhs.p,expected_p)
    compare('mean_number_conservation',[np.trace(x) for x in drho],np.zeros(2))
    compare('tangent_number_conservation',[np.trace(x) for x in expected_delta],np.zeros(2))

    # Literal occupation-bit CAR, independent of production Wick contraction.
    def annihilator(mode):
        a=np.zeros((2**size,2**size),complex)
        for bits in range(2**size):
            if bits & (1<<mode):
                parity=(bits & ((1<<mode)-1)).bit_count()
                a[bits^(1<<mode),bits]=(-1)**parity
        return a
    car=tuple(annihilator(i) for i in range(size))
    def lift(a):
        return sum(a[i,j]*car[i].conj().T@car[j] for i in range(size) for j in range(size))
    fock_weights=[]
    for f in occupations:
        fock_weights.append(np.array([np.prod([f[j] if bits&(1<<j) else 1-f[j] for j in range(size)]) for bits in range(2**size)]))
    weights=np.kron(*fock_weights)
    fock_density=np.diag(weights)
    identity_fock=np.eye(2**size)
    def collective(matrices):
        return np.sqrt(ncells)*weight*(np.kron(lift(matrices[0]),identity_fock)+np.kron(identity_fock,lift(matrices[1])))
    amat=(hermitian(.4)+1j*hermitian(.1),hermitian(.2)+1j*hermitian(.07))
    bmat=(hermitian(.3)+1j*hermitian(.15),hermitian(.4)+1j*hermitian(.12))
    ai=np.array([.3+.1j,-.2,.15j,.12])
    bi=np.array([-.1j,.25,.17,-.08+.02j])
    def probe(matrices,ion,name,constant=0j):
        return PreparedProbe(tuple(lambda x,a=a:a@x for a in matrices),
                             tuple(lambda x,a=a:a.conj().T@x for a in matrices),
                             ion[:ndof],ion[ndof:],constant,name,'synthetic')
    a=probe(amat,ai,'A',.12+.03j)
    b=probe(bmat,bi,'B',-.07j)
    widths=np.array([.7,1.1])
    sigma=np.block([[np.diag(widths),.5j*np.eye(ndof)],[-.5j*np.eye(ndof),np.diag(1/(4*widths))]])
    ions=IonicPreparation(sigma,'product full-line pure Gaussian with these exact u variances')
    product=FactorCrossPreparation.product(model,'explicit product of independent quasifree spins and full-line ions')
    prep=PreparedCovariance(model,electronic_authority='independent Bernoulli modes, not fixed-N',ionic=ions,cross=product)
    ea,eb=collective(amat),collective(bmat)
    def expectation(op):
        return np.trace(fock_density@op)
    ce=expectation(ea@eb)-expectation(ea)*expectation(eb)
    reverse=expectation(eb@ea)-expectation(eb)*expectation(ea)
    sab,sba=prep.source(b,ordering='AB'),prep.source(b,ordering='BA')
    for ordering,source in [('AB',sab),('BA',sba)]:
        for spin,represented in enumerate(reconstruct(state,source.electronic)):
            direct=np.zeros((size,size),complex)
            for i in range(size):
                for j in range(size):
                    coordinate=car[j].conj().T@car[i]
                    xop=np.sqrt(ncells)*(np.kron(coordinate,identity_fock) if spin==0
                                         else np.kron(identity_fock,coordinate))
                    product_op=xop@eb if ordering=='AB' else eb@xop
                    direct[i,j]=expectation(product_op)-expectation(xop)*expectation(eb)
            compare('literal_CAR_full_source_'+ordering+'_'+str(spin),represented,direct)
            mutation('saved_band_source_'+ordering+'_'+str(spin),
                     np.pad(direct[:rank,:rank],((0,1),(0,1))),direct)
    rab=ordered_readout(model,state,sab,a)
    rba=ordered_readout(model,state,sba,a)
    compare('CAR_product_AB_electronic',rab['electronic'],ce)
    compare('CAR_product_BA_electronic',rba['electronic'],reverse)
    compare('canonical_ions_AB',rab['ionic'],ai@sigma@bi)
    compare('canonical_ions_BA',rba['ionic'],bi@sigma@ai)
    compare('explicit_product_total_AB',rab['total'],ce+ai@sigma@bi)
    compare('explicit_product_total_BA',rba['total'],reverse+bi@sigma@ai)
    mutation('wrong_collective_normalization',ce/ncells,rab['electronic'])
    mutation('conjugating_ordered_right',expectation(ea@eb.conj().T)-expectation(ea)*expectation(eb.conj().T),ce)
    badj=ordered_readout(model,state,prep.source(b.adjoint()),a)
    compare('explicit_adjoint_Fock',badj['electronic'],expectation(ea@eb.conj().T)-expectation(ea)*expectation(eb.conj().T))
    expected_mean=.12+.03j+expectation(ea)/np.sqrt(ncells)+ai[:ndof]@state.u+ai[ndof:]@state.p
    compare('uncentered_mean_distinct_from_covariance',probe_mean(model,state,a),expected_mean)
    force_probe=total_force_probe(model,np.array([1.,-.3]),identity='full total force')
    compare('total_force_probe_mean',probe_mean(model,state,force_probe),np.array([1.,-.3])@force)

    # Population/equal-occupation covariance survives although commutator vanishes.
    population=np.diag([1.,0.,0.]).astype(complex)
    coherence=np.array([[0,1,0],[1,0,0],[0,0,0]],complex)
    for label,op in [('population',population),('degenerate_coherence',coherence)]:
        pp=probe((op,np.zeros((size,size))),np.zeros(2*ndof),label)
        vv=ordered_readout(model,state,prep.source(pp),pp)['electronic']
        ff=collective((op,np.zeros((size,size))))
        compare('CAR_'+label,vv,expectation(ff@ff)-expectation(ff)**2)
        mutation('dropping_'+label,0j,vv)

    # A constructive nonzero-cross joint state: Bernoulli occupation-conditioned
    # displacement of an infinite-line canonical Gaussian. Electronic marginal
    # remains exactly the same independent quasifree state.
    displacement=np.array([.18,-.09,.12,.07])
    occupied=np.array([((bits//(2**size)) & 1) for bits in range(2**(2*size))])
    centered=occupied-occupations[0][0]
    conditional_means=centered[:,None]*displacement[None,:]
    classical=conditional_means.T@(weights[:,None]*conditional_means)
    joint_ions=IonicPreparation(sigma+classical,'occupation-conditioned displaced full-line Gaussian mixture')
    cross_coeff=[]
    for i,s in enumerate(sectors):
        coeff=np.zeros((2*ndof,rank,rank),complex)
        if i==0:
            inv=np.linalg.inv(transforms[i])
            raw=inv@np.diag([occupations[0][0]*(1-occupations[0][0]),0])@inv.conj().T
            coeff=np.sqrt(ncells)*displacement[:,None,None]*raw[None,:,:]
        cross_coeff.append(coeff)
    joint=FactorCrossPreparation(tuple(cross_coeff),'explicit Bernoulli-conditioned displaced Gaussian mixture',
                                'same independent quasifree electronic marginals; cross exactly in raw factor span')
    joint_prep=PreparedCovariance(model,electronic_authority='same independent modes',ionic=joint_ions,cross=joint)
    joint_result=ordered_readout(model,state,joint_prep.source(b),a)
    emean_a=np.diag(ea)-expectation(ea)
    emean_b=np.diag(eb)-expectation(eb)
    direct_cross=np.sum(weights*(emean_a*(conditional_means@bi)+emean_b*(conditional_means@ai)))
    compare('constructive_joint_cross_Fock_Gaussian',joint_result['cross'],direct_cross)
    compare('constructive_joint_total',joint_result['total'],ce+ai@(sigma+classical)@bi+direct_cross)
    joint_ba=ordered_readout(model,state,joint_prep.source(b,ordering='BA'),a)
    compare('constructive_joint_reverse',joint_ba['total'],reverse+bi@(sigma+classical)@ai+direct_cross)
    partial=PreparedCovariance(model,electronic_authority='same marginals',ionic=ions)
    partial_result=ordered_readout(model,state,partial.source(b),a)
    if partial_result['total'] is not None or partial_result['cross'] is not None or partial_result['complete']:
        raise AssertionError('missing cross silently completed')
    CHECKS.append(dict(name='missing_cross_stays_partial',passed=True))
    try:
        PreparedCovariance(model,electronic_authority='x',ionic=ions,cross={'unsupported':'joint density'})
    except ValueError:
        CHECKS.append(dict(name='unsupported_joint_representation_rejected',passed=True))
    else:
        raise AssertionError('unsupported joint representation accepted')

    # Raw invertible factor-coordinate change: same physical mean and source.
    change=np.array([[1.1,.13j],[.08,.93]],complex)
    rebased_sectors=tuple(sector(i,t@change) for i,t in enumerate(transforms))
    rebased_model=PreparedFactorDynamics(reference,rebased_sectors)
    rebased_state=rebased_model.initial_mean(state.u,state.p)
    rebased_prep=PreparedCovariance(rebased_model,electronic_authority='same physical marginals',ionic=ions,
                                   cross=FactorCrossPreparation.product(rebased_model,'explicit same product'))
    rebased_source=rebased_prep.source(b)
    for i,(x,y) in enumerate(zip(reconstruct(state,sab.electronic),reconstruct(rebased_state,rebased_source.electronic))):
        compare('raw_frame_covariance_invariance_'+str(i),x,y)
    compare('raw_frame_mean_invariance',rebased_model.evaluate_mean(rebased_state).energy,ev.energy)

    # Instantaneous moving expected commutator: Jdot=L J+J L^T.
    # Hermitian complete coordinates retain every population and complement.
    basis=[]
    for i in range(size):
        e=np.zeros((size,size),complex); e[i,i]=1; basis.append(e)
    for i in range(size):
        for j in range(i+1,size):
            e=np.zeros((size,size),complex);e[i,j]=e[j,i]=1/np.sqrt(2);basis.append(e)
            e=np.zeros((size,size),complex);e[i,j]=-1j/np.sqrt(2);e[j,i]=1j/np.sqrt(2);basis.append(e)
    length=2*ndof+2*len(basis)
    def pack(delta,u,p):
        return np.concatenate((u,p,np.array([np.trace(e@x) for x in delta for e in basis])))
    def unpack(vector):
        out=tuple(sum(vector[2*ndof+i*len(basis)+j]*e for j,e in enumerate(basis)) for i in range(2))
        return out,vector[:ndof],vector[ndof:2*ndof]
    generator=np.column_stack([pack(*original_tangent(*unpack(col))) for col in np.eye(length)])
    comm=np.zeros((length,length),complex)
    comm[:2*ndof,:2*ndof]=sigma-sigma.T
    commdot=np.zeros_like(comm)
    for spin in range(2):
        start=2*ndof+spin*len(basis)
        for i,e in enumerate(basis):
            for j,f in enumerate(basis):
                comm[start+i,start+j]=ncells*np.trace(rhos[spin]@(e@f-f@e))
                commdot[start+i,start+j]=ncells*np.trace(drho[spin]@(e@f-f@e))
    compare('instantaneous_moving_expected_commutator',generator@comm+comm@generator.T,commdot)
    mutation('fixed_reference_commutator',np.zeros_like(commdot),commdot)
    for column in range(length):
        ion=np.zeros(2*ndof,complex)
        ops=[np.zeros((size,size),complex),np.zeros((size,size),complex)]
        if column<2*ndof:
            ion[column]=1
        else:
            spin,j=divmod(column-2*ndof,len(basis));ops[spin]=basis[j]/weight
        pr=probe(tuple(ops),ion,'complete-coordinate-'+str(column))
        forward,backward=prep.source(pr),prep.source(pr,ordering='BA')
        def combine(source):
            parts=(source.electronic,source.ionic,source.cross)
            return model.tangent(tuple(sum(x.ket[i] for x in parts) for i in range(2)),
                                 tuple(sum(x.bra[i] for x in parts) for i in range(2)),
                                 sum(x.u for x in parts),sum(x.p for x in parts))
        ft,bt=combine(forward),combine(backward)
        ct=model.tangent(tuple(x-y for x,y in zip(ft.ket,bt.ket)),tuple(x-y for x,y in zip(ft.bra,bt.bra)),ft.u-bt.u,ft.p-bt.p)
        compare('CAR_canonical_column_'+str(column),pack(reconstruct(state,ct),ct.u,ct.p),comm[:,column])
        dot=model.tangent_rhs(ev,ct)
        compare('moving_factor_commutator_column_'+str(column),
                pack(reconstruct_dot(state,ev.rhs,ct,dot),dot.u,dot.p),generator@comm[:,column])


def bindings():
    paths=[ROOT/'src/substrate_framework'/name for name in (
        'prepared_reference.py','prepared_factor_dynamics.py','prepared_covariance_source.py',
        'physical_factor_preparation.py','physical_bilinear_preparation.py')]
    paths += [Path(__file__).resolve()]
    paths += [Path(__file__).with_name(name) for name in (
        'implementation-registration.json','reference-preparation-grounding.json',
        'action-dynamics-grounding.json','additional-source-access-001.json')]
    return [dict(path=str(p),bytes=str(len(raw)),sha256=hashlib.sha256(raw).hexdigest())
            for p in paths for raw in [p.read_bytes()]]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result',type=Path,required=True)
    args=parser.parse_args()
    try:
        symbolic_scalar_first()
        finite_blocks()
    except Exception as exc:
        FAILURES.append(dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()))
    result=dict(passed=not FAILURES,checks=CHECKS,failures=FAILURES,original_source_bindings=bindings(),
                scope='new finite synthetic mathematical primitives only; no time step, material state or PBE recalibration',
                comparison_allowance=dict(rtol=1e-10,atol=1e-12),
                execution_order='independent SymPy scalar FIRST; mean; tangent; CAR/canonical sources; moving commutator')
    args.result.parent.mkdir(parents=True,exist_ok=True)
    # Exclusive creation preserves a first failed claim rather than overwriting it.
    with args.result.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    return 0 if result['passed'] else 1


if __name__=='__main__':
    raise SystemExit(main())

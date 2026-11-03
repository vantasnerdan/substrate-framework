"""Executable finite-dose oracle. No material source, exporter, or cache needed.

Run from repository root:
PYTHONPATH=src:proposals/P254-ni-collective-response/attempts/0029 \
  python proposals/P254-ni-collective-response/attempts/0029/oracle_reduced_electronic_dynamics.py
"""
from copy import copy
from dataclasses import replace
import json

import numpy as np
from scipy.integrate import solve_ivp
import sympy as sp

from substrate_framework.material_clamped_array import ClampedArray
from reduced_coefficients import QCoefficients, ReducedCoefficients
from reduced_electronic_dynamics import ActiveBandDynamics, BaselineState, DeltaState, compare_observables
from reduced_finite_reference import DenseFiniteReference


def symbolic_oracle():
    # Exact noncommuting 2x2 finite commutator, not a scalar tangent surrogate.
    matrices = [sp.Matrix(2,2,sp.symbols(f'{name}0:4')) for name in ('H','B','D','V')]
    h,b,d,v = matrices
    comm = lambda a,c:a*c-c*a
    residual = comm(h+v,b+d)-comm(h,b)-comm(h,d)-comm(v,b)-comm(v,d)
    assert all(sp.expand(value) == 0 for value in residual)
    assert any(sp.expand(value) != 0 for value in comm(v,d))
    u,du,r,dr,k,g = sp.symbols('u du r dr k g')
    scalar = lambda x,y:g*x*y+k*y*y/2
    exact = scalar(u+du,r+dr)-scalar(u,r)
    split = g*(u*dr+du*r+du*dr)+k*r*dr+k*dr*dr/2
    assert sp.expand(exact-split) == 0
    # Bilinear complex q/-q pairing uses transpose, never conjugation.
    left,right,kl = sp.symbols('left right kl')
    assert sp.diff(left*kl*right,left) == kl*right
    fourier = sp.Matrix([[1,1],[1,-1]])/sp.sqrt(2)
    assert fourier.T*fourier == sp.eye(2)
    p,m,z,a = sp.symbols('p m z a',nonzero=True)
    momenta = fourier*sp.Matrix([p,0])
    assert sp.simplify(sum(value**2/(2*m) for value in momenta)-p*p/(2*m)) == 0
    canonical_scalar = (p-sp.sqrt(2)*z*a)**2/(2*m)
    assert sp.simplify(sp.diff(canonical_scalar,p)-(p-sp.sqrt(2)*z*a)/m) == 0
    # Self-inverse Nyquist q=-1/2: G=0 alone is not reality closed.
    nyquist = -sp.Rational(1,2)
    involution = lambda integer_g:-integer_g-nyquist-nyquist
    closed = {sp.Integer(0),involution(sp.Integer(0))}
    assert closed == {sp.Integer(0),sp.Integer(1)}
    assert all(involution(involution(value)) == value for value in closed)
    pairing = sp.Matrix([[0,1],[1,0]])
    features = sp.Matrix(sp.symbols('n0:2'))
    assert sp.expand((features.T*pairing*features)[0]/2-features[0]*features[1]) == 0
    return {'exact_commutator':True,'finite_screening_scalar':True,
            'orthonormal_impulse':True,'canonical_field_velocity':True,
            'self_inverse_feature_closure':True}


def fixture():
    rng = np.random.default_rng(739410)
    mesh = (3,1,1)
    n,ns,m,d = 3,6,2,3
    indices = np.indices(mesh).reshape(3,-1).T
    target = np.asarray([np.asarray([spin*n+(k+q)%n for spin in range(2) for k in range(n)])
                         for q in range(n)],dtype=int)
    minus = np.array([0,2,1])
    def hermitian(size,scale):
        raw = rng.normal(size=(size,size))+1j*rng.normal(size=(size,size))
        return scale*(raw+raw.conj().T)/2
    def edges(channels,scale):
        value = np.zeros((n,ns,channels,m,m),complex)
        for q in range(n):
            for slot in range(ns):
                inverse = (int(minus[q]),int(target[q,slot]))
                if (q,slot) > inverse:
                    continue
                for channel in range(channels):
                    raw = scale*(rng.normal(size=(m,m))+1j*rng.normal(size=(m,m)))
                    if (q,slot) == inverse:
                        raw = (raw+raw.conj().T)/2
                    value[q,slot,channel] = raw
                    value[inverse[0],inverse[1],channel] = raw.conj().T
        return value
    ionic = edges(d,.17)
    # Ragged closed families with a nontrivial Gamma inversion and a permuted
    # partner list: no fixed feature dimension or identical-index pairing works.
    gs = [np.array([[-1,0,0],[0,0,0],[1,0,0]]),
          np.array([[0,0,0],[1,0,0]]),np.array([[-1,0,0],[0,0,0]])]
    partners = [np.array([2,1,0,5,4,3]),np.array([1,0,3,2]),np.array([1,0,3,2])]
    profiles = [np.zeros((ns,2*len(g),m,m),complex) for g in gs]
    for q in range(n):
        ng = len(gs[q])
        for slot in range(ns):
            spin = slot//n
            for r in range(spin*ng,(spin+1)*ng):
                inverse = (int(minus[q]),int(target[q,slot]),int(partners[q][r]))
                if (q,slot,r) > inverse:
                    continue
                raw = .13*(rng.normal(size=(m,m))+1j*rng.normal(size=(m,m)))
                if (q,slot,r) == inverse:
                    raw = (raw+raw.conj().T)/2
                profiles[q][slot,r] = raw
                profiles[inverse[0]][inverse[1],inverse[2]] = raw.conj().T
    field = np.asarray([[hermitian(m,.12) for _ in range(3)] for _ in range(ns)])
    h0 = np.asarray([np.diag([-.4,.6])+hermitian(m,.21) for _ in range(ns)])
    rho0 = np.asarray([.45*np.eye(m)+hermitian(m,.08) for _ in range(ns)])
    lattice = np.zeros(mesh+(d,d),complex)
    lattice[0,0,0] = np.diag([.7,.9,1.2])
    lattice[1,0,0] = lattice[2,0,0] = -.06*np.eye(d)
    mechanical = ClampedArray(lattice,np.array([2.,3.,4.]),np.array([1.3]),
        np.zeros((d,d)),np.diag([.08,.07,.06]),np.full((d,3),.011),
        np.array([.013,-.019,.007]),np.array([.012,-.003,.002]),
        np.array([[1.,0.,0.]]),1.7,np.eye(2*d)*.5,np.eye(2)*.5,np.array([.2,.3,.4]))
    records = []
    for q in range(n):
        inverse = minus[q]
        reverse = np.asarray([ionic[inverse,target[q,s]].transpose(0,2,1) for s in range(ns)])
        density = np.asarray([profiles[inverse][target[q,s],partners[q]].transpose(0,2,1) for s in range(ns)])
        kernel = np.kron(np.array([[.31,.073],[.073,.23]]),np.eye(len(gs[q])))*(1+.1*(q != 0))
        records.append(QCoefficients(ionic[q],reverse,density,np.moveaxis(profiles[q],1,-1),
                                     kernel,partners[q],gs[q],np.repeat(['up','down'],len(gs[q]))))
    c = ReducedCoefficients(mesh=mesh,k_indices=indices,q_indices=indices,
        slot_spin=np.repeat(['up','down'],n),slot_k=np.tile(np.arange(n),2),target=target,minus=minus,
        h0=h0,rho0=rho0,field=field,field_reverse=field.transpose(0,1,3,2),
        mechanical=mechanical,q_provider=lambda q:records[q],field_polarizations=np.array([[1.,0.,0.]]))
    baseline = BaselineState(rho0+np.asarray([hermitian(m,.019) for _ in range(ns)]),
                             rng.normal(size=c.mechanical_sizes[0])*.037+0j)
    delta_rho = np.asarray([hermitian(n*m,.027) for _ in range(2)])
    x = [rng.normal(size=size)*.07+0j for size in c.mechanical_sizes]
    x[1] = .07*(rng.normal(size=2*d)+1j*rng.normal(size=2*d))
    x[2] = x[1].conj()
    return c,records,baseline,DeltaState(delta_rho,tuple(x))


def close(left,right,label,atol=3e-11):
    error = float(np.max(np.abs(np.asarray(left)-np.asarray(right))))
    if error > atol:
        raise AssertionError(f'{label}: max absolute error {error:.6g} > {atol}')
    return error


def run_oracle():
    symbolic = symbolic_oracle()
    c,records,b,delta = fixture()
    model = ActiveBandDynamics(c)
    reference = DenseFiniteReference(c)
    xb = tuple(b.x if q == 0 else np.zeros(size,complex) for q,size in enumerate(c.mechanical_sizes))
    rb = reference.diagonal(b.rho)
    rk = rb+delta.rho
    xk = tuple(base+kick for base,kick in zip(xb,delta.x))
    dense_b = reference.rhs(rb,xb)
    dense_k = reference.rhs(rk,xk)
    expected = DeltaState(dense_k[0]-dense_b[0],tuple(a-bb for a,bb in zip(dense_k[1],dense_b[1])))
    observed = model.delta_rhs_state(b,delta)
    rhs_error = close(model.pack_delta(observed),model.pack_delta(expected),'finite nonlinear RHS')
    baseline_dot = model.unpack_baseline(model.baseline_rhs(0.,model.pack_baseline(b)))
    close(baseline_dot.rho,reference.diagonal_blocks(dense_b[0]),'nonstationary baseline density')
    close(baseline_dot.x,dense_b[1][0],'affine baseline mechanics/current')
    if np.linalg.norm(baseline_dot.rho) < 1e-4:
        raise AssertionError('fixture lost actual baseline drift')
    initial_drift = model.unpack_baseline(model.baseline_rhs(0.,model.pack_baseline(model.initial_baseline())))
    if np.linalg.norm(initial_drift.rho) < 1e-4:
        raise AssertionError('actual compressed preparation was silently made stationary')
    eb,ek = reference.energy(rb,xb),reference.energy(rk,xk)
    split = model.energy_difference(b,delta).as_dict()
    energy_errors = {key:close(split[key],ek[key]-eb[key],f'{key} exact energy difference') for key in eb}
    close(split['total'],sum(ek.values())-sum(eb.values()),'total exact scalar difference')
    close(model.canonical_velocity(delta),tuple(value[:c.d] for value in reference.rhs(reference.rho0,delta.x)[1]),
          'global-field canonical velocity')
    impulse = model.initial_delta()
    close(model.energy_difference(model.initial_baseline(),impulse).total,
          1/(2*c.mechanical.masses[0]),'unit local canonical impulse')
    close(np.sum(np.abs(np.asarray([x[c.d] for x in impulse.x]))**2),1.,'unit Fourier Parseval')
    # Scalar conservation is an independent local derivative check of the dense
    # Hamiltonian/feedback, not a source-text or bare-nonthrow test.
    step = 2e-6
    ep = reference.energy(rk+step*dense_k[0],tuple(x+step*v for x,v in zip(xk,dense_k[1])))
    em = reference.energy(rk-step*dense_k[0],tuple(x-step*v for x,v in zip(xk,dense_k[1])))
    energy_slope = abs((sum(ep.values())-sum(em.values()))/(2*step))
    if energy_slope > 2e-8:
        raise AssertionError(f'independent scalar gradient inconsistency: {energy_slope}')
    # Each mutation changes actual public numeric behavior in the finite model.
    mutations = {}
    def mutation(name,mutated):
        changed = ActiveBandDynamics(mutated).delta_rhs_state(b,delta)
        gap = float(np.linalg.norm(model.pack_delta(changed)-model.pack_delta(expected)))
        if gap < 1e-5:
            raise AssertionError(f'finite fixture failed to expose {name}: {gap}')
        mutations[name] = gap
    def bank(values,**kwargs):
        return replace(c,q_provider=lambda q:values[q],**kwargs)
    mutation('wrong_q_target_fold',bank([records[c.minus[q]] for q in range(c.nq)]))
    mutation('lost_inverse_sqrt_N',bank([replace(value,ionic=value.ionic*c.sqrt_n,
        ionic_reverse=value.ionic_reverse*c.sqrt_n,density=value.density*c.sqrt_n,
        potential=value.potential*c.sqrt_n) for value in records]))
    mutation('duplicate_explicit_spin_feedback',bank([replace(value,ionic_reverse=2*value.ionic_reverse)
        for value in records],field_reverse=2*c.field_reverse))
    mutation('density_conjugation',bank([replace(value,density=value.density.conj()) for value in records]))
    mutation('potential_transpose',bank([replace(value,potential=value.potential.swapaxes(1,2)) for value in records]))
    mutation('reverse_q_dagger',bank([replace(value,ionic_reverse=value.ionic_reverse.conj()) for value in records]))
    mutation('duplicate_screening',bank([replace(value,kernel=2*value.kernel) for value in records]))
    wrong_pairing = bank([replace(value,feature_partner=np.arange(len(value.kernel))) for value in records])
    pairing_energy = ActiveBandDynamics(wrong_pairing).energy_difference(b,delta).screening
    pairing_gap = float(abs(pairing_energy-split['screening']))
    if pairing_gap < 1e-5:
        raise AssertionError('ragged fixture failed to expose missing physical feature inversion')
    mutations['missing_feature_partner_energy'] = pairing_gap
    doubled = copy(c.mechanical)
    doubled.h = [h.copy() for h in c.mechanical.h]
    for h in doubled.h:
        h[:c.d,:c.d] *= 2
    doubled.generators = tuple(j@h for j,h in zip(doubled.j,doubled.h))
    mutation('duplicate_ionic_stiffness',bank(records,mechanical=doubled))
    dh = model.delta_hamiltonian(delta)
    nonlinear = -1j*(dh@delta.rho-delta.rho@dh)
    missing = model.pack_delta(DeltaState(observed.rho-nonlinear,observed.x))
    mutations['missing_nonlinear_delta_commutator'] = float(np.linalg.norm(missing-model.pack_delta(expected)))
    if mutations['missing_nonlinear_delta_commutator'] < 1e-5:
        raise AssertionError('finite pulse did not expose omitted nonlinear commutator')
    mutations['lost_baseline_nonstationarity'] = float(np.linalg.norm(baseline_dot.rho))
    # Exercise real SciPy integration, dense baseline callback and finite energy
    # on a small nonstationary, screened finite pulse, against independent full
    # kicked-minus-baseline dense integrations.
    interval = (0.,.14)
    baseline = model.integrate_baseline(interval,b,rtol=2e-11,atol=2e-13,max_step=.025)
    finite = model.integrate_finite(interval,delta,baseline=baseline,rtol=2e-11,atol=2e-13,max_step=.025)
    dense_base = solve_ivp(reference.packed_rhs,interval,reference.pack(rb,xb),
                           rtol=2e-11,atol=2e-13,method='DOP853',dense_output=True,max_step=.025)
    dense_kick = solve_ivp(reference.packed_rhs,interval,reference.pack(rk,xk),
                           rtol=2e-11,atol=2e-13,method='DOP853',dense_output=True,max_step=.025)
    if not dense_base.success or not dense_kick.success:
        raise AssertionError('independent dense integration failed')
    integration_error = 0.
    for time in (0.,.037,.091,.14):
        base_rho,base_x = reference.unpack(dense_base.sol(time))
        kick_rho,kick_x = reference.unpack(dense_kick.sol(time))
        exact = DeltaState(kick_rho-base_rho,tuple(a-bb for a,bb in zip(kick_x,base_x)))
        integration_error = max(integration_error,close(model.pack_delta(finite.state(time)),model.pack_delta(exact),
                                                       'finite trajectory',atol=2e-9))
        rho_at,h_at = baseline(time)
        close(reference.diagonal(rho_at),base_rho,'dense baseline access',atol=2e-9)
        close(reference.diagonal(h_at),reference.hamiltonian(base_rho,base_x),'nonstationary baseline Hamiltonian',atol=2e-9)
        close(finite.energy(time).total,split['total'],'finite excitation energy conservation',atol=2e-9)
    coarse = model.integrate_finite(interval,delta,baseline_initial=b,rtol=2e-7,atol=2e-9,max_step=.07)
    comparison = compare_observables(coarse,finite,[0.,.07,.14],lambda solution,t:np.asarray([
        solution.state(t).x[0][0],solution.energy(t).total]))
    if np.max(comparison['absolute_error']) > 2e-7:
        raise AssertionError('coarse/fine physical observable disagreement')
    return dict(symbolic=symbolic,rhs_max_absolute_error=rhs_error,energy_split_errors=energy_errors,
                independent_scalar_derivative=energy_slope,finite_trajectory_max_error=integration_error,
                mutation_gaps=mutations,coarse_fine_max_observable_error=float(np.max(comparison['absolute_error'])))


if __name__ == '__main__':
    print(json.dumps(run_oracle(),indent=2,sort_keys=True))

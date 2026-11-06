#!/usr/bin/env python3
"""Independent finite-Fock/SymPy/SciPy oracle; ROOT runs, never material proof.

Dense PW matrices and complete finite fundamental matrices exist ONLY here.
The integrated nonlinear physical mean supplies a moving nonstationary baseline.
The comparator differentiates full rectangular density matrices, not factor legs.
Both nonzero times are pulled back through an independently integrated tangent;
finite electron Fock and finite oscillator vacuum expectations check each origin.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace as NS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT/'src'))

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag
import sympy as sp

from substrate_framework.prepared_bloch_response import ResponseEdge
from substrate_framework.material_joint_field_scalar import JointFieldCoefficients
from substrate_framework.joint_bloch_response import (
    JointSector, JointStage, CartesianFieldActions, JointBlochResponse)
from substrate_framework.joint_ordered_sources import (
    JointProbe, JointProductPreparation, linear_combination)
from substrate_framework.joint_response_observables import (
    readout, source_readout, TwoSidedJointCovariance, total_force_probe,
    local_momentum_kick, spatial_pair)

CHECKS = []


def check(name, actual, expected, *, mutation=False, tolerance=2e-8):
    a, b = np.broadcast_arrays(np.asarray(actual), np.asarray(expected))
    error = float(np.max(np.abs(a-b), initial=0))
    close = a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all() and np.allclose(a, b, atol=tolerance, rtol=tolerance)
    passed = not close if mutation else close
    CHECKS.append(dict(name=name, passed=bool(passed), mutation=mutation, max_abs=error))
    if not passed:
        raise AssertionError(name+': '+str(error))


def reject(name, fn):
    try:
        fn()
    except (ValueError, TypeError):
        CHECKS.append(dict(name=name, passed=True, rejected=True))
    else:
        raise AssertionError('not rejected: '+name)


def symbolic():
    u, p, a, pi, m, mf, n, z, c, dc, cf, t, ell, g = sp.symbols('u p a pi m mf n z c dc cf t ell g', positive=True)
    energy = (c+dc)*u*u/2 + (p-sp.sqrt(n)*z*a)**2/(2*m) + sp.sqrt(n)*t*u*a + n*cf*a*a/2 + pi*pi/(2*mf) + g*u + sp.sqrt(n)*ell*a
    expected = [ (p-sp.sqrt(n)*z*a)/m, -(c+dc)*u-g-sp.sqrt(n)*t*a,
                 pi/mf, -sp.sqrt(n)*(ell+t*u-z*p/m)-n*(cf+z*z/m)*a ]
    actual = [sp.diff(energy,p), -sp.diff(energy,u), sp.diff(energy,pi), -sp.diff(energy,a)]
    for label, x, y in zip(('minimal_coupling','DeltaC_mixed','Maxwell_velocity','current_contact_feedback'), actual, expected):
        if sp.simplify(x-y) != 0:
            raise AssertionError('symbolic '+label)
        CHECKS.append(dict(name='symbolic_'+label, passed=True, exact=True))
    check('symbolic_missing_sqrtN_sensitive', float((actual[0]-(p-z*a)/m).subs({p:1,z:2,a:.3,n:3,m:2})), 0, mutation=True)


class MatrixAction:
    def __init__(self, matrix):
        self.matrix = matrix
    def apply(self, x):
        return self.matrix @ x
    def adjoint_apply(self, x):
        return self.matrix.conj().T @ x


class DenseCoupling:
    """Independent tiny rectangular action adapter; NOT a material Bloch grid."""
    def __init__(self, g, v):
        self.g, self.v = g, v
        self.bridge = NS(local_apply=lambda f,x: f[0]*v @ x,
                         local_adjoint_apply=lambda f,x: f[0].conjugate()*v.conj().T @ x)
    def prepare(self, u):
        return MatrixAction(np.einsum('a,aij->ij', u.ravel(), self.g))
    def moments(self, l, r, *, weight):
        x = l @ r.conj().T
        return np.array([[weight*np.vdot(self.v, x)]]), weight*np.einsum('aij,ij->a', self.g.conj(), x)


class FiniteSystem:
    def __init__(self):
        self.rng = np.random.default_rng(2252819)
        self.sizes, self.ranks = (2,3,2), (2,2,1)
        self.n, self.w, self.dof = 3, 1/3, 3
        self.mass, self.z, self.e, self.mf = np.array([2.,3.,4.]), 1.2, np.array([[.8,.6,0.]]), 2.7
        self.charge = self.z*self.e.T
        self.g = [np.stack([self.hermitian(n,.09) for _ in range(3)]) for n in self.sizes]
        self.v = [self.hermitian(n,.12) for n in self.sizes]
        self.h0 = [self.hermitian(n,.3) for n in self.sizes]
        self.ell = [np.stack([self.hermitian(n,.11) for _ in range(3)]) for n in self.sizes]
        self.c = np.array([[.9,.03,.02],[.03,1.2,.04],[.02,.04,1.4]],complex)
        self.dc = np.array([[.07,.02,0],[.02,.09,.01],[0,.01,.11]])
        self.contact = np.array([[.7,.04,.02],[.04,.8,.01],[.02,.01,.6]])
        self.mixed = np.array([[.12,.03,-.02],[-.04,.09,.01],[.02,.05,.08]])
        self.k = .37
        self.frames, self.ds = [], []
        for j, (n, rank) in enumerate(zip(self.sizes,self.ranks)):
            q = np.linalg.qr(self.random((n,rank)))[0]
            transform = np.eye(rank,dtype=complex)
            if rank == 2:
                transform[0,1] = .23+.31j
            occ = np.linspace(.21,.68,rank)
            if j == 1:
                occ[0] = 0
            inverse = np.linalg.inv(transform)
            self.frames.append(q @ transform)
            self.ds.append(inverse @ np.diag(occ) @ inverse.conj().T)
        self.frame_sizes = [p.size for p in self.frames]
        self.mean0 = np.concatenate([*(p.ravel() for p in self.frames), [.07,-.03,.04], [.03,.02,-.04], [.08], [.05]]).astype(complex)
        self.mean = solve_ivp(self.mean_rhs, (0,.31), self.mean0, method='DOP853', atol=2e-12, rtol=2e-12, dense_output=True)
        if not self.mean.success:
            raise RuntimeError(self.mean.message)
        self.wi = self.oscillator_w(self.mass*np.array([.8,1.2,.9]))
        self.wf = self.oscillator_w(np.array([1.25]))

    def random(self, shape):
        return self.rng.normal(size=shape)+1j*self.rng.normal(size=shape)
    def hermitian(self, n, scale):
        x = self.random((n,n)); return scale*(x+x.conj().T)
    @staticmethod
    def oscillator_w(stiffness):
        n = len(stiffness)
        return np.block([[np.diag(.5/stiffness), .5j*np.eye(n)],[-.5j*np.eye(n), np.diag(.5*stiffness)]])
    def unpack_mean(self, vector):
        frames, start = [], 0
        for n, rank in zip(self.sizes,self.ranks):
            frames.append(vector[start:start+n*rank].reshape(n,rank)); start += n*rank
        return frames, vector[start:start+3], vector[start+3:start+6], vector[start+6:start+7], vector[start+7:start+8]
    def mean_values(self, vector):
        frames,u,p,a,pi = self.unpack_mean(vector)
        rho = [ps @ d @ ps.conj().T for ps,d in zip(frames,self.ds)]
        density = self.w*sum(np.trace(v @ x) for v,x in zip(self.v,rho)).real
        hs = [h+self.k*(density-.04)*v+np.einsum('a,aij->ij',u.real,g)
              +np.einsum('a,aij->ij',(self.e.T @ a).real,l)
              for h,v,g,l in zip(self.h0,self.v,self.g,self.ell)]
        return frames,u,p,a,pi,rho,hs
    def mean_rhs(self, time, vector):
        frames,u,p,a,pi,rho,hs = self.mean_values(vector)
        grad = self.w*sum(np.einsum('aij,ji->a',g,x) for g,x in zip(self.g,rho))
        ell = self.w*sum(np.einsum('aij,ji->a',l,x) for l,x in zip(self.ell,rho))
        velocity = (p-self.charge @ a)/self.mass
        force = -(self.c+self.dc) @ u-grad-self.mixed @ self.e.T @ a
        pidot = -self.n*(self.e @ ell+self.e @ self.contact @ self.e.T @ a
                 +self.e @ self.mixed.T @ u-self.charge.T @ velocity)
        return np.concatenate([*((-1j*h @ ps).ravel() for h,ps in zip(hs,frames)),velocity,force,pi/self.mf,pidot])
    def stage(self, time):
        frames,*_,hs = self.mean_values(self.mean.sol(time))
        return JointStage(tuple(frames), tuple((lambda x,h=h: h @ x) for h in hs),
                          tuple((lambda x,h=h: h.conj().T @ x) for h in hs), float(time), 'finite-moving-physical-stage')

    def response(self, gamma):
        targets = list(range(3)) if gamma else [1,2,0]
        couplings = [DenseCoupling(self.g[s],self.v[s]) if gamma else
                     DenseCoupling(.1*self.random((3,self.sizes[t],self.sizes[s])),.1*self.random((self.sizes[t],self.sizes[s])))
                     for s,t in enumerate(targets)]
        cq = self.c.copy()
        if not gamma:
            cq[0,1] += .08j; cq[1,0] -= .08j
        coefficients = JointFieldCoefficients(self.dc,self.contact,self.mixed,'finite fixed physical scalar coefficients')
        return JointBlochResponse([JointSector(p.shape,d,self.w,0,'finite-basis-'+str(i)) for i,(p,d) in enumerate(zip(self.frames,self.ds))],
            [ResponseEdge(s,t,c) for s,(t,c) in enumerate(zip(targets,couplings))], ncells=self.n,
            q_fractional=[0,0,0 if gamma else 1/3], q_bohr_inv=[0,0,0 if gamma else .2],
            density_shape=(1,1), kernel=lambda n,**kw:self.k*n, common_curvature=cq,
            coefficients=coefficients,masses=self.mass,charges=[self.z],polarizations=self.e,field_mass=self.mf,
            field_vertices=[CartesianFieldActions(tuple(MatrixAction(x) for x in l)) for l in self.ell],
            scalar_identity='finite joint-second-order scalar',coverage_identity='three explicit finite slots, no material coverage claim')


class DenseTangent:
    """Independent full-density tangent and actual two-sided finite propagator."""
    def __init__(self, system, response, *, bare=False):
        self.s, self.r, self.bare = system, response, bare
        self.shapes = [(system.sizes[e.target], system.sizes[e.source]) for e in response.edges]
        self.sizes = [int(np.prod(x)) for x in self.shapes]
        self.ne = sum(self.sizes)
        self.size = self.ne+6+2*response.nfield
        self.solution = None
    def unpack(self, x):
        blocks, start = [], 0
        for shape,size in zip(self.shapes,self.sizes):
            blocks.append(x[start:start+size].reshape(shape)); start += size
        return blocks, x[start:start+3], x[start+3:start+6], x[start+6:start+6+self.r.nfield], x[start+6+self.r.nfield:]
    def pack(self, blocks, u, p, a, pi):
        return np.concatenate([*(x.ravel() for x in blocks),u,p,a,pi])
    def physical(self, stage, tangent):
        xs = [l @ stage.frames[e.source].conj().T+stage.frames[e.target] @ b
              for e,l,b in zip(self.r.edges,tangent.ket,tangent.bra)]
        return self.pack(xs,tangent.u,tangent.p,tangent.coordinate,tangent.momentum)
    def rhs(self, time, vector):
        s, r = self.s, self.r
        x,u,p,a,pi = self.unpack(vector)
        frames,_,_,_,_,rho,hs = s.mean_values(s.mean.sol(time))
        density = sum(s.w*np.vdot(e.coupling.v,b) for e,b in zip(r.edges,x))
        grad = sum(s.w*np.einsum('aij,ij->a',e.coupling.g.conj(),b) for e,b in zip(r.edges,x))
        dx = []
        for e,b in zip(r.edges,x):
            drive = s.k*density*e.coupling.v+np.einsum('a,aij->ij',u,e.coupling.g)
            if r.gamma:
                drive += np.einsum('a,aij->ij',np.sqrt(s.n)*s.e.T @ a,s.ell[e.source])
            commutator = hs[e.target] @ b-b @ hs[e.source]
            if not self.bare:
                commutator += drive @ rho[e.source]-rho[e.target] @ drive
            dx.append(-1j*commutator)
        if self.bare:
            return self.pack(dx,np.zeros(3),np.zeros(3),np.zeros(r.nfield),np.zeros(r.nfield))
        curvature = s.c+s.dc
        if not r.gamma:
            curvature = curvature+np.array([[0,.08j,0],[-.08j,0,0],[0,0,0]])
        udot, pdot = p/s.mass, -curvature @ u-grad
        adot, pidot = np.zeros(0), np.zeros(0)
        if r.gamma:
            ell = s.w*sum(np.einsum('aij,ji->a',l,b) for l,b in zip(s.ell,x))
            udot = (p-np.sqrt(s.n)*s.charge @ a)/s.mass
            pdot -= np.sqrt(s.n)*s.mixed @ s.e.T @ a
            adot = pi/s.mf
            pidot = -np.sqrt(s.n)*(s.e @ ell+s.e @ s.mixed.T @ u-s.charge.T @ (p/s.mass))
            pidot -= s.n*(s.e @ s.contact @ s.e.T+s.charge.T @ (s.charge/s.mass[:,None])) @ a
        return self.pack(dx,udot,pdot,adot,pidot)
    def transition(self, time):
        if self.solution is None:
            eye = np.eye(self.size,dtype=complex)
            def derivative(t, flat):
                generator = np.column_stack([self.rhs(t,col) for col in eye.T])
                return (generator @ flat.reshape(self.size,self.size)).ravel()
            self.solution = solve_ivp(derivative,(0,.31),eye.ravel(),method='DOP853',rtol=3e-11,atol=3e-12,dense_output=True)
            if not self.solution.success:
                raise RuntimeError(self.solution.message)
        return self.solution.sol(time).reshape(self.size,self.size)
    def coefficient(self, probe):
        matrices = [a(np.eye(shape[1])) for a,shape in zip(probe.actions,self.shapes)]
        return self.pack([self.s.w*a for a in matrices],probe.u,probe.p,probe.coordinate,probe.momentum)
    def probe(self, coefficient, identity):
        matrices,u,p,a,pi = self.unpack(coefficient)
        actions = [MatrixAction(b/self.s.w) for b in matrices]
        return JointProbe(tuple(x.apply for x in actions),tuple(x.adjoint_apply for x in actions),u,p,a,pi,identity,'finite observable')
    def pullback(self, time, probe):
        return self.probe(self.transition(time).conj().T @ self.coefficient(probe),probe.identity+' pulled to original t0')
    def covariance(self, ordering):
        rho = [p @ d @ p.conj().T for p,d in zip(self.s.frames,self.s.ds)]
        blocks = []
        for e in self.r.edges:
            rs,rt = rho[e.source],rho[e.target]
            blocks.append(self.s.n*np.kron(np.eye(len(rt))-rt,rs.T) if ordering == 'AB' else
                          self.s.n*np.kron(rt,(np.eye(len(rs))-rs).T))
        wi = self.s.wi if ordering == 'AB' else self.s.wi.T
        wf = (self.s.wf if ordering == 'AB' else self.s.wf.T) if self.r.gamma else np.zeros((0,0))
        return block_diag(*blocks,wi,wf)


def fermion_annihilators(n):
    result = []
    for orbital in range(n):
        op = np.zeros((2**n,2**n),complex)
        for state in range(2**n):
            if state & (1 << orbital):
                sign = (-1)**((state & ((1 << orbital)-1)).bit_count())
                op[state ^ (1 << orbital), state] = sign
        result.append(op)
    return result


class FockOracle:
    """Exact finite electronic Fock product; oscillator vacuum expectations.

    Product origin factorization avoids allocating the full tensor-product
    Hilbert matrix, without discarding its connected covariance contributions.
    """
    def __init__(self, system):
        self.s = system
        rho = block_diag(*[p @ d @ p.conj().T for p,d in zip(system.frames,system.ds)])
        self.ann = fermion_annihilators(len(rho))
        self.bilinears = np.array([[a.conj().T @ b for b in self.ann] for a in self.ann])
        occupations,modes = np.linalg.eigh(rho)
        identity = np.eye(len(self.ann[0]))
        state = identity.copy()
        for j,occupation in enumerate(occupations):
            annihilator = sum(modes[i,j].conjugate()*op for i,op in enumerate(self.ann))
            number = annihilator.conj().T @ annihilator
            state = state @ ((1-occupation)*identity+(2*occupation-1)*number)
        self.state = state
        check('Fock_state_trace',np.trace(state),1)
        self.offsets = np.cumsum((0,)+system.sizes)
    def observable(self, response, probe):
        matrix = np.zeros((sum(self.s.sizes),)*2,complex)
        for e,action in zip(response.edges,probe.actions):
            block = action(np.eye(self.s.sizes[e.source])).conj().T
            matrix[self.offsets[e.source]:self.offsets[e.source+1],self.offsets[e.target]:self.offsets[e.target+1]] += np.sqrt(self.s.n)*self.s.w*block
        return np.einsum('ij,ijab->ab',matrix,self.bilinears)
    @staticmethod
    def oscillator_pair(left, right, stiffness, ordering):
        ladder = np.diag(np.sqrt(np.arange(1,3)),1).astype(complex)
        total = 0j
        n = len(stiffness)
        for j,k in enumerate(stiffness):
            q = (ladder+ladder.conj().T)/np.sqrt(2*k)
            p = -1j*np.sqrt(k/2)*(ladder-ladder.conj().T)
            a = left[j].conjugate()*q+left[n+j].conjugate()*p
            b = right[j]*q+right[n+j]*p
            total += (a @ b if ordering == 'AB' else b @ a)[0,0]
        return total
    def parts(self, response, left, right, ordering):
        a, b = self.observable(response,left), self.observable(response,right).conj().T
        electronic = np.trace(self.state @ (a @ b if ordering == 'AB' else b @ a))
        electronic -= np.trace(self.state @ a)*np.trace(self.state @ b)
        ionic = self.oscillator_pair(np.r_[left.u,left.p],np.r_[right.u,right.p],self.s.mass*np.array([.8,1.2,.9]),ordering)
        field = self.oscillator_pair(np.r_[left.coordinate,left.momentum],np.r_[right.coordinate,right.momentum],np.array([1.25]) if response.gamma else np.zeros(0),ordering)
        return dict(electronic=electronic,ionic=ionic,field=field)


def pack_legs(tangent):
    return np.concatenate([*(x.ravel() for x in tangent.ket),*(x.ravel() for x in tangent.bra),tangent.u,tangent.p,tangent.coordinate,tangent.momentum])


def unpack_legs(response, vector):
    zero, start = response.zero_tangent(), 0
    groups = []
    for templates in (zero.ket,zero.bra):
        values = []
        for x in templates:
            values.append(vector[start:start+x.size].reshape(x.shape)); start += x.size
        groups.append(tuple(values))
    return response.tangent(*groups,vector[start:start+3],vector[start+3:start+6],
                            vector[start+6:start+6+response.nfield],vector[start+6+response.nfield:])


def leg_propagation(system, response, initial, time):
    sol = solve_ivp(lambda t,x: pack_legs(response.rhs(system.stage(t),unpack_legs(response,x))),
                    (0,time),pack_legs(initial),method='DOP853',rtol=2e-11,atol=2e-12)
    if not sol.success:
        raise RuntimeError(sol.message)
    return unpack_legs(response,sol.y[:,-1])


def exercise(system, fock, gamma):
    r = system.response(gamma)
    dense = DenseTangent(system,r)
    tag = 'Gamma' if gamma else 'nonGamma'
    preparation = JointProductPreparation(r,system.stage(0),ionic_ordered=system.wi,
        ionic_minus_ordered=system.wi,field_ordered=system.wf if gamma else np.zeros((0,0)),
        electronic_identity='original finite raw singular complex D',ionic_identity='three oscillator product vacuum W',
        field_identity='one global oscillator vacuum W, Gamma only',product_identity='explicit original electron x ionic x field tensor product')
    left = dense.probe(system.random((dense.size,)),tag+' random A')
    right = dense.probe(system.random((dense.size,)),tag+' random B')
    zero = r.zero_tangent()
    tangent = unpack_legs(r,.07*system.random(pack_legs(zero).shape))
    time = .173
    stage, rhs = system.stage(time), r.rhs(system.stage(time),tangent)
    reconstructed = dense.physical(stage,rhs)
    xs,du,dp,da,dpi = dense.unpack(reconstructed.copy())
    for i,(e,l,b) in enumerate(zip(r.edges,tangent.ket,tangent.bra)):
        ps,pt = stage.frames[e.source],stage.frames[e.target]
        xs[i] += l @ (-1j*stage.h_actions[e.source](ps)).conj().T+(-1j*stage.h_actions[e.target](pt)) @ b
    expected = dense.rhs(time,dense.physical(stage,tangent))
    check(tag+'_moving_full_density_RHS',dense.pack(xs,du,dp,da,dpi),expected)
    check(tag+'_omitted_moving_connection',reconstructed,expected,mutation=True)
    check(tag+'_independent_bra_not_dagger',tangent.bra[0],tangent.ket[0].conj().T if tangent.bra[0].shape == tangent.ket[0].conj().T.shape else np.zeros_like(tangent.bra[0]),mutation=True)
    check(tag+'_physical_readout',readout(r,stage,tangent,left)['total'],np.vdot(dense.coefficient(left),dense.physical(stage,tangent)))
    evolved = leg_propagation(system,r,tangent,.27)
    check(tag+'_actual_leg_trajectory',dense.physical(system.stage(.27),evolved),dense.transition(.27) @ dense.physical(system.stage(0),tangent))
    _,_,_,mean_field,_,_,current_h = system.mean_values(system.mean.sol(time))
    fieldless_h = [h-np.einsum('a,aij->ij',system.e.T @ mean_field,l)
                   for h,l in zip(current_h,system.ell)]
    fieldless_stage = JointStage(stage.frames,tuple((lambda x,h=h:h @ x) for h in fieldless_h),
        tuple((lambda x,h=h:h.conj().T @ x) for h in fieldless_h),time,'wrong omitted current uniform mean field')
    check(tag+'_omitted_uniform_field_in_endpoint_H',pack_legs(r.rhs(fieldless_stage,tangent)),pack_legs(rhs),mutation=True)
    for sign,label in ((-1,'omitted'),(1,'doubled')):
        wrong = copy.copy(r); wrong.curvature = r.curvature+sign*system.dc
        check(tag+'_DeltaC_'+label,wrong.rhs(stage,tangent).p,rhs.p,mutation=True)
    if gamma:
        wrong = copy.copy(r); wrong.sqrt_n = 1
        check('Gamma_wrong_sqrtN',pack_legs(wrong.rhs(stage,tangent)),pack_legs(rhs),mutation=True)
        wrong = copy.copy(r); wrong.mixed = np.zeros_like(r.mixed)
        check('Gamma_missing_mixed',pack_legs(wrong.rhs(stage,tangent)),pack_legs(rhs),mutation=True)
        wrong = copy.copy(r); wrong.contact = np.zeros_like(r.contact)
        check('Gamma_missing_contact',wrong.rhs(stage,tangent).momentum,rhs.momentum,mutation=True)
        wrong = copy.copy(r); wrong.field_vertices = tuple(CartesianFieldActions(tuple(MatrixAction(np.zeros_like(x)) for x in l)) for l in system.ell)
        check('Gamma_missing_current_vertices',pack_legs(wrong.rhs(stage,tangent)),pack_legs(rhs),mutation=True)
        a = np.array([.2+.3j,-.1+.4j,.2-.5j])
        vertex = r.field_vertices[0]
        matrix = np.einsum('a,aij->ij',a,system.ell[0])
        vectors = system.random((system.sizes[0],2))
        check('complex_field_forward',vertex.apply(a,vectors),matrix @ vectors)
        check('complex_field_independent_adjoint',vertex.apply(a,vectors,adjoint=True),matrix.conj().T @ vectors)
        check('complex_field_wrong_adjoint',vertex.apply(a,vectors),matrix.conj().T @ vectors,mutation=True)
        check('Maxwell_current_feedback',2*r.current_variation(stage,tangent)['total'],rhs.momentum)
    else:
        reject('nonGamma_rejects_extra_photon',lambda:r.tangent(zero.ket,zero.bra,zero.u,zero.p,[0],[0]))
        check('nonGamma_complex_curvature_preserved',r.curvature[0,1].imag,.08)
    force = total_force_probe(r,np.array([.3+.2j,.1,-.2j]),identity='full joint total force')
    check(tag+'_total_force_readout',readout(r,stage,tangent,force)['total'],np.vdot(np.array([.3+.2j,.1,-.2j]),rhs.p))
    bmatrices = [action(np.eye(shape[1])) for action,shape in zip(right.actions,dense.shapes)]
    wrong_right = JointProbe(right.actions,tuple((lambda x,b=b:b.T @ x) for b in bmatrices),
        right.u,right.p,right.coordinate,right.momentum,'wrong transpose instead of independent adjoint',right.units)
    check(tag+'_wrong_source_adjoint',
        source_readout(r,system.stage(0),preparation.source(wrong_right,ordering='BA'),left)['total'],
        source_readout(r,system.stage(0),preparation.source(right,ordering='BA'),left)['total'],mutation=True)
    for ordering in ('AB','BA'):
        src = preparation.source(right,ordering=ordering)
        complete = linear_combination(*[(1,x) for x in src.parts.values()])
        check(tag+'_'+ordering+'_source_matrix',dense.physical(system.stage(0),complete),dense.covariance(ordering) @ dense.coefficient(right))
        parts = source_readout(r,system.stage(0),src,left)
        exact = fock.parts(r,left,right,ordering)
        for origin in ('electronic','ionic','field'):
            check(tag+'_'+ordering+'_Fock_'+origin,parts[origin],exact[origin])
        check(tag+'_'+ordering+'_missing_N_source',parts['electronic']/system.n,exact['electronic'],mutation=True)
        check(tag+'_'+ordering+'_wrong_ionic_mean_normalization',parts['ionic']/system.n,exact['ionic'],mutation=True)
        for origin, initial in src.parts.items():
            if not gamma and origin == 'field':
                continue
            final = leg_propagation(system,r,initial,.11)
            check(tag+'_'+ordering+'_propagated_'+origin,dense.physical(system.stage(.11),final),dense.transition(.11) @ dense.physical(system.stage(0),initial))
    ab,ba = preparation.source(right,ordering='AB'),preparation.source(right,ordering='BA')
    check(tag+'_coherent_impulse_distinct',dense.physical(system.stage(0),preparation.coherent_impulse(right)),
          -1j*(dense.covariance('AB')-dense.covariance('BA')) @ dense.coefficient(right))
    check(tag+'_noise_is_not_coherent',source_readout(r,system.stage(0),ab,left)['total'],readout(r,system.stage(0),preparation.coherent_impulse(right),left)['total'],mutation=True)
    covariance = TwoSidedJointCovariance(preparation,dense.pullback,dynamics='screened-joint',propagation_identity='independent nonlinear moving mean and full density DOP853 fundamental solution')
    for t,s in ((.27,.11),(.11,.27),(.173,.173),(.27,0.)):
        pair = covariance.pair(left,right,t,s)
        for ordering in ('AB','BA'):
            direct = np.vdot(dense.coefficient(left),dense.transition(t) @ dense.covariance(ordering) @ dense.transition(s).conj().T @ dense.coefficient(right))
            check(f'{tag}_{ordering}_two_sided_{t}_{s}',pair[ordering]['total'],direct)
            pulled_a,pulled_b = dense.pullback(t,left),dense.pullback(s,right)
            exact = fock.parts(r,pulled_a,pulled_b,ordering)
            for origin in exact:
                check(f'{tag}_{ordering}_two_sided_Fock_{origin}_{t}_{s}',pair[ordering][origin],exact[origin])
        check(f'{tag}_commutator_{t}_{s}',pair['commutator'],pair['AB']['total']-pair['BA']['total'])
    actual = covariance.pair(left,right,.27,.11)['AB']['total']
    check(tag+'_wrong_second_time_zero',covariance.pair(left,right,.27,0)['AB']['total'],actual,mutation=True)
    check(tag+'_wrong_stationarity',covariance.pair(left,right,.16,0)['AB']['total'],actual,mutation=True)
    equal = dense.transition(.173) @ dense.covariance('AB') @ dense.transition(.173).conj().T
    check(tag+'_generated_electron_ion_cross',equal[:dense.ne,dense.ne:dense.ne+6],0,mutation=True)
    if gamma:
        check('Gamma_generated_ion_field_cross',equal[dense.ne:dense.ne+6,dense.ne+6:],0,mutation=True)
    reset = block_diag(equal[:dense.ne,:dense.ne],equal[dense.ne:dense.ne+6,dense.ne:dense.ne+6],equal[dense.ne+6:,dense.ne+6:])
    check(tag+'_wrong_product_reset',np.vdot(dense.coefficient(left),reset @ dense.coefficient(right)),covariance.pair(left,right,.173,.173)['AB']['total'],mutation=True)
    bare = DenseTangent(system,r,bare=True)
    pure_e = dense.coefficient(left).copy(); pure_e[dense.ne:] = 0
    electronic_probe = dense.probe(pure_e,'electronic-only noise probe')
    bare_covariance = TwoSidedJointCovariance(preparation,bare.pullback,dynamics='bare-electronic',propagation_identity='independent bare endpoint-H transport; no screened feedback')
    check(tag+'_bare_not_screened',bare_covariance.pair(electronic_probe,electronic_probe,.27,.11)['AB']['total'],covariance.pair(electronic_probe,electronic_probe,.27,.11)['AB']['total'],mutation=True)
    reject(tag+'_no_later_product_repreparation',lambda:JointProductPreparation(r,system.stage(.1),ionic_ordered=system.wi,ionic_minus_ordered=system.wi,field_ordered=system.wf if gamma else np.zeros((0,0)),electronic_identity='original',ionic_identity='original',field_identity='original',product_identity='wrong reset'))
    kick = local_momentum_kick(r,[0,0,1],[1,0,0],.4)
    check(tag+'_canonical_local_kick',kick.p[0],.4*np.exp(-2j*np.pi*r.q_fractional[2])/np.sqrt(3))


def bind(path):
    data = path.read_bytes()
    return dict(path=str(path.relative_to(ROOT)),bytes=len(data),sha256=hashlib.sha256(data).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('preserve existing receipt: '+str(args.output))
    registration = json.loads((HERE/'joint-response-registration.json').read_text())
    for item in registration['sources']+registration['data']:
        if bind(ROOT/item['path']) != item:
            raise ValueError('registered source/data bytes changed: '+item['path'])
    symbolic()
    system = FiniteSystem()
    check('moving_nonstationary_baseline',system.mean.sol(.27),system.mean0,mutation=True)
    check('original_complex_D_retained',system.ds[0].imag,0,mutation=True)
    check('original_singular_D_retained',np.linalg.det(system.ds[1]),0)
    fock = FockOracle(system)
    exercise(system,fock,True)
    exercise(system,fock,False)
    q = np.array([[0,0,0],[0,0,1/3],[0,0,2/3]])
    values = np.array([.4+.2j,.1-.7j,-.2+.3j])
    expected = sum(np.exp(2j*np.pi*x[2])*v/3 for x,v in zip(q,values))
    check('q_minus_q_spatial_normalization',spatial_pair(q,values,[0,0,1],[0,0,0],ncells=3),expected)
    check('wrong_spatial_dagger',spatial_pair(q,values.conj(),[0,0,1],[0,0,0],ncells=3),expected,mutation=True)
    check('wrong_spatial_N',3*expected,expected,mutation=True)
    result = dict(schema='P254-joint-response-oracle-v1',status='PASS',checks=CHECKS,
        source_bindings=[bind(ROOT/p) for p in registration['implementation']],registration=bind(HERE/'joint-response-registration.json'),
        data_bindings=registration['data'],command=registration['command'],
        scope='Finite mathematical oracle only; no material trajectory, convergence, equilibrium, damping or FDT claim',
        numerical_methods='SymPy, finite-Fock, oscillator vacuum, independent SciPy DOP853 moving mean and full-density fundamental matrices',
        times=[0,.11,.173,.27],q_fractional=[[0,0,0],[0,0,1/3]],endpoint_sizes=list(system.sizes),ranks=list(system.ranks))
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(json.dumps(dict(status=result['status'],checks=len(CHECKS),output=str(args.output))))


if __name__ == '__main__':
    main()

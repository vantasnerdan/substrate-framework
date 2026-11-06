#!/usr/bin/env python3
"""Independent finite dense/SymPy oracle, not a material trajectory.

Uses A27's finite byte fixture solely for checked Bloch geometry. Dense PW
matrices exist ONLY here. Literal Fourier sums and original density matrices
supply comparators; production retains its independent full-row factor legs.
Root runs this new oracle; no prior oracle or material calculation is replayed.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace as NS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path[:0] = [str(ROOT/'src'), str(HERE.parent/'0027')]

import numpy as np
from scipy.linalg import expm
import sympy as sp
from bloch_oracle_fixture import fixture
from substrate_framework.qe_bloch_grid import build_bloch_grid, UNITARY_POLICY
from substrate_framework.qe_bloch_edge import BlochEdge
from substrate_framework.qe_bloch_finiteq import BlochFiniteQBridge
from substrate_framework.material_reference_state import prepare_source
from substrate_framework.physical_factor_preparation import FactorFunctional, FactorPreparation
from substrate_framework.prepared_reference import PreparedReference, PreparedSector
from substrate_framework.prepared_factor_dynamics import PreparedFactorDynamics
from substrate_framework.prepared_bloch_response import PreparedBlochResponse, ResponseEdge

CHECKS = []


def compare(name, actual, expected, *, mutation=False):
    a, b = np.asarray(actual), np.asarray(expected)
    residual = float(np.max(np.abs(a-b), initial=0))
    close = a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all() and np.allclose(a, b, atol=3e-11, rtol=3e-10)
    passed = not close if mutation else close
    CHECKS.append(dict(name=name, passed=bool(passed), mutation=mutation, max_abs=residual))
    if not passed:
        raise AssertionError(name + ': ' + str(residual))


def rejects(name, fn):
    try:
        fn()
    except (ValueError, TypeError):
        CHECKS.append(dict(name=name, passed=True, rejection=True))
    else:
        raise AssertionError('not rejected: '+name)


def symbolic_scalar():
    # Treat +q and -q coordinates as independent complex analytic variables;
    # the reality constraint is applied only after differentiating the scalar.
    n, nm, u, um, p, pm, k, c, m, g, gm, dv = sp.symbols('n nm u um p pm k c m g gm dv')
    pair_energy = dv*nm*k*n + um*g + gm*u + um*c*u + pm*p/m
    for name, delta in [('Kq',sp.diff(pair_energy,nm)-dv*k*n),
                        ('ionic_sign',-sp.diff(pair_energy,um)+g+c*u),
                        ('velocity',sp.diff(pair_energy,pm)-p/m)]:
        if sp.simplify(delta) != 0:
            raise AssertionError(name)
        CHECKS.append(dict(name='symbolic_'+name, passed=True, exact=True))
    # Gamma/self-inverse real scalar has its own 1/2, not twice the RHS.
    real_energy = dv*k*n*n/2+u*g+c*u*u/2+p*p/(2*m)
    if sp.simplify(-sp.diff(real_energy,u)+g+c*u) != 0:
        raise AssertionError('self_inverse_count')
    CHECKS.append(dict(name='symbolic_self_inverse_no_double',passed=True,exact=True))
    ht, hs, l, r, ps, pt, sh, th, ds, dt, v = [sp.MatrixSymbol(x,2,2) for x in ('Ht','Hs','L','R','Ps','Pt','PsH','PtH','Ds','Dt','V')]
    # Independent ket/bra symbols avoid assuming Hermiticity of the tangent.
    x = l*sh+pt*r
    reconstructed = -sp.I*(ht*l+v*ps*ds)*sh + l*(sp.I*sh*hs) - sp.I*ht*pt*r + pt*sp.I*(r*hs+dt*th*v)
    expected = -sp.I*(ht*x-x*hs+v*ps*ds*sh-pt*dt*th*v)
    if sp.expand(reconstructed-expected) != sp.ZeroMatrix(2,2):
        raise AssertionError('moving symbolic product rule')
    CHECKS.append(dict(name='symbolic_moving_product_rule',passed=True,exact=True))


def finite_case():
    rng = np.random.default_rng(2252804)
    def complex_random(shape):
        return rng.normal(size=shape)+1j*rng.normal(size=shape)
    def hermitian(n, scale):
        a=complex_random((n,n)); return scale*(a+a.conj().T)
    with TemporaryDirectory(prefix='p254-response-oracle-') as directory:
        inputs, bindings, _, _ = fixture(directory)
        # A separate finite 5^3 collocation scalar; no cutoff convergence claim.
        inputs.reference=replace(inputs.reference,fft_grid=(5,5,5),fft_smooth=(5,5,5),fft_box=(5,5,5))
        grid=build_bloch_grid(inputs.reference,bindings,symmetry_policy=UNITARY_POLICY)
        frames=[grid.load_frame(slot) for slot in grid.slots]
        lookup={(f.ispin,f.basis.slot.mesh_index):i for i,f in enumerate(frames)}
        nr, ndof, volume = 125, 6, abs(np.linalg.det(inputs.reference.cell_bohr))
        dv=volume/nr
        coords=np.array(list(np.ndindex(5,5,5)),float)[:,::-1]/5
        # np.ndindex is last-axis-fastest; reversed coordinates are F-order.
        waves=[np.exp(2j*np.pi*coords@f.basis.miller) for f in frames]
        hs=[hermitian(f.coefficients.shape[0],.025) for f in frames]
        gs=[np.stack([hermitian(len(h),.008) for _ in range(ndof)]) for h in hs]
        spin_kernel=np.array([[.8,.17],[.17,.6]])
        curvature=np.diag(np.linspace(.8,1.4,ndof))+.03*np.ones((ndof,ndof))
        masses=np.linspace(2.,5.,ndof)
        class MeanCoupling:
            def __init__(self,i):
                self.i=i; self.npw=len(hs[i]); self.nr=nr; self.ndof=ndof
                self.grid_measure_bohr3=dv; self.cell_volume_bohr3=volume
            def moments(self,l,r,coefficient=None):
                x=l@r.conj().T if coefficient is None else l@coefficient@r.conj().T
                w=waves[self.i]
                return np.einsum('ri,ij,rj->r',w,x,w.conj())/volume, np.einsum('aij,ji->a',gs[self.i],x)
            def potential_action(self,v,x):
                w=waves[self.i]; return w.conj().T@(v[:,None]*(w@x))/nr
            def ionic_action(self,u,x):
                return np.einsum('a,aij->ij',u,gs[self.i])@x
        sectors=[]
        for i,f in enumerate(frames):
            # Deliberately unequal ranks across endpoints, singular D on some.
            rank=1+i%2
            psi=f.coefficients[:,:rank]
            populations=np.linspace(.12,.67,rank)
            if i%7==0: populations[0]=0.
            state=prepare_source(psi,populations,source_identity='finite-'+str(i),basis_identity=f.basis.identity,authority='finite oracle declared probabilities, not saved MV labels')
            functional=FactorFunctional(state.gram,{}, {}, {},'finite-'+str(i),f.basis.identity,len(psi))
            prep=FactorPreparation(functional,state.density_coefficients,'finite raw nonorthogonal state, empty complement')
            sectors.append(PreparedSector(prep,psi,MeanCoupling(i),lambda x,h=hs[i]:h@x,f.ispin-1,f.basis.slot.weight,f.basis.identity,'literal finite scalar operators'))
        ref=PreparedReference(nbar=np.full((nr,2),.004),ionic_force=np.linspace(.01,.06,ndof),
            held_density=np.full((nr,2),.001),held_gradient=np.linspace(-.01,.01,ndof),curvature=curvature,masses=masses,
            support_force=np.full(ndof,.003),kernel=lambda n:n@spin_kernel.T,ncells=64,
            expansion_authority='finite collocation scalar',ionic_authority='finite linear force and fixed C',coverage_authority='all128 finite slots',
            support_authority='explicit constant support',coordinate_identity='six finite coordinates',state_authority='declared independent quasifree states')
        model=PreparedFactorDynamics(ref,sectors)
        mean0=model.initial_mean(np.linspace(.02,.07,ndof),np.linspace(-.03,.04,ndof))
        # Distinct admissible current state, reached by unitary drift without
        # reinterpreting saved energies as the current H or freezing the mean.
        mean1=model.mean(tuple(expm(-.17j*h)@p for h,p in zip(hs,mean0.frames)),mean0.u+.04,mean0.p-.01)
        ev0,ev1=model.evaluate_mean(mean0),model.evaluate_mean(mean1)
        for qz in (.25,0.,.5):
            q=np.array([0.,0.,qz]); shift=int(round(qz*4))
            kq=spin_kernel.astype(complex)
            if qz==.25: kq[0,1]+=.09j; kq[1,0]-=.09j
            cq=curvature.astype(complex)
            if qz==.25: cq[0,1]+=.07j; cq[1,0]-=.07j
            links=[]; gm=[]; phases=[]
            class EdgeCoupling:
                def __init__(self,edge,g):
                    self.edge=edge; self.inputs=inputs; self.ndof=ndof; self.g=g
                    self.bridge=BlochFiniteQBridge(inputs,edge,block_size=2)
                def prepare(self,u):
                    matrix=np.einsum('a,aij->ij',u.ravel(),self.g)
                    return NS(apply=lambda x:matrix@x,adjoint_apply=lambda x:matrix.conj().T@x)
                def moments(self,l,r,coefficient=None,*,weight):
                    x=l@r.conj().T if coefficient is None else l@coefficient@r.conj().T
                    return self.bridge.density_from_dyads(l,r,coefficient,weight=weight),weight*np.einsum('aij,ij->a',self.g.conj(),x)
            for i,f in enumerate(frames):
                index=list(f.basis.slot.mesh_index); index[2]=(index[2]+shift)%4
                j=lookup[f.ispin,tuple(index)]
                edge=BlochEdge.create(inputs.reference,f.basis,frames[j].basis,q)
                g=gs[i] if qz==0 else .025*complex_random((ndof,len(hs[j]),len(hs[i])))
                gm.append(g); phases.append(np.exp(-2j*np.pi*coords@edge.reciprocal_fold))
                links.append(ResponseEdge(i,j,EdgeCoupling(edge,g)))
            response=PreparedBlochResponse(model,links,kernel=lambda n,*,q_bohr_inv,k=kq:n@k.T,
                curvature=cq,authority_text='finite collocation q scalar with exact directed transpose')
            tangent=response.tangent([.04*complex_random((len(hs[e.target]),sectors[e.source].initial_frame.shape[1])) for e in links],
                [.04*complex_random((sectors[e.target].initial_frame.shape[1],len(hs[e.source]))) for e in links],
                .04*complex_random((ndof,)),.03*complex_random((ndof,)))
            def reconstruct(ev,tan):
                return [l@ev.state.frames[e.source].conj().T+ev.state.frames[e.target]@b for e,l,b in zip(links,tan.ket,tan.bra)]
            def derivative(ev,tan,dot,*,connection=True):
                out=[]
                for e,l,b,dl,db in zip(links,tan.ket,tan.bra,dot.ket,dot.bra):
                    value=dl@ev.state.frames[e.source].conj().T+ev.state.frames[e.target]@db
                    if connection: value+=l@ev.rhs.frames[e.source].conj().T+ev.rhs.frames[e.target]@b
                    out.append(value)
                return out
            def dense(ev,tan):
                rho=[p@s.d@p.conj().T for p,s in zip(ev.state.frames,sectors)]
                nmean=np.zeros((nr,2),complex)
                for i,(w,x,s) in enumerate(zip(waves,rho,sectors)):
                    nmean[:,s.density_component]+=s.weight*np.einsum('ri,ij,rj->r',w,x,w.conj())/volume
                vmean=(nmean+ref.held_density-ref.nbar)@spin_kernel.T
                h=[a+w.conj().T@(vmean[:,s.density_component,None]*w)/nr+np.einsum('a,aij->ij',ev.state.u,g) for a,w,s,g in zip(hs,waves,sectors,gs)]
                x=reconstruct(ev,tan); n=np.zeros((nr,2),complex); grad=np.zeros(ndof,complex)
                for e,xx,phase,g in zip(links,x,phases,gm):
                    s=sectors[e.source]; wt,ws=waves[e.target],waves[e.source]
                    n[:,s.density_component]+=s.weight*phase.conj()*np.einsum('ri,ij,rj->r',wt,xx,ws.conj())/volume
                    grad+=s.weight*np.einsum('aij,ij->a',g.conj(),xx)
                field=n@kq.T
                drives=[waves[e.target].conj().T@((phase*field[:,sectors[e.source].density_component])[:,None]*waves[e.source])/nr+np.einsum('a,aij->ij',tan.u,g) for e,phase,g in zip(links,phases,gm)]
                dx=[-1j*(h[e.target]@xx-xx@h[e.source]+v@rho[e.source]-rho[e.target]@v) for e,xx,v in zip(links,x,drives)]
                return x,n,grad,h,rho,drives,dx
            for stage,ev in enumerate((ev0,ev1)):
                label=f'q{qz}_stage{stage}'
                x,n,g,h,rho,drives,dx=dense(ev,tangent)
                dn,dg=response.moments(ev,tangent)
                compare(label+'_density',dn,n); compare(label+'_gradient',dg,g)
                dot=response.rhs(ev,tangent); result=derivative(ev,tangent,dot)
                for i,(a,b) in enumerate(zip(result,dx)): compare(label+'_edge'+str(i),a,b)
                compare(label+'_velocity',dot.u,tangent.p/masses)
                compare(label+'_force',dot.p,-cq@tangent.u-g)
                compare(label+'_wrong_force_sign',-cq@tangent.u+g,dot.p,mutation=True)
                compare(label+'_wrong_weight',64*dn,n,mutation=True)
                compare(label+'_Gamma_bilinear_gradient',g.conj(),g,mutation=True)
                no_connection=derivative(ev,tangent,dot,connection=False)
                compare(label+'_frozen_frame',np.concatenate([a.ravel() for a in no_connection]),np.concatenate([a.ravel() for a in dx]),mutation=True)
                if stage==1:
                    stale=response.rhs(ev0,tangent)
                    compare(label+'_stale_current_H_and_frames',np.concatenate([a.ravel() for a in derivative(ev,tangent,stale)]),np.concatenate([a.ravel() for a in dx]),mutation=True)
                # Ordered source and observable, against original dense rho.
                probes=[complex_random(xx.shape)*.03 for xx in x]
                actions=[lambda z,a=a:a@z for a in probes]
                adjoints=[lambda z,a=a:a.conj().T@z for a in probes]
                for ordering in ('AB','BA'):
                    source=response.vertex_source(ev,actions,adjoints,ordering=ordering)
                    represented=reconstruct(ev,source)
                    for i,(e,b,xx) in enumerate(zip(links,probes,represented)):
                        direct=((np.eye(len(rho[e.target])))-rho[e.target])@b@rho[e.source] if ordering=='AB' else rho[e.target]@b@(np.eye(len(rho[e.source]))-rho[e.source])
                        compare(label+'_'+ordering+'_source'+str(i),xx,direct)
                    expected=sum(sectors[e.source].weight*np.vdot(a,xx) for e,a,xx in zip(links,probes,represented))
                    compare(label+'_'+ordering+'_readout',response.readout(ev,source,actions,adjoints)['total'],expected)
                compare(label+'_readout_complex',response.readout(ev,tangent,actions,adjoints,ionic_u=tangent.p)['total'],sum(sectors[e.source].weight*np.vdot(a,xx) for e,a,xx in zip(links,probes,x))+np.vdot(tangent.p,tangent.u))
                # Mutations on a nonempty edge with nontrivial fold and matching
                # ranks expose wrong population placement without shape failures.
                chosen=next(i for i,e in enumerate(links) if np.linalg.norm(rho[e.source])>0 and np.linalg.norm(rho[e.target])>0 and (qz==0 or np.any(e.coupling.edge.reciprocal_fold)))
                i=chosen; e=links[i]; xx=x[i]; v=drives[i]
                compare(label+'_lost_source_population',-1j*(h[e.target]@xx-xx@h[e.source]+v-rho[e.target]@v),dx[i],mutation=True)
                compare(label+'_nonHermitian_source_conjugation',-1j*(h[e.target]@xx.conj()-xx.conj()@h[e.source]+v@rho[e.source]-rho[e.target]@v),dx[i],mutation=True)
                compare(label+'_wrong_drive_adjoint',-1j*(h[e.target]@xx-xx@h[e.source]+v.conj()@rho[e.source]-rho[e.target]@v.conj()),dx[i],mutation=True)
                compare(label+'_wrong_target_population',-1j*(h[e.target]@xx-xx@h[e.source]+v@rho[e.source]-v),dx[i],mutation=True)
                if qz:
                    wrong=waves[e.target].conj().T@((phases[i].conj()*(n@kq.T)[:,sectors[e.source].density_component])[:,None]*waves[e.source])/nr+np.einsum('a,aij->ij',tangent.u,gm[i])
                    compare(label+'_wrong_fold',wrong,v,mutation=True)
                    rejects(label+'_swapped_endpoints',lambda:PreparedBlochResponse(model,[replace(links[0],source=links[0].target,target=links[0].source)]+links[1:],kernel=response.kernel,curvature=cq,authority_text='mutation'))
                    bad=replace(e.coupling.edge,q_fractional=np.array([0.,0.,-qz]))
                    rejects(label+'_wrong_q',lambda:bad.validate(inputs.reference))
                if qz==0:
                    ordinary=model.tangent(tangent.ket,tangent.bra,tangent.u,tangent.p)
                    ordinary_rhs=model.tangent_rhs(ev,ordinary)
                    for i in range(len(links)):
                        compare(label+'_Gamma_ket'+str(i),dot.ket[i],ordinary_rhs.ket[i])
                        compare(label+'_Gamma_bra'+str(i),dot.bra[i],ordinary_rhs.bra[i])
                    compare(label+'_Gamma_force',dot.p,ordinary_rhs.p)
            rejects(f'q{qz}_missing_source',lambda:PreparedBlochResponse(model,links[:-1],kernel=response.kernel,curvature=cq,authority_text='mutation'))
            rejects(f'q{qz}_duplicate_target',lambda:PreparedBlochResponse(model,[replace(links[0],target=links[1].target)]+links[1:],kernel=response.kernel,curvature=cq,authority_text='mutation'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result',type=Path,required=True)
    args=parser.parse_args()
    # Fail closed on changed original authorities before consuming the fixture.
    from substrate_framework.qe_saved_native import read_bound_bytes
    registration=json.loads((HERE/'finiteq-response-registration.json').read_text())
    for receipt in registration['source_bindings']+registration['additional_original_bytes']:
        read_bound_bytes(receipt)
    symbolic_scalar()
    finite_case()
    source=Path(__file__).read_bytes()
    report=dict(scope='independent finite scalar algebra/operator oracle; not material dynamics',passed=True,
        checks=CHECKS,source=dict(path=str(Path(__file__).resolve()),bytes=str(len(source)),sha256=hashlib.sha256(source).hexdigest()),
        limitations=['No actual material trajectory or H0 stationarity claim','Finite synthetic Fourier scalar, no cutoff/thermodynamic convergence','No material populations inferred from MV labels'])
    with args.result.open('x') as stream: json.dump(report,stream,indent=2,allow_nan=False); stream.write('\n')
    print(json.dumps(dict(passed=True,checks=len(CHECKS),result=str(args.result))))


if __name__=='__main__':
    main()

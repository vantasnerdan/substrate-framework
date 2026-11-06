"""Small independent A32 Hamiltonian/CAR oracle; no material inputs or actions.

Run explicitly with Python after both A32 implementation slices have settled.
"""
from types import SimpleNamespace
import numpy as np
from scipy.linalg import expm
from hybrid_gamma_model import HybridGammaModel


def close(actual, expected, label, tol=2e-9):
    error = float(np.max(np.abs(actual-expected), initial=0))
    if error > tol:
        raise AssertionError(f'{label}: {error} > {tol}')
    return error


def fixture():
    rng = np.random.default_rng(32017)
    def raw(shape):
        return rng.normal(size=shape)+1j*rng.normal(size=shape)
    def herm(shape):
        v = raw(shape)
        return (v+v.swapaxes(-1,-2).conj())/2
    ns,m,b,d,f,n = 3,2,4,2,1,4
    nx = 2*d+2*f
    frame = raw((b,3))/4
    gram = frame.conj().T@frame
    density = np.diag([.22,.37,.51])
    rho = frame@density@frame.conj().T
    q,_ = np.linalg.qr(raw((b,b)))
    q = q[:,:m]
    r0 = np.asarray([q.conj().T@rho@q, .23*np.eye(m), .61*np.eye(m)],complex)
    ionic = herm((ns,d,m,m))/7
    fields = herm((ns,f,m,m))/8
    features = herm((ns,2,m,m))/6
    block = SimpleNamespace(ionic=ionic,ionic_reverse=ionic.swapaxes(-1,-2),
        density=features.swapaxes(-1,-2),potential=np.moveaxis(features,1,-1),
        kernel=np.asarray([[.31,.07],[.07,.28]]))
    j = np.zeros((nx,nx)); j[:d,d:2*d]=np.eye(d); j[d:2*d,:d]=-np.eye(d)
    j[2*d,2*d+1]=1; j[2*d+1,2*d]=-1
    hmech = np.diag(np.linspace(.7,1.2,nx))
    w = .4*np.eye(nx)+.5j*j
    mechanical = SimpleNamespace(generators=[j@hmech],j=[j],linear=np.arange(nx)/70,w0=[w],h=[hmech])
    c = SimpleNamespace(n=n,sqrt_n=np.sqrt(n),ns=ns,m=m,d=d,f=f,rho0=r0,
        h0=herm((ns,m,m))/5,field_modes=fields,
        field_reverse_modes=fields.swapaxes(-1,-2),mechanical=mechanical,
        mechanical_sizes=(nx,),q=lambda index:block)
    hb = herm((b,b))/3+.003*raw((b,b))
    # Deliberately non-Hermitian raw legs: projection must use BOTH supplied legs.
    forces = raw((2,b,b))/9
    adjoints = forces.swapaxes(-1,-2).conj()+raw((2,b,b))/30
    bank = dict(projected_frame=frame,gram=gram,density_coefficients=density,
        projected_density=rho,projected_active_coordinates=q,hamiltonian=hb,
        hamiltonian_adjoint=hb.conj().T,
        force_mode0=forces[0],force_mode1=forces[1],
        adjoint_force_mode0=adjoints[0],adjoint_force_mode1=adjoints[1])
    directions = raw((d,2))/3
    initial_x = raw((nx,)).real/20
    model = HybridGammaModel(c,initial_x,bank,directions)
    return model,c,bank,directions,rng


def independent_physical(model,c,bank,u,time,state):
    """Literal physical Gamma Hamiltonian and force; no model math helpers."""
    nx=model.nx; m=c.m; b=len(bank['hamiltonian']); root=np.sqrt(c.n)
    x=state[:nx]; other=state[nx:nx+(c.ns-1)*m*m].reshape(c.ns-1,m,m)
    ri=state[-b*b:].reshape(b,b)
    # This eigensystem is independently reconstructed from literal raw H.
    eps,v=np.linalg.eig(bank['hamiltonian']); vi=np.linalg.inv(v)
    a=expm(-1j*bank['hamiltonian']*time)@v
    ai=vi@expm(1j*bank['hamiltonian']*time)
    rb=a@ri@ai; q=bank['projected_active_coordinates']
    r=np.concatenate(((q.conj().T@rb@q)[None],other))
    data=c.q(0)
    def active(xx,rr):
        density=np.einsum('srab,sab->r',data.density,rr-c.rho0)/c.n
        return (c.h0+np.einsum('siab,i->sab',data.ionic,xx[:c.d])/root
            +np.einsum('sabr,r->sab',data.potential,data.kernel@density)
            +np.einsum('sfab,f->sab',c.field_modes,xx[2*c.d:2*c.d+c.f]))
    rinitial=c.rho0.copy(); rinitial[0]=q.conj().T@bank['projected_density']@q
    hi=active(model.initial[:nx],rinitial)[0]
    ha=active(x,r)
    g=data.ionic[0]
    fold=-np.einsum('im,iab->mab',u,g)
    df=np.asarray([bank['force_mode0'],bank['force_mode1']])-q@fold@q.conj().T
    da=np.asarray([bank['adjoint_force_mode0'],bank['adjoint_force_mode1']])-q@fold.swapaxes(-1,-2).conj()@q.conj().T
    extra=-.5*(np.einsum('im,i,mab->ab',u.conj(),x[:c.d],df)
                 +np.einsum('im,i,mab->ab',u,x[:c.d],da))/root
    hb=bank['hamiltonian']+q@(ha[0]-hi)@q.conj().T+extra
    grad=np.zeros(nx,complex)
    grad[:c.d]=np.einsum('siab,sab->i',data.ionic_reverse,r-c.rho0)/root
    grad[2*c.d:2*c.d+c.f]=np.einsum('sfab,sab->f',c.field_reverse_modes,r-c.rho0)
    drb=rb-bank['projected_density']
    grad[:c.d]-=.5*(u.conj()@np.einsum('mab,ba->m',df,drb)
                         +u@np.einsum('mab,ba->m',da,drb))/root
    dx=c.mechanical.generators[0]@x+c.mechanical.j[0]@(c.mechanical.linear+grad)
    rd=-1j*(ha[1:]@other-other@ha[1:])
    rbd=-1j*(hb@rb-rb@hb)
    # Remove only literal H0, not current H or a diagonalized density.
    rid=ai@(rbd+1j*(bank['hamiltonian']@rb-rb@bank['hamiltonian']))@a
    return np.concatenate((dx,rd.ravel(),rid.ravel())),(hb,rb,df,da)


def finite_hamiltonian(model,c,bank,u,state,time):
    """Scalar finite Hamiltonian, differentiated independently below."""
    b=len(bank['hamiltonian']); nx=model.nx; root=np.sqrt(c.n)
    eps,v=np.linalg.eig(bank['hamiltonian']); vi=np.linalg.inv(v)
    a=v@np.diag(np.exp(-1j*eps*time)); ai=np.diag(np.exp(1j*eps*time))@vi
    x=state[:nx]; other=state[nx:-b*b].reshape(c.ns-1,c.m,c.m)
    rb=a@state[-b*b:].reshape(b,b)@ai; q=bank['projected_active_coordinates']
    active=np.concatenate(((q.conj().T@rb@q)[None],other))
    delta=active-c.rho0; data=c.q(0)
    density=np.einsum('srab,sab->r',data.density,delta)/root
    energy=.5*x@c.mechanical.h[0]@x+c.mechanical.linear@x
    energy+=np.einsum('sab,sba->',c.h0[1:],other)+np.trace(bank['hamiltonian']@rb)
    ionic=np.einsum('siab,i->sab',data.ionic,x[:c.d])/root
    field=np.einsum('sfab,f->sab',c.field_modes,x[2*c.d:2*c.d+c.f])
    energy+=np.einsum('sab,sba->',ionic+field,delta)+.5*density@data.kernel@density
    initial_x=model.initial[:nx]
    initial_active=c.rho0.copy(); initial_active[0]=q.conj().T@bank['projected_density']@q
    initial_density=np.einsum('srab,sab->r',data.density,initial_active-c.rho0)/root
    initial_controls=(np.einsum('iab,i->ab',data.ionic[0],initial_x[:c.d])/root
        +np.einsum('abr,r->ab',data.potential[0],data.kernel@initial_density)/root
        +np.einsum('fab,f->ab',c.field_modes[0],initial_x[2*c.d:2*c.d+c.f]))
    energy-=np.trace((q@initial_controls@q.conj().T)@rb)
    old=-np.einsum('im,iab->mab',u,data.ionic[0])
    forward=np.asarray([bank['force_mode0'],bank['force_mode1']])-q@old@q.conj().T
    adjoint=np.asarray([bank['adjoint_force_mode0'],bank['adjoint_force_mode1']])-q@old.swapaxes(-1,-2).conj()@q.conj().T
    extra=-.5*(np.einsum('im,i,mab->ab',u.conj(),x[:c.d],forward)
                 +np.einsum('im,i,mab->ab',u,x[:c.d],adjoint))/root
    return energy+np.trace(extra@(rb-bank['projected_density']))


def hamiltonian_gradient(model,c,bank,u,state,time):
    """Full coordinate gradient; no analytic force or commutator reused."""
    step=2e-6
    gradient=np.empty(model.size,complex)
    for coordinate in range(model.size):
        direction=np.zeros(model.size); direction[coordinate]=step
        gradient[coordinate]=(finite_hamiltonian(model,c,bank,u,state+direction,time)
            -finite_hamiltonian(model,c,bank,u,state-direction,time))/(2*step)
    return gradient


def hamiltonian_and_derivative():
    model,c,bank,u,rng=fixture()
    state=model.initial.copy()+.01*(rng.normal(size=model.size)+1j*rng.normal(size=model.size))
    t=.43
    expected,physical=independent_physical(model,c,bank,u,t,state)
    close(model.mean_rhs(t,state),expected,'literal global mean / interaction picture')
    gradient=hamiltonian_gradient(model,c,bank,u,state,t)
    b=len(bank['hamiltonian']); nx=model.nx
    other=state[nx:-b*b].reshape(c.ns-1,c.m,c.m)
    rb_i=state[-b*b:].reshape(b,b)
    h_other=gradient[nx:-b*b].reshape(c.ns-1,c.m,c.m).swapaxes(-1,-2)
    eps,_=np.linalg.eig(bank['hamiltonian'])
    h_interaction=gradient[-b*b:].reshape(b,b).T-np.diag(eps)
    generated=np.concatenate((c.mechanical.j[0]@gradient[:nx],
        (-1j*(h_other@other-other@h_other)).ravel(),
        (-1j*(h_interaction@rb_i-rb_i@h_interaction)).ravel()))
    close(model.mean_rhs(t,state),generated,'independent scalar Hamiltonian gradient',3e-8)
    values=rng.normal(size=(model.size,3))+1j*rng.normal(size=(model.size,3))
    step=2e-6
    derivatives=np.column_stack([(independent_physical(model,c,bank,u,t,state+step*v)[0]
        -independent_physical(model,c,bank,u,t,state-step*v)[0])/(2*step) for v in values.T])
    close(model.tangent(t,state,values),derivatives,'canonical exact tangent',2e-8)
    dual=rng.normal(size=(model.size,2))+1j*rng.normal(size=(model.size,2))
    close(dual.T@model.tangent(t,state,values),model.transpose(t,state,dual).T@values,
          'nonconjugating batched transpose',2e-9)
    close(model.tangent(t,state,values[:,0]),model.tangent(t,state,values)[:,0],'vector tangent')
    close(model.transpose(t,state,dual[:,0]),model.transpose(t,state,dual)[:,0],'vector transpose')
    # The same numerical predicate must reject each plausible changed model.
    def reject(actual,correct,label,tolerance=1e-8):
        error=float(np.max(np.abs(actual-correct),initial=0))
        if error<=tolerance:
            raise AssertionError(label+' mutation escaped the owning predicate')
        return error
    witnesses={}
    _,rb,df,da=physical
    dr=rb-bank['projected_density']
    force=-.5*(u.conj()@np.einsum('mab,ba->m',df,dr)+u@np.einsum('mab,ba->m',da,dr))/np.sqrt(c.n)
    omitted=expected.copy()
    mechanical_force=np.zeros(model.nx,complex); mechanical_force[:c.d]=force
    omitted[:model.nx]-=c.mechanical.j[0]@mechanical_force
    witnesses['omitted_backreaction']=reject(omitted,expected,'omitted backreaction')
    data=c.q(0)
    doubled=SimpleNamespace(**vars(data))
    doubled.density=data.density.copy(); doubled.density[0]*=2
    doubled.ionic_reverse=data.ionic_reverse.copy(); doubled.ionic_reverse[0]*=2
    cc=SimpleNamespace(**vars(c)); cc.q=lambda index:doubled
    cc.field_reverse_modes=c.field_reverse_modes.copy(); cc.field_reverse_modes[0]*=2
    duplicated=HybridGammaModel(cc,model.initial[:model.nx],bank,u)
    witnesses['duplicated_source']=reject(duplicated.mean_rhs(t,state),expected,'duplicated source contribution')
    cc=SimpleNamespace(**vars(c)); cc.n=4*c.n
    wrong_root=HybridGammaModel(cc,model.initial[:model.nx],bank,u)
    witnesses['wrong_rootN']=reject(wrong_root.mean_rhs(t,state),expected,'wrong rootN')
    hermitian_dual=model.transpose(t,state,dual.conj()).conj()
    witnesses['conjugated_transpose']=reject(hermitian_dual.T@values,
        dual.T@model.tangent(t,state,values),'conjugated transpose')
    altered=dict(bank)
    for mode in range(2):
        average=(bank[f'force_mode{mode}']+bank[f'adjoint_force_mode{mode}'])/2
        altered[f'force_mode{mode}']=average
        altered[f'adjoint_force_mode{mode}']=average
    projected=HybridGammaModel(c,model.initial[:model.nx],altered,u)
    witnesses['extra_quadrature']=reject(projected.mean_rhs(t,state),expected,'extra quadrature/leg averaging')
    altered=dict(bank); altered['adjoint_force_mode0']=bank['adjoint_force_mode0']+.13j*np.eye(len(rb))
    mutation=HybridGammaModel(c,model.initial[:model.nx],altered,u)
    witnesses['discarded_adjoint']=reject(mutation.mean_rhs(t,state),expected,'discarded supplied adjoint leg')
    witnesses['adjoint_tangent']=reject(mutation.tangent(t,state,values),
        model.tangent(t,state,values),'discarded adjoint tangent')
    return dict(state_size=model.size,mean_error=float(np.max(np.abs(model.mean_rhs(t,state)-expected))),
                hamiltonian_error=float(np.max(np.abs(model.mean_rhs(t,state)-generated))),
                mutation_separations=witnesses)


def annihilators(n):
    out=[]
    for j in range(n):
        a=np.zeros((2**n,2**n),complex)
        for bits in range(2**n):
            if bits&(1<<j): a[bits^(1<<j),bits]=(-1)**((bits&((1<<j)-1)).bit_count())
        out.append(a)
    return out


def fock_covariance(rho,left,right):
    """Full finite Fock trace, independent of matrix Wick contraction."""
    vals,v=np.linalg.eigh(rho)
    a=annihilators(len(rho))
    # d_k=sum_i conj(V_ik)c_i; preparation product in d eigenmodes.
    modes=[sum(v[i,k].conjugate()*a[i] for i in range(len(rho))) for k in range(len(rho))]
    identity=np.eye(2**len(rho)); preparation=identity.copy()
    for occupation,d in zip(vals,modes):
        number=d.conj().T@d
        preparation=preparation@(occupation*number+(1-occupation)*(identity-number))
    def second(o):
        return sum(o[i,j]*(a[i].conj().T@a[j]) for i in range(len(rho)) for j in range(len(rho)))
    l,r=second(left),second(right)
    return np.trace(preparation@l@r)-np.trace(preparation@l)*np.trace(preparation@r)


def original_car_and_fock():
    model,c,bank,u,rng=fixture()
    left=rng.normal(size=(model.size,2))+1j*rng.normal(size=(model.size,2))
    right=rng.normal(size=(model.size,2))+1j*rng.normal(size=(model.size,2))
    expected=left[:model.nx].T@c.mechanical.w0[0]@right[:model.nx]
    b=len(bank['hamiltonian']); eps,v=np.linalg.eig(bank['hamiltonian']); vi=np.linalg.inv(v)
    frame=bank['projected_frame']; density=bank['density_coefficients']
    bathleft=left[-b*b:].reshape(b,b,2); bathright=right[-b*b:].reshape(b,b,2)
    for p in range(2):
        for r in range(2):
            ol=v@bathleft[:,:,p].T@vi; ort=v@bathright[:,:,r].T@vi
            # Original raw frame/D CAR includes finite-basis complement.
            x=frame.conj().T@ol@ort@frame
            y=frame.conj().T@ol@frame; z=frame.conj().T@ort@frame
            car=np.trace(density@x-density@y@density@z)
            close(car,fock_covariance(bank['projected_density'],ol,ort),'complex bath CAR/Fock')
            expected[p,r]+=car
            for slot in range(1,c.ns):
                start=model.nx+(slot-1)*c.m*c.m
                oa=left[start:start+c.m*c.m,p].reshape(c.m,c.m).T
                ob=right[start:start+c.m*c.m,r].reshape(c.m,c.m).T
                expected[p,r]+=fock_covariance(c.rho0[slot],oa,ob)
    close(model.initial_ordered(left,right),expected,'same original total Gamma ordered covariance')
    close(model.initial_ordered(left[:,0],right[:,0]),expected[:1,:1],'vector initial covariance')
    close(model.probe(0)@model.initial,model.initial[0],'original coordinate probe')
    if model.size != model.nx+(c.ns-1)*c.m*c.m+b*b:
        raise AssertionError('old source0 coordinate duplicated')
    # With raw non-Hermitian H, V is NOT a unitary observable transform.
    actual,adjoint=model.preparation_observables(left)
    for p in range(2):
        physical=v@bathleft[:,:,p].T@vi
        close(actual[p],physical,'physical source0 observable')
        close(adjoint[p],physical.conj().T,'true physical adjoint')
    wrong_adjoint=v@bathleft[:,:,0].conj()@vi
    if np.max(np.abs(wrong_adjoint-adjoint[0]))<1e-7:
        raise AssertionError('nonunitary physical-adjoint mutation escaped')
    if np.max(np.abs(expected-expected/c.n))<1e-7:
        raise AssertionError('extra Fourier/cell factor mutation escaped')
    return dict(covariance_max_abs=float(np.max(np.abs(expected))))


def main():
    import json
    print(json.dumps(dict(schema='P254-A32-hybrid-model-oracle-v1',
        hamiltonian=hamiltonian_and_derivative(),covariance=original_car_and_fock(),
        scope='independent synthetic finite Hamiltonian/global-coordinate derivative and exact finite Fock/CAR; no material or physical calibration'),sort_keys=True))


if __name__=='__main__':
    main()

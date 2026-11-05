"""Independent small directed-commutator, Schur response and Fock oracle.

Only A42's new insertion, density interpolation/integration and same-source
substitution are tested. This is not a material numerical result or old suite.
"""
from pathlib import Path
import json
import sys
from types import SimpleNamespace
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / name) for name in ('0038', '0031', '0030', '0029')]
import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad_vec
from canonical_density_filter import density_insertion, build_canonical_density_filter
from density_integration import integrate_density_residual
from density_joint_readout import compose_same_source, substitution_readout
from reduced_electronic_correlations import ActiveBandCorrelations


def fixture():
    rng = np.random.default_rng(4206)
    z = lambda shape: rng.normal(size=shape)+1j*rng.normal(size=shape)
    ns, m, r, n = 2, 2, 2, 9
    j = np.array([[0, 1], [-1, 0]], complex)
    mechanics = SimpleNamespace(generators=[np.array([[0, 1], [-.7, -.05]], complex)]*2,
        j=[j, j], w0=[np.eye(2), np.eye(2)])
    rho = z((ns, m, m))*.08
    h = z((ns, m, m))*.04
    block = SimpleNamespace(density=z((ns, r, m, m))*.1,
        potential=z((ns, m, m, r))*.1, kernel=np.array([[.8+.1j,.2],[-.3j,.6]]),
        ionic=z((ns, 1, m, m))*.1, ionic_reverse=z((ns, 1, m, m))*.1,
        feature_g=np.array([[0,0,0],[1,0,0]]), feature_spin=np.array(['up','down']))
    c = SimpleNamespace(n=n, nq=2, ns=ns, m=m, d=1, f=0, mechanical=mechanics,
        target=np.array([[0,1],[1,0]]), rho0=rho, q=lambda q: block,
        minus=np.array([0,1]))
    c.nk=1;c.nspin=2;c.mesh=(1,1,1);c.mechanical_sizes=(2,2)
    c.sqrt_n=3.;c.h0=h;c.field_modes=np.empty((2,0,2,2),complex)
    corr = ActiveBandCorrelations(c, lambda t: (rho, h), rtol=2e-11, atol=2e-13)
    return corr, block, rho, h, z


def direct_density_drive(corr, block, rho, external, *, wrong_target=False):
    """Forward external n -> dsigma by literal entrywise matrix products."""
    drive = np.zeros((2, 2, 2, 1), complex)
    potential = block.kernel @ external
    for s in range(2):
        dh = sum(block.potential[s,:,:,r]*potential[r,0] for r in range(2))/3
        rt = rho[s if wrong_target else 1-s]
        drive[s,:,:,0] = -1j*(dh @ rho[s]-rt @ dh)
    return drive


def annihilation(mode, modes=3):
    a = np.zeros((2**modes, 2**modes), complex)
    for bits in range(2**modes):
        if bits & (1 << mode):
            target = bits ^ (1 << mode)
            a[target,bits] = (-1)**((bits & ((1 << mode)-1)).bit_count())
    return a


def physical_projection_oracle(z):
    """Independent explicit grid/PW sum tests the full-PW profile/action join."""
    from density_projection import BorrowedBathFrame,FourierDensityProjector
    from substrate_framework.qe_material_grid import PhysicalRhoFFT
    shape=(5,1,1);omega=7.
    points=np.arange(5)/5
    waves=np.exp(2j*np.pi*points[:,None]*np.arange(3)[None,:])
    fft=PhysicalRhoFFT.create(shape,np.array([[0,1,2],[0,0,0],[0,0,0]]),np.eye(3),omega)
    def action(field,vectors):
        return waves.conj().T@(field[:,None]*(waves@vectors))/5
    coupling=SimpleNamespace(factors=SimpleNamespace(fft=fft),potential_action=action)
    psi=z((3,2));inverse=np.array([[1,.2j],[.1,1]],complex)
    extra=z((3,1));frame=BorrowedBathFrame(psi,inverse,extra)
    forward,adjoint=FourierDensityProjector(coupling,frame,frame,ncells=9,block_size=1).project([1,0,0])
    b=np.concatenate((psi@inverse,extra),axis=1)
    profile_conjugate=np.exp(-2j*np.pi*points)/np.sqrt(omega)
    direct=np.zeros((3,3),complex)
    for i in range(3):
        for j in range(3):
            for point in range(5):
                left=sum(waves[point,k]*b[k,i] for k in range(3))
                right=sum(waves[point,k]*b[k,j] for k in range(3))
                direct[i,j]+=left.conjugate()*profile_conjugate[point]*right/15
    return float(np.max(abs(forward-direct))),float(np.max(abs(adjoint-direct.conj().T)))


def durable_reader_oracle(corr,z):
    """Real persisted numeric reader, not an in-memory producer-shaped stand-in."""
    import tempfile
    from durable_active_baseline import SavedMovingBaseline
    from substrate_framework.streamed_material_state import file_binding
    c=corr.c
    base=np.r_[c.rho0.ravel(),np.zeros(2,complex)]
    direction=z(base.shape)*.02
    times=np.array([0.,.31,.7])
    seventh=lambda t:np.polynomial.chebyshev.chebval(2*t/.7-1,[0]*7+[1])
    polynomial=lambda t:base+(seventh(t)+1)*direction
    coefficients=[]
    for left,right in zip(times[:-1],times[1:]):
        x=np.cos(np.pi*np.arange(8)/7)
        samples=np.stack([polynomial(left+(right-left)*(node+1)/2) for node in x])
        coefficients.append(np.polynomial.chebyshev.chebfit(x,samples,7))
    with tempfile.TemporaryDirectory() as directory:
        root=Path(directory);bindings={}
        for name,array in dict(boundary_times=times,boundary_states=np.stack([polynomial(t) for t in times]),
                               chebyshev=np.stack(coefficients)).items():
            path=root/(name+'.npy')
            with path.open('xb') as stream:np.save(stream,array,allow_pickle=False)
            bindings[name]=file_binding(path)
        receipt=dict(schema='P254-A42-durable-moving-active-baseline-v1',status='complete',arrays=bindings,
            coefficient_specification=dict(mesh=[1,1,1],bands=2,indexing=dict(slot_spin=[1,2])))
        path=root/'receipt.json';path.write_text(json.dumps(receipt))
        saved=SavedMovingBaseline(path,c)
        errors=[]
        for time in [0.,.123,.31,.511,.7]:
            state=saved.state(time)
            errors.append(float(np.max(abs(np.r_[state.rho.ravel(),state.x]-polynomial(time)))))
        # Product of moving degree7 rho and synthetic degree7 dual, restricted
        # on common intervals by the actual capture hook; high degrees survive.
        from canonical_density_filter import _DensityCapture
        moving=ActiveBandCorrelations(c,saved,rtol=2e-11,atol=2e-13)
        capture=_DensityCapture(moving,1,times)
        dual0=z((10,1));dual1=z((10,1))
        dense=lambda t:dual0+seventh(t)*dual1
        capture.segment(0.,.7,dense,dual0.shape)
        for left,right,coefficient in capture.segments:
            time=left+.39*(right-left)
            evaluated=np.polynomial.chebyshev.chebval(2*(time-left)/(right-left)-1,coefficient)
            errors.append(float(np.max(abs(evaluated-density_insertion(moving,1,time,dense(time))))))
        degree14=max(float(np.max(abs(coefficient[14]))) for _,_,coefficient in capture.segments)
    return max(errors),degree14


def run():
    corr, block, rho, h, z = fixture()
    dual = z((10,1)); external = z((2,1))
    computed = density_insertion(corr,1,.3,dual)
    zs = dual[2:].reshape(2,2,2,1)
    direct = np.sum(zs*direct_density_drive(corr,block,rho,external))
    pairing_error = float(abs(np.sum(computed*external)-direct))
    feature = np.linalg.solve(block.kernel.T,computed)
    mutations = dict(missing_feedback=float(abs(direct)),
        kernel_dagger=float(abs(np.sum((block.kernel.conj().T@feature)*external)-direct)),
        wrong_source_target=float(abs(np.sum(zs*direct_density_drive(corr,block,rho,external,wrong_target=True))-direct)),
        conjugating_pairing=float(abs(np.vdot(computed,external)-direct)),
        missing_second_rootN=float(abs(3*direct-direct)),
        duplicated_spin=float(abs(2*direct-direct)),
        missing_volume_profile=float(abs(np.sqrt(7)*direct-direct)))
    filt = build_canonical_density_filter(corr,1,.7,(0,),baseline_knots=[0.,.31,.7])
    # Full augmented forward Schur response, not the backward implementation.
    eye = np.eye(10,dtype=complex)
    l = np.column_stack([corr.tangent(1,0,eye[:,i]) for i in range(10)])
    b = np.zeros((10,2),complex)
    for r in range(2):
        e=np.zeros((2,1),complex);e[r]=1
        b[2:,r]=direct_density_drive(corr,block,rho,e).ravel()
    augmented=np.zeros((12,12),complex);augmented[:10,:10]=l;augmented[:10,10:]=b
    schur=expm(.7*augmented)[0,10:]
    integrated=quad_vec(lambda t:filt.evaluate_density(t)[:,0],0,.7,
                        points=filt.density_boundary_times)[0]
    schur_error=float(np.max(abs(integrated-schur)))
    # Actual spectral integrator on asymmetric complex-gap data, compare direct expm.
    hb=np.array([[.2+.01j,.13,.02j],[.04,-.1-.02j,.07],[.01j,.02,.4]],complex)
    energies,v=np.linalg.eig(hb);vi=np.linalg.inv(v)
    from filtered_quadratic_bath import _chebyshev_to_legendre
    conversion=_chebyshev_to_legendre(8)
    saved=dict(eigenvalues=energies,right_eigenvectors=v,inverse_right_eigenvectors=vi,
        energy_gaps=energies[:,None]-energies[None,:],chebyshev_to_legendre=conversion)
    vertices=[z((3,3))*.2 for r in range(2)]
    projector=np.diag([1.,0.,0.])
    provide=lambda r:(vertices[r],vertices[r].conj().T)
    spectral=integrate_density_residual(saved,provide,range(2),projector,
        filt.density_segment_bounds,filt.density_chebyshev)
    def integrand(t):
        coeff=filt.evaluate_density(t)[:,0]
        operator=sum(coeff[r]*(vertices[r]-projector@vertices[r]@projector) for r in range(2))
        return expm(1j*hb*t)@operator@expm(-1j*hb*t)
    exact=quad_vec(integrand,0,.7,points=filt.density_boundary_times)[0]
    spectral_error=float(np.max(abs(exact-spectral['operators'][0])))
    # Pure one-particle Gaussian in three modes. Direct 8-dimensional Fock state,
    # no Wick formula reused in this independent contraction.
    occupied=np.array([1,1j,.7-.2j]);occupied/=np.linalg.norm(occupied)
    state=np.zeros(8,complex)
    for i in range(3):state[1<<i]=occupied[i]
    a=[annihilation(i) for i in range(3)]
    lift=lambda o:sum(o[i,j]*(a[i].conj().T@a[j]) for i in range(3) for j in range(3))
    prep=z((3,3))*.1;ionic=z((3,3))*.1;field=z((3,3))*.1
    density=spectral['operators'][0]
    matrices,adjoints=compose_same_source(prep,prep.conj().T,ionic,ionic.conj().T,
        field,field.conj().T,density,density.conj().T)
    rho_f=np.outer(occupied,occupied.conj())
    wick=lambda o,p:np.trace(rho_f@o@p-rho_f@o@rho_f@p)
    fock=lambda o,p:(state.conj()@lift(o)@lift(p)@state
        -(state.conj()@lift(o)@state)*(state.conj()@lift(p)@state))
    moments=dict(covariance=[dict(left=x,right=y,value=wick(o,p))
        for x,o in matrices.items() for y,p in matrices.items()])
    sub=substitution_readout(moments,2.3+wick(prep,prep),wick(prep,prep))
    fock_error=max(float(abs(wick(o,p)-fock(o,p))) for o in matrices.values() for p in matrices.values())
    replacement_error=float(abs(sub['conditional_matched_active_total']-(2.3+fock(matrices['matched'],matrices['matched']))))
    mutations['cross_order_loss']=float(abs(sub['forbidden_independent_variance_addition_error']))
    mutations['bad_represented_replacement']=float(abs(wick(prep,prep)))
    # Exact symbolic changed algebra: the Chebyshev product can reach degree14.
    # Truncating to the old force degree7 is exposed independently of the ODE.
    import sympy as sp
    x=sp.symbols('x')
    product=sp.expand(sp.chebyshevt(7,x)**2)
    product_defect=sp.expand(product-(sp.chebyshevt(14,x)+1)/2)
    mutations['truncated_back_to_degree7']=float(abs((product-sp.Rational(1,2)).subs(x,sp.Rational(4,5))))
    # Force nonzero high orders through the ACTUAL higher-order integrator.
    high=np.zeros((1,15,2,1),complex);high[0,14,0,0]=1
    high_result=integrate_density_residual(saved,provide,range(2),projector,[[0,.7]],high)
    residual=vertices[0]-projector@vertices[0]@projector
    high_direct=quad_vec(lambda t:float(sp.chebyshevt(14,x).subs(x,2*t/.7-1))*
        (expm(1j*hb*t)@residual@expm(-1j*hb*t)),0,.7)[0]
    high_error=float(np.max(abs(high_result['operators'][0]-high_direct)))
    physical_error,physical_adjoint_error=physical_projection_oracle(z)
    # Directed unequal endpoint spectra: no assumption k+q has source's H.
    ht=np.array([[.17+.02j,.04],[-.03j,-.2]],complex)
    et,vt=np.linalg.eig(ht);vti=np.linalg.inv(vt)
    target=dict(eigenvalues=et,right_eigenvectors=vt,inverse_right_eigenvectors=vti,
        energy_gaps=et[:,None]-et[None,:],chebyshev_to_legendre=conversion)
    rectangular=z((3,2))*.1;pt=np.diag([1.,0.])
    rc=np.zeros((1,15,1,1),complex);rc[0,14,0,0]=1
    rr=integrate_density_residual(saved,lambda r:(rectangular,rectangular.conj().T),
        [0],projector,[[0,.7]],rc,target_saved=target,target_projector=pt)
    delta=rectangular-projector@rectangular@pt
    rd=quad_vec(lambda t:float(sp.chebyshevt(14,x).subs(x,2*t/.7-1))*
        (expm(1j*hb*t)@delta@expm(-1j*ht*t)),0,.7)[0]
    rectangular_error=float(np.max(abs(rr['operators'][0]-rd)))
    rda=quad_vec(lambda t:float(sp.chebyshevt(14,x).subs(x,2*t/.7-1))*
        (expm(1j*ht.conj().T*t)@delta.conj().T@expm(-1j*hb.conj().T*t)),0,.7)[0]
    rectangular_adjoint_error=float(np.max(abs(rr['adjoint_operators'][0]-rda)))
    durable_error,moving_high=durable_reader_oracle(corr,z)
    mutations['actual_moving_degree14_truncation']=moving_high
    errors=dict(directed_pairing=pairing_error,schur=schur_error,
        spectral=spectral_error,fock=fock_error,represented_replacement=replacement_error,
        cross_reconstruction=float(abs(sub['raw_all16_reconstruction'])),
        degree14_moments=high_error,symbolic_product=float(abs(product_defect)),
        fullPW_projection=physical_error,fullPW_projection_adjoint=physical_adjoint_error,
        directed_rectangular_spectral=rectangular_error,
        directed_rectangular_adjoint=rectangular_adjoint_error,durable_reader_and_product=durable_error)
    predicates={name:bool(error<2e-8) for name,error in errors.items()}
    predicates.update({name+'_exposed':bool(error>1e-9) for name,error in mutations.items()})
    return dict(schema='P254-A42-focused-density-oracle-v1',errors=errors,mutations=mutations,
        predicates=predicates,passed=bool(all(predicates.values())),
        product_scope='exact represented rho*z on union knots; no degree7 truncation',
        scope='new encoded density insertion/spectral residual/same-source substitution; not material closure')


if __name__=='__main__':
    result=run()
    print(json.dumps(result,indent=2,allow_nan=False))
    raise SystemExit(0 if result['passed'] else 1)

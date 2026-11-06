"""Actual A53 density/Gamma covariance -> common mean, independent U, Gamma S."""
import argparse,json,os,resource,time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from mean_feedback_functional import HERE,SOURCE,CommonMeanFunctional,CompactCubic,binding,save,authenticate
from material_electron_ion_inputs import load_material
from general_moving_reference import SavedMovingReference
from mean_feedback_oracle import boundary_oracle
from mean_feedback_reference import MeanFeedbackReference


def _node_value(t,times,values):
    if not times[0]<=t<=times[-1]:raise ValueError('forcing outside actual input nodes')
    i=min(int(np.searchsorted(times,t,side='right')-1),len(times)-2)
    z=(t-times[i])/(times[i+1]-times[i])
    return (1-z)*values[i]+z*values[i+1]


def integrate(model,times,rhos,gradients,covariance_forces,*,rtol,atol):
    """Genuine simultaneous x/U/S integration, separate at forcing breakpoints."""
    c=model.c;nx=model.nx;nu=c.ns*c.m*c.m
    initial_x=model.base.initial_baseline().x
    u0=np.broadcast_to(np.eye(c.m),(c.ns,c.m,c.m)).astype(complex).copy()
    s0=np.eye(nx,dtype=complex)
    y=np.concatenate((initial_x,u0.ravel(),s0.ravel()))
    saved=[y.copy()];work=[]
    def rhs(t,y):
        x=y[:nx];u=y[nx:nx+nu].reshape(c.ns,c.m,c.m);s=y[nx+nu:].reshape(nx,nx)
        rho=_node_value(t,times,rhos);gradient=_node_value(t,times,gradients);force=_node_value(t,times,covariance_forces)
        xdot=model.mean_rhs(x,gradient,force)
        udot=-1j*(model.hamiltonian(rho,x)@u)
        sdot=(c.mechanical.j[0]@model.curvature(x))@s
        return np.concatenate((xdot,udot.ravel(),sdot.ravel()))
    for left,right in zip(times[:-1],times[1:]):
        start=time.monotonic()
        sol=solve_ivp(rhs,(left,right),y,method='DOP853',rtol=rtol,atol=atol,t_eval=[right],dense_output=False)
        if not sol.success:raise RuntimeError(sol.message)
        y=sol.y[:,-1];saved.append(y.copy())
        work.append(dict(left=float(left),right=float(right),nfev=sol.nfev,seconds=time.monotonic()-start))
    saved=np.asarray(saved)
    return saved[:,:nx],saved[:,nx:nx+nu].reshape(len(times),c.ns,c.m,c.m),saved[:,nx+nu:].reshape(len(times),nx,nx),work


def produce_mean_feedback(dataset,baseline_catalog,output,*,rtol=2e-8,atol=2e-10):
    started=time.monotonic();output=Path(output).resolve();output.mkdir(parents=True,exist_ok=False)
    input_binding=binding(dataset);catalog_binding=binding(baseline_catalog)
    catalog=json.loads(Path(baseline_catalog).read_bytes())
    c,manifest,coefficient_binding=load_material()
    baseline=SavedMovingReference(c,authenticate(catalog['reference_bundle']),coefficient_binding)
    try:
        with np.load(dataset,allow_pickle=False) as inp:
            times=inp['times'];rho=-1j*np.asarray([inp['electron_less'][i,i] for i in range(len(times))])
            # theta(0)=1: preserve the raw retarded diagonal in greater=less+ret.
            wg=1j*np.asarray([inp['gamma_less'][i,i]+inp['gamma_ret'][i,i] for i in range(len(times))])
            source_rho=-1j*inp['source_electron_less0']
            source_wg=1j*(inp['source_gamma_less'][0,0]+inp['source_gamma_ret'][0,0])
        if (times.tolist()!=[0.,.5,1.] or c.n!=64 or rho.shape!=(3,128,12,12) or wg.shape!=(3,114,114)
                or not np.array_equal(rho[0],c.rho0) or not np.array_equal(source_rho,c.rho0)):
            raise ValueError('actual fine original density/preparation/time binding differs')
        source=CompactCubic.from_source(SOURCE);model=CommonMeanFunctional(c,source)
        ix=baseline.position_indices(0,'positions')
        initial_covariance_error=float(np.max(np.abs(wg[0]-c.mechanical.w0[0][np.ix_(ix,ix)])))
        if initial_covariance_error>2e-12:raise ValueError('dressed initial ordered covariance changed')
        oracle=boundary_oracle(model,rho[1],wg[1,:c.d,:c.d],baseline)
        save(output/'boundary-oracle.json',oracle)
        print(json.dumps(dict(stage='A57_NEW_BOUNDARY_ORACLE',findings=oracle)),flush=True)
        force=np.zeros((len(times),c.d),complex);qcost=[]
        for k in range(c.nq):
            start=time.monotonic();qs=source.indices[[0,k,int(source.minus[k])]]
            for a in range(len(times)):
                w=wg[a,:c.d,:c.d] if k==0 else 1j*baseline.boson(k,a,a)[0]
                force[a]-=.5*source.contract(qs,covariance=w)
            entry=dict(internal_q=k,seconds=time.monotonic()-start);qcost.append(entry)
            save(output/'node-force-checkpoint.json',dict(completed_internal_q=list(range(k+1)),timings=qcost,forces_real=force.real.tolist(),forces_imag=force.imag.tolist(),complete=k==63))
            print(json.dumps(dict(stage='A57_NODE_FORCE',**entry)),flush=True)
        gradients=np.asarray([model.electronic_gradient(r) for r in rho])
        mean,u,s,work=integrate(model,times,rho,gradients,force,rtol=rtol,atol=atol)
        archive=output/'mean-feedback.npz'
        np.savez(archive,times=times,mean_gamma=mean,U=u,S_gamma=s,rho0=c.rho0,W0_gamma=c.mechanical.w0[0],node_rho=rho,node_ordered_gamma_covariance=wg,node_cubic_covariance_force=force,node_electronic_gradient=gradients)
        j=c.mechanical.j[0]
        defects=dict(initial_mean_error=float(np.max(np.abs(mean[0]))),initial_U_error=float(np.max(np.abs(u[0]-np.eye(c.m)))),initial_S_error=float(np.max(np.abs(s[0]-np.eye(model.nx)))),
            initial_ordered_covariance_error=initial_covariance_error,
            raw_rho_hermiticity=[float(np.linalg.norm(r-r.conj().swapaxes(-1,-2))) for r in rho],
            raw_rho_particle_trace=[float(np.trace(r,axis1=-2,axis2=-1).real.sum()) for r in rho],
            U_unitarity=[float(np.linalg.norm(v@v.conj().swapaxes(-1,-2)-np.eye(c.m))) for v in u],
            Gamma_S_symplectic=[float(np.linalg.norm(v@j@v.T-j)) for v in s],
            Gamma_hessian_transpose=[float(np.linalg.norm(model.curvature(x)-model.curvature(x).T)) for x in mean],
            feedback_mean_norm=[float(np.linalg.norm(x)) for x in mean],
            cubic_covariance_force_norm=[float(np.linalg.norm(f)) for f in force],
            density_change_norm=[float(np.linalg.norm(r-c.rho0)) for r in rho],
            dressed_vs_source_gamma_covariance_norm=float(np.linalg.norm(wg[-1]-1j*baseline.ordered(0,2,2)[0])),
            forcing_midpoint_nonlinearity_rho=float(np.linalg.norm(rho[1]-(rho[0]+rho[2])/2)),
            forcing_midpoint_nonlinearity_cubic=float(np.linalg.norm(force[1]-(force[0]+force[2])/2)),
            covariance_source_initial_error=float(np.max(np.abs(source_wg-wg[0]))))
        source_paths=[HERE/n for n in ('mean_feedback_functional.py','mean_feedback_reference.py','produce_mean_feedback.py','mean_feedback_oracle.py')]
        source_paths.extend([HERE.parent/'0029/reduced_electronic_dynamics.py',HERE.parent/'0047/cubic_feedback.py',HERE.parent/'0051/cubic_contour.py',HERE.parent/'0052/general_moving_reference.py',HERE.parent/'0041/material_electron_ion_inputs.py'])
        spec=dict(schema='p254-common-mean-feedback-reference-v1',status='complete',archive=binding(archive),coefficient_manifest=coefficient_binding,input_dataset=input_binding,baseline_catalog=catalog_binding,
            physical_history_bindings=[catalog['reference_bundle'],source.source_binding,source.source_record['arrays']],source_bindings=[binding(p) for p in source_paths],
            times_inverse_ry=times.tolist(),method=dict(solver='DOP853',rtol=rtol,atol=atol,segments=work),oracle=binding(output/'boundary-oracle.json'),defects=defects,node_force_timings=qcost,
            approximation='Order1 interpolation ONLY of actual nodewise raw density/electronic gradient and all64 contracted cubic covariance force. No U/S or physical-kernel interpolation. Temporal forcing error unresolved; midpoint nonlinearity is a diagnostic, not an enclosure.',
            functional='A29 common scalar with original contacts/passive complement/mechanical source/active K plus A47 E3=<Phi uuu>/6. E3 homogeneous gradient=Phi000 mu^2/2+sumq Phi0q-q Wq/2; Gamma Hessian=Horiginal+Phi000 mu. No field cubic or XC counterterm. Fixed-rho phase Hessian; electronic bubble belongs to common Pi, not stacked contact.',
            reference_source_kind='bare-mechanical-with-original-contacts-plus-common-cubic-mean-curvature-at-held-dressed-rho-and-covariance; no-coupled-electronic-bubble-in-S',
            phase_layout='Original full228 (u111,p111,a3,Pi3); other-q coherent means exactly zero by homogeneous preparation/translation-invariant covariance.',
            bath_scope='Non-Gamma raw original prescribed Gaussian bath and physical calibration frozen; not closed all-q energy conservation or full selfconsistent iteration.',
            elapsed_seconds=time.monotonic()-started,maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            limits=['Single bounded feedback update on actual [0,.5,1] saved material input, not converged selfconsistency.','No clipping, projection, refill, Hermitian repair, mode dropping or Ufromrho.','No temporal accuracy, continuum calibration, XC doublecount cure, lifetime, rate, heat or T30 claim.'])
        save(output/'result.json',spec)
        ref=MeanFeedbackReference(c,output/'result.json',coefficient_binding)
        try:
            eg,el=ref.electron(2,1);bg,bl=ref.ordered(0,2,1,'positions');phase_g,phase_l=ref.ordered(0,2,1,'phase')
            probe=np.linspace(.1,.9,c.d+c.f)*(1+.37j)
            oldg,oldl=baseline.electron(2,1);oldbg,_=baseline.ordered(0,2,1)
            consumer=dict(electron_greater_norm=float(np.linalg.norm(eg)),electron_lesser_norm=float(np.linalg.norm(el)),fullGamma_ordered_greater_norm=float(np.linalg.norm(bg)),fullGamma_ordered_lesser_norm=float(np.linalg.norm(bl)),fullGamma_ordered_action_norm=float(np.linalg.norm(bg@probe)),fullGamma_phase_shape=list(phase_g.shape),
                electron_feedback_sensitivity=float(np.linalg.norm(el-oldl)),Gamma_feedback_sensitivity=float(np.linalg.norm(bg-oldbg)),
                Gamma_transition_feedback_sensitivity=float(np.linalg.norm(ref.S(0,1.)-baseline.S(0,1.))),
                electron_shape=list(el.shape),Gamma_positions_shape=list(bg.shape),u_integration_calls=ref.u_integration_calls,Gamma_integration_calls=ref.gamma_integration_calls,archive=binding(archive),result=binding(output/'result.json'))
            # Exercise original paired-q bath orientation without new propagation laws.
            ng,nl=ref.boson(1,2,1);og,ol=baseline.boson(1,2,1)
            consumer['original_nonGamma_bath_max_error']=float(max(np.max(np.abs(ng-og)),np.max(np.abs(nl-ol))))
            if min(consumer['electron_feedback_sensitivity'],consumer['Gamma_feedback_sensitivity'])<=0:raise AssertionError('feedback not exercised')
            save(output/'saved-consumer.json',consumer)
        finally:ref.close()
        result=dict(result=binding(output/'result.json'),archive=binding(archive),oracle=binding(output/'boundary-oracle.json'),consumer=binding(output/'saved-consumer.json'),findings=defects,consumer_findings=consumer,elapsed_seconds=time.monotonic()-started)
        save(output/'handoff.json',result)
        print(json.dumps(dict(stage='A57_ACTUAL_COMPLETE',**result)),flush=True)
        return result
    finally:baseline.close()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dataset',required=True);p.add_argument('--catalog',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    produce_mean_feedback(a.dataset,a.catalog,a.output)

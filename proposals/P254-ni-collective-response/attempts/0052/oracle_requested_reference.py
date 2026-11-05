"""One new finite noncommuting-H / paired-q arbitrary Gaussian field oracle."""
import sys
from types import SimpleNamespace
import numpy as np
from scipy.linalg import expm
from general_moving_reference import produce_u,GaussianTransitions,atomic_json,file_binding

def main(output):
    sx=np.array([[0,1],[1,0]],complex);sz=np.diag([1.,-1.]);omega=.83
    h0=.61*sx+.22*sz
    def h(t):
        r=expm(-1j*omega*t*sz/2);return r@h0@r.conj().T
    times=np.array([0.,.375,1.,1.625,2.]);u,work=produce_u(h,times,np.eye(2))
    exact=np.asarray([expm(-1j*omega*t*sz/2)@expm(-1j*(h0-omega*sz/2)*t) for t in times])
    defect=float(np.max(np.abs(u-exact)))
    wrong=float(np.max(np.abs(u[-1]-expm(-1j*times[-1]*h0))))
    rng=np.random.default_rng(52);d=2;f=1
    def J(n):return np.block([[np.zeros((n,n)),np.eye(n)],[-np.eye(n),np.zeros((n,n))]])
    # Gamma order is (u2,p2,a1,Pi1), not all-positions then momenta.
    jg=np.zeros((6,6));jg[:4,:4]=J(2);jg[4:,4:]=J(1)
    r=rng.normal(size=(6,6));hg=(r+r.T)*.08;ag=jg@hg
    r=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));hq=(r+r.conj().T)*.09
    aq=J(2)@hq;am=J(2)@hq.T
    vg=rng.normal(size=(6,6));vg=vg@vg.T+np.eye(6)
    vq=rng.normal(size=(4,4));vq=vq@vq.T+np.eye(4)
    wg=vg+.5j*jg;wq=vq+.5j*J(2)
    c=SimpleNamespace(d=d,f=f,minus=np.array([0,2,1]),mechanical=SimpleNamespace(generators=(ag,aq,am),w0=(wg,wq,wq)))
    provider=GaussianTransitions(c,times,cache_bytes=1024**2)
    records=[]
    for q in (0,1,2):
        qm=int(c.minus[q]);t,s=times[3],times[1]
        full=expm(t*c.mechanical.generators[q])@c.mechanical.w0[q]@expm(s*c.mechanical.generators[qm]).T
        for selection in ('phase','positions','ionic')+ (('field',) if q==0 else ()):
            ix=provider.position_indices(q,selection)
            value=provider.connected(q,t,s,selection)
            error=float(np.max(np.abs(value-full[np.ix_(ix,ix)])))
            assert error<2e-12
            records.append(dict(q=q,selection=selection,shape=list(value.shape),maxabs_error=error))
    calls=provider.expm_calls
    provider.connected(1,times[3],times[1]);assert provider.expm_calls==calls
    gg,gl=provider.ordered(1,3,1,'phase')
    reverse=expm(times[1]*am)@wq@expm(times[3]*aq).T
    reverse_error=float(np.max(np.abs(gl+1j*reverse.T)))
    dagger_gap=float(np.linalg.norm(expm(times[3]*aq)@wq@expm(times[1]*am).T-expm(times[3]*aq)@wq@expm(times[1]*am).conj().T))
    assert defect<2e-8 and wrong>.1 and reverse_error<2e-12 and dagger_gap>.1
    atomic_json(output,dict(schema='p254-a52-new-finite-provider-oracle-v1',status='passed',U_exact_rotating_frame_maxabs=defect,wrong_frozen_H_gap=wrong,paired_q_checks=records,reverse_transpose_maxabs=reverse_error,wrong_dagger_gap=dagger_gap,cache_expm_calls=calls,cache_repeated_request_extra_expm=provider.expm_calls-calls,fractional_times=times.tolist(),integration_work=work,scope='new generalized provider only; no old proof replay or global accuracy bound'))
    print(file_binding(output),flush=True)
if __name__=='__main__':main(sys.argv[1])

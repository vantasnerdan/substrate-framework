"""Time-ordered U and original constant bare Gaussian S; never density factors."""
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from moving_reference_inputs import atomic_json, file_binding, authenticate


def integrate_u(hamiltonian, left, right, initial, *, rtol=2e-8, atol=2e-10):
    shape=initial.shape
    result=solve_ivp(lambda t,y: (-1j*(hamiltonian(t)@y.reshape(shape))).ravel(),
        (float(left),float(right)),initial.ravel(),method='DOP853',rtol=rtol,atol=atol)
    if not result.success:
        raise RuntimeError(result.message)
    return result.y[:,-1].reshape(shape), dict(nfev=int(result.nfev),
        accepted_steps=len(result.t)-1, last_time=float(result.t[-1]))


class AffineSavedHamiltonian:
    """Ephemeral exact-degree interpolation of ACTUAL baseline_hamiltonian.

    Saved state is degree7 on each admitted A42 segment and its existing h is
    affine in rho/x. Eight Lobatto evaluations therefore determine its degree7
    h on that same segment, up to floating arithmetic. This is not a new mean,
    basis or physical Hamiltonian. No derived h bank is persisted.
    """
    def __init__(self, baseline):
        self.baseline=baseline
        self.times=baseline.times
        nodes=np.cos(np.pi*np.arange(8)/7)
        inverse=np.linalg.inv(np.polynomial.chebyshev.chebvander(nodes,7))
        self.coefficients=[]
        started=perf_counter()
        for left,right in zip(self.times[:-1],self.times[1:]):
            ts=left+(right-left)*(nodes+1)/2
            ts[0],ts[-1]=right,left
            samples=np.asarray([baseline(float(t))[1] for t in ts])
            self.coefficients.append(np.einsum('ij,jsab->isab',inverse,samples))
        self.construction_seconds=perf_counter()-started

    def __call__(self,time):
        index=min(len(self.times)-2,max(0,int(np.searchsorted(self.times,time,side='right'))-1))
        left,right=self.times[index:index+2]
        return np.polynomial.chebyshev.chebval(2*(time-left)/(right-left)-1,self.coefficients[index])


def checkpoint_u(directory, time, u, bindings, work):
    directory=Path(directory)
    path=directory/('U-boundary-%02d.npy'%time)
    with path.open('xb') as stream:
        np.save(stream,u,allow_pickle=False)
        stream.flush()
        import os
        os.fsync(stream.fileno())
    record=dict(schema='p254-a45-u-time-boundary-v1',time_inverse_ry=float(time),
        U=file_binding(path),history_bindings=bindings,work=work,
        recovery='actual U(t,0) at named boundary; DOP853 adaptive step/error history NOT preserved')
    receipt=path.with_suffix('.json');atomic_json(receipt,record)
    restored=np.load(authenticate(record['U']),allow_pickle=False)
    if not np.array_equal(restored,u):
        raise ValueError('persisted U boundary differs')
    return restored,file_binding(receipt)


def bare_sequence(generator,times):
    return np.asarray([np.eye(len(generator),dtype=complex) if t==0 else expm(t*generator)
                       for t in times])

"""Exact requested-time moving mean reference, with bounded raw Gaussian cache.

Producer integrates the existing saved Chebyshev mean Hamiltonian, never rho.
Saved reader consumes U bytes only. Fractions are admitted by produce_u, not
interpolated from U samples; the saved reader rejects times absent from its bank.
Coordinates are original (u,p[,a,Pi]); no instability or reciprocity repair.
"""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/n) for n in ('0045','0042','0041','0030','0029')]
from collections import OrderedDict
import json
import numpy as np
from scipy.linalg import expm
from moving_reference_inputs import inputs,atomic_json,file_binding,authenticate
from moving_reference_propagation import AffineSavedHamiltonian,integrate_u
from material_electron_ion_inputs import FrozenReference

SETTINGS=dict(method='DOP853',rtol=2e-8,atol=2e-10)

def requested_times(times,left=0.,right=30.):
    value=np.asarray(times,dtype=float)
    if value.ndim!=1 or not len(value) or not np.isfinite(value).all() or np.any(value<left) or np.any(value>right):
        raise ValueError('finite requested times within saved mean domain required')
    return value

def produce_u(hamiltonian,times,initial,*,rtol=2e-8,atol=2e-10):
    """Integrate genuinely at each distinct requested time; preserve caller order.

    Fractions use the same DOP853 Hamiltonian evaluation as integer times.
    Initial is U(0), not density. No new mean trajectory is computed.
    """
    times=requested_times(times)
    u=np.asarray(initial,dtype=complex).copy();previous=0.;values={0.:u.copy()};work=[]
    for t in np.unique(times):
        if t==0.:continue
        u,entry=integrate_u(hamiltonian,previous,float(t),u,rtol=rtol,atol=atol)
        work.append(dict(left=previous,right=float(t),**entry));values[float(t)]=u.copy();previous=float(t)
    return np.asarray([values[float(t)] for t in times]),work

class GaussianTransitions:
    """Direct expm per q/time, LRU byte bound; no per-vertex repeated expm.

    Cache limit counts retained S payload. transitions retains at most the current
    q/-q history pair (at31nodes <50MiB), in addition to this <=128MiB LRU.
    There is no persistent S bank or all64q history allocation.
    """
    def __init__(self,c,times,cache_bytes=128*1024**2):
        self.c,self.times=c,np.asarray(times,float)
        self.cache_bytes=int(cache_bytes);self.cache=OrderedDict();self.cache_size=0
        self.expm_calls=0;self.cache_hits=0;self._q=None;self._transitions=None
    def S(self,q,time):
        q=int(q);t=float(time)
        if not 0<=q<len(self.c.mechanical.generators) or not np.isfinite(t) or not 0<=t<=30:
            raise ValueError('original q and time in [0,30] required')
        key=(q,t)
        if key in self.cache:
            self.cache_hits+=1;self.cache.move_to_end(key);return self.cache[key]
        a=self.c.mechanical.generators[q]
        if t==0.:value=np.eye(len(a),dtype=complex)
        else:value=expm(t*a);self.expm_calls+=1
        value.flags.writeable=False
        while self.cache and self.cache_size+value.nbytes>self.cache_bytes:
            _,old=self.cache.popitem(last=False);self.cache_size-=old.nbytes
        if value.nbytes<=self.cache_bytes:
            self.cache[key]=value;self.cache_size+=value.nbytes
        return value
    def transitions(self,q):
        q=int(q)
        if self._q!=q:
            qm=int(self.c.minus[q])
            sq=np.asarray([self.S(q,t) for t in self.times])
            sm=sq if q==qm else np.asarray([self.S(qm,t) for t in self.times])
            self._q,self._transitions=q,(sq,sm)
        return self._transitions
    def position_indices(self,q,selection='ionic'):
        q=int(q);d=self.c.d;f=self.c.f if q==0 else 0
        if selection=='ionic':return np.arange(d)
        if selection=='positions':return np.r_[np.arange(d),np.arange(2*d,2*d+f)]
        if selection=='field':
            if q!=0:raise ValueError('global field exists at Gamma only')
            return np.arange(2*d,2*d+f)
        if selection=='phase':return np.arange(len(self.c.mechanical.generators[q]))
        raise ValueError('selection must be ionic/positions/field/phase')
    def connected(self,q,t,s,selection='positions'):
        qm=int(self.c.minus[q]);ix=self.position_indices(q,selection)
        # Select rows BEFORE multiplication; no unneeded full covariance copy.
        left=self.S(q,t)[ix];right=self.S(qm,s)[ix]
        return left@self.c.mechanical.w0[q]@right.T
    def ordered(self,q,a,b,selection='positions'):
        t,s=self.times[a],self.times[b];qm=int(self.c.minus[q])
        return (-1j*self.connected(q,t,s,selection),
                -1j*self.connected(qm,s,t,selection).T)
    def boson(self,q,a,b):return self.ordered(q,a,b,'ionic')

class SavedMovingReference(GaussianTransitions):
    """A41/A50 duck-compatible electron/boson/transitions provider.

    Reader performs zero U integration. Requested_times can select/reorder exact
    saved times, including fractional times produced by produce_u. Missing times
    fail closed; use the producer for additional genuine integrated nodes.
    """
    label='actual-requested-time-saved-U-original-bare-Gaussian-reference'
    def __init__(self,c,bundle,coefficient_binding,times=None,cache_bytes=128*1024**2):
        record_binding=file_binding(bundle);spec=json.loads(Path(bundle).read_bytes())
        if spec['schema']!='p254-a52-requested-time-reference-v1' or spec['coefficient_manifest']!=coefficient_binding:
            raise ValueError('reference/material binding mismatch')
        for entry in spec['physical_history_bindings']+spec['source_bindings']:
            authenticate(entry)
        saved_times=np.asarray(spec['times_inverse_ry'],float)
        all_u=np.load(authenticate(spec['U']),mmap_mode='r',allow_pickle=False)
        if all_u.shape!=(len(saved_times),c.ns,c.m,c.m) or not np.isfinite(all_u).all():
            raise ValueError('complete finite explicit-sector U required')
        selected=saved_times if times is None else requested_times(times)
        indices=[]
        for t in selected:
            found=np.flatnonzero(saved_times==t)
            if not len(found):raise ValueError('requested time absent: run genuine producer; no U interpolation')
            indices.append(int(found[0]))
        self.u=all_u if times is None else all_u[indices]
        zero=np.flatnonzero(selected==0.)
        if len(zero) and not np.array_equal(self.u[zero[0]],np.broadcast_to(np.eye(c.m),(c.ns,c.m,c.m))):
            raise ValueError('literal U0 changed')
        self.bindings=[record_binding,spec['U'],*spec['physical_history_bindings'],*spec['source_bindings']]
        self.spec=spec;self.u_integration_calls=0
        super().__init__(c,selected,cache_bytes)
    electron=FrozenReference.electron
    def close(self):
        self.cache.clear();self.cache_size=0;self._transitions=None

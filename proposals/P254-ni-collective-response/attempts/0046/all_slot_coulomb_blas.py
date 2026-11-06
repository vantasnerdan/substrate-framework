"""Exact local12^4 BLAS regrouping, static vertices shared across all16 pairs.

Only band-local tensors [externalbatch,64h,12,12,12,12] are formed.
No1536^4 tensor. Physical XX=DD/XD=DX license is the already-exposed
original source identity; no new Green symmetry, screening or approximation.
"""
import numpy as np
from all_slot_coulomb_contractions import KEYS,CHANNELS,zeros,validate_physical


def direct_q(factors,q,forward,reverse):
    c = factors.c
    qm = int(c.minus[q]); p = c.feature_partner[q][:len(factors.weights[q])]
    r = len(p); m = c.m; slots = np.arange(c.ns)
    a = c.target[qm,slots]
    y = c.density[qm][:,p].conj() @ forward[:,None]
    z = c.density[q][a].conj() @ reverse[a,None]
    yy = y.transpose(1,0,2,3).reshape(r,-1)
    zz = z.swapaxes(-1,-2).transpose(1,0,2,3).reshape(r,-1)
    w = factors.weights[q]/c.n
    kernel = (yy @ zz.T)*w[:,None]*w[None,:]
    b = c.target[qm,slots]
    x = c.density[q][b].conj() @ forward[b,None]
    right = c.density[qm][:,p].conj().reshape(c.ns,r,m*m)
    dressed = (kernel[None] @ right).reshape(c.ns,r,m,m)
    return x.transpose(0,2,1,3).reshape(c.ns,m,r*m) @ dressed.reshape(c.ns,r*m,m)


def local_exchange(left,right,fb,fc,reverse):
    """Unreduced band contraction L_ibac Fb_bd Fc_ce R_fa Rvertex_dfej.

    Four skinny BLAS matrix products, all ordinary/bilinear transposes.
    Both static local vertices are reused by the caller across time pairs.
    """
    u,m = left.shape[:2]
    t = left.transpose(0,1,3,4,2).reshape(u,m**3,m) @ fb
    t = t.reshape(u,m,m,m,m).transpose(0,1,2,4,3).reshape(u,m**3,m) @ fc
    t = t.reshape(u,m,m,m,m).transpose(0,1,3,4,2).reshape(u,m**3,m) @ reverse.swapaxes(-1,-2)
    t = t.reshape(u,m,m,m,m).transpose(0,1,2,4,3).reshape(u,m,m**3)
    return -(t @ right.reshape(u,m**3,m))


def static_exchange(factors,q,external):
    c = factors.c
    m,nh = c.m,c.nq
    qm = int(c.minus[q]); p = c.feature_partner[q][:len(factors.weights[q])]
    i = np.asarray(external)
    b = c.target[qm,i]
    cc = c.target[c.minus[:,None],i[None,:]].T
    aa = c.target[qm,cc]
    ng = len(p)
    x = c.density[q][b].conj()*(factors.weights[q]/c.n)[None,:,None,None]
    y = c.density[qm][cc][:,:,p].conj()
    left = (x.transpose(0,2,3,1)[:,None].reshape(len(i),1,m*m,ng)
            @ y.reshape(len(i),nh,ng,m*m)).reshape(len(i)*nh,m,m,m,m)
    right = np.empty((len(i),nh,m,m,m,m),complex)
    for h in range(nh):
        hm = int(c.minus[h]); ph = c.feature_partner[h][:len(factors.weights[h])]
        z = c.density[h][aa[:,h]].conj()*(factors.weights[h]/c.n)[None,:,None,None]
        w = c.density[hm][i][:,ph].conj()
        right[:,h] = (z.transpose(0,2,3,1).reshape(len(i),m*m,-1)
            @w.reshape(len(i),len(ph),m*m)).reshape(len(i),m,m,m,m)
    return left,right.reshape(len(i)*nh,m,m,m,m),np.repeat(b,nh),cc.ravel(),aa.ravel()


def all_time_q(factors,q,electrons,times,*,external_batch=4):
    c,nt = factors.c,len(times)
    delta = zeros(c,nt)
    for ti in range(nt):
        for si in range(nt):
            f,r = electrons[ti,si],electrons[si,ti]
            delta['greater_direct'][:,ti,si] = direct_q(factors,q,f[0],r[1])
            delta['lesser_direct'][:,ti,si] = direct_q(factors,q,f[1],r[0])
    for start in range(0,c.ns,external_batch):
        i = np.arange(start,min(start+external_batch,c.ns))
        left,right,b,cc,aa = static_exchange(factors,q,i)
        for ti in range(nt):
            for si in range(nt):
                f,r = electrons[ti,si],electrons[si,ti]
                for order,ff,rr in (('greater',f[0],r[1]),('lesser',f[1],r[0])):
                    value = local_exchange(left,right,ff[b],ff[cc],rr[aa])
                    delta[order+'_exchange'][i,ti,si] = value.reshape(len(i),c.nq,c.m,c.m).sum(axis=1)
    return delta


def evaluate_time_grid(factors,electrons,times,checkpoint,*,external_batch=4,after_commit=None):
    validate_physical(factors)
    c,nt = factors.c,len(times)
    arrays,done = checkpoint.restore()
    if arrays is None:
        arrays = zeros(c,nt)
    for q in range(done,c.nq):
        delta = all_time_q(factors,q,electrons,times,external_batch=external_batch)
        for order in ('greater','lesser'):
            for channel in ('direct','exchange'):
                arrays[order+'_'+channel] += delta[order+'_'+channel]
            arrays[order+'_total'][:] = arrays[order+'_direct']+arrays[order+'_exchange']
        mask = np.greater_equal(np.asarray(times)[:,None],np.asarray(times)[None,:])
        for channel in CHANNELS:
            arrays['retarded_'+channel][:] = (arrays['greater_'+channel]-arrays['lesser_'+channel])*mask[None,:,:,None,None]
        del delta
        checkpoint.commit(arrays,q+1)
        if after_commit is not None:
            after_commit(q+1)
    saved,completed = checkpoint.restore()
    if completed != c.nq:
        raise ValueError('optimized full numeric prefix incomplete')
    return saved,dict(new_q=c.nq-done,reused_q=done,completed_q=completed,
        external_slots=c.ns,ordered_time_pairs=nt*nt,
        local_static_vertex_pairs=(c.nq-done)*c.ns*c.nq,
        exchange_ordered_local_contractions=2*nt*nt*(c.nq-done)*c.ns*c.nq)


def bounds(batch=4):
    return dict(local_static_vertex_bytes=2*batch*64*12**4*16,
        local_band_transform_bytes_upper=3*batch*64*12**4*16,
        global_four_index_tensor=False,exchange_complex_multiply_adds_upper=32*64*128*64*4*12**5,
        static_vertex_complex_multiply_adds_upper=2*64*128*64*46*12**4,
        final_numeric_bytes=42467328,output_cap_bytes=268435456)

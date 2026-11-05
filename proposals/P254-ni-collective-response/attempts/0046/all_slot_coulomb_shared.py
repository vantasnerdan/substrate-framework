"""Exact shared-(q,c,timepair) bare exchange superoperator.

For b=i-q,c=i-h,a=c-q, dress X_q(b)=w_q L_q(b) F_b/N
and Y_q(c)=R_a L_-q(c) F_c ONCE, shared across every external i.
K_h(i)[d,j;f,e]=sum_s w_h L_h(a)[s,d,f] L_-h(i)[p(s),e,j]/N
is a static local144x144 map, reused across all32 greater/lesser pairs.
Exchange=-sum_{h,r,d} X_r[i,d]*(K Y_r)[d,j]. Ordinary complex
bilinear products only; no reverse inference or extra physical permutation.
"""
import numpy as np
from all_slot_coulomb_blas import direct_q
from all_slot_coulomb_contractions import CHANNELS,zeros,validate_physical


def exchange_superoperator(matrix,x,y):
    """matrix[u,d*j,f*e], x/y[u,time,r,band,band]. No conjugation."""
    u,nt,nr,m,_ = x.shape
    dressed = (y.reshape(u,nt*nr,m*m) @ matrix.swapaxes(-1,-2)).reshape(u,nt,nr*m,m)
    external = x.transpose(0,1,3,2,4).reshape(u,nt,m,nr*m)
    return -(external @ dressed)


def static_superoperators(factors,q,external):
    c,m = factors.c,factors.c.m
    i = np.asarray(external)
    qm = int(c.minus[q])
    cc = c.target[c.minus[:,None],i[None,:]].T
    aa = c.target[qm,cc]
    b = np.repeat(c.target[qm,i],c.nq)
    matrix = np.empty((len(i),c.nq,m*m,m*m),complex)
    for h in range(c.nq):
        hm = int(c.minus[h]); p = c.feature_partner[h][:len(factors.weights[h])]
        z = c.density[h][aa[:,h]].conj()*(factors.weights[h]/c.n)[None,:,None,None]
        w = c.density[hm][i][:,p].conj()
        raw = (z.transpose(0,2,3,1).reshape(len(i),m*m,len(p))
            @w.reshape(len(i),len(p),m*m)).reshape(len(i),m,m,m,m)
        matrix[:,h] = raw.transpose(0,1,4,2,3).reshape(len(i),m*m,m*m)
    return matrix.reshape(len(i)*c.nq,m*m,m*m),b,cc.ravel()


def shared_operands(factors,q,electrons,times):
    c,nt = factors.c,len(times)
    qm = int(c.minus[q]); p = c.feature_partner[q][:len(factors.weights[q])]
    aa = c.target[qm,np.arange(c.ns)]
    raw = c.density[qm][:,p].conj()
    lq = c.density[q].conj()
    x = np.empty((2*nt*nt,c.ns,len(p),c.m,c.m),complex)
    y = np.empty_like(x)
    labels = []
    for ti in range(nt):
        for si in range(nt):
            f,r = electrons[ti,si],electrons[si,ti]
            for order,ff,rr in (('greater',f[0],r[1]),('lesser',f[1],r[0])):
                k = len(labels)
                x[k] = (lq @ ff[:,None])*(factors.weights[q]/c.n)[None,:,None,None]
                y[k] = rr[aa,None] @ raw @ ff[:,None]
                labels.append((order,ti,si))
    return x,y,labels


def all_time_q(factors,q,electrons,times,*,external_batch=4,time_batch=8):
    c,nt = factors.c,len(times)
    delta = zeros(c,nt)
    for ti in range(nt):
        for si in range(nt):
            f,r = electrons[ti,si],electrons[si,ti]
            delta['greater_direct'][:,ti,si] = direct_q(factors,q,f[0],r[1])
            delta['lesser_direct'][:,ti,si] = direct_q(factors,q,f[1],r[0])
    x,y,labels = shared_operands(factors,q,electrons,times)
    for start in range(0,c.ns,external_batch):
        i = np.arange(start,min(start+external_batch,c.ns))
        matrix,b,cc = static_superoperators(factors,q,i)
        for begin in range(0,len(labels),time_batch):
            end = min(begin+time_batch,len(labels))
            # Gathered batches are finite; time8/external4 <=declared4GiB.
            xb = x[begin:end,b].transpose(1,0,2,3,4)
            yb = y[begin:end,cc].transpose(1,0,2,3,4)
            value = exchange_superoperator(matrix,xb,yb)
            value = value.reshape(len(i),c.nq,end-begin,c.m,c.m).sum(axis=1)
            for k,(order,ti,si) in enumerate(labels[begin:end]):
                delta[order+'_exchange'][i,ti,si] = value[:,k]
    return delta


def evaluate_time_grid(factors,electrons,times,checkpoint,*,external_batch=4,
                       time_batch=8,after_commit=None):
    validate_physical(factors)
    c,nt = factors.c,len(times)
    arrays,done = checkpoint.restore()
    if arrays is None:
        arrays = zeros(c,nt)
    owned_q = checkpoint.owned_q
    if checkpoint.global_nq != c.nq:
        raise ValueError('partition inner-h global momentum count differs')
    for prefix_index in range(done,len(owned_q)):
        q = owned_q[prefix_index]
        delta = all_time_q(factors,q,electrons,times,external_batch=external_batch,time_batch=time_batch)
        for order in ('greater','lesser'):
            for channel in ('direct','exchange'):
                arrays[order+'_'+channel] += delta[order+'_'+channel]
            arrays[order+'_total'][:] = arrays[order+'_direct']+arrays[order+'_exchange']
        mask = np.greater_equal(np.asarray(times)[:,None],np.asarray(times)[None,:])
        for channel in CHANNELS:
            arrays['retarded_'+channel][:] = (arrays['greater_'+channel]-arrays['lesser_'+channel])*mask[None,:,:,None,None]
        del delta
        checkpoint.commit(arrays,prefix_index+1)
        if after_commit is not None:
            after_commit(prefix_index+1)
    saved,completed = checkpoint.restore()
    if completed != len(owned_q):
        raise ValueError('shared-superoperator numeric prefix incomplete')
    return saved,dict(new_q=len(owned_q)-done,reused_q=done,completed_q=completed,
        owned_outer_q=owned_q,exchange_inner_h_count=c.nq,
        external_slots=c.ns,ordered_time_pairs=nt*nt,
        shared_q_c_time_operands=2*nt*nt*(len(owned_q)-done)*c.ns,
        static_superoperators=(len(owned_q)-done)*c.ns*c.nq,
        ordered_local_contractions=2*nt*nt*(len(owned_q)-done)*c.ns*c.nq)


def bounds(external_batch=4,time_batch=8):
    return dict(shared_x_y_bytes_upper=2*32*128*46*12*12*16,
        static_superoperator_bytes=external_batch*64*12**4*16,
        time_local_numeric_batches_bytes_upper=4*external_batch*64*time_batch*46*12*12*16,
        density_bytes_upper=868220928,all_numeric_accumulators_bytes_upper=3*42467328,
        exchange_complex_multiply_adds_upper=32*64*128*64*46*(12**4+12**3),
        global_four_index_allocated=False,final_numeric_bytes=42467328,output_cap_bytes=268435456)

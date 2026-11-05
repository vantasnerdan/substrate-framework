"""All-slot bare second Born; physical pair symmetry, NOT generic vertex symmetry.

Let L_q(b)=conj(D_q(b)), p_q the original opposite-feature permutation.
V_ia;bc=sum_r w_q[r]/N L_q(b)[r,i,b] L_-q(c)[p_q(r),a,c].
Commutation of the two scalar factors, p_-q p_q=I and
w_-q[p_q(r)]=w_q[r] imply V_ia;bc=V_ai;cb, without any
adjoint relation between density edges or any Green-leg symmetry.
In the COMPLETE b,c sums, swapping b<->c,d<->e then applying this
pair identity to the right vertex proves XX=DD and XD=DX. Thus the
original .5(DD+XX-DX-XD) is DD-DX. Both channels remain separate.
This license applies only to these physical factor banks; arbitrary generic
vertices must use A44's unreduced polynomial. No JVP/transpose API is exposed.
"""
import numpy as np

ORDERS = ('greater', 'lesser', 'retarded')
CHANNELS = ('direct', 'exchange', 'total')
KEYS = tuple(o+'_'+c for o in ORDERS for c in CHANNELS)


def validate_physical(factors):
    c = factors.c
    slots = np.arange(c.ns)
    for q in range(c.nq):
        qm = int(c.minus[q])
        p = np.asarray(c.feature_partner[q][:len(factors.weights[q])])
        pm = np.asarray(c.feature_partner[qm][:len(factors.weights[qm])])
        if (int(c.minus[qm]) != q or not np.array_equal(pm[p], np.arange(len(p)))
                or not np.allclose(factors.weights[q], factors.weights[qm][p], rtol=2e-14, atol=0)
                or not np.array_equal(c.target[qm, c.target[q]], slots)
                or not np.array_equal(c.slot_spin[c.target[q]], c.slot_spin)):
            raise ValueError('physical pair-commutation prerequisites differ')
    # Transfer group must commute: a=(i-h)-q=(i-q)-h.
    for q in range(c.nq):
        for h in range(c.nq):
            if not np.array_equal(c.target[q,c.target[h]], c.target[h,c.target[q]]):
                raise ValueError('noncommuting momentum translations')


def zeros(c, nt):
    return {k:np.zeros((c.ns,nt,nt,c.m,c.m), complex) for k in KEYS}


def ordered_q(factors, q, forward, reverse, *, external_batch=4):
    """One complete outer q, ALL external slots, both physical channels.

    Forward and reverse are arbitrary independent complex block matrices.
    Direct bubble is formed once, not once per external. Exchange is a
    feature-pair matrix chain batched over external slots; no 4-index vertex.
    """
    c = factors.c
    qm = int(c.minus[q])
    p = c.feature_partner[q][:len(factors.weights[q])]
    wq = factors.weights[q]/c.n
    # a=c-q for every internal slot c, including BOTH explicit spin sectors.
    slots = np.arange(c.ns)
    a = c.target[qm,slots]
    # L_-q(c) F_c L_q(a) R_a, traced in band space, shared by all i.
    y = c.density[qm][:,p].conj() @ forward[:,None]
    z = c.density[q][a].conj() @ reverse[a,None]
    bubble = np.einsum('urac,utca->rt', y, z, optimize=False)
    kernel = bubble*wq[:,None]*wq[None,:]
    direct = np.zeros((c.ns,c.m,c.m),complex)
    exchange = np.zeros_like(direct)
    for start in range(0,c.ns,external_batch):
        i = slots[start:start+external_batch]
        b = c.target[qm,i]
        x = c.density[q][b].conj() @ forward[b,None]
        right = c.density[qm][i][:,p].conj()
        # Sum features before band multiplication (small batch*ng*m*m).
        dressed = np.einsum('rt,utdj->urdj',kernel,right,optimize=False)
        direct[i] = np.einsum('urib,urbj->uij',x,dressed,optimize=False)
        x *= wq[None,:,None,None]
        for h in range(c.nq):
            hm = int(c.minus[h])
            ph = c.feature_partner[h][:len(factors.weights[h])]
            cc = c.target[hm,i]
            aa = c.target[qm,cc]
            z = c.density[h][aa].conj() @ reverse[aa,None]
            z *= (factors.weights[h]/c.n)[None,:,None,None]
            y = c.density[qm][cc][:,p].conj() @ forward[cc,None]
            right = c.density[hm][i][:,ph].conj()
            # Only large intermediate: [batch,ng_q,ng_h,m,m].
            chain = np.einsum('urib,usbf->ursif',x,z,optimize=False)
            reduced = np.einsum('ursif,urfc->usic',chain,y,optimize=False)
            exchange[i] -= np.einsum('usic,uscj->uij',reduced,right,optimize=False)
    return direct,exchange


def evaluate_time_grid(factors, electrons, times, checkpoint, *, external_batch=4,
                       after_commit=None):
    """q-major prefix; one q commits all external batches and all ordered pairs.

    electrons[(t,s)] owns (> , <) independently of electrons[(s,t)].
    Recovery discards at most ONE all-external q, never persists q archives.
    """
    c,nt = factors.c,len(times)
    if not 1 <= external_batch <= 8:
        raise ValueError('external_batch must be1..8 for the declared RAM bound')
    validate_physical(factors)
    arrays,done = checkpoint.restore()
    if arrays is None:
        arrays = zeros(c,nt)
    new_q = 0
    for q in range(done,c.nq):
        for ti,t in enumerate(times):
            for si,s in enumerate(times):
                forward,backward = electrons[ti,si],electrons[si,ti]
                for order,f,r in (('greater',forward[0],backward[1]),
                                  ('lesser',forward[1],backward[0])):
                    d,x = ordered_q(factors,q,f,r,external_batch=external_batch)
                    arrays[order+'_direct'][:,ti,si] += d
                    arrays[order+'_exchange'][:,ti,si] += x
        for order in ('greater','lesser'):
            arrays[order+'_total'][:] = arrays[order+'_direct']+arrays[order+'_exchange']
        mask = np.greater_equal(np.asarray(times)[:,None],np.asarray(times)[None,:])
        for channel in CHANNELS:
            arrays['retarded_'+channel][:] = (arrays['greater_'+channel]-arrays['lesser_'+channel])*mask[None,:,:,None,None]
        checkpoint.commit(arrays,q+1)
        new_q += 1
        if after_commit is not None:
            after_commit(q+1)
    # Genuine saved numerical consumer; no echo of unsaved RAM.
    saved,completed = checkpoint.restore()
    if completed != c.nq:
        raise ValueError('complete numeric prefix not retained')
    return saved,dict(new_q=new_q,reused_q=done,completed_q=completed,
        external_slots=c.ns,ordered_time_pairs=nt*nt,
        reduced_direct_bubbles=2*nt*nt*new_q,
        exchange_q_h_slot_terms=2*nt*nt*new_q*c.nq*c.ns)


def source_bounds(ns=128,nq=64,m=12,nt=4,features=46,batch=4):
    payload = ns*nt*nt*9*m*m*16
    return dict(final_numeric_bytes=payload,rotating_numeric_generations=2,
        transient_numeric_archives_including_final_and_stage_upper=4,
        output_bytes_upper=4*payload+8*1024**2,output_cap_bytes=256*1024**2,
        density_resident_bytes_upper=nq*ns*features*m*m*16,
        electron_legs_resident_bytes=nt*nt*2*ns*m*m*16,
        contraction_large_intermediate_bytes=batch*features**2*m*m*16,
        live_numeric_accumulators_bytes_upper=3*payload,
        conservative_algorithm_ram_bytes_upper=2*1024**3,
        direct_bubbles=2*nt*nt*nq,
        exchange_q_h_slot_terms=2*nt*nt*nq*nq*ns,
        complex_multiply_adds_upper=2*nt*nt*(nq*ns*(4*features*m**3+features**2*m**2)
            +nq*nq*ns*(2*features**2*m**3+4*features*m**3+features*m**2)),
        interruption_loss='at most current all-external q including all16 ordered pairs',
        ownership='ascending outerq prefix; every committed q owns all128 externals/all64 innerh/all16 pairs',
        physical_symmetry='XX=DD and XD=DX after complete dummy-index exchange; pair-commuting bare density vertices only',
        full_four_index_allocation=False,jvp_transpose_exposed=False,
        measured_material_speed=None,material_run_authorized=False)

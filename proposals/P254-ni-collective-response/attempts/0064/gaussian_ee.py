"""Exact bare-Gaussian EE time separation, with all internal modes retained.

A46 left L[i,b,a,c], right V[d,f,e,j] and independent routed U legs:
 Fb=Bt Bs, Fc=Ct Cs, R=As At, with Bt=Ub(t) Cb,
 Bs=Ub(s)^dagger, Ct=Uc(t) Cc, Cs=Uc(s)^dagger,
 As=Ua(s) Da, At=Ua(t)^dagger. No relation between L and V is used.
 LT[i,B,A,C]=sum_bac L[i,b,a,c] Bt[b,B] At[A,a] Ct[c,C].
 RS[B,A,C,j]=sum_dfe Bs[B,d] As[f,A] Cs[C,e] V[d,f,e,j].
 Exchange is minus the ordinary bilinear LT @ RS over ALL m^3 modes.
 Raw C,D=(empty,rho) or (rho,empty); the three Green phases are -i,+i.

Direct bubble uses Tr(Y Bt Bs Z As At)=Tr((At Y Bt)(Bs Z As)).
Feature/time Gram factors are blocked over internal slots; exchange factors
are blocked over inner h. No all-slot/time three-mode bank is allocated.
Only genuine Gaussian saved U is supported, not a dressed two-time bank.
"""
import numpy as np
from all_slot_coulomb_blas import static_exchange
from all_slot_coulomb_contractions import validate_physical, KEYS
from streamed_ee import Block, _pair_tiles


def _mode(tensor, matrix, axis):
    """Batched old->new mode multiplication, ordinary product (not adjoint)."""
    u,m=matrix.shape[:2]
    moved=np.moveaxis(tensor,axis,-1)
    value=(moved.reshape(u,-1,m) @ matrix).reshape(moved.shape)
    return np.moveaxis(value,-1,axis)


def _left_bank(left,u,b,cc,aa,occupations):
    nt,nlocal,m=len(u),len(b),u.shape[-1]
    bank=np.empty((2,nt,nlocal,m,m**3),complex)
    for order,C in enumerate(occupations):
        for t in range(nt):
            v=_mode(left,u[t,b]@C[b],2)
            v=_mode(v,u[t,aa].conj(),3)  # At.T, not At or an inferred reverse G
            v=_mode(v,u[t,cc]@C[cc],4)
            bank[order,t]=v.reshape(nlocal,m,m**3)
    return bank


def _right_bank(right,u,b,cc,aa,occupations):
    nt,nlocal,m=len(u),len(b),u.shape[-1]
    bank=np.empty((2,nt,nlocal,m**3,m),complex)
    for order,D in enumerate(occupations[::-1]):
        for s in range(nt):
            v=_mode(right,u[s,b].conj(),1)  # Bs.T
            v=_mode(v,u[s,aa]@D[aa],2)
            v=_mode(v,u[s,cc].conj(),3)  # Cs.T
            bank[order,s]=v.reshape(nlocal,m**3,m)
    return bank


def _direct_kernels(factors,q,u,occupations,internal_batch):
    """[order,t,s,r,v], raw two-Green phases multiply to +1."""
    c=factors.c;nt=len(u);qm=int(c.minus[q]);nr=len(factors.weights[q])
    p=c.feature_partner[q][:nr];m=c.m
    kernels=np.zeros((2,nt,nt,nr,nr),complex)
    peak=0
    for start in range(0,c.ns,internal_batch):
        cc=np.arange(start,min(start+internal_batch,c.ns));aa=c.target[qm,cc]
        ly=c.density[qm][cc][:,p].conj();lz=c.density[q][aa].conj()
        for order,C in enumerate(occupations):
            D=occupations[1-order]
            left=np.empty((nt,nr,len(cc),m,m),complex);right=np.empty_like(left)
            for t in range(nt):
                at=u[t,aa].conj().swapaxes(-1,-2)
                bs=u[t,cc].conj().swapaxes(-1,-2)
                left[t]=(at[:,None]@ly@(u[t,cc]@C[cc])[:,None]).transpose(1,0,2,3)
                right[t]=(bs[:,None]@lz@(u[t,aa]@D[aa])[:,None]).swapaxes(-1,-2).transpose(1,0,2,3)
            peak=max(peak,left.nbytes+right.nbytes)
            gram=left.reshape(nt*nr,-1)@right.reshape(nt*nr,-1).T
            kernels[order]+=gram.reshape(nt,nr,nt,nr).transpose(0,2,1,3)
    w=factors.weights[q]/c.n
    kernels*=w[None,None,None,:,None]*w[None,None,None,None,:]
    return kernels,peak


def stream_gaussian_ee(factors,reference,times,consume,*,owned_q=None,
                       owned_pairs=None,pair_batch=4,inner_batch=8,
                       internal_batch=8,source_bindings=()):
    """A60 Block/consumer ABI: additive q, all slots and nine channels.

    reference must be a genuine saved bare-Gaussian U provider with original
    c.rho0 and matching exact times. It is NOT a generic electron callback.
    Resident temporal exchange payload <=4*nt*inner_batch*m^4*16; direct
    kernels <=2*nt^2*ng_q^2*16. No rank truncation or stationarity is used.
    consume must persist/reduce synchronously, as for A60.
    """
    c=factors.c;times=np.asarray(times,float);u=reference.u
    if times.ndim!=1 or not len(times) or not np.isfinite(times).all():
        raise ValueError('nonempty finite exact saved grid required')
    if not np.array_equal(times,reference.times) or reference.c is not c:
        raise ValueError('same original material and exact saved U grid required')
    if u.shape!=(len(times),c.ns,c.m,c.m):
        raise ValueError('all explicit saved U slots and bands required')
    if any(type(v) is not int or v<=0 for v in (pair_batch,inner_batch,internal_batch)):
        raise ValueError('positive integer block sizes required')
    qs=tuple(range(c.nq)) if owned_q is None else tuple(owned_q)
    if len(set(qs))!=len(qs) or any(type(q) is not int or not 0<=q<c.nq for q in qs):
        raise ValueError('unique owned outer q required')
    pairs=tuple((a,b) for tile in _pair_tiles(len(times),pair_batch,owned_pairs) for a,b in tile)
    if len(set(pairs))!=len(pairs):raise ValueError('duplicate owned time pairs')
    validate_physical(factors)
    occupations=(np.eye(c.m)-c.rho0,c.rho0)
    manifest=dict(schema='p254-a64-gaussian-ee-v1',times=times.tolist(),
        owned_outer_q=qs,owned_pairs=pairs,external_slots=c.ns,bands=c.m,
        exchange_inner_h=c.nq,channels=KEYS,additive_axis='outer_q',
        internal_modes=c.m**3,rank_truncation=False,blocks=0,
        source_bindings=list(source_bindings),temporal_payload_peak=0,
        direct_kernel_payload_peak=0,static_payload_peak=0,
        scope='bare Gaussian first-map only; no dressed-G reconstruction',
        full_physical_scope=len(qs)==c.nq and len(pairs)==len(times)**2)
    for q in qs:
        kernels,peak=_direct_kernels(factors,q,u,occupations,internal_batch)
        manifest['direct_kernel_payload_peak']=max(manifest['direct_kernel_payload_peak'],kernels.nbytes)
        manifest['temporal_payload_peak']=max(manifest['temporal_payload_peak'],peak)
        qm=int(c.minus[q]);p=c.feature_partner[q][:len(factors.weights[q])]
        for i in range(c.ns):
            left,right,b,cc,aa=static_exchange(factors,q,np.array([i]))
            manifest['static_payload_peak']=max(manifest['static_payload_peak'],left.nbytes+right.nbytes)
            exchange=np.zeros((2,len(pairs),c.m,c.m),complex)
            for start in range(0,c.nq,inner_batch):
                sl=slice(start,min(start+inner_batch,c.nq))
                lt=_left_bank(left[sl],u,b[sl],cc[sl],aa[sl],occupations)
                rs=_right_bank(right[sl],u,b[sl],cc[sl],aa[sl],occupations)
                manifest['temporal_payload_peak']=max(manifest['temporal_payload_peak'],lt.nbytes+rs.nbytes)
                for k,(t,s) in enumerate(pairs):
                    exchange[0,k]+=1j*(lt[0,t]@rs[0,s]).sum(axis=0)
                    exchange[1,k]+=-1j*(lt[1,t]@rs[1,s]).sum(axis=0)
                del lt,rs
            bi=int(c.target[qm,i]);lq=c.density[q][bi].conj()
            rvertex=c.density[qm][i][p].conj().reshape(len(p),-1)
            for tile in _pair_tiles(len(times),pair_batch,pairs):
                values={};offset=pairs.index(tile[0]);indices=slice(offset,offset+len(tile))
                for order,name in enumerate(('greater','lesser')):
                    phase=(-1j,1j)[order]
                    direct=np.empty((1,len(tile),c.m,c.m),complex)
                    for k,(t,s) in enumerate(tile):
                        f=phase*(u[t,bi]@occupations[order][bi]@u[s,bi].conj().T)
                        x=lq@f
                        dressed=(kernels[order,t,s]@rvertex).reshape(len(p),c.m,c.m)
                        direct[0,k]=x.transpose(1,0,2).reshape(c.m,-1)@dressed.reshape(-1,c.m)
                    values[name+'_direct']=direct
                    values[name+'_exchange']=exchange[order,indices][None].copy()
                    values[name+'_total']=direct+values[name+'_exchange']
                mask=np.asarray([times[t]>=times[s] for t,s in tile])[None,:,None,None]
                for channel in ('direct','exchange','total'):
                    values['retarded_'+channel]=(values['greater_'+channel]-values['lesser_'+channel])*mask
                consume(Block(q,(i,),tile,values));manifest['blocks']+=1
            del left,right,exchange
    return manifest

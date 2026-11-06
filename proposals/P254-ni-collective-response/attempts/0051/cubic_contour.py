"""Leading cubic ionic contour callback, arbitrary external q, full internal mesh.

D greater/lesser = -i <delta u delta u>; Pi=-i<delta F delta F>.
V=Phi uuu/6 and F=-Phi uu/2. The complete ordered internal-q sum
combines the two Wick orders into -i/2 Phi W W Phi (A43 ROUTED_TERMS).
Canonical 1/sqrt(N) is inside each source vertex, never applied again.
Only one independently assembled forward/reverse vertex pair is resident.
"""
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'0047'))
from cubic_feedback import CompactCubic


def vertex(source,qs,piece='total'):
    """Vectorized association of the qualified compact source, not a q bank."""
    q=source._qs(qs,piece); n=source.nat
    out=np.zeros((n,n,n,3,3,3),complex)
    atoms=np.arange(n)
    if piece in ('total','electronic_local','electronic_nonlocal'):
        onsite=(source.local+source.nonlocal_term if piece=='total' else
                source.local if piece=='electronic_local' else source.nonlocal_term)
        out[atoms,atoms,atoms]+=onsite/np.sqrt(source.n)
    if piece in ('total','ionic_real'):
        phases=np.exp(2j*np.pi*(source.lattice@(q/source.mesh).T))
        for mask in range(8):
            labels=tuple(source.i if mask&(1<<leg) else source.j for leg in range(3))
            coefficient=np.prod(np.stack([np.ones(len(source.i)) if mask&(1<<leg)
                else -phases[:,leg] for leg in range(3)]),axis=0)
            np.add.at(out,labels,coefficient[:,None,None,None]*source.real)
    if piece in ('total','ionic_reciprocal'):
        i,j=np.indices((n,n))
        for mask,coefficient,pair in source._reciprocal(q):
            labels=tuple(i if mask&(1<<leg) else j for leg in range(3))
            np.add.at(out,labels,coefficient[:,:,None,None,None]*pair)
    return np.ascontiguousarray(out.transpose(0,3,1,4,2,5).reshape(source.d,source.d,source.d))


def ordered_matrix(left,right,first_d,second_d):
    """One routed Wick reduction, independent raw D legs, no conjugation."""
    d=left.shape[1]
    first_w=1j*np.asarray(first_d);second_w=1j*np.asarray(second_d)
    if first_w.shape!=(d,d) or second_w.shape!=(d,d):
        raise ValueError('complete ordered connected Cartesian D legs required')
    # L_iab W_ac W_bd, followed by the independent right_jcd.
    temp=(left.reshape(-1,d)@second_w).reshape(left.shape[0],d,d)
    transformed=np.einsum('iab,ac->icb',temp,first_w,optimize=True)
    return -.5j*(transformed.reshape(left.shape[0],-1)@right.reshape(right.shape[0],-1).T)


def ordered_action(left,right,first_d,second_d,value,transpose=False):
    """Matrix-free action or TRUE bilinear transpose, not Hilbert adjoint."""
    if transpose:
        return ordered_action(right,left,np.asarray(first_d).T,np.asarray(second_d).T,value)
    cot=np.asarray(value)
    if cot.shape!=(right.shape[0],):raise ValueError('full coordinate vector required')
    dual=np.einsum('j,jcd->cd',cot,right,optimize=True)
    covariance=(1j*np.asarray(first_d))@dual@(1j*np.asarray(second_d)).T
    return -.5j*np.einsum('iab,ab->i',left,covariance,optimize=True)


class CubicContour:
    """Stateless full-matrix/action callback; caller owns independent D histories.

    legs(k,l) returns (Dg_k,Dg_l,Dl_k,Dl_l) at the requested ordered time
    pair, each [111,111]. All64 calls are compulsory, including onsite terms.
    theta(0)=1. No coherent-mean-generated quadratic vertices are invented.
    Mean/JVP/transpose remain the existing CompactCubic APIs.
    """
    def __init__(self,source):
        self.source=source;self.d=source.d;self.n=source.n

    @classmethod
    def from_source(cls,path):return cls(CompactCubic.from_source(path))

    def pair(self,external_q,internal_q,piece='total'):
        s=self.source;q=int(external_q);k=int(internal_q)
        if not 0<=q<s.n or not 0<=k<s.n:raise IndexError('original mesh q outside range')
        l=s.qindex(s.indices[q]-s.indices[k])
        left=vertex(s,s.indices[[int(s.minus[q]),k,l]],piece)
        right=vertex(s,s.indices[[q,int(s.minus[k]),int(s.minus[l])]],piece)
        return l,left,right

    def evaluate(self,external_q,legs,theta,value=None,transpose=False):
        if theta not in (0,1):raise ValueError('explicit causal theta=0/1 required')
        shape=(self.d,self.d) if value is None else (self.d,)
        out={key:np.zeros(shape,complex) for key in ('greater','lesser')}
        for k in range(self.n):
            l,left,right=self.pair(external_q,k)
            dg1,dg2,dl1,dl2=legs(k,l)
            for key,a,b in (('greater',dg1,dg2),('lesser',dl1,dl2)):
                out[key]+=(ordered_matrix(left,right,a,b) if value is None else
                           ordered_action(left,right,a,b,value,transpose))
        out['retarded']=theta*(out['greater']-out['lesser'])
        return out

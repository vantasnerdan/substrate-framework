"""Exact represented degree14 density insertion on common baseline/dual knots.

n_ext adds to D:sigma/sqrt(N), BEFORE K and second 1/sqrt(N). Its dual
is K.T P.T C_rho.T z_sigma/sqrt(N). Physical source density vertices carry
the OTHER 1/sqrt(N), no quadrature weight, spin2 or inverse Fourier factor.
"""
from dataclasses import dataclass
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent/n) for n in ('0030','0029')]
import numpy as np
from canonical_force_filter import CanonicalForceFilter, build_canonical_force_filter, chebyshev_lobatto_transform
from reduced_correlation_transpose import commutator_transpose
from substrate_framework.streamed_material_state import authenticate


def density_insertion(correlations,q,time,dual,block=None):
    """Same-q nonconjugating external density dual; actual source/target rho.

    This is ActiveBandCorrelations._transpose's existing density_dual, not
    an independent stochastic bath. h drops out of commutator's dh dual.
    """
    c=correlations.c
    block=c.q(q) if block is None else block
    _,zs,_=correlations._split(q,dual)
    rho,h=correlations._baseline(time)
    target=c.target[q]
    _,dh_dual=commutator_transpose(h[target],h,rho[target],rho,zs)
    feature=np.einsum('sabr,sabp->rp',block.potential,dh_dual,optimize=True)/correlations.root_n
    return block.kernel.T@feature


@dataclass(frozen=True)
class CanonicalDensityFilter:
    force_filter: object
    density_segment_bounds: np.ndarray
    density_chebyshev: np.ndarray
    density_boundary_times: np.ndarray
    density_boundary: np.ndarray
    metadata: dict

    @property
    def arrays(self):
        return dict(**self.force_filter.arrays,density_segment_bounds=self.density_segment_bounds,
            density_chebyshev=self.density_chebyshev,density_boundary_times=self.density_boundary_times,
            density_boundary=self.density_boundary)

    def evaluate_density(self,time):
        time=float(time)
        if not np.isfinite(time) or not 0<=time<=self.force_filter.final_time:
            raise ValueError('density insertion time outside interval')
        right=int(np.searchsorted(self.density_boundary_times,time,side='left'))
        if right<len(self.density_boundary_times) and self.density_boundary_times[right]==time:
            return self.density_boundary[right].copy()
        left,end=self.density_segment_bounds[right-1]
        return np.polynomial.chebyshev.chebval(2*(time-left)/(end-left)-1,
                                              self.density_chebyshev[right-1])


def load_canonical_density_filter(packet):
    """Consume the actual completed numeric filter; never repeat its dual solve."""
    metadata=packet['filter']
    if packet['status']!='complete' or metadata['schema']!='P254-A42-canonical-density-forcing-v1':
        raise ValueError('complete saved canonical density filter required')
    force_metadata=metadata['force_filter']
    with np.load(authenticate(packet['archive']),allow_pickle=False) as archive:
        force=CanonicalForceFilter(
            final_time=force_metadata['final_time_inverse_ry'],metadata=force_metadata,
            **{name:archive[name] for name in (
                'initial_dual','segment_bounds','ionic_chebyshev','field_chebyshev',
                'boundary_times','ionic_boundary','field_boundary')})
        return CanonicalDensityFilter(force,
            *(archive[name] for name in (
                'density_segment_bounds','density_chebyshev',
                'density_boundary_times','density_boundary')),metadata)


class _DensityCapture:
    def __init__(self,correlations,q,baseline_knots):
        self.correlations,self.q=correlations,q
        self.block=correlations.c.q(q)
        self.baseline_knots=np.asarray(baseline_knots,float)
        self.nodes,self.transform=chebyshev_lobatto_transform(14)
        self.segments=[];self.boundaries={}

    def insertion(self,time,dual):
        value=density_insertion(self.correlations,self.q,time,dual,self.block)
        if not np.isfinite(value).all():raise FloatingPointError('nonfinite density dual')
        return value

    def boundary(self,time,dual):
        self.boundaries[float(time)]=self.insertion(time,dual)

    def segment(self,left,right,dense,shape):
        # No full-state history bank: the current public interpolant dies at
        # return. Split BEFORE multiplying represented degree7 rho and z.
        cuts=np.r_[left,self.baseline_knots[(self.baseline_knots>left)&(self.baseline_knots<right)],right]
        for a,b in zip(cuts[:-1],cuts[1:]):
            values=[]
            for j,node in enumerate(self.nodes):
                t=b if j==0 else a if j==14 else a+(b-a)*(node+1)/2
                values.append(self.insertion(t,dense(t).reshape(shape)))
            coefficients=np.tensordot(self.transform,np.asarray(values),axes=(1,0))
            self.segments.append((a,b,coefficients))
            self.boundaries[float(a)]=values[-1]
            self.boundaries.setdefault(float(b),values[0])


def build_canonical_density_filter(correlations,q,final_time,coordinates,*,baseline_knots=None):
    """ONE NEW backward solve, exact degree14 product of represented polynomials.

    baseline_knots MUST cover the supplied degree<=7 baseline representation.
    For SavedMovingBaseline they default to its public .times. A non-polynomial
    arbitrary callable is unsupported, not silently sampled at degree7/14.
    On every union interval rho and z are degree7; fixed P/K and the commutator
    insertion are bilinear, hence degree<=14. Fifteen Lobatto samples recover
    that product without truncation or new tolerance/accuracy gate. Baseline
    and dual ODE approximation errors remain, not certified by this identity.
    """
    if baseline_knots is None:
        baseline_knots=getattr(correlations.baseline,'times',None)
    if baseline_knots is None:raise ValueError('declared degree7 baseline knots required')
    knots=np.asarray(baseline_knots,float)
    if knots.ndim!=1 or len(knots)<2 or knots[0]!=0 or knots[-1]<final_time or np.any(np.diff(knots)<=0):
        raise ValueError('ascending baseline knots must cover [0,terminal]')
    capture=_DensityCapture(correlations,q,knots)
    force=build_canonical_force_filter(correlations,q,final_time,coordinates,capture=capture)
    segments=sorted(capture.segments,key=lambda row:row[0])
    times=np.array(sorted(capture.boundaries))
    r,p=capture.boundaries[times[0]].shape
    bounds=np.asarray([(a,b) for a,b,_ in segments],float).reshape(-1,2)
    coefficient=np.stack([row[2] for row in segments]) if segments else np.empty((0,15,r,p),complex)
    metadata=dict(schema='P254-A42-canonical-density-forcing-v1',force_filter=force.metadata,
        q_index=int(q),feature_count=r,feature_g=np.asarray(capture.block.feature_g).tolist(),
        feature_spin=np.asarray(capture.block.feature_spin).tolist(),polynomial_degree=14,
        polynomial='exact degree<=14 product of represented degree7 baseline rho and backward dual on union knots',
        pairing='nonconjugating same-q; K.T, no dagger/no -q swap',
        computation='NEW backward solve using explicitly reconstructed durable baseline; no mean replay here',
        extra_interpolation_approximation=False,full_state_history_saved=False,
        error_scope='original ODE approximation errors retained; no new tolerance gate or enclosure')
    return CanonicalDensityFilter(force,bounds,coefficient,times,
        np.stack([capture.boundaries[t] for t in times]),metadata)

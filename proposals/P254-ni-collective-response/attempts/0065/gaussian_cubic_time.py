"""Exact all-phase-rank cubic Pi for the bare Gaussian first map only.

W^g_q(t,s)=P S_q(t) W0_q S_-q(s)^T P^T;
W^l_q(t,s)=P S_q(t) W0_-q^T S_-q(s)^T P^T.
For X_q=P S_q W0 and Y_q=P S_-q, form
F_iuv(t)=sum_ab left_iab X_k,au(t) X_l,bv(t),
H_juv(s)=sum_cd right_jcd Y_k,cu(s) Y_l,dv(s).
Pi_ij=-i/2 sum_uv F_iuv H_juv. Transposes are bilinear, never
Hilbert adjoints. Every phase coordinate is retained, even singular modes.
The first phase axis is chunked; no factor bank grows with nt squared.
"""
from pathlib import Path
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / '0052'))
from general_moving_reference import GaussianTransitions


def transform_chunk(vertex, first, second, start, stop):
    """[external,u,v] with only u chunked; no full-rank intermediate."""
    partial = np.einsum('iab,au->iub', vertex, first[:,start:stop], optimize=True)
    return (partial.reshape(-1, vertex.shape[2]) @ second).reshape(
        vertex.shape[0], stop-start, second.shape[1])


def stationarity_defect(reference, q):
    """Raw necessary/sufficient generator identity, without declaring tolerance zero."""
    q = int(q)
    qm = int(reference.c.minus[q])
    mechanical = reference.c.mechanical
    a, b, w = mechanical.generators[q], mechanical.generators[qm], mechanical.w0[q]
    return a @ w + w @ b.T


class GaussianCubicTimeSeparated:
    """Reuse each time's transformed vertex across all ordered time pairs.

    The caller supplies owning A51 independent vertices, q partner and a running
    greater/lesser sum. Call accumulate_q for each declared internal q, including
    both partner routes where the owning ordered sum includes both. Nothing here
    supplies missing q routes or embeds the ionic block into field coordinates.
    """
    def __init__(self, reference, chunk=16):
        if not isinstance(reference, GaussianTransitions):
            raise TypeError('only original bare Gaussian factors are supported')
        if int(chunk) != chunk or chunk < 1:
            raise ValueError('positive integer phase chunk required')
        self.reference = reference
        self.times = np.asarray(reference.times)
        self.chunk = int(chunk)

    def factors(self, q):
        ref = self.reference
        q = int(q)
        qm = int(ref.c.minus[q])
        ix = ref.position_indices(q, 'ionic')
        im = ref.position_indices(qm, 'ionic')
        sq = np.asarray([ref.S(q,t)[ix] for t in self.times])
        sm = np.asarray([ref.S(qm,t)[im] for t in self.times])
        return (sq @ ref.c.mechanical.w0[q],
                sq @ ref.c.mechanical.w0[qm].T, sm)

    def zeros(self, external_coordinates):
        shape = (len(self.times), len(self.times), external_coordinates, external_coordinates)
        return {key: np.zeros(shape, complex) for key in ('greater','lesser')}

    def accumulate_q(self, left, right, k, partner, out):
        started = time.monotonic()
        if left.shape != right.shape or left.ndim != 3 or left.shape[1:] != (self.reference.c.d,)*2:
            raise ValueError('complete independent ionic vertices required')
        nt, d = len(self.times), left.shape[0]
        if any(out[key].shape != (nt,nt,d,d) for key in ('greater','lesser')):
            raise ValueError('complete ordered time output axes required')
        first, second = self.factors(k), self.factors(partner)
        factors_seconds = time.monotonic() - started
        transform_seconds = assembly_seconds = 0.
        rank1, rank2 = first[2].shape[2], second[2].shape[2]
        peak_feature_bytes = 0
        for start in range(0,rank1,self.chunk):
            stop = min(start+self.chunk,rank1)
            tick = time.monotonic()
            # Construct directly into the bank, avoiding a list plus stacked copy.
            bank_shape = (nt,d,(stop-start)*rank2)
            rb = np.empty(bank_shape,complex)
            for b in range(nt):
                rb[b] = transform_chunk(right,first[2][b],second[2][b],start,stop).reshape(d,-1)
            transform_seconds += time.monotonic()-tick
            for component,key in enumerate(('greater','lesser')):
                tick = time.monotonic()
                lb = np.empty(bank_shape,complex)
                for a in range(nt):
                    lb[a] = transform_chunk(left,first[component][a],second[component][a],start,stop).reshape(d,-1)
                transform_seconds += time.monotonic()-tick
                peak_feature_bytes = max(peak_feature_bytes,rb.nbytes+lb.nbytes)
                tick = time.monotonic()
                for a in range(nt):
                    for b in range(nt):
                        out[key][a,b] += (-.5j)*(lb[a] @ rb[b].T)
                assembly_seconds += time.monotonic()-tick
                del lb
            del rb
        return dict(internal_q=int(k),partner_q=int(partner),ranks=[rank1,rank2],
                    discarded_phase_modes=0,chunk=self.chunk,factor_seconds=factors_seconds,
                    transform_seconds=transform_seconds,assembly_seconds=assembly_seconds,
                    total_seconds=time.monotonic()-started,peak_feature_bank_bytes=peak_feature_bytes)

    def retarded(self, out):
        theta = self.times[:,None] >= self.times[None,:]
        return theta[:,:,None,None]*(out['greater']-out['lesser'])

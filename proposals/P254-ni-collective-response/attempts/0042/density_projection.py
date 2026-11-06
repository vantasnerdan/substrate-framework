"""Exact finite-grid Fourier density vertices in borrowed physical PW frames.

Multiplication by one Fourier character is an integer shift of the admitted
FFT slots, including its actual grid aliases and directed reciprocal fold.
This changes representation, not the field, cutoff, volume or state.
"""
import numpy as np
from scipy.linalg.blas import zgemm


class BorrowedBathFrame:
    """B=[psi S^-1/2,E]; form its active columns once, borrow the paid extra."""
    def __init__(self, psi, inverse_sqrt_gram, extra_basis):
        self.rank = psi.shape[1]
        self.size = self.rank + extra_basis.shape[1]
        if inverse_sqrt_gram.shape != (self.rank, self.rank) or extra_basis.shape[0] != len(psi):
            raise ValueError('paid frame coordinate shapes differ')
        self.active = zgemm(1.,np.asfortranarray(psi),
                           np.asfortranarray(inverse_sqrt_gram))
        self.extra = np.asfortranarray(extra_basis)

    def columns(self, start, stop):
        if stop <= self.rank:
            return self.active[:,start:stop]
        if start >= self.rank:
            return self.extra[:,start-self.rank:stop-self.rank]
        return np.asfortranarray(np.concatenate(
            (self.active[:,start:self.rank],self.extra[:,:stop-self.rank]),axis=1))

    def left_project(self, action):
        # BLAS's actual dagger flag avoids conjugating/copying the large banks.
        active = zgemm(1.,self.active,action,trans_a=2)
        if not self.extra.shape[1]:
            return active
        return np.concatenate(
            (active,zgemm(1.,self.extra,action,trans_a=2)),axis=0)


class FourierDensityProjector:
    """Independent O_st=A_ts[phi_G]^dagger/sqrt(N) and opposite action.

    phi_G=exp(2pi i G.x)/sqrt(Omega). The bridge multiplier is
    exp(-2pi i fold.x)*phi_G. No active embedding, spin/slot weight, mass,
    additional volume factor or continuum/unaliased accuracy is assumed.
    """
    def __init__(self, bridge, source_frame, target_frame, *, ncells, block_size=128):
        if type(ncells) is not int or ncells<1 or type(block_size) is not int or block_size<1:
            raise ValueError('positive exact N and PW column block size required')
        if hasattr(bridge,'source_fft'):
            self.source_fft,self.target_fft = bridge.source_fft,bridge.target_fft
            self.fold = np.asarray(bridge.edge.reciprocal_fold)
        else:
            self.source_fft = self.target_fft = bridge.factors.fft
            self.fold = np.zeros(3,dtype=int)
        if (self.source_fft.physical!=self.target_fft.physical
                or self.source_fft.omega!=self.target_fft.omega):
            raise ValueError('same physical cell/grid required for Fourier multiplication')
        self.source,self.target = source_frame,target_frame
        if (len(source_frame.active)!=len(self.source_fft.nl)
                or len(target_frame.active)!=len(self.target_fft.nl)):
            raise ValueError('borrowed frame does not retain its full PW endpoint')
        self.block_size = block_size
        self.normalization = np.sqrt(ncells)*np.sqrt(self.source_fft.omega)
        self.source_lookup = self._lookup(self.source_fft)
        self.target_lookup = (self.source_lookup if self.target_fft is self.source_fft
                              else self._lookup(self.target_fft))
        self.forward_workspace = np.empty(
            (len(self.source_fft.nl),min(block_size,target_frame.size)),complex,order='F')
        self.adjoint_workspace = np.empty(
            (len(self.target_fft.nl),min(block_size,source_frame.size)),complex,order='F')

    @staticmethod
    def _lookup(fft):
        lookup = np.full(fft.nnr,-1,dtype=np.intp)
        lookup[fft.nl] = np.arange(len(fft.nl))
        return lookup

    def _leg(self, destination, incoming, fft, lookup, shift, workspace):
        shape = fft.physical
        wanted = fft.miller-shift[:,None]
        slots = np.ravel_multi_index(
            tuple(wanted[i]%shape[i] for i in range(3)),shape,order='F')
        incoming_rows = lookup[slots]
        rows = np.flatnonzero(incoming_rows>=0)
        incoming_rows = incoming_rows[rows]
        result = np.empty((destination.size,incoming.size),complex,order='F')
        for start in range(0,incoming.size,self.block_size):
            stop = min(start+self.block_size,incoming.size)
            action = workspace[:,:stop-start]
            action.fill(0)
            action[rows] = incoming.columns(start,stop)[incoming_rows]
            result[:,start:stop] = destination.left_project(action)
        result /= self.normalization
        return result

    def project(self, g):
        g = np.asarray(g)
        if g.shape!=(3,) or g.dtype.kind not in 'iu':
            raise ValueError('integer Fourier profile G required')
        forward = self._leg(self.source,self.target,self.source_fft,
            self.target_lookup,self.fold-g,self.forward_workspace)
        adjoint = self._leg(self.target,self.source,self.target_fft,
            self.source_lookup,g-self.fold,self.adjoint_workspace)
        return forward,adjoint

"""Full-PW Fourier density vertices on borrowed source/target bath frames.

For directed sigma_ts, D_r:sigma_ts=Tr(O_st sigma_ts), where
O_st=A_ts[phi_r]^dagger, phi=exp(2pi i G.r)/sqrt(Omega).
This is the existing Bloch bridge action, not an active-band zero embedding.
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / '0029'))
import numpy as np
from reduced_projection_algebra import profile


class BorrowedBathFrame:
    """B=[psi S^-1/2,E] without storing/rebuilding B; paid E and S borrowed."""
    def __init__(self, psi, inverse_sqrt_gram, extra_basis):
        self.psi, self.inverse, self.extra = psi, inverse_sqrt_gram, extra_basis
        self.rank = psi.shape[1]
        self.size = self.rank + extra_basis.shape[1]
        if inverse_sqrt_gram.shape != (self.rank, self.rank) or extra_basis.shape[0] != len(psi):
            raise ValueError('paid frame coordinate shapes differ')

    def columns(self, start, stop):
        if stop <= self.rank:
            return self.psi @ self.inverse[:, start:stop]
        if start >= self.rank:
            return self.extra[:, start-self.rank:stop-self.rank]
        return np.concatenate((self.psi @ self.inverse[:, start:self.rank],
                               self.extra[:, :stop-self.rank]), axis=1)

    def left_project(self, action):
        return np.concatenate((self.inverse.conj().T @ (self.psi.conj().T @ action),
                               self.extra.conj().T @ action), axis=0)


def project_density_feature(bridge, source_frame, target_frame, g, *, ncells,
                            block_size=8):
    """Return O_st/sqrt(N), O_ts†/sqrt(N), independently applied.

    Works for every actual directed q bridge; source/target can have distinct
    PW and finite-bath dimensions. Caller selects the explicit spin feature.
    No quadrature weight, spin multiplier or final inverse Fourier 1/N occurs.
    One PW column block and two projected matrices live at a time.
    """
    if type(ncells) is not int or ncells < 1 or block_size < 1:
        raise ValueError('positive exact N and column block size required')
    if hasattr(bridge, 'potential_action'):
        # Existing physical q=0 dyad coupling is a multiplication action.
        fft = bridge.factors.fft
        field = profile(fft.physical, g, fft.omega)
        apply_density = lambda vectors: bridge.potential_action(field.conj(), vectors)
        apply_adjoint = lambda vectors: bridge.potential_action(field, vectors)
    else:
        fft = bridge.source_fft if hasattr(bridge, 'source_fft') else bridge.orbital_fft
        field = profile(fft.physical, g, fft.omega)
        apply_density = lambda vectors: bridge.local_adjoint_apply(field, vectors)
        apply_adjoint = lambda vectors: bridge.local_apply(field, vectors)
    forward = np.empty((source_frame.size, target_frame.size), complex)
    adjoint = np.empty((target_frame.size, source_frame.size), complex)
    for start in range(0, target_frame.size, block_size):
        stop = min(start + block_size, target_frame.size)
        forward[:, start:stop] = source_frame.left_project(
            apply_density(target_frame.columns(start, stop)))
    for start in range(0, source_frame.size, block_size):
        stop = min(start + block_size, source_frame.size)
        adjoint[:, start:stop] = target_frame.left_project(
            apply_adjoint(source_frame.columns(start, stop)))
    forward /= np.sqrt(ncells)
    adjoint /= np.sqrt(ncells)
    return forward, adjoint

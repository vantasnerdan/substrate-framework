"""Rectangular finite-q full-PW multiplication and directed dyad density.

For X_ts=L_t C R_s^dagger the envelope is exp(+2 pi i ell.x) times
U_t L C (U_s R)^*/Omega. The local action instead uses exp(-2 pi i ell.x).
Physical pairing is Omega/NFFT; no volume enters number multiplication.
No saved-band closure, implicit spin degeneracy or k weight is supplied.
"""
import numpy as np

from .qe_bloch_edge import BlochEdge
from .qe_material_grid import PhysicalRhoFFT


class BlochFiniteQBridge:
    def __init__(self, inputs, edge, *, block_size=8):
        if not isinstance(edge, BlochEdge):
            raise ValueError('checked directed BlochEdge required')
        edge.validate(inputs.reference)
        if type(block_size) is not int or not 1 <= block_size <= 16:
            raise ValueError('bounded column block_size in 1..16 required')
        ref = inputs.reference
        if ref.fft_grid != ref.fft_smooth:
            raise ValueError('same physical smooth/dense grid required')
        self.edge, self.block_size = edge, block_size
        self.ispin = edge.ispin
        volume = abs(float(np.linalg.det(ref.cell_bohr)))
        self.source_fft = PhysicalRhoFFT.create(ref.fft_smooth, edge.source_basis.miller,
                                               ref.reciprocal_bohr_inv, volume)
        self.target_fft = PhysicalRhoFFT.create(ref.fft_smooth, edge.target_basis.miller,
                                               ref.reciprocal_bohr_inv, volume)
        self.density_fft = PhysicalRhoFFT.create(ref.fft_grid, inputs.density.miller,
                                                ref.reciprocal_bohr_inv, volume,
                                                edge.q_cart_bohr_inv)
        shape = self.source_fft.physical
        # Separable construction avoids a stored three-coordinate real-space grid.
        phase = np.ones(shape, dtype=np.complex128)
        for axis in range(3):
            axis_shape = [1, 1, 1]
            axis_shape[axis] = shape[axis]
            phase *= np.exp(-2j*np.pi*edge.reciprocal_fold[axis]
                            * np.arange(shape[axis])/shape[axis]).reshape(axis_shape)
        self.multiplication_phase = phase.reshape(-1, order='F')
        self.multiplication_phase.setflags(write=False)
        self.grid_measure_bohr3 = volume/self.source_fft.nnr

    @staticmethod
    def _matrix(value, rows):
        value = np.asarray(value)
        if value.ndim != 2 or value.shape[0] != rows or not np.isfinite(value).all():
            raise ValueError('finite full-PW columns with the correct endpoint rows required')
        return value

    def _field(self, value):
        value = np.asarray(value)
        if value.shape != (self.source_fft.nnr,) or not np.isfinite(value).all():
            raise ValueError('finite F-order physical scalar envelope required')
        return value

    def local_apply(self, field, vectors):
        """A_ts[v] on [nsource,ncolumn]; v is Ry or a number profile."""
        field = self._field(field)
        vectors = self._matrix(vectors, self.edge.source_dimension)
        out = np.empty((self.edge.target_dimension, vectors.shape[1]), dtype=np.complex128)
        multiplier = self.multiplication_phase*field
        for start in range(0, vectors.shape[1], self.block_size):
            sl = slice(start, start+self.block_size)
            out[:, sl] = self.target_fft.to_g(
                multiplier[:, None]*self.source_fft.from_g(vectors[:, sl]))
        return out

    def local_adjoint_apply(self, field, vectors):
        """Independent exact A_ts[v]^dagger; no rebuilt/truncated -q field."""
        field = self._field(field)
        vectors = self._matrix(vectors, self.edge.target_dimension)
        out = np.empty((self.edge.source_dimension, vectors.shape[1]), dtype=np.complex128)
        multiplier = (self.multiplication_phase*field).conj()
        for start in range(0, vectors.shape[1], self.block_size):
            sl = slice(start, start+self.block_size)
            out[:, sl] = self.source_fft.to_g(
                multiplier[:, None]*self.target_fft.from_g(vectors[:, sl]))
        return out

    def iter_dyad_blocks(self, left, right, coefficient=None):
        """Yield bounded factors for L_t C R_s^dagger, absorbing C on smaller leg."""
        left = self._matrix(left, self.edge.target_dimension)
        right = self._matrix(right, self.edge.source_dimension)
        nl, nr = left.shape[1], right.shape[1]
        if coefficient is None:
            if nl != nr:
                raise ValueError('identity coefficient requires equal column counts')
        else:
            coefficient = np.asarray(coefficient)
            if coefficient.shape != (nl, nr) or not np.isfinite(coefficient).all():
                raise ValueError('finite rectangular coefficient[nleft,nright] required')
        for start in range(0, min(nl, nr), self.block_size):
            stop = min(start+self.block_size, min(nl, nr))
            if coefficient is None:
                yield left[:, start:stop], right[:, start:stop]
            elif nr <= nl:
                yield left @ coefficient[:, start:stop], right[:, start:stop]
            else:
                yield left[:, start:stop], right @ coefficient[start:stop, :].conj().T

    @staticmethod
    def _weight(weight):
        if not np.isscalar(weight) or not np.isreal(weight) or not np.isfinite(weight):
            raise ValueError('explicit finite real sector weight required')
        return float(np.real(weight))

    def density_from_dyads(self, left, right, coefficient=None, *, weight):
        """Complex [nr,2] bohr^-3 density; weight is supplied exactly once."""
        weight = self._weight(weight)
        density = np.zeros((self.source_fft.nnr, 2), dtype=np.complex128)
        channel = density[:, self.ispin-1]
        for lblock, rblock in self.iter_dyad_blocks(left, right, coefficient):
            lwave = self.target_fft.from_g(lblock)
            rwave = self.source_fft.from_g(rblock)
            channel += np.sum(lwave*rwave.conj(), axis=1)
        channel *= weight*self.multiplication_phase.conj()/self.source_fft.omega
        return density

    def density_gradient(self, density):
        """i(q+G) on the declared density support, never orbital i(k+G)."""
        density = np.asarray(density)
        if density.shape != (self.source_fft.nnr, 2) or not np.isfinite(density).all():
            raise ValueError('finite explicit-spin physical density required')
        return self.density_fft.gradient(density)

    def density_pair(self, density, field):
        """Omega/N sum conj(n_q) v_q; equals weight Tr(X_ts^dagger A_ts[v])."""
        density = np.asarray(density)
        field = self._field(field)
        if density.shape != (self.source_fft.nnr, 2) or not np.isfinite(density).all():
            raise ValueError('finite explicit-spin physical density required')
        return complex(self.grid_measure_bohr3*np.vdot(density[:, self.ispin-1], field))

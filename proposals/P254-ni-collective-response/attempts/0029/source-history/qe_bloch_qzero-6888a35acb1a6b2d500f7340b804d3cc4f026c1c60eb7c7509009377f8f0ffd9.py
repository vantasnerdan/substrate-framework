"""Same-k full-PW periodic number/density bridge; density gradients use G.

Orbital FFTs reconstruct the periodic envelope; k cancels in q=0 dyads.
Density is bohr^-3, while integrated number multiplication is dimensionless.
No implicit spin multiplicity, star weight, or preparation is supplied.
"""
import numpy as np
from .qe_material_grid import PhysicalRhoFFT


class BlochQZeroBridge:
    def __init__(self, inputs, basis):
        basis.validate(inputs.reference)
        ref = inputs.reference
        volume = abs(np.linalg.det(ref.cell_bohr))
        self.basis = basis
        self.ispin = basis.slot.ispin
        self.orbital_fft = PhysicalRhoFFT.create(ref.fft_smooth, basis.miller,
                                                basis.reciprocal_bohr_inv, volume)
        self.density_fft = PhysicalRhoFFT.create(ref.fft_grid, inputs.density.miller,
                                                ref.reciprocal_bohr_inv, volume)
        if ref.fft_grid != ref.fft_smooth:
            raise ValueError('same physical dense/smooth grid required')

    def _block(self, vectors):
        value = np.asarray(vectors, dtype=np.complex128)
        if (value.ndim != 2 or value.shape[0] != self.basis.miller.shape[1]
                or not np.isfinite(value).all()):
            raise ValueError('finite full PW block [npw,ncolumn] required')
        return value

    def local_apply(self, field, vectors):
        """Periodic potential or integrated-number profile on all PW rows."""
        field = np.asarray(field)
        if field.shape != (self.orbital_fft.nnr,) or not np.isfinite(field).all():
            raise ValueError('finite F-order physical scalar field required')
        x = self._block(vectors)
        result = np.empty(x.shape, dtype=np.complex128)
        for start in range(0, x.shape[1], 8):
            block = x[:, start:start+8]
            result[:, start:start+8] = self.orbital_fft.to_g(field[:, None]*self.orbital_fft.from_g(block))
        return result

    def local_adjoint_apply(self, field, vectors):
        return self.local_apply(np.asarray(field).conj(), vectors)

    def density_from_dyads(self, left, right, coefficients, *, weight):
        """Density of weight*sum_j coefficients[j]|left_j><right_j|, explicit spin.

        Arbitrary complex dyads are allowed. The caller declares the mode weight;
        no normalization/probability interpretation is inferred from coefficients.
        """
        left, right = self._block(left), self._block(right)
        c = np.asarray(coefficients, dtype=np.complex128)
        if left.shape != right.shape or c.shape != (left.shape[1],) or not np.isfinite(c).all() or not np.isfinite(weight):
            raise ValueError('matched dyad columns, finite coefficients/weight required')
        result = np.zeros((self.orbital_fft.nnr, 2), dtype=np.complex128)
        for start in range(0, left.shape[1], 8):
            sl = slice(start, start+8)
            l, r = self.orbital_fft.from_g(left[:, sl]), self.orbital_fft.from_g(right[:, sl])
            result[:, self.ispin-1] += np.sum(l*r.conj()*c[sl], axis=1)*weight/self.orbital_fft.omega
        return result

    def density_gradient(self, density):
        """iG gradient on saved density support, NOT orbital i(k+G)."""
        density = np.asarray(density)
        if density.shape != (self.density_fft.nnr, 2) or not np.isfinite(density).all():
            raise ValueError('finite explicit up/down density required')
        return self.density_fft.gradient(density)

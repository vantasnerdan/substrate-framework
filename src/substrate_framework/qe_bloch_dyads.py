"""Unweighted q=0 density/ionic-gradient maps in a checked Bloch space.

The admitted Gamma implementation's block contractions and trace reductions
are k-independent algorithms: orbital envelopes supply density, while the
placed projectors and their derivatives use k+G. This adapter replaces only
its Gamma-specific identity admission with BlochIonicAction's checked basis.
No old reader/source contract changes and no derivative algorithm is copied.

Ordinary NC, no NLCC; positive dH/dR in Ry/bohr, density in bohr^-3.
Spin/full-grid measure and physical preparation belong to the caller. This
is not a rectangular finite-q map or a stationarity assertion.
"""
from types import SimpleNamespace
import numpy as np
from scipy.linalg.blas import zgemm

from .qe_bloch_ionic_action import BlochIonicAction
from .qe_material_dyads import GammaDyadCoupling
from .qe_material_gradient import GammaGradientReduction


class BlochDyadCoupling(GammaDyadCoupling):
    """Reuse the existing DyadCoupling API on all rows of one Bloch basis.

    The inherited methods do not assume k=0. The gradient's geometry view
    borrows the checked ionic action's arrays; in particular its p contains
    k+G, whereas its g and local Fourier support contain density G only.
    """

    def __init__(self, inputs, basis, *, factors=None, block_size=8):
        if type(block_size) is not int or not 1 <= block_size <= 16:
            raise ValueError("block_size must be an integer in 1..16")
        self.ionic = BlochIonicAction(inputs, basis, factors=factors)
        self.factors = self.ionic.factors
        self.block_size = block_size
        fft, ref = self.factors.fft, inputs.reference
        self.npw = fft.miller.shape[1]
        self.nr = fft.nnr
        self.ndof = 3 * ref.nat
        self.cell_volume_bohr3 = fft.omega
        self.grid_measure_bohr3 = fft.omega / fft.nnr
        geometry = SimpleNamespace(
            p=self.ionic.p, g=self.ionic.g, shape=fft.physical,
            species=ref.atom_names, tables=inputs.tables, vloc=self.ionic.vloc,
            positions=ref.positions_bohr,
            local_indices=tuple(inputs.density.miller[i] % fft.physical[i] for i in range(3)),
        )
        self._gradient = GammaGradientReduction(self.factors, geometry)

    def prepared_moments(self, psi, coefficient):
        """Same raw dyad, with wide BLAS and only the FFTs column-blocked.

        One full-row weighted frame trades bounded extra memory for fewer
        skinny GEMMs. No rank truncation, Hermitization or changed pairing.
        """
        psi = self._matrix(psi, "raw frame")
        coefficient = np.asarray(coefficient)
        if coefficient.shape != (psi.shape[1], psi.shape[1]) or not np.isfinite(coefficient).all():
            raise ValueError("finite raw-frame coefficient matrix required")
        weighted = zgemm(1., psi, coefficient)
        gradient = self._gradient.nonlocal_block(weighted, psi)
        density = np.zeros(self.nr, dtype=np.complex128)
        fft = self.factors.fft
        for start in range(0, psi.shape[1], self.block_size):
            stop = min(start + self.block_size, psi.shape[1])
            left = fft.from_g(weighted[:, start:stop])
            right = fft.from_g(psi[:, start:stop])
            density += np.sum(left * right.conj(), axis=1) / self.cell_volume_bohr3
        gradient += self._gradient.local(density)
        return density, gradient

"""Isolated full-PW physical response to +v + g*u, with both q branches.

For exp(-i*z*t), the density-matrix equation is
    z*delta_rho - [h0, delta_rho] = [delta_h, rho0].
Expanding the two terms of the right commutator in occupied source columns
produces the plus ket and minus bra below. In particular the second ket uses
conj(z), then its bra is conjugated; neither branch has an extra minus sign.
The returned ionic observable is delta<g_-q>, not force or positive chi to -xB.
No thermal population derivative, screening, or physical damping is added.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .pw_response_data import PWResponseChannel, PWResponseModel, PWResponseValue


from .pw_transition_space import PWTransitionSpace, finite_array


class PWIndependentResponse:
    """Reusable retarded response at one energy z (Ry), strictly Im(z)>0.

    Source occupations already include k quadrature and their explicit spin.
    Full target eigenvectors have PW rows and eigenstate columns. Target FFT
    slots are zero-based storage offsets, independently mapped for each branch.
    Models are treated as immutable for the lifetime of this response object.
    """

    def __init__(self, model: PWResponseModel, z: complex):
        self.model = model
        self.z = complex(z)
        if not np.isfinite(self.z) or self.z.imag <= 0:
            raise ValueError("z must be finite with strictly positive imaginary part")
        self._space = PWTransitionSpace(model)

    def _ket(self, channel: PWResponseChannel, branch: str, potential: NDArray,
             displacement: NDArray, coordinates: tuple[NDArray, ...]) -> NDArray:
        vectors = getattr(channel, f"{branch}_vectors")
        energies = getattr(channel, f"{branch}_energies")
        drive = np.einsum("atn,a->tn", getattr(channel, f"g_{branch}"), displacement,
                          dtype=np.complex128)
        drive += self._space.to_pw(potential[:, None] * channel.source_grid, coordinates)
        shift = self.z if branch == "plus" else -self.z.conjugate()
        denominator = channel.source_energies[None, :] + shift - energies[:, None]
        return vectors @ ((vectors.conj().T @ drive) / denominator)

    def apply(self, local_potential, displacement) -> PWResponseValue:
        """Return physical spin density (bohr^-3) and ionic gradient (Ry/bohr).

        Complex amplitudes enter linearly: minus-branch conjugation belongs to
        its ket construction and is undone in the observable bra contraction.
        Padding is not physical volume and contributes no density integral.
        """
        fft = self.model.fft
        potential = finite_array(local_potential, (fft.nnr, 2), "local potential")
        u = finite_array(displacement, (6,), "displacement")
        density = np.zeros((fft.nnr, 2), dtype=np.complex128)
        gradient = np.zeros(6, dtype=np.complex128)
        for channel, (plus_map, minus_map) in zip(self.model.channels, self._space.maps):
            v = potential[:, channel.spin]
            weights = channel.weighted_occupations
            plus = self._ket(channel, "plus", v, u, plus_map)
            density[:, channel.spin] += np.einsum(
                "rn,rn,n->r", channel.source_grid.conj(),
                self._space.to_grid(plus, plus_map), weights) / fft.omega
            gradient += np.einsum("atn,tn,n->a", channel.g_plus.conj(), plus, weights)
            minus = self._ket(channel, "minus", v.conj(), u.conj(), minus_map)
            density[:, channel.spin] += np.einsum(
                "rn,rn,n->r", self._space.to_grid(minus, minus_map).conj(),
                channel.source_grid, weights) / fft.omega
            gradient += np.einsum("tn,atn,n->a", minus.conj(), channel.g_minus, weights)
        return PWResponseValue(density, gradient)

    def ionic_matrix(self) -> NDArray[np.complex128]:
        """Six-column physical electronic curvature; no Cfrozen or screening."""
        potential = np.zeros((self.model.fft.nnr, 2), dtype=np.complex128)
        unit = np.eye(6, dtype=np.complex128)
        return np.column_stack([self.apply(potential, u).ionic_gradient for u in unit])

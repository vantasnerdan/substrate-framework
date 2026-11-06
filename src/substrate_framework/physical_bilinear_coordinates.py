"""Complete finite Gamma Hermitian CAR coordinates; no orbit/gap quotient.

For N complete uniform BvK cells X_ka=sqrt(N)*(dGamma(T_ka)-mean).
The density-matrix fluctuation is sum_a T_a X_ka; physical collective
observables use weight 1/N. Ionic u,p are the normalized Gamma collective
coordinates (sqrt(N) times uniform per-cell displacement/momentum), [u,p]=i
in Ry time, hbar=1. These are NOT arbitrary IBZ-weighted independent modes.
General q requires A(q)^dagger=A(-q), actual reciprocal/gauge maps, one
independent real pair for non-self-inverse q, and only one copy at Gamma or
self-inverse boundary q. This implementation covers Gamma transfer only.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


def hermitian_basis(dimension):
    """Trace-orthonormal d^2 matrices: populations, real then imaginary pairs."""
    if not isinstance(dimension, int) or dimension < 1:
        raise ValueError("positive integer dimension required")
    basis = []
    for i in range(dimension):
        matrix = np.zeros((dimension, dimension), complex)
        matrix[i, i] = 1
        basis.append(matrix)
    for i in range(dimension):
        for j in range(i + 1, dimension):
            real = np.zeros((dimension, dimension), complex)
            real[i, j] = real[j, i] = 1 / np.sqrt(2)
            imag = np.zeros_like(real)
            imag[i, j], imag[j, i] = -1j / np.sqrt(2), 1j / np.sqrt(2)
            basis.extend((real, imag))
    return np.asarray(basis)


@dataclass(frozen=True)
class UniformModeIndex:
    """One explicit spin/k channel per full-grid point; no implicit spin factor.

k_indices are integer triples on grid_shape. Every spin has the entire grid
exactly once. Finite-system mode is grid_shape=(1,1,1), one channel per spin.
Labels assert a supplied complete physical grid, not a symmetry reconstruction.
"""
    grid_shape: tuple[int, int, int]
    k_indices: tuple[tuple[int, int, int], ...]

    @property
    def ncells(self):
        return int(np.prod(self.grid_shape))

    def validate(self, channels):
        if (len(self.grid_shape) != 3 or any(
                not isinstance(n, int) or n <= 0 for n in self.grid_shape)):
            raise ValueError("three positive integer BvK dimensions required")
        if len(self.k_indices) != len(channels):
            raise ValueError("one explicit k index per channel required")
        expected = set(np.ndindex(self.grid_shape))
        by_spin = {}
        for channel, k in zip(channels, self.k_indices):
            k = tuple(k)
            if k not in expected or k in by_spin.setdefault(channel.spin, set()):
                raise ValueError("invalid or duplicate spin/k mode")
            by_spin[channel.spin].add(k)
            if channel.weight != 1.0 / self.ncells:
                raise ValueError("only complete uniform weights 1/Ncells are supported")
        if not by_spin or any(points != expected for points in by_spin.values()):
            raise ValueError("each independent spin must contain the complete BvK grid")


class PhysicalBilinearCoordinates:
    """Finite coordinate adapter, [u,p,X_channel]; allocates d^2 basis matrices.

Use matrix-valued probes in physical_observable_transfer for selected
contractions without constructing this finite adapter or a full covariance.
"""
    def __init__(self, channels, nion, modes):
        self.channels = tuple(channels)
        modes.validate(self.channels)
        if not isinstance(nion, int) or nion < 1:
            raise ValueError("positive ionic dimension required")
        self.modes, self.nion = modes, nion
        self.bases = tuple(hermitian_basis(len(c.occupations)) for c in self.channels)
        cursor = 2 * nion
        self.slices = []
        for basis in self.bases:
            self.slices.append(slice(cursor, cursor + len(basis)))
            cursor += len(basis)
        self.size = cursor

    def unpack(self, state):
        state = np.asarray(state, complex)
        if state.shape != (self.size,):
            raise ValueError("incorrect complete-coordinate state shape")
        matrices = tuple(np.einsum("a,aij->ij", state[s], basis)
                         for s, basis in zip(self.slices, self.bases))
        return state[:self.nion], state[self.nion:2*self.nion], matrices

    def pack(self, u, p, matrices):
        if len(matrices) != len(self.channels):
            raise ValueError("one original-basis matrix per channel required")
        parts = [np.asarray(u), np.asarray(p)]
        if any(v.shape != (self.nion,) for v in parts):
            raise ValueError("incorrect ionic state shape")
        for matrix, basis in zip(matrices, self.bases):
            if np.shape(matrix) != basis.shape[1:]:
                raise ValueError("incorrect electronic matrix shape")
            parts.append(np.einsum("aij,ji->a", basis, matrix))
        return np.concatenate(parts).astype(complex)

    def expected_commutator(self):
        """Return <[Z_a,Z_b]>, NOT divided by i; CAR state expectation only."""
        result = np.zeros((self.size, self.size), complex)
        n = self.nion
        result[:n, n:2*n] = 1j * np.eye(n)
        result[n:2*n, :n] = -1j * np.eye(n)
        for channel, basis, section in zip(self.channels, self.bases, self.slices):
            result[section, section] = self.modes.ncells * np.asarray([
                [np.trace(channel.rho0 @ (a @ b - b @ a)) for b in basis]
                for a in basis])
        return result

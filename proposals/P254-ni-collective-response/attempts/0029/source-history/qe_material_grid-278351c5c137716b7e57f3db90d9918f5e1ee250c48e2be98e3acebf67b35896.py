"""Physical serial QE Fourier mapping, independent of fixture storage guards.

Miller columns and Cartesian reciprocal columns are explicit; all flattened
fields are Fortran ordered. Forward transforms are normalized, inverse
transforms reconstruct amplitudes. No spin or cell-volume factor enters FFTs.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np


@dataclass(frozen=True)
class PhysicalRhoFFT:
    physical: tuple[int, int, int]
    miller: np.ndarray
    reciprocal_bohr_inv: np.ndarray
    q_bohr_inv: np.ndarray
    omega: float
    nl: np.ndarray
    wavevectors: np.ndarray
    nnr: int

    @classmethod
    def create(cls, physical, miller, reciprocal_bohr_inv, omega, q_bohr_inv=(0., 0., 0.)):
        shape = tuple(physical)
        if len(shape) != 3 or any(int(n) != n or n < 1 for n in shape):
            raise ValueError("three positive physical grid dimensions required")
        shape = tuple(int(n) for n in shape)
        miller = np.asarray(miller)
        reciprocal = np.asarray(reciprocal_bohr_inv, dtype=float)
        q = np.asarray(q_bohr_inv, dtype=float)
        if miller.ndim != 2 or miller.shape[0] != 3 or miller.dtype.kind not in "iu":
            raise ValueError("integer Miller columns (3,ng) required")
        if reciprocal.shape != (3, 3) or q.shape != (3,) or not np.isfinite(reciprocal).all() or not np.isfinite(q).all():
            raise ValueError("finite physical reciprocal columns and q required")
        if not np.isfinite(omega) or omega <= 0:
            raise ValueError("positive physical cell volume required")
        slots = np.ravel_multi_index(tuple(miller[i] % shape[i] for i in range(3)), shape, order="F")
        if np.unique(slots).size != slots.size:
            raise ValueError("Miller map aliases on the physical FFT grid")
        return cls(shape, miller, reciprocal, q, float(omega), slots,
                   reciprocal @ miller + q[:, None], int(np.prod(shape)))

    def at_q(self, q_bohr_inv):
        q = np.asarray(q_bohr_inv, dtype=float)
        if q.shape != (3,) or not np.isfinite(q).all():
            raise ValueError("finite Cartesian physical q required")
        return replace(self, q_bohr_inv=q, wavevectors=self.reciprocal_bohr_inv @ self.miller + q[:, None])

    def physical_view(self, field):
        field = np.asarray(field)
        if field.shape[0] != self.nnr:
            raise ValueError("first dimension must equal physical grid size")
        return field.reshape((*self.physical, *field.shape[1:]), order="F")

    def to_g(self, field):
        field = np.asarray(field)
        spectrum = np.fft.fftn(self.physical_view(field), axes=(0, 1, 2), norm="forward")
        return spectrum.reshape((self.nnr, *field.shape[1:]), order="F")[self.nl]

    def from_g(self, coefficients):
        coefficients = np.asarray(coefficients)
        if coefficients.shape[0] != self.nl.size:
            raise ValueError("coefficient rows must match the explicit Miller map")
        spectrum = np.zeros((self.nnr, *coefficients.shape[1:]), dtype=np.complex128)
        spectrum[self.nl] = coefficients
        field = np.fft.ifftn(self.physical_view(spectrum), axes=(0, 1, 2), norm="forward")
        return field.reshape(spectrum.shape, order="F")

    def project(self, field):
        return self.from_g(self.to_g(field))

    def gradient(self, field):
        coefficients = self.to_g(field)
        return np.stack([self.from_g(1j * k[:, None] * coefficients) for k in self.wavevectors])

    def divergence(self, vector):
        vector = np.asarray(vector)
        if vector.ndim != 3 or vector.shape[:2] != (3, self.nnr):
            raise ValueError("vector field must have shape (3,nr,nchannels)")
        return self.from_g(sum(1j * self.wavevectors[i, :, None] * self.to_g(vector[i]) for i in range(3)))

    def pair(self, density, potential):
        if np.shape(density) != (self.nnr, 2) or np.shape(potential) != (self.nnr, 2):
            raise ValueError("pairing requires explicit up/down fields")
        return complex(self.omega / self.nnr * np.vdot(density, potential))

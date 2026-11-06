"""Shared full-PW channel validation and Native FFT geometry.

Extracted from the settled independent response: no frequency or time law lives
here. Models and their arrays are immutable for the lifetime of a consumer.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from .pw_response_data import PWResponseChannel, PWResponseModel


def finite_array(value, shape: tuple[int, ...], name: str) -> NDArray:
    array = np.asarray(value)
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be finite with shape {shape}")
    return array


class PWTransitionSpace:
    """Complete target spaces with independently mapped plus/minus branches."""

    def __init__(self, model: PWResponseModel):
        self.model = model
        fft = model.fft
        if not np.isfinite(fft.omega) or fft.omega <= 0:
            raise ValueError("FFT volume must be finite and positive")
        if len(fft.physical) != 3 or len(fft.storage) != 3 or any(
                p < 1 or s < p for p, s in zip(fft.physical, fft.storage)):
            raise ValueError("FFT physical/storage dimensions are inconsistent")
        if fft.nnr != int(np.prod(fft.storage)):
            raise ValueError("FFT nnr must equal the storage size")
        finite_array(model.q_cart_2pi_alat, (3,), "q")
        self.maps = tuple(self._validate_channel(channel) for channel in model.channels)

    def _coordinates(self, slots, count: int, name: str) -> tuple[NDArray, ...]:
        slots = np.asarray(slots)
        fft = self.model.fft
        if slots.shape != (count,) or slots.dtype.kind not in "iu":
            raise ValueError(f"{name} must contain {count} integer FFT offsets")
        if np.any(slots < 0) or np.any(slots >= fft.nnr) or np.unique(slots).size != count:
            raise ValueError(f"{name} contains invalid or duplicated FFT offsets")
        coordinates = np.unravel_index(slots, fft.storage, order="F")
        if any(np.any(c >= n) for c, n in zip(coordinates, fft.physical)):
            raise ValueError(f"{name} maps into FFT padding")
        return coordinates

    def _validate_channel(self, channel: PWResponseChannel):
        if channel.spin not in (0, 1):
            raise ValueError("channel spin must be zero-based 0 or 1")
        energies = np.asarray(channel.source_energies)
        if energies.ndim != 1 or energies.dtype.kind not in "fiu":
            raise ValueError("source energies must be a real vector")
        count = energies.size
        finite_array(energies, (count,), "source energies")
        weights = finite_array(channel.weighted_occupations, (count,), "occupations")
        if weights.dtype.kind not in "fiu" or np.any(weights < 0):
            raise ValueError("weighted occupations must be real and nonnegative")
        finite_array(channel.source_grid, (self.model.fft.nnr, count), "source grid")
        maps = []
        for branch in ("plus", "minus"):
            target = np.asarray(getattr(channel, f"{branch}_energies"))
            if target.ndim != 1 or target.dtype.kind not in "fiu":
                raise ValueError(f"{branch} target energies must be a real vector")
            size = target.size
            finite_array(target, (size,), f"{branch} energies")
            finite_array(getattr(channel, f"{branch}_vectors"), (size, size),
                         f"{branch} full target eigenvectors")
            finite_array(getattr(channel, f"g_{branch}"), (6, size, count),
                         f"{branch} ionic operators")
            maps.append(self._coordinates(getattr(channel, f"{branch}_slots"),
                                          size, f"{branch} slots"))
        return tuple(maps)

    def to_pw(self, grid: NDArray, coordinates: tuple[NDArray, ...]) -> NDArray:
        spectrum = np.fft.fftn(self.model.fft.physical_view(grid),
                               axes=(0, 1, 2), norm="forward")
        return spectrum[coordinates]

    def to_grid(self, coefficients: NDArray, coordinates: tuple[NDArray, ...]) -> NDArray:
        fft = self.model.fft
        spectrum = np.zeros((*fft.physical, coefficients.shape[1]), dtype=np.complex128)
        spectrum[coordinates] = coefficients
        physical = np.fft.ifftn(spectrum, axes=(0, 1, 2), norm="forward")
        return fft._embed(physical)

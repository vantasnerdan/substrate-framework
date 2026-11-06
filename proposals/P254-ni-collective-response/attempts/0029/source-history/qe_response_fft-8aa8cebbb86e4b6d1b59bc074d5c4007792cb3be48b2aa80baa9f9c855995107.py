"""Serial QE Rho FFT actions on actual native G/nl maps and padded storage.

This deliberately supports only the registered dense=smooth 24^3 pilot. It does
not silently replace differing grids, distributed slabs, gamma-packed spectra or
nonperiodic Coulomb operators with a general-grid approximation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .qe_response_bundle import QEResponseBundle


@dataclass(frozen=True)
class NativeRhoFFT:
    physical: tuple[int, int, int]
    storage: tuple[int, int, int]
    nnr: int
    nl: NDArray[Any]
    wavevectors: NDArray[Any]
    omega: float

    @classmethod
    def from_bundle(cls, bundle: QEResponseBundle) -> NativeRhoFFT:
        bg = bundle.background
        meta = bg.metadata
        physical = tuple(meta[f"nr{i}"] for i in (1, 2, 3))
        storage = tuple(meta[f"nr{i}x"] for i in (1, 2, 3))
        for key in ("nnr", "ngm", "nspin", *(f"nr{i}" for i in (1, 2, 3)),
                    *(f"nr{i}x" for i in (1, 2, 3))):
            if type(meta[key]) is not int or meta[key] < 1:
                raise ValueError(f"background {key} must be a positive native integer")
        if physical != (24, 24, 24) or any(s < p for s, p in zip(storage, physical)):
            raise ValueError("kernel scope is physical 24^3 with explicit serial FFT padding")
        if meta["nnr"] != int(np.prod(storage)):
            raise ValueError("distributed/nonrectangular FFT storage is outside this pilot")
        for key in ("nnr", *(f"nr{i}" for i in (1, 2, 3)),
                    *(f"nr{i}x" for i in (1, 2, 3))):
            if meta[f"smooth_{key}"] != meta[key]:
                raise ValueError("registered kernel requires identical actual dense/smooth grids")
        if meta["nspin"] != 2 or meta["functional"] != "PBE":
            raise ValueError("registered kernel requires collinear two-spin PBE")
        for key in ("gamma_only", "do_comp_mt", "do_cutoff_2D", "noncolin", "nlcc_any"):
            if meta[key] is not False:
                raise ValueError(f"unsupported native flag {key}: expected false")
        # lgamma identifies q=0; gamma_only instead changes PW/FFT storage.
        for key in ("lgamma", "lrpa", "lnolr"):
            if type(meta[key]) is not bool:
                raise ValueError(f"native flag {key} must be a JSON boolean")
        for key in ("omega", "alat", "tpiba", "e2"):
            if not isinstance(meta[key], (int, float)) or not np.isfinite(meta[key]) or meta[key] <= 0:
                raise ValueError(f"invalid native scalar {key}")
        if not np.isclose(meta["tpiba"] * meta["alat"], 2 * np.pi, rtol=1e-13, atol=0):
            raise ValueError("native alat/tpiba units disagree")
        g = bg.array("g", (3, meta["ngm"]), "float64")
        nl = bg.array("nl", (meta["ngm"],), "int32")
        if bg.descriptors["nl"].index_base != 1:
            raise ValueError("nl must explicitly declare native one-based indices")
        slots = nl.astype(np.int64) - 1
        if np.any(slots < 0) or np.any(slots >= meta["nnr"]) or np.unique(slots).size != slots.size:
            raise ValueError("invalid or duplicated native nl indices")
        coords = np.unravel_index(slots, storage, order="F")
        if any(np.any(c >= n) for c, n in zip(coords, physical)):
            raise ValueError("native reciprocal map points into padding")
        q = bg.array("q_cart_2pi_alat", (3,), "float64")
        return cls(physical, storage, meta["nnr"], slots,
                   meta["tpiba"] * (g + q[:, None]), float(meta["omega"]))

    def physical_view(self, field: NDArray[Any]) -> NDArray[Any]:
        if field.shape[0] != self.nnr:
            raise ValueError("density/potential first dimension must equal native nnr")
        shaped = field.reshape((*self.storage, *field.shape[1:]), order="F")
        return shaped[:self.physical[0], :self.physical[1], :self.physical[2], ...]

    def _embed(self, physical: NDArray[Any]) -> NDArray[Any]:
        padded = np.zeros((*self.storage, *physical.shape[3:]), dtype=np.complex128)
        padded[:self.physical[0], :self.physical[1], :self.physical[2], ...] = physical
        return padded.reshape((self.nnr, *physical.shape[3:]), order="F")

    def to_g(self, field: NDArray[Any]) -> NDArray[Any]:
        # QE fwfft is normalized; invfft reconstructs coefficient amplitudes.
        spectrum = np.fft.fftn(self.physical_view(field), axes=(0, 1, 2), norm="forward")
        return self._embed(spectrum)[self.nl]

    def from_g(self, coefficients: NDArray[Any]) -> NDArray[Any]:
        if coefficients.shape[0] != self.nl.size:
            raise ValueError("spectral coefficients must have native ngm rows")
        padded = np.zeros((self.nnr, *coefficients.shape[1:]), dtype=np.complex128)
        padded[self.nl] = coefficients
        physical = np.fft.ifftn(self.physical_view(padded), axes=(0, 1, 2), norm="forward")
        return self._embed(physical)

    def project(self, field: NDArray[Any]) -> NDArray[Any]:
        return self.from_g(self.to_g(field))

    def gradient(self, field: NDArray[Any]) -> NDArray[Any]:
        coefficients = self.to_g(field)
        return np.stack([self.from_g(1j * k[:, None] * coefficients)
                         for k in self.wavevectors], axis=0)

    def divergence(self, vector: NDArray[Any]) -> NDArray[Any]:
        if vector.shape != (3, self.nnr, 2):
            raise ValueError("spin vector field must have shape (3,nnr,2)")
        coefficients = sum(1j * self.wavevectors[i, :, None] * self.to_g(vector[i])
                           for i in range(3))
        return self.from_g(coefficients)

    def pair(self, density: NDArray[Any], potential: NDArray[Any]) -> complex:
        if density.shape != (self.nnr, 2) or potential.shape != density.shape:
            raise ValueError("pairing requires up/down arrays of shape (nnr,2)")
        # Spin channels are explicit: there is no additional spin factor two.
        return complex(self.omega / np.prod(self.physical) *
                       np.vdot(self.physical_view(density), self.physical_view(potential)))

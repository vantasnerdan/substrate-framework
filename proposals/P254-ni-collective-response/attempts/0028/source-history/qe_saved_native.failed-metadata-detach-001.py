"""Bound QE 7.5 serial saved records, without gauge or normalization repair.

Supported ABI: little-endian int32/logical32, float64/complex128, positive
four-byte Fortran record markers. Full-complex scalar orbitals only. WFC k/B
are already bohr^-1 (PW/pw_restart_new.f90:1027-1031); density B is likewise
physical (PW/io_rho_xml.f90:63-66). ngw is a GLOBAL G-index bound, not igwx.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
import struct
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np

MAX_SOURCE_BYTES = 1 << 30
ABI = "little-endian:marker4:int32:logical32:float64:complex128"


def read_bound_bytes(binding: Mapping[str, Any], *, max_bytes: int = MAX_SOURCE_BYTES):
    """Verify the actual immutable decoding buffer, with a bounded file read.

    Follows qe_frozen_curvature.verify_buffer's exact count/SHA256 contract,
    without importing its unrelated phonon XML dependency. Counts stay integers.
    """
    size = int(binding["bytes"])
    if str(size) != str(binding["bytes"]) or not 0 < size <= max_bytes:
        raise ValueError("invalid or unsupported source byte count")
    with Path(binding["path"]).open("rb") as stream:
        data = stream.read(size + 1)
    digest = hashlib.sha256(data).hexdigest()
    if len(data) != size or digest != binding["sha256"]:
        raise ValueError(f"source byte count/SHA256 mismatch: {binding['path']}")
    source = MappingProxyType(dict(path=str(binding["path"]), bytes=str(size), sha256=digest))
    return data, source


class _Records:
    def __init__(self, data: bytes):
        self.data = data
        self.offset = 0

    def take(self, expected: int) -> memoryview:
        start = self.offset
        if start + 8 + expected > len(self.data):
            raise ValueError("truncated QE saved record")
        size = struct.unpack_from("<i", self.data, start)[0]
        if size != expected:
            raise ValueError(f"unsupported ABI or malformed record: {size} != {expected}")
        stop = start + 4 + size
        if struct.unpack_from("<i", self.data, stop)[0] != size:
            raise ValueError("malformed Fortran trailing record marker")
        self.offset = stop + 4
        return memoryview(self.data)[start + 4:stop]

    def finish(self):
        if self.offset != len(self.data):
            raise ValueError("unexpected trailing saved records")


def _full_complex(logical: int) -> bool:
    if logical not in (0, 1):
        raise ValueError("unsupported logical32 ABI")
    if logical:
        raise ValueError("gamma-compressed records are unsupported; no implicit expansion")
    return False


def _array(record, dtype, shape):
    result = np.frombuffer(record, dtype=dtype).reshape(shape, order="F")
    if result.dtype.kind in "fc" and not np.isfinite(result).all():
        raise ValueError("nonfinite saved numeric record")
    return result


@dataclass(frozen=True)
class SavedWavefunctions:
    ik: int                         # original spin-local one-based index
    ispin: int                      # explicit original one-based spin
    gamma_only: bool
    scalef: float
    ngw: int                        # maximum GLOBAL G index, not row count
    igwx: int                       # compact physical PW coefficient count
    npol: int
    nbnd: int
    k_cart_bohr_inv: np.ndarray      # (3,)
    reciprocal_bohr_inv: np.ndarray # (3,3), reciprocal vectors are columns
    miller: np.ndarray               # (3,igwx), original int32
    coefficients: np.ndarray         # (igwx,nbnd), all saved bands, complex128
    source: Mapping[str, str]
    abi: str = ABI


def load_wavefunctions(binding: Mapping[str, Any]) -> SavedWavefunctions:
    data, source = read_bound_bytes(binding)
    records = _Records(data)
    ik, kx, ky, kz, spin, logical, scale = struct.unpack("<i3diid", records.take(44))
    gamma = _full_complex(logical)
    ngw, count, npol, nbnd = struct.unpack("<4i", records.take(16))
    if ik < 1 or spin not in (1, 2) or not np.isfinite([kx, ky, kz, scale]).all():
        raise ValueError("invalid scalar WFC k/spin/header")
    if npol != 1:
        raise ValueError("spinor WFC unsupported")
    if scale != 1.0:
        raise ValueError("nonunit WFC scalef unsupported; no normalization repair")
    if count < 1 or nbnd < 1 or ngw < count:
        raise ValueError("invalid WFC ngw/igwx/band dimensions")
    expected_bytes = 52 + 24 + 80 + 8 + 12 * count + nbnd * (8 + 16 * count)
    if len(data) != expected_bytes:
        raise ValueError("WFC dimensions disagree with bounded record byte count")
    # Long-lived bases retain this small metadata, not the whole WFC byte buffer.
    reciprocal = _array(records.take(72), "<f8", (3, 3)).copy(order="F")
    miller = _array(records.take(12 * count), "<i4", (3, count)).copy(order="F")
    reciprocal.setflags(write=False)
    miller.setflags(write=False)
    coefficients = np.empty((count, nbnd), dtype=np.complex128, order="F")
    for band in range(nbnd):
        coefficients[:, band] = _array(records.take(16 * count), "<c16", (count,))
    records.finish()
    coefficients.setflags(write=False)
    k = np.array([kx, ky, kz]); k.setflags(write=False)
    return SavedWavefunctions(ik, spin, gamma, scale, ngw, count, npol, nbnd,
                             k, reciprocal, miller, coefficients, source)


@dataclass(frozen=True)
class SavedDensity:
    gamma_only: bool
    ngm: int
    nspin: int
    reciprocal_bohr_inv: np.ndarray # (3,3), columns, already physical
    miller: np.ndarray               # (3,ngm), original int32
    components: np.ndarray           # (ngm,nspin), Fourier density bohr^-3
    component_names: tuple[str, ...] # total then magnetization; NOT up/down
    source: Mapping[str, str]
    abi: str = ABI


def load_density(binding: Mapping[str, Any]) -> SavedDensity:
    data, source = read_bound_bytes(binding)
    records = _Records(data)
    logical, count, nspin = struct.unpack("<3i", records.take(12))
    gamma = _full_complex(logical)
    if nspin not in (1, 2):
        raise ValueError("noncollinear density components unsupported")
    if count < 1 or len(data) != 20 + 80 + 8 + 12 * count + nspin * (8 + 16 * count):
        raise ValueError("density dimensions disagree with bounded record byte count")
    reciprocal = _array(records.take(72), "<f8", (3, 3))
    miller = _array(records.take(12 * count), "<i4", (3, count))
    components = np.empty((count, nspin), dtype=np.complex128, order="F")
    for spin in range(nspin):
        components[:, spin] = _array(records.take(16 * count), "<c16", (count,))
    records.finish()
    components.setflags(write=False)
    names = ("total", "magnetization")[:nspin]
    return SavedDensity(gamma, count, nspin, reciprocal, miller, components, names, source)

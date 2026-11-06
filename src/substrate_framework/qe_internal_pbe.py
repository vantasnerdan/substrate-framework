"""Byte-bound original QE 7.5 internal PBE, collinear ordinary NC without NLCC.

Inputs use physical bohr units: rho(total,magnetization)[nr,2] and spin
up/down gradients[3,nr,2]. Full XC is local_xc_ry - div(gradient_flux_ry_bohr).
The five coefficient arrays retain Native LR spin axes, thresholds, finite
differences and Ry conversion; no further factor of two is appropriate.
Noncollinear, NLCC, Libxc, hybrids and meta-GGA are not supported here.
"""
from __future__ import annotations

import ctypes
from dataclasses import dataclass, fields
import fcntl
import os
import threading
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .qe_saved_native import read_bound_bytes

_CALL_LOCK = threading.Lock()  # XClib has mutable threshold/functional globals.
_DOUBLE = ctypes.POINTER(ctypes.c_double)


@dataclass(frozen=True)
class NativePBEResult:
    """Caller-owned, Fortran contiguous float64 arrays in Native axis order."""

    local_xc_ry: np.ndarray
    gradient_flux_ry_bohr: np.ndarray
    dmuxc: np.ndarray
    dvxc_rr: np.ndarray
    dvxc_sr: np.ndarray
    dvxc_ss: np.ndarray
    dvxc_s: np.ndarray

    @classmethod
    def allocate(cls, nr: int) -> NativePBEResult:
        return cls(*(np.empty(shape, dtype=np.float64, order="F")
                     for shape in _shapes(nr)))


def _shapes(nr):
    return ((nr, 2), (3, nr, 2)) + ((nr, 2, 2),) * 5


def _output_buffer(array, shape):
    if (not isinstance(array, np.ndarray) or array.dtype != np.float64
            or array.shape != shape or not array.flags.f_contiguous
            or not array.flags.aligned or not array.flags.writeable):
        raise ValueError(f"output/workspace must be writable aligned F-float64 {shape}")
    return array


def _input_buffer(value, shape):
    array = np.asarray(value)
    if array.shape != shape or array.dtype.kind not in "fiu":
        raise ValueError(f"expected real numeric input with shape {shape}")
    array = np.require(array, dtype=np.float64, requirements=["F", "A"])
    if not np.isfinite(array).all():
        raise ValueError("nonfinite PBE input")
    return array


def _sealed_library(data):
    """dlopen the verified buffer, never reopen an unbound mutable source path."""
    fd = os.memfd_create("qe-internal-pbe", os.MFD_CLOEXEC | os.MFD_ALLOW_SEALING)
    try:
        with os.fdopen(os.dup(fd), "wb") as stream:
            stream.write(data)
        fcntl.fcntl(fd, fcntl.F_ADD_SEALS, fcntl.F_SEAL_WRITE | fcntl.F_SEAL_GROW
                    | fcntl.F_SEAL_SHRINK | fcntl.F_SEAL_SEAL)
        library = ctypes.CDLL(f"/proc/self/fd/{fd}", mode=os.RTLD_LOCAL | os.RTLD_NOW)
    except BaseException:
        os.close(fd)
        raise
    return fd, library


class NativePBE:
    """Load a trusted Root-built wrapper bound by path/bytes/SHA256.

    This Linux-only loader pins a sealed copy of the verified ELF bytes. The
    build receipt supplies provenance for original archive and linked runtime
    dependencies; a content hash alone is not an authorization to execute code.
    Calls are serialized; input buffers must not be mutated concurrently.
    Repeated calls may reuse ``out`` and ``workspace[nr,19]``. The adapter
    allocates no compiled scratch; original XClib routines retain their own
    allocations. The FD stays live with the CDLL to prevent dlopen path reuse.
    """

    def __init__(self, bound_library: Mapping[str, Any]):
        data, self.source = read_bound_bytes(bound_library, max_bytes=64 << 20)
        if not data.startswith(b"\x7fELF"):
            raise ValueError("expected shared ELF wrapper, not a static archive")
        dependencies = bound_library.get("dependencies")
        if not isinstance(dependencies, list) or not dependencies:
            raise ValueError("bound_library requires the build receipt's runtime dependencies")
        self.dependencies = tuple(read_bound_bytes(item, max_bytes=256 << 20)[1]
                                  for item in dependencies)
        fd, self._library = _sealed_library(data)
        # Tie descriptor lifetime to CDLL, including ctypes function references.
        import weakref
        self._fd_finalizer = weakref.finalize(self._library, os.close, fd)
        mapped = {str(Path(line.split(maxsplit=5)[5].strip()).resolve())
                  for line in Path("/proc/self/maps").read_text().splitlines()
                  if len(line.split(maxsplit=5)) == 6 and line.split(maxsplit=5)[5].startswith("/")}
        if any(str(Path(item["path"]).resolve()) not in mapped for item in self.dependencies):
            raise RuntimeError("loaded runtime closure differs from the bound build receipt")
        abi = self._library.ni_pbe_abi
        abi.argtypes = []
        abi.restype = ctypes.c_int
        if abi() != 1:
            raise ValueError("unsupported Native PBE wrapper ABI")
        self._evaluate = self._library.ni_pbe_evaluate
        self._evaluate.argtypes = [ctypes.c_int] + [_DOUBLE] * 10 + [ctypes.POINTER(ctypes.c_int)]
        self._evaluate.restype = None

    def evaluate(self, rho_total_magnetization, grho_updown, *, out=None, workspace=None):
        """Return all registered fields without clipping, shifts or Hermitization.

        Default outputs are fresh. Optional output and scratch buffers must be
        disjoint from each other and inputs; returned arrays alias ``out``.
        Native's own near-zero/negative-density conventions remain untouched.
        """
        shape = np.shape(rho_total_magnetization)
        if len(shape) != 2 or shape[1] != 2 or not 0 < shape[0] <= (2**31-1)//8:
            raise ValueError("expected nonempty rho[nr,2] within Native int32 indexing")
        nr = shape[0]
        rho = _input_buffer(rho_total_magnetization, (nr, 2))
        gradient = _input_buffer(grho_updown, (3, nr, 2))
        if out is None:
            out = NativePBEResult.allocate(nr)
        if not isinstance(out, NativePBEResult):
            raise TypeError("out must be NativePBEResult")
        outputs = [_output_buffer(getattr(out, field.name), expected)
                   for field, expected in zip(fields(out), _shapes(nr))]
        if workspace is None:
            workspace = np.empty((nr, 19), dtype=np.float64, order="F")
        workspace = _output_buffer(workspace, (nr, 19))
        arrays = [rho, gradient, *outputs, workspace]
        for i in range(2, len(arrays)):
            if any(np.shares_memory(arrays[i], prior) for prior in arrays[:i]):
                raise ValueError("PBE outputs/workspace must not overlap other buffers")
        status = ctypes.c_int(-1)
        with _CALL_LOCK:
            self._evaluate(nr, *(a.ctypes.data_as(_DOUBLE) for a in arrays), ctypes.byref(status))
        if status.value:
            raise RuntimeError(f"Native PBE wrapper rejected call: status {status.value}")
        if any(not np.isfinite(a).all() for a in outputs):
            raise FloatingPointError("Native PBE returned nonfinite fields; no repair applied")
        return out

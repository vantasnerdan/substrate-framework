"""Exact Schwinger--Hopf algebra for a supplied pair of Euler KKS modes.

The functions in this module assume that two positive, equally normalized
canonical modes have already been constructed.  They do not establish that
Euler supplies an invariant mode doublet, a physical analyzer, quantization,
or a measurement law.
"""

from __future__ import annotations

from collections.abc import Sequence
from math import pi, sqrt

import numpy as np


def canonical_doublet(
    coordinates: Sequence[float], *, kks_scale: float
) -> np.ndarray:
    """Return ``z_a=sqrt(B/2)(q_a+i p_a)`` for two canonical pairs."""

    if kks_scale <= 0:
        raise ValueError("kks_scale must be strictly positive")
    values = np.asarray(coordinates, dtype=float)
    if values.shape != (4,) or not np.all(np.isfinite(values)):
        raise ValueError("coordinates must be four finite values (q1,p1,q2,p2)")
    q1, p1, q2, p2 = values
    factor = sqrt(kks_scale / 2.0)
    return factor * np.array([q1 + 1j * p1, q2 + 1j * p2])


def total_action(z: Sequence[complex]) -> float:
    """Return ``I=(|z1|^2+|z2|^2)/2``."""

    vector = _doublet(z)
    return float(np.vdot(vector, vector).real / 2.0)


def stokes_vector(z: Sequence[complex]) -> np.ndarray:
    """Return ``S_i=z^dagger sigma_i z/2`` in the x, y, z convention."""

    z1, z2 = _doublet(z)
    cross = np.conjugate(z1) * z2
    return np.array(
        [cross.real, cross.imag, (abs(z1) ** 2 - abs(z2) ** 2) / 2.0],
        dtype=float,
    )


def hopf_identity_residual(z: Sequence[complex]) -> float:
    """Return ``|S|^2-I^2``, which vanishes for every supplied doublet."""

    stokes = stokes_vector(z)
    action = total_action(z)
    return float(np.dot(stokes, stokes) - action**2)


def reduced_kks_area(action: float) -> float:
    """Return the area ``4*pi*I`` of the reduced radius-I KKS sphere."""

    if action <= 0:
        raise ValueError("action must be strictly positive")
    return 4.0 * pi * action


def apply_mode_mixing(z: Sequence[complex], mixing: Sequence[Sequence[complex]]) -> np.ndarray:
    """Apply a supplied two-mode unitary after checking its exact domain."""

    vector = _doublet(z)
    matrix = np.asarray(mixing, dtype=complex)
    if matrix.shape != (2, 2) or not np.all(np.isfinite(matrix)):
        raise ValueError("mixing must be a finite 2x2 matrix")
    if not np.allclose(matrix.conj().T @ matrix, np.eye(2), atol=1e-12, rtol=1e-12):
        raise ValueError("mixing must be unitary")
    return matrix @ vector


def _doublet(z: Sequence[complex]) -> np.ndarray:
    vector = np.asarray(z, dtype=complex)
    if vector.shape != (2,) or not np.all(np.isfinite(vector)):
        raise ValueError("z must be a finite complex doublet")
    return vector

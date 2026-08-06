"""Localized conserved 3+1 sources built from sine-Gordon channel stress."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .dimensional_sine_gordon import (
    DimensionalSineGordonCoefficients,
    DimensionalSineGordonStress,
    dimensional_sine_gordon_stress,
)


@dataclass(frozen=True)
class CanonicalFiberStress:
    """A localized 3+1 distribution built from one exact 1+1 stress."""

    sine_gordon_stress: DimensionalSineGordonStress
    transverse_density: sp.Expr
    stress_energy: sp.ImmutableMatrix
    divergence: sp.ImmutableMatrix


def canonical_sine_gordon_fiber_stress(
    field: Any,
    longitudinal_coordinate: sp.Symbol,
    time: sp.Symbol,
    transverse_y: sp.Symbol,
    transverse_z: sp.Symbol,
    coefficients: DimensionalSineGordonCoefficients,
) -> CanonicalFiberStress:
    r"""Embed a sine-Gordon channel as an isolated conserved 3+1 source.

    The channel lies along ``x`` and has transverse density
    ``delta(y)*delta(z)``.  In orthonormal coordinates ``(c*t,x,y,z)``, the
    only nonzero block is the accepted physical 1+1 stress.  Its four-
    divergence is therefore
    ``delta(y)delta(z)*(u_t*R/c,-u_x*R,0,0)`` and vanishes exactly on shell.
    """

    for coordinate, name in (
        (longitudinal_coordinate, "longitudinal_coordinate"),
        (time, "time"),
        (transverse_y, "transverse_y"),
        (transverse_z, "transverse_z"),
    ):
        if not isinstance(coordinate, sp.Symbol):
            raise ValueError(f"{name} must be a SymPy Symbol; got {coordinate!r}")
    stress = dimensional_sine_gordon_stress(
        field,
        longitudinal_coordinate,
        time,
        coefficients,
    )
    transverse_density = sp.DiracDelta(transverse_y) * sp.DiracDelta(transverse_z)
    embedded = sp.zeros(4)
    embedded[:2, :2] = stress.contravariant
    stress_energy = sp.ImmutableMatrix(embedded * transverse_density)
    divergence = sp.ImmutableMatrix(
        [
            stress.divergence[0] * transverse_density,
            stress.divergence[1] * transverse_density,
            0,
            0,
        ]
    )
    return CanonicalFiberStress(
        sine_gordon_stress=stress,
        transverse_density=transverse_density,
        stress_energy=stress_energy,
        divergence=divergence,
    )

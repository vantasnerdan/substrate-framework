"""Exact induced geometry for narrowly named Lorentz-group orbits.

Authority status: ``unit_timelike_vector_orbit_metric`` implements accepted
claim C-LOR-001 in release v0.160.0. The volume helper remains conditional
unpromoted infrastructure linked to open goal issue #28. The current public
surface describes only the Lorentz orbit of one future-directed unit timelike
vector, namely ``H^2`` or ``H^3``. It does not describe a timelike two-plane or
tube worldsheet orbit and makes no claim about a Lorentz-invariant tube
ensemble.
"""

from __future__ import annotations

from typing import Sequence

import sympy as sp


def _induced_metric(
    embedding: sp.Matrix,
    parameters: Sequence[sp.Symbol],
    signature: Sequence[int],
) -> sp.Matrix:
    """Pull back a constant ambient metric to the supplied parameters."""

    dimension = len(parameters)
    ambient_metric = sp.diag(*signature)
    return sp.Matrix(
        dimension,
        dimension,
        lambda i, j: sp.trigsimp(
            (
                embedding.diff(parameters[i]).T
                * ambient_metric
                * embedding.diff(parameters[j])
            )[0]
        ),
    )


def unit_timelike_vector_orbit_metric(spacetime_dimension: int = 4) -> sp.Matrix:
    """Return the induced metric on the future unit-vector orbit.

    ``spacetime_dimension=3`` gives ``H^2`` and the previous default of four
    gives ``H^3``. Both are induced directly from the standard hyperboloid in
    mostly-plus ambient signature.
    """

    if spacetime_dimension not in (3, 4):
        raise ValueError("spacetime_dimension must be 3 or 4")
    eta, theta, phi = sp.symbols("eta theta phi", positive=True, real=True)
    if spacetime_dimension == 3:
        unit_vector = sp.Matrix(
            [
                sp.cosh(eta),
                sp.sinh(eta) * sp.cos(theta),
                sp.sinh(eta) * sp.sin(theta),
            ]
        )
        return _induced_metric(unit_vector, [eta, theta], [-1, 1, 1])
    unit_vector = sp.Matrix(
        [
            sp.cosh(eta),
            sp.sinh(eta) * sp.sin(theta) * sp.cos(phi),
            sp.sinh(eta) * sp.sin(theta) * sp.sin(phi),
            sp.sinh(eta) * sp.cos(theta),
        ]
    )
    return _induced_metric(unit_vector, [eta, theta, phi], [-1, 1, 1, 1])


def unit_timelike_vector_orbit_volume(spacetime_dimension: int = 4) -> sp.Expr:
    """Return the infinite invariant volume of ``H^2`` or ``H^3``.

    The radial factors are ``sinh(eta)`` for ``H^2`` and ``sinh(eta)^2`` for
    ``H^3`` (with their standard angular measures). These results apply only
    to a unit timelike vector, not to a tube worldsheet orbit.
    """

    if spacetime_dimension not in (3, 4):
        raise ValueError("spacetime_dimension must be 3 or 4")
    radial_cutoff = sp.Symbol("radial_cutoff", positive=True, real=True)
    if spacetime_dimension == 3:
        radial_volume = sp.limit(
            sp.cosh(radial_cutoff) - 1,
            radial_cutoff,
            sp.oo,
        )
        return 2 * sp.pi * radial_volume
    radial_volume = sp.limit(
        sp.sinh(2 * radial_cutoff) / 4 - radial_cutoff / 2,
        radial_cutoff,
        sp.oo,
    )
    theta, phi = sp.symbols("theta phi", positive=True)
    angular_volume = sp.integrate(
        sp.sin(theta),
        (theta, 0, sp.pi),
        (phi, 0, 2 * sp.pi),
    )
    return radial_volume * angular_volume

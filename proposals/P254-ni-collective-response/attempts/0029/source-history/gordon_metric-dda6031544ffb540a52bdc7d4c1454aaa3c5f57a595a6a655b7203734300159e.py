"""Exact Gordon effective metrics in the framework's mostly-plus convention.

For ``eta = diag(-1, 1, 1, 1)`` and a unit timelike medium velocity
``eta(u, u) = -1``, the signature-consistent Gordon inverse metric is

``g^ab = eta^ab + (1 - n**2) u^a u^b``.

Its inverse is ``g_ab = eta_ab + (1 - 1/n**2) u_a u_b``.  These formulas are
the mostly-plus translation of the commonly quoted mostly-minus convention;
mixing the two signs creates a spurious pole at ``n=sqrt(2)``.

The metric is an effective geometry for wave propagation in a declared medium.
This module does not infer an Einstein source, physical gravity, a material
action, or a coupled medium solution from its curvature.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


MINKOWSKI_MOSTLY_PLUS = sp.diag(-1, 1, 1, 1)


def _exact_real(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.has(sp.Float):
        raise ValueError(f"{name} must be exact")
    if result.is_real is not True:
        raise ValueError(f"{name} must be declared real")
    return sp.simplify(result)


def _positive_exact(value: Any, name: str) -> sp.Expr:
    result = _exact_real(value, name)
    if result.is_positive is not True:
        raise ValueError(f"{name} must be declared positive")
    return result


@dataclass(frozen=True)
class GordonMetric:
    """Exact covariant and contravariant effective metric data."""

    refractive_index: sp.Expr
    four_velocity_up: sp.Matrix
    four_velocity_down: sp.Matrix
    contravariant: sp.Matrix
    covariant: sp.Matrix
    contravariant_determinant: sp.Expr
    covariant_determinant: sp.Expr
    rest_phase_speed: sp.Expr


@dataclass(frozen=True)
class TransverseProfileEinstein:
    """Einstein tensor for constant z velocity and ``n=n(x)``.

    Coordinates are ordered ``(t, x, y, z)``.  The result follows by deriving
    the rest-frame curvature and applying the constant Lorentz transformation;
    ``x`` is transverse to the boost and is unchanged.
    """

    metric: GordonMetric
    coordinate: sp.Symbol
    velocity: sp.Expr
    gamma_squared: sp.Expr
    curvature_kernel: sp.Expr
    einstein_covariant: sp.Matrix
    ricci_scalar: sp.Expr


def gordon_metric_mostly_plus(
    refractive_index: Any,
    four_velocity_up: Any,
) -> GordonMetric:
    """Return the exact Gordon metric for a unit timelike four-velocity.

    ``refractive_index`` must be exact, real, and provably positive.  The four
    vector must have four exact real entries and satisfy ``eta(u,u)=-1``.
    """

    index = _positive_exact(refractive_index, "refractive_index")
    velocity = sp.Matrix(four_velocity_up)
    if velocity.shape not in {(4, 1), (1, 4)}:
        raise ValueError("four_velocity_up must contain exactly four entries")
    velocity = velocity.reshape(4, 1)
    for component in velocity:
        _exact_real(component, "four_velocity_up")
    velocity = velocity.applyfunc(sp.simplify)

    norm = sp.simplify((velocity.T * MINKOWSKI_MOSTLY_PLUS * velocity)[0])
    if norm != -1:
        raise ValueError("four_velocity_up must have mostly-plus norm -1")

    velocity_down = MINKOWSKI_MOSTLY_PLUS * velocity
    contravariant = (
        MINKOWSKI_MOSTLY_PLUS + (1 - index**2) * velocity * velocity.T
    ).applyfunc(sp.simplify)
    covariant = (
        MINKOWSKI_MOSTLY_PLUS
        + (1 - index**-2) * velocity_down * velocity_down.T
    ).applyfunc(sp.simplify)
    inverse_residual = (contravariant * covariant - sp.eye(4)).applyfunc(
        sp.simplify
    )
    if inverse_residual != sp.zeros(4):
        raise AssertionError("closed Gordon forms are not mutual inverses")

    determinant_up = sp.simplify(contravariant.det())
    determinant_down = sp.simplify(covariant.det())
    if sp.simplify(determinant_up + index**2) != 0:
        raise AssertionError("unexpected contravariant Gordon determinant")
    if sp.simplify(determinant_down + index**-2) != 0:
        raise AssertionError("unexpected covariant Gordon determinant")

    return GordonMetric(
        refractive_index=index,
        four_velocity_up=velocity,
        four_velocity_down=velocity_down,
        contravariant=contravariant,
        covariant=covariant,
        contravariant_determinant=determinant_up,
        covariant_determinant=determinant_down,
        rest_phase_speed=sp.simplify(1 / index),
    )


def transverse_profile_einstein(
    refractive_index: Any,
    coordinate: sp.Symbol,
    velocity: Any,
) -> TransverseProfileEinstein:
    """Return exact curvature for a uniform z flow and transverse index.

    The declared four-velocity is ``gamma*(1,0,0,v)`` with ``|v|<1`` and the
    index depends only on ``coordinate=x``.  Define
    ``K=(n*n.diff(x,2)-2*n.diff(x)**2)/n**2``.  In the medium rest frame the
    only nonzero covariant Einstein components are ``G_yy=G_zz=-K``.  A
    constant z boost gives the returned laboratory-frame tensor.
    """

    if not isinstance(coordinate, sp.Symbol) or coordinate.is_real is not True:
        raise ValueError("coordinate must be a real SymPy Symbol")
    index = _positive_exact(refractive_index, "refractive_index")
    speed = _exact_real(velocity, "velocity")
    speed_squared = sp.simplify(speed**2)
    subluminal_margin = sp.simplify(1 - speed_squared)
    if subluminal_margin.is_positive is not True:
        raise ValueError("velocity must satisfy |velocity| < 1 exactly")

    gamma_squared = sp.simplify(1 / subluminal_margin)
    gamma = sp.sqrt(gamma_squared)
    four_velocity = sp.Matrix([gamma, 0, 0, gamma * speed])
    metric = gordon_metric_mostly_plus(index, four_velocity)
    kernel = sp.simplify(
        (index * sp.diff(index, coordinate, 2) - 2 * sp.diff(index, coordinate) ** 2)
        / index**2
    )
    einstein = sp.zeros(4)
    einstein[0, 0] = sp.simplify(-gamma_squared * speed_squared * kernel)
    einstein[0, 3] = sp.simplify(gamma_squared * speed * kernel)
    einstein[3, 0] = einstein[0, 3]
    einstein[2, 2] = sp.simplify(-kernel)
    einstein[3, 3] = sp.simplify(-gamma_squared * kernel)

    return TransverseProfileEinstein(
        metric=metric,
        coordinate=coordinate,
        velocity=speed,
        gamma_squared=gamma_squared,
        curvature_kernel=kernel,
        einstein_covariant=einstein,
        ricci_scalar=sp.simplify(2 * kernel),
    )

"""Conditional rational-map geometry and sphere cubature.

The APIs in this module act on separately declared maps of the unit Riemann
sphere.  They do not select a minimizing map, derive a Skyrme action or radial
profile, identify a physical baryon or nucleus, or supply a reaction yield.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import math
from typing import Any

import numpy as np
from numpy.typing import NDArray
import sympy as sp

ComplexArray = NDArray[np.complex128]


@dataclass(frozen=True)
class RationalMapDegreeEvidence:
    """Exact coprimality reduction and algebraic degree of ``p(z)/q(z)``."""

    numerator: sp.Expr
    denominator: sp.Expr
    common_factor: sp.Expr
    reduced_numerator: sp.Expr
    reduced_denominator: sp.Expr
    numerator_degree: int
    denominator_degree: int
    common_factor_degree: int
    degree: int

    @property
    def is_coprime(self) -> bool:
        """Return whether the supplied polynomial representation was reduced."""

        return self.common_factor_degree == 0


@dataclass(frozen=True)
class RationalMapSphereEvidence:
    """Resolution-bounded Gauss-Legendre evidence for one declared map."""

    numerator_coefficients: ComplexArray
    denominator_coefficients: ComplexArray
    declared_degree: int
    polar_order: int
    azimuthal_order: int
    normalized_area: float
    angular_integral: float
    maximum_conformal_jacobian: float
    minimum_homogeneous_norm_squared: float
    coordinate: str = "u=cos(theta), phi"
    method: str = "tensor_gauss_legendre"
    precision: str = "IEEE-754 binary64/complex128"

    @property
    def degree_area_relative_error(self) -> float:
        """Return the relative pullback-area error against the declared degree."""

        if self.declared_degree == 0:
            return abs(self.normalized_area)
        return abs(self.normalized_area - self.declared_degree) / self.declared_degree


def _exact_polynomial(coefficients: Sequence[Any], variable: sp.Symbol, name: str) -> sp.Poly:
    values = tuple(sp.sympify(value) for value in coefficients)
    if not values:
        raise ValueError(f"{name} must contain at least one coefficient")
    if any(value.has(sp.Float) for value in values):
        raise ValueError(f"{name} coefficients must be exact")
    polynomial = sp.Poly.from_list(list(values), variable, extension=True)
    if polynomial.is_zero:
        raise ValueError(f"{name} must not be the zero polynomial")
    return polynomial


def exact_rational_map_degree(
    numerator_coefficients: Sequence[Any],
    denominator_coefficients: Sequence[Any],
) -> RationalMapDegreeEvidence:
    """Reduce exact descending coefficients and return ``max(deg p, deg q)``.

    A common nonconstant factor changes the apparent polynomial degrees but not
    the rational map.  The returned evidence keeps both surfaces visible.
    """

    coordinate = sp.Symbol("z")
    numerator = _exact_polynomial(numerator_coefficients, coordinate, "numerator")
    denominator = _exact_polynomial(
        denominator_coefficients,
        coordinate,
        "denominator",
    )
    common = sp.gcd(numerator, denominator).monic()
    reduced_numerator = sp.exquo(numerator, common)
    reduced_denominator = sp.exquo(denominator, common)
    numerator_degree = int(numerator.degree())
    denominator_degree = int(denominator.degree())
    common_degree = int(common.degree())
    degree = max(int(reduced_numerator.degree()), int(reduced_denominator.degree()))
    return RationalMapDegreeEvidence(
        numerator=numerator.as_expr(),
        denominator=denominator.as_expr(),
        common_factor=common.as_expr(),
        reduced_numerator=reduced_numerator.as_expr(),
        reduced_denominator=reduced_denominator.as_expr(),
        numerator_degree=numerator_degree,
        denominator_degree=denominator_degree,
        common_factor_degree=common_degree,
        degree=degree,
    )


def rational_map_angular_lower_bound(degree: int) -> sp.Integer:
    r"""Return the exact Cauchy bound ``I[R] >= B**2``.

    The theorem additionally requires an orientation-preserving holomorphic
    degree-``B`` map and the normalized round-sphere measure.  This function
    returns the bound, not a verdict that a supplied map saturates it.
    """

    if isinstance(degree, bool) or int(degree) != degree or degree < 1:
        raise ValueError("degree must be a positive integer")
    return sp.Integer(int(degree) ** 2)


def axial_rational_map_angular_integral(degree: int) -> sp.Expr:
    r"""Return the exact angular functional for ``R(z)=z**B``.

    For integer ``B>=1``, beta integration gives

    ``I_B = B**3/3 * (1 + Gamma(2-1/B)*Gamma(2+1/B))``.

    The expression includes ``I_1=1`` and for ``B=2`` reduces to
    ``8/3 + pi``.  It evaluates one axial map and is not a global-minimizer
    theorem for arbitrary degree.
    """

    if isinstance(degree, bool) or int(degree) != degree or degree < 1:
        raise ValueError("degree must be a positive integer")
    value = sp.Integer(int(degree))
    inverse = sp.Rational(1, int(degree))
    result = value**3 * (
        1 + sp.gamma(2 - inverse) * sp.gamma(2 + inverse)
    ) / 3
    return sp.simplify(sp.expand_func(result))


def _numeric_polynomial(coefficients: Sequence[complex], name: str) -> ComplexArray:
    values = np.asarray(tuple(coefficients), dtype=np.complex128)
    if values.ndim != 1 or values.size == 0:
        raise ValueError(f"{name} must be a nonempty one-dimensional sequence")
    if not np.all(np.isfinite(values.real)) or not np.all(np.isfinite(values.imag)):
        raise ValueError(f"{name} must contain only finite coefficients")
    nonzero = np.flatnonzero(values != 0.0)
    if nonzero.size == 0:
        raise ValueError(f"{name} must not be the zero polynomial")
    return values[int(nonzero[0]) :]


def _quadrature_order(value: int, name: str) -> int:
    if isinstance(value, bool) or int(value) != value or value < 4:
        raise ValueError(f"{name} must be an integer at least four")
    return int(value)


def rational_map_sphere_integrals(
    numerator_coefficients: Sequence[complex],
    denominator_coefficients: Sequence[complex],
    *,
    declared_degree: int,
    polar_order: int = 32,
    azimuthal_order: int = 64,
) -> RationalMapSphereEvidence:
    r"""Evaluate normalized pullback area and ``I[R]`` on the round sphere.

    Coefficients are ordered from highest to lowest power.  The stable
    homogeneous formula

    ``J=(1+|z|**2)**2 * |p'*q-p*q'|**2 / (|p|**2+|q|**2)**2``

    avoids division by the rational-map denominator.  Tensor Gauss-Legendre
    nodes integrate in ``u=cos(theta)`` and ``phi`` without sampling either
    stereographic pole.  The caller must establish exact coprimality and degree
    separately; ``declared_degree`` is used only for the reported area error.
    """

    numerator = _numeric_polynomial(numerator_coefficients, "numerator")
    denominator = _numeric_polynomial(denominator_coefficients, "denominator")
    if (
        isinstance(declared_degree, bool)
        or int(declared_degree) != declared_degree
        or declared_degree < 0
    ):
        raise ValueError("declared_degree must be a nonnegative integer")
    n_polar = _quadrature_order(polar_order, "polar_order")
    n_azimuthal = _quadrature_order(azimuthal_order, "azimuthal_order")

    coefficient_scale = max(float(np.max(np.abs(numerator))), float(np.max(np.abs(denominator))))
    numerator = numerator / coefficient_scale
    denominator = denominator / coefficient_scale

    u, u_weights = np.polynomial.legendre.leggauss(n_polar)
    phi_base, phi_base_weights = np.polynomial.legendre.leggauss(n_azimuthal)
    phi = math.pi * (phi_base + 1.0)
    phi_weights = math.pi * phi_base_weights
    radius = np.sqrt((1.0 - u) / (1.0 + u))
    coordinate = radius[:, None] * np.exp(1j * phi[None, :])

    numerator_value = np.polyval(numerator, coordinate)
    denominator_value = np.polyval(denominator, coordinate)
    numerator_prime = np.polyval(np.polyder(numerator), coordinate)
    denominator_prime = np.polyval(np.polyder(denominator), coordinate)
    wronskian = numerator_prime * denominator_value - numerator_value * denominator_prime
    homogeneous_norm = np.abs(numerator_value) ** 2 + np.abs(denominator_value) ** 2
    minimum_norm = float(np.min(homogeneous_norm))
    if not math.isfinite(minimum_norm) or minimum_norm <= np.finfo(np.float64).tiny:
        raise ValueError("numerator and denominator are numerically simultaneous zeros")

    with np.errstate(over="raise", invalid="raise", divide="raise"):
        conformal_jacobian = (
            (1.0 + np.abs(coordinate) ** 2) ** 2
            * np.abs(wronskian) ** 2
            / homogeneous_norm**2
        )
    if not np.all(np.isfinite(conformal_jacobian)):
        raise ValueError("conformal Jacobian became non-finite")
    weights = u_weights[:, None] * phi_weights[None, :]
    normalized_area = float(np.sum(weights * conformal_jacobian) / (4.0 * math.pi))
    angular_integral = float(
        np.sum(weights * conformal_jacobian**2) / (4.0 * math.pi)
    )
    if normalized_area < 0.0 or angular_integral < 0.0:
        raise AssertionError("nonnegative sphere integrals became negative")
    return RationalMapSphereEvidence(
        numerator_coefficients=numerator.copy(),
        denominator_coefficients=denominator.copy(),
        declared_degree=int(declared_degree),
        polar_order=n_polar,
        azimuthal_order=n_azimuthal,
        normalized_area=normalized_area,
        angular_integral=angular_integral,
        maximum_conformal_jacobian=float(np.max(conformal_jacobian)),
        minimum_homogeneous_norm_squared=minimum_norm,
    )


def rotate_rational_map_about_axis(
    numerator_coefficients: Sequence[complex],
    denominator_coefficients: Sequence[complex],
    *,
    domain_angle: float = 0.0,
    target_angle: float = 0.0,
) -> tuple[ComplexArray, ComplexArray]:
    r"""Return coefficients for ``exp(i*b)*R(exp(i*a)*z)``.

    These are orientation-preserving rotations about the stereographic axis.
    Their sphere degree and angular integral must be unchanged.  General
    Möbius rotations are not inferred from this restricted helper.
    """

    numerator = _numeric_polynomial(numerator_coefficients, "numerator")
    denominator = _numeric_polynomial(denominator_coefficients, "denominator")
    alpha = float(domain_angle)
    beta = float(target_angle)
    if not math.isfinite(alpha) or not math.isfinite(beta):
        raise ValueError("rotation angles must be finite")

    def rotate_domain(coefficients: ComplexArray) -> ComplexArray:
        powers = np.arange(coefficients.size - 1, -1, -1)
        return coefficients * np.exp(1j * alpha * powers)

    return (
        np.exp(1j * beta) * rotate_domain(numerator),
        rotate_domain(denominator),
    )

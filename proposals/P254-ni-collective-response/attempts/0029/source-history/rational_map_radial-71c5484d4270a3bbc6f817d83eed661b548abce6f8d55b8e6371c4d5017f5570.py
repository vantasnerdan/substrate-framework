"""Conditional rational-map radial profiles and energy evidence.

The exact APIs in this module start from a separately declared reduced radial
functional with positive degree ``B`` and angular integral ``I``.  The numeric
API constructs one stationary shooting branch.  Neither surface establishes a
global rational-map minimum, a solution of the full three-dimensional field
theory, or an identification with a physical particle or nucleus.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import brentq
import sympy as sp

from .numerics import (
    IVPEvidence,
    SolverTolerances,
    solve_ivp_evidence,
    trapezoid_integral,
)

FloatArray = NDArray[np.float64]


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or int(value) != value or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _positive_float(value: Any, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return result


def _positive_symbolic(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.is_number and result.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return result


def _sample_vector(values: ArrayLike, name: str) -> FloatArray:
    result = np.asarray(values, dtype=np.float64)
    if result.ndim != 1 or result.size < 4:
        raise ValueError(
            f"{name} must be a one-dimensional array with at least four entries"
        )
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def rational_map_radial_energy_density(
    field: Any,
    radial_derivative: Any,
    radius: Any,
    degree: Any,
    angular_integral: Any,
) -> sp.Expr:
    r"""Return the declared dimensionless radial density for ``r>0``.

    The full energy is ``4*pi*integral L dr`` with

    ``L=r**2*f'**2 + 2*B*sin(f)**2*(1+f'**2) + I*sin(f)**4/r**2``.

    The declaration is conditional on the rational-map reduction and its
    normalization; this helper does not derive that physical action.
    """

    f = sp.sympify(field)
    fp = sp.sympify(radial_derivative)
    r = _positive_symbolic(radius, "radius")
    b = _positive_symbolic(degree, "degree")
    angular = _positive_symbolic(angular_integral, "angular_integral")
    sine_squared = sp.sin(f) ** 2
    return sp.simplify(
        r**2 * fp**2
        + 2 * b * sine_squared * (1 + fp**2)
        + angular * sine_squared**2 / r**2
    )


def rational_map_radial_euler_lagrange_residual(
    profile: Any,
    radius: sp.Symbol,
    degree: Any,
    angular_integral: Any,
) -> sp.Expr:
    r"""Return one-half of ``d/dr(dL/df')-dL/df`` exactly."""

    if not isinstance(radius, sp.Symbol):
        raise ValueError("radius must be a SymPy symbol")
    f = sp.sympify(profile)
    fp = sp.diff(f, radius)
    fpp = sp.diff(f, radius, 2)
    b = _positive_symbolic(degree, "degree")
    angular = _positive_symbolic(angular_integral, "angular_integral")
    return sp.simplify(
        (radius**2 + 2 * b * sp.sin(f) ** 2) * fpp
        + 2 * radius * fp
        + b * sp.sin(2 * f) * (fp**2 - 1)
        - angular * sp.sin(2 * f) * sp.sin(f) ** 2 / radius**2
    )


@dataclass(frozen=True)
class RadialEndpointExponents:
    """Exact regular-origin and massless-tail powers for declared ``B``."""

    degree: sp.Expr
    origin_power: sp.Expr
    tail_power: sp.Expr


def rational_map_radial_endpoint_exponents(degree: Any) -> RadialEndpointExponents:
    r"""Return ``sigma=(sqrt(1+8B)-1)/2`` and ``p=sigma+1``.

    The regular branch has ``f=pi-A*r**sigma+...`` and the decaying massless
    branch has ``f=C*r**(-p)+...``.  At ``B=1`` nonlinear terms enter the
    leading coefficient equation but do not change the powers ``sigma=1`` and
    ``p=2``.
    """

    b = _positive_symbolic(degree, "degree")
    origin = sp.simplify((sp.sqrt(1 + 8 * b) - 1) / 2)
    tail = sp.simplify((sp.sqrt(1 + 8 * b) + 1) / 2)
    return RadialEndpointExponents(
        degree=b,
        origin_power=origin,
        tail_power=tail,
    )


def regular_origin_boundary_residual(
    radius: Any,
    field: Any,
    radial_derivative: Any,
    origin_power: Any,
) -> Any:
    r"""Return ``r*f' + sigma*(pi-f)`` for regular-origin data."""

    return (
        radius * radial_derivative
        + origin_power * (np.pi - field)
    )


def massless_tail_boundary_residual(
    radius: Any,
    field: Any,
    radial_derivative: Any,
    tail_power: Any,
) -> Any:
    r"""Return ``r*f' + p*f`` for a decaying ``r**(-p)`` tail."""

    return radius * radial_derivative + tail_power * field


def _radial_rhs_unchecked(
    radius: float,
    state: ArrayLike,
    degree: float,
    angular_integral: float,
) -> FloatArray:
    values = np.asarray(state, dtype=np.float64)
    field, derivative = values
    sine = np.sin(field)
    sine_twice = np.sin(2.0 * field)
    coefficient = radius**2 + 2.0 * degree * sine**2
    second = (
        -2.0 * radius * derivative
        - degree * sine_twice * (derivative**2 - 1.0)
        + angular_integral * sine_twice * sine**2 / radius**2
    ) / coefficient
    return np.asarray((derivative, second), dtype=np.float64)


def rational_map_radial_rhs(
    radius: float,
    state: ArrayLike,
    degree: int,
    angular_integral: Any,
) -> FloatArray:
    """Return the first-order stationary radial system for ``r>0``."""

    r = _positive_float(radius, "radius")
    b = float(_positive_integer(degree, "degree"))
    angular = _positive_float(angular_integral, "angular_integral")
    values = np.asarray(state, dtype=np.float64)
    if values.shape != (2,) or not np.all(np.isfinite(values)):
        raise ValueError("state must contain finite field and derivative entries")
    return _radial_rhs_unchecked(r, values, b, angular)


def rational_map_radial_energy_components(
    radius: ArrayLike,
    field: ArrayLike,
    radial_derivative: ArrayLike,
    degree: int,
    angular_integral: Any,
) -> tuple[float, float]:
    r"""Return sampled ``E2`` and ``E4`` using the shared quadrature API."""

    r = _sample_vector(radius, "radius")
    f = _sample_vector(field, "field")
    fp = _sample_vector(radial_derivative, "radial_derivative")
    if not (r.shape == f.shape == fp.shape):
        raise ValueError("radius, field, and radial_derivative must have equal shapes")
    if np.any(r <= 0.0) or np.any(np.diff(r) <= 0.0):
        raise ValueError("radius must be positive and strictly increasing")
    b = float(_positive_integer(degree, "degree"))
    angular = _positive_float(angular_integral, "angular_integral")
    sine_squared = np.sin(f) ** 2
    two_density = r**2 * fp**2 + 2.0 * b * sine_squared
    four_density = (
        2.0 * b * sine_squared * fp**2
        + angular * sine_squared**2 / r**2
    )
    return (
        4.0 * np.pi * trapezoid_integral(two_density, r),
        4.0 * np.pi * trapezoid_integral(four_density, r),
    )


@dataclass(frozen=True)
class RationalMapRadialProfileEvidence:
    """Resolution-bounded shooting evidence for one declared radial branch."""

    degree: int
    angular_integral: float
    radius: FloatArray
    field: FloatArray
    radial_derivative: FloatArray
    origin_amplitude: float
    origin_power: float
    tail_power: float
    inner_boundary_residual: float
    outer_boundary_residual: float
    domain_two_derivative_energy: float
    domain_four_derivative_energy: float
    origin_two_derivative_estimate: float
    origin_four_derivative_estimate: float
    tail_two_derivative_estimate: float
    tail_four_derivative_estimate: float
    method: str
    function_evaluations: int

    @property
    def two_derivative_energy(self) -> float:
        """Return the sampled-domain value plus leading endpoint estimates."""

        return (
            self.domain_two_derivative_energy
            + self.origin_two_derivative_estimate
            + self.tail_two_derivative_estimate
        )

    @property
    def four_derivative_energy(self) -> float:
        """Return the sampled-domain value plus leading endpoint estimates."""

        return (
            self.domain_four_derivative_energy
            + self.origin_four_derivative_estimate
            + self.tail_four_derivative_estimate
        )

    @property
    def energy_coefficient(self) -> float:
        """Return the conventional conditional coefficient ``E/(12*pi**2)``."""

        return (
            self.two_derivative_energy + self.four_derivative_energy
        ) / (12.0 * np.pi**2)

    @property
    def per_degree_energy_coefficient(self) -> float:
        """Return the conditional energy coefficient divided by declared degree."""

        return self.energy_coefficient / self.degree

    @property
    def virial_relative_imbalance(self) -> float:
        """Return ``|E2-E4|/(E2+E4)`` for the corrected energy evidence."""

        two = self.two_derivative_energy
        four = self.four_derivative_energy
        return abs(two - four) / (two + four)


def solve_rational_map_radial_profile(
    degree: int,
    angular_integral: Any,
    *,
    inner_radius: float = 1.0e-4,
    outer_radius: float = 32.0,
    sample_points: int = 3201,
    amplitude_bracket: tuple[float, float] = (0.5, 4.0),
    tolerances: SolverTolerances = SolverTolerances(
        rtol=1.0e-10,
        atol=1.0e-12,
        max_step=0.02,
    ),
) -> RationalMapRadialProfileEvidence:
    r"""Shoot one branch between regular-origin and massless-tail conditions.

    Initial data use ``f(eps)=pi-A*eps**sigma`` and its derivative.  Brent's
    method selects ``A`` only by the declared outer residual
    ``R*f'(R)+p*f(R)`` inside the caller-visible amplitude bracket.  Leading
    asymptotic estimates account for the omitted ``[0,eps]`` and ``[R,infty)``
    energy intervals; they remain estimates rather than exact tail integrals.
    """

    b_int = _positive_integer(degree, "degree")
    b = float(b_int)
    angular = _positive_float(angular_integral, "angular_integral")
    inner = _positive_float(inner_radius, "inner_radius")
    outer = _positive_float(outer_radius, "outer_radius")
    if outer <= inner:
        raise ValueError("outer_radius must exceed inner_radius")
    if isinstance(sample_points, bool) or int(sample_points) != sample_points:
        raise ValueError("sample_points must be an integer")
    points = int(sample_points)
    if points < 101:
        raise ValueError("sample_points must be at least 101")
    if len(amplitude_bracket) != 2:
        raise ValueError("amplitude_bracket must contain two values")
    low = _positive_float(amplitude_bracket[0], "amplitude_bracket lower bound")
    high = _positive_float(amplitude_bracket[1], "amplitude_bracket upper bound")
    if high <= low:
        raise ValueError("amplitude_bracket must be strictly increasing")

    endpoint_exponents = rational_map_radial_endpoint_exponents(b_int)
    sigma = float(endpoint_exponents.origin_power)
    tail_power = float(endpoint_exponents.tail_power)
    endpoint_times = np.asarray((inner, outer), dtype=np.float64)

    def integrate(amplitude: float, sample_times: FloatArray) -> IVPEvidence:
        # Evolve g=pi-f rather than f near the origin.  For B>1 the
        # load-bearing perturbation A*eps**sigma can be smaller than the local
        # error allowed on f~=pi; the vacuum-complement equation has the same
        # first-order form and retains that small signal directly.
        initial = np.asarray(
            (
                amplitude * inner**sigma,
                amplitude * sigma * inner ** (sigma - 1.0),
            ),
            dtype=np.float64,
        )
        return solve_ivp_evidence(
            lambda radius, state: _radial_rhs_unchecked(
                radius,
                state,
                b,
                angular,
            ),
            (inner, outer),
            initial,
            sample_times=sample_times,
            tolerances=tolerances,
            method="DOP853",
        )

    def tail_residual(amplitude: float) -> float:
        endpoint = integrate(amplitude, endpoint_times)
        return float(
            massless_tail_boundary_residual(
                outer,
                np.pi - endpoint.state[0, -1],
                -endpoint.state[1, -1],
                tail_power,
            )
        )

    low_residual = tail_residual(low)
    high_residual = tail_residual(high)
    if low_residual * high_residual >= 0.0:
        raise ValueError(
            "amplitude_bracket does not bracket the asymptotic residual root"
        )
    amplitude = float(
        brentq(
            tail_residual,
            low,
            high,
            xtol=np.finfo(np.float64).eps,
            rtol=4.0 * np.finfo(np.float64).eps,
        )
    )
    radius = np.linspace(inner, outer, points, dtype=np.float64)
    solution = integrate(amplitude, radius)
    field = np.pi - solution.state[0]
    derivative = -solution.state[1]
    inner_residual = float(
        regular_origin_boundary_residual(
            inner,
            field[0],
            derivative[0],
            sigma,
        )
    )
    outer_residual = float(
        massless_tail_boundary_residual(
            outer,
            field[-1],
            derivative[-1],
            tail_power,
        )
    )
    domain_two, domain_four = rational_map_radial_energy_components(
        radius,
        field,
        derivative,
        b_int,
        angular,
    )

    origin_two = (
        4.0
        * np.pi
        * amplitude**2
        * (sigma**2 + 2.0 * b)
        * inner ** (2.0 * sigma + 1.0)
        / (2.0 * sigma + 1.0)
    )
    origin_four = (
        4.0
        * np.pi
        * amplitude**4
        * (2.0 * b * sigma**2 + angular)
        * inner ** (4.0 * sigma - 1.0)
        / (4.0 * sigma - 1.0)
    )
    tail_amplitude = float(field[-1] * outer**tail_power)
    tail_two = (
        4.0
        * np.pi
        * tail_amplitude**2
        * (tail_power**2 + 2.0 * b)
        * outer ** (1.0 - 2.0 * tail_power)
        / (2.0 * tail_power - 1.0)
    )
    tail_four = (
        4.0
        * np.pi
        * tail_amplitude**4
        * (2.0 * b * tail_power**2 + angular)
        * outer ** (-4.0 * tail_power - 1.0)
        / (4.0 * tail_power + 1.0)
    )
    return RationalMapRadialProfileEvidence(
        degree=b_int,
        angular_integral=angular,
        radius=radius,
        field=field,
        radial_derivative=derivative,
        origin_amplitude=amplitude,
        origin_power=sigma,
        tail_power=tail_power,
        inner_boundary_residual=inner_residual,
        outer_boundary_residual=outer_residual,
        domain_two_derivative_energy=domain_two,
        domain_four_derivative_energy=domain_four,
        origin_two_derivative_estimate=float(origin_two),
        origin_four_derivative_estimate=float(origin_four),
        tail_two_derivative_estimate=float(tail_two),
        tail_four_derivative_estimate=float(tail_four),
        method=(
            "DOP853 vacuum-complement amplitude shooting with asymptotic "
            "Robin tail"
        ),
        function_evaluations=solution.function_evaluations,
    )

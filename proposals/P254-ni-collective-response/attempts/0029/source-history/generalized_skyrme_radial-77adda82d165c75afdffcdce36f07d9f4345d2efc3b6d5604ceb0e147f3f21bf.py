"""Conditional rational-map profiles for a declared L2+L4+L6+L0 model.

The exact APIs start from independently supplied dimensionless coefficients
and rational-map data.  The numerical API constructs one stationary branch on
a truncated radial domain with asymptotic Robin data.  It does not establish a
half-line existence or uniqueness theorem, a global minimum, a full
three-dimensional solution, or a physical particle or nucleus.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.interpolate import make_interp_spline
from scipy.special import kve
import sympy as sp

from .numerics import BVPEvidence, solve_bvp_evidence, trapezoid_integral
from .rational_map_radial import rational_map_radial_endpoint_exponents


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


def _nonnegative_float(value: Any, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result < 0.0:
        raise ValueError(f"{name} must be nonnegative and finite")
    return result


def _positive_symbolic(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.is_number and result.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return result


def _nonnegative_symbolic(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.is_number and result.is_nonnegative is not True:
        raise ValueError(f"{name} must be nonnegative")
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


def generalized_skyrme_radial_energy_density(
    field: Any,
    radial_derivative: Any,
    radius: Any,
    degree: Any,
    angular_integral: Any,
    sextic_coefficient: Any,
    potential_coefficient: Any,
) -> sp.Expr:
    r"""Return the declared dimensionless radial density for ``r>0``.

    The full conditional energy is ``4*pi*integral L dr``.  The extra terms
    are ``c6*I*sin(f)**4*f'**2/r**2`` and
    ``c0*r**2*(1-cos(f))``.
    """

    f = sp.sympify(field)
    fp = sp.sympify(radial_derivative)
    r = _positive_symbolic(radius, "radius")
    b = _positive_symbolic(degree, "degree")
    angular = _positive_symbolic(angular_integral, "angular_integral")
    c6 = _nonnegative_symbolic(sextic_coefficient, "sextic_coefficient")
    c0 = _nonnegative_symbolic(potential_coefficient, "potential_coefficient")
    sine_squared = sp.sin(f) ** 2
    return sp.simplify(
        r**2 * fp**2
        + 2 * b * sine_squared * (1 + fp**2)
        + angular * sine_squared**2 / r**2
        + c6 * angular * sine_squared**2 * fp**2 / r**2
        + c0 * r**2 * (1 - sp.cos(f))
    )


def generalized_skyrme_radial_euler_lagrange_residual(
    profile: Any,
    radius: sp.Symbol,
    degree: Any,
    angular_integral: Any,
    sextic_coefficient: Any,
    potential_coefficient: Any,
) -> sp.Expr:
    r"""Return one-half of ``d/dr(dL/df')-dL/df`` exactly."""

    if not isinstance(radius, sp.Symbol):
        raise ValueError("radius must be a SymPy symbol")
    f = sp.sympify(profile)
    fp = sp.diff(f, radius)
    fpp = sp.diff(f, radius, 2)
    b = _positive_symbolic(degree, "degree")
    angular = _positive_symbolic(angular_integral, "angular_integral")
    c6 = _nonnegative_symbolic(sextic_coefficient, "sextic_coefficient")
    c0 = _nonnegative_symbolic(potential_coefficient, "potential_coefficient")
    sine = sp.sin(f)
    sine_squared = sine**2
    sine_twice = sp.sin(2 * f)
    kinetic = radius**2 + 2 * b * sine_squared + c6 * angular * sine**4 / radius**2
    return sp.simplify(
        kinetic * fpp
        + (2 * radius - 2 * c6 * angular * sine**4 / radius**3) * fp
        + b * sine_twice * (fp**2 - 1)
        + c6 * angular * sine_squared * sine_twice * fp**2 / radius**2
        - angular * sine_squared * sine_twice / radius**2
        - c0 * radius**2 * sine / 2
    )


def generalized_skyrme_reduced_coefficients(
    lambda_bps: Any,
    potential_scale: Any,
    skyrme_coupling: Any,
    decay_scale: Any,
) -> tuple[sp.Expr, sp.Expr]:
    r"""Return ``(c6,c0)`` in the accepted BPS lambda convention.

    With ``lambda_A=pi**2*lambda_BPS``, the same sextic coefficient can be
    written ``lambda_A**2*e**4*F**2/(8*pi**4)``.  All arguments remain
    supplied; this function derives no physical identification or value.
    """

    lam = _positive_symbolic(lambda_bps, "lambda_bps")
    mu = _positive_symbolic(potential_scale, "potential_scale")
    coupling = _positive_symbolic(skyrme_coupling, "skyrme_coupling")
    scale = _positive_symbolic(decay_scale, "decay_scale")
    return (
        sp.simplify(lam**2 * coupling**4 * scale**2 / 8),
        sp.simplify(32 * mu**2 / (coupling**2 * scale**4)),
    )


def generalized_skyrme_scaling_residual(
    two_derivative_energy: Any,
    four_derivative_energy: Any,
    sextic_energy: Any,
    potential_energy: Any,
) -> sp.Expr:
    r"""Return the stationary Derrick residual ``E2-E4-3E6+3E0``."""

    return sp.simplify(
        sp.sympify(two_derivative_energy)
        - sp.sympify(four_derivative_energy)
        - 3 * sp.sympify(sextic_energy)
        + 3 * sp.sympify(potential_energy)
    )


@dataclass(frozen=True)
class GeneralizedSkyrmeEndpointData:
    """Regular-origin and linear-tail data for the declared radial model."""

    origin_power: float
    tail_power: float
    tail_mass: float
    bessel_order: float


def generalized_skyrme_endpoint_data(
    degree: int,
    potential_coefficient: Any,
) -> GeneralizedSkyrmeEndpointData:
    """Return endpoint powers and the potential-induced tail mass."""

    b = _positive_integer(degree, "degree")
    c0 = _nonnegative_float(potential_coefficient, "potential_coefficient")
    powers = rational_map_radial_endpoint_exponents(b)
    return GeneralizedSkyrmeEndpointData(
        origin_power=float(powers.origin_power),
        tail_power=float(powers.tail_power),
        tail_mass=float(np.sqrt(c0 / 2.0)),
        bessel_order=float(np.sqrt(0.25 + 2.0 * b)),
    )


def generalized_skyrme_tail_robin_coefficient(
    radius: Any,
    degree: int,
    potential_coefficient: Any,
) -> float:
    r"""Return ``q`` for the linear asymptotic condition ``f'+q*f=0``."""

    r = _positive_float(radius, "radius")
    data = generalized_skyrme_endpoint_data(degree, potential_coefficient)
    if data.tail_mass == 0.0:
        return data.tail_power / r
    z = data.tail_mass * r
    ratio = (
        kve(data.bessel_order - 1.0, z)
        + kve(data.bessel_order + 1.0, z)
    ) / kve(data.bessel_order, z)
    return float(1.0 / (2.0 * r) + 0.5 * data.tail_mass * ratio)


def generalized_skyrme_radial_rhs(
    radius: ArrayLike,
    state: ArrayLike,
    degree: int,
    angular_integral: Any,
    sextic_coefficient: Any,
    potential_coefficient: Any,
) -> FloatArray:
    """Return the first-order stationary radial system for positive radii."""

    r = np.asarray(radius, dtype=np.float64)
    if np.any(~np.isfinite(r)) or np.any(r <= 0.0):
        raise ValueError("radius must contain positive finite values")
    y = np.asarray(state, dtype=np.float64)
    if y.shape[0] != 2 or y.shape[-1] != r.size or np.any(~np.isfinite(y)):
        raise ValueError("state must have shape (2, radius size) and be finite")
    b = float(_positive_integer(degree, "degree"))
    angular = _positive_float(angular_integral, "angular_integral")
    c6 = _nonnegative_float(sextic_coefficient, "sextic_coefficient")
    c0 = _nonnegative_float(potential_coefficient, "potential_coefficient")
    field, derivative = y
    sine = np.sin(field)
    sine_squared = sine**2
    sine_twice = np.sin(2.0 * field)
    kinetic = r**2 + 2.0 * b * sine_squared + c6 * angular * sine**4 / r**2
    numerator = (
        -2.0 * r * derivative
        - b * sine_twice * (derivative**2 - 1.0)
        + angular * sine_squared * sine_twice / r**2
        - c6 * angular * sine_squared * sine_twice * derivative**2 / r**2
        + 2.0 * c6 * angular * sine**4 * derivative / r**3
        + 0.5 * c0 * r**2 * sine
    )
    return np.vstack((derivative, numerator / kinetic))


def generalized_skyrme_energy_components(
    radius: ArrayLike,
    field: ArrayLike,
    radial_derivative: ArrayLike,
    degree: int,
    angular_integral: Any,
    sextic_coefficient: Any,
    potential_coefficient: Any,
) -> tuple[float, float, float, float]:
    """Return sampled ``(E2,E4,E6,E0)`` using shared quadrature."""

    r = _sample_vector(radius, "radius")
    f = _sample_vector(field, "field")
    fp = _sample_vector(radial_derivative, "radial_derivative")
    if not (r.shape == f.shape == fp.shape):
        raise ValueError("radius, field, and radial_derivative must have equal shapes")
    if np.any(r <= 0.0) or np.any(np.diff(r) <= 0.0):
        raise ValueError("radius must be positive and strictly increasing")
    b = float(_positive_integer(degree, "degree"))
    angular = _positive_float(angular_integral, "angular_integral")
    c6 = _nonnegative_float(sextic_coefficient, "sextic_coefficient")
    c0 = _nonnegative_float(potential_coefficient, "potential_coefficient")
    sine_squared = np.sin(f) ** 2
    densities = (
        r**2 * fp**2 + 2.0 * b * sine_squared,
        2.0 * b * sine_squared * fp**2 + angular * sine_squared**2 / r**2,
        c6 * angular * sine_squared**2 * fp**2 / r**2,
        c0 * r**2 * (1.0 - np.cos(f)),
    )
    return tuple(
        4.0 * np.pi * trapezoid_integral(density, r) for density in densities
    )  # type: ignore[return-value]


@dataclass(frozen=True)
class GeneralizedSkyrmeRadialProfileEvidence:
    """Resolution-bounded evidence for one declared stationary branch."""

    degree: int
    angular_integral: float
    sextic_coefficient: float
    potential_coefficient: float
    radius: FloatArray
    field: FloatArray
    radial_derivative: FloatArray
    two_derivative_energy: float
    four_derivative_energy: float
    sextic_energy: float
    potential_energy: float
    origin_power: float
    tail_robin_coefficient: float
    inner_boundary_residual: float
    outer_boundary_residual: float
    solver_nodes: int
    solver_iterations: int
    max_rms_residual: float
    continuation_steps: int
    method: str

    @property
    def energy_coefficient(self) -> float:
        """Return the finite-domain coefficient ``E/(12*pi**2)``."""

        total = (
            self.two_derivative_energy
            + self.four_derivative_energy
            + self.sextic_energy
            + self.potential_energy
        )
        return total / (12.0 * np.pi**2)

    @property
    def virial_relative_residual(self) -> float:
        """Return the scale-stationarity residual relative to total energy."""

        residual = (
            self.two_derivative_energy
            - self.four_derivative_energy
            - 3.0 * self.sextic_energy
            + 3.0 * self.potential_energy
        )
        total = (
            self.two_derivative_energy
            + self.four_derivative_energy
            + self.sextic_energy
            + self.potential_energy
        )
        return abs(residual) / total


def solve_generalized_skyrme_radial_profile(
    degree: int,
    angular_integral: Any,
    sextic_coefficient: Any,
    potential_coefficient: Any,
    *,
    inner_radius: float = 1.0e-4,
    outer_radius: float = 20.0,
    initial_points: int = 501,
    sample_points: int = 4001,
    continuation_steps: int = 10,
    tolerance: float = 1.0e-6,
    max_nodes: int = 100_000,
) -> GeneralizedSkyrmeRadialProfileEvidence:
    r"""Construct one continuation BVP branch with asymptotic Robin data.

    The continuation starts at the C-RPROF-001 ``c6=c0=0`` equation and ramps
    both supplied coefficients linearly.  Every stage passes through
    :func:`solve_bvp_evidence`, which rejects unsuccessful or nonfinite solves.
    Energy is sampled independently of the adaptive solver mesh.
    """

    b = _positive_integer(degree, "degree")
    angular = _positive_float(angular_integral, "angular_integral")
    c6 = _nonnegative_float(sextic_coefficient, "sextic_coefficient")
    c0 = _nonnegative_float(potential_coefficient, "potential_coefficient")
    inner = _positive_float(inner_radius, "inner_radius")
    outer = _positive_float(outer_radius, "outer_radius")
    if outer <= inner:
        raise ValueError("outer_radius must exceed inner_radius")
    initial_n = _positive_integer(initial_points, "initial_points")
    sample_n = _positive_integer(sample_points, "sample_points")
    steps = _positive_integer(continuation_steps, "continuation_steps")
    if initial_n < 101:
        raise ValueError("initial_points must be at least 101")
    if sample_n < 101:
        raise ValueError("sample_points must be at least 101")
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be positive and finite")
    if max_nodes < initial_n:
        raise ValueError("max_nodes cannot be smaller than initial_points")

    endpoint = generalized_skyrme_endpoint_data(b, c0)
    coordinate = np.linspace(inner, outer, initial_n, dtype=np.float64)
    width = 1.5 + 0.8 * b
    field_guess = np.pi * np.exp(-(coordinate / width) ** 1.2)
    state = np.vstack((field_guess, np.gradient(field_guess, coordinate)))
    final: BVPEvidence | None = None

    for fraction in np.linspace(0.0, 1.0, steps + 1):
        stage_c6 = c6 * float(fraction)
        stage_c0 = c0 * float(fraction)
        stage_tail = generalized_skyrme_tail_robin_coefficient(
            outer,
            b,
            stage_c0,
        )

        def equations(radius: FloatArray, values: FloatArray) -> FloatArray:
            return generalized_skyrme_radial_rhs(
                radius,
                values,
                b,
                angular,
                stage_c6,
                stage_c0,
            )

        def boundary(left: FloatArray, right: FloatArray) -> FloatArray:
            return np.asarray(
                (
                    inner * left[1]
                    + endpoint.origin_power * (np.pi - left[0]),
                    right[1] + stage_tail * right[0],
                ),
                dtype=np.float64,
            )

        final = solve_bvp_evidence(
            equations,
            boundary,
            coordinate,
            state,
            tolerance=tolerance,
            max_nodes=max_nodes,
        )
        coordinate = final.coordinate
        state = final.state

    if final is None:  # pragma: no cover - positive steps make this unreachable
        raise RuntimeError("continuation produced no solve")
    radius = np.linspace(inner, outer, sample_n, dtype=np.float64)
    sampled = make_interp_spline(final.coordinate, final.state, axis=1)(radius)
    field = np.asarray(sampled[0], dtype=np.float64)
    derivative = np.asarray(sampled[1], dtype=np.float64)
    components = generalized_skyrme_energy_components(
        radius,
        field,
        derivative,
        b,
        angular,
        c6,
        c0,
    )
    tail = generalized_skyrme_tail_robin_coefficient(outer, b, c0)
    inner_residual = float(
        inner * final.state[1, 0]
        + endpoint.origin_power * (np.pi - final.state[0, 0])
    )
    outer_residual = float(final.state[1, -1] + tail * final.state[0, -1])
    return GeneralizedSkyrmeRadialProfileEvidence(
        degree=b,
        angular_integral=angular,
        sextic_coefficient=c6,
        potential_coefficient=c0,
        radius=radius,
        field=field,
        radial_derivative=derivative,
        two_derivative_energy=components[0],
        four_derivative_energy=components[1],
        sextic_energy=components[2],
        potential_energy=components[3],
        origin_power=endpoint.origin_power,
        tail_robin_coefficient=tail,
        inner_boundary_residual=inner_residual,
        outer_boundary_residual=outer_residual,
        solver_nodes=final.coordinate.size,
        solver_iterations=final.iterations,
        max_rms_residual=final.max_rms_residual,
        continuation_steps=steps,
        method="adaptive collocation continuation with asymptotic Robin data",
    )


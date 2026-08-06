"""Exact conditional reduction of a dimensional cosine field model.

The declared density is

``lambda*u_t**2/2 - tension*u_x**2/2 - mu*(1-cos(u))``

for a dimensionless real field and physical space and time coordinates. This
module derives its coefficient ratios and the pullback of the accepted
normalized sine-Gordon breather. It does not derive the density from a
material, select any coefficient, or remove the common coefficient scale.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

import sympy as sp

from .exact_symbolic import exact_real as _real_exact
from .exact_symbolic import positive_exact as _positive_exact
from .sine_gordon import (
    breather_action,
    breather_energy,
    breather_field,
    breather_inverse_width,
)


def _nonnegative_exact(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True or expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be provably nonnegative and real")
    return expression


def _normalized_frequency(value: Any) -> sp.Expr:
    frequency = sp.sympify(value)
    if frequency.has(sp.Float):
        raise ValueError("frequency must be exact rather than floating")
    if frequency.is_number and (
        frequency.is_real is not True
        or not bool(frequency > 0)
        or not bool(frequency < 1)
    ):
        raise ValueError("frequency must satisfy 0 < frequency < 1")
    if frequency.is_real is False:
        raise ValueError("frequency must be real")
    return frequency


@dataclass(frozen=True)
class DimensionalSineGordonCoefficients:
    """Positive coefficients of the declared dimensional density."""

    inertia: sp.Expr
    gradient: sp.Expr
    onsite: sp.Expr


@dataclass(frozen=True)
class DimensionalSineGordonScales:
    """Ratios and dimensional normalization scales of the declared model."""

    signal_speed: sp.Expr
    gap_frequency: sp.Expr
    length: sp.Expr
    energy: sp.Expr
    action: sp.Expr


@dataclass(frozen=True)
class DimensionalBreatherObservables:
    """Exact physical-coordinate observables of a pulled-back breather."""

    angular_frequency: sp.Expr
    period: sp.Expr
    inverse_width: sp.Expr
    profile_length: sp.Expr
    energy: sp.Expr
    action: sp.Expr


@dataclass(frozen=True)
class DimensionalSineGordonStress:
    """Physical Noether stress in orthonormal coordinates ``(c*t, x)``."""

    energy_density: sp.Expr
    energy_flux: sp.Expr
    longitudinal_pressure: sp.Expr
    contravariant: sp.ImmutableMatrix
    divergence: sp.ImmutableMatrix
    field_equation_residual: sp.Expr


TailBehavior = Literal["evanescent", "threshold", "oscillatory"]


@dataclass(frozen=True)
class DimensionalSineGordonLinearSpectrum:
    """Exact positive-frequency branch of the linearized physical spectrum."""

    wavenumber: sp.Expr
    angular_frequency_squared: sp.Expr
    angular_frequency: sp.Expr
    group_velocity: sp.Expr
    phase_velocity: sp.Expr | None


@dataclass(frozen=True)
class DimensionalSineGordonTimeHarmonicLedger:
    """Exterior and whole-line classification at one real angular frequency."""

    angular_frequency: sp.Expr
    spatial_coefficient: sp.Expr
    behavior: TailBehavior
    spatial_rate: sp.Expr
    right_half_line_l2_dimension: int
    left_half_line_l2_dimension: int
    whole_line_l2_dimension: int


def _validated_coefficients(
    coefficients: DimensionalSineGordonCoefficients,
) -> DimensionalSineGordonCoefficients:
    """Validate a coefficient record even when its dataclass was built directly."""

    if not isinstance(coefficients, DimensionalSineGordonCoefficients):
        raise TypeError("coefficients must be a DimensionalSineGordonCoefficients record")
    return dimensional_sine_gordon_coefficients(
        coefficients.inertia,
        coefficients.gradient,
        coefficients.onsite,
    )


def dimensional_sine_gordon_coefficients(
    inertia: Any,
    gradient: Any,
    onsite: Any,
) -> DimensionalSineGordonCoefficients:
    """Return the exact positive coefficients in canonical field order."""

    return DimensionalSineGordonCoefficients(
        inertia=_positive_exact(inertia, "inertia"),
        gradient=_positive_exact(gradient, "gradient"),
        onsite=_positive_exact(onsite, "onsite"),
    )


def dimensional_sine_gordon_scales(
    coefficients: DimensionalSineGordonCoefficients,
) -> DimensionalSineGordonScales:
    r"""Return ``c``, ``omega_0``, ``ell``, and energy/action scales.

    For coefficients ``(lambda,T,mu)``,
    ``c=sqrt(T/lambda)``, ``omega_0=sqrt(mu/lambda)``,
    ``ell=sqrt(T/mu)``, ``E_scale=sqrt(T*mu)``, and
    ``J_scale=sqrt(lambda*T)``.
    """

    lam = _positive_exact(coefficients.inertia, "inertia")
    tension = _positive_exact(coefficients.gradient, "gradient")
    mu = _positive_exact(coefficients.onsite, "onsite")
    return DimensionalSineGordonScales(
        signal_speed=sp.sqrt(tension / lam),
        gap_frequency=sp.sqrt(mu / lam),
        length=sp.sqrt(tension / mu),
        energy=sp.sqrt(tension * mu),
        action=sp.sqrt(lam * tension),
    )


def dimensional_sine_gordon_coefficients_from_speed_gap(
    inertia_scale: Any,
    signal_speed: Any,
    gap_frequency: Any,
) -> DimensionalSineGordonCoefficients:
    r"""Return the coefficient ray at supplied ``lambda``, ``c``, ``omega_0``.

    The inverse family is
    ``(lambda,T,mu)=(lambda,lambda*c**2,lambda*omega_0**2)``.
    Thus ``c`` and ``omega_0`` do not determine the common positive scale.
    """

    lam = _positive_exact(inertia_scale, "inertia_scale")
    speed = _positive_exact(signal_speed, "signal_speed")
    gap = _positive_exact(gap_frequency, "gap_frequency")
    return dimensional_sine_gordon_coefficients(
        lam,
        lam * speed**2,
        lam * gap**2,
    )


def rescale_dimensional_sine_gordon_coefficients(
    coefficients: DimensionalSineGordonCoefficients,
    multiplier: Any,
) -> DimensionalSineGordonCoefficients:
    """Apply a common positive multiplier to all three coefficients."""

    factor = _positive_exact(multiplier, "multiplier")
    return dimensional_sine_gordon_coefficients(
        factor * coefficients.inertia,
        factor * coefficients.gradient,
        factor * coefficients.onsite,
    )


def dimensional_sine_gordon_lagrangian_density(
    field: Any,
    coordinate: sp.Symbol,
    time: sp.Symbol,
    coefficients: DimensionalSineGordonCoefficients,
) -> sp.Expr:
    """Return the declared dimensional cosine Lagrangian density."""

    expression = sp.sympify(field)
    validated = _validated_coefficients(coefficients)
    return (
        validated.inertia * sp.diff(expression, time) ** 2 / 2
        - validated.gradient * sp.diff(expression, coordinate) ** 2 / 2
        - validated.onsite * (1 - sp.cos(expression))
    )


def dimensional_sine_gordon_hamiltonian_density(
    field: Any,
    coordinate: sp.Symbol,
    time: sp.Symbol,
    coefficients: DimensionalSineGordonCoefficients,
) -> sp.Expr:
    """Return the positive Noether energy density for constant coefficients."""

    expression = sp.sympify(field)
    validated = _validated_coefficients(coefficients)
    return (
        validated.inertia * sp.diff(expression, time) ** 2 / 2
        + validated.gradient * sp.diff(expression, coordinate) ** 2 / 2
        + validated.onsite * (1 - sp.cos(expression))
    )


def dimensional_sine_gordon_stress(
    field: Any,
    coordinate: sp.Symbol,
    time: sp.Symbol,
    coefficients: DimensionalSineGordonCoefficients,
) -> DimensionalSineGordonStress:
    r"""Return the exact physical stress and its off-shell divergence.

    Indices refer to the orthonormal coordinates ``(x**0,x**1)=(c*t,x)``.
    With ``R=lambda*u_tt-T*u_xx+mu*sin(u)``, the exact identity is

    ``partial_mu T**(mu nu) = (u_t*R/c, -u_x*R)``.

    Consequently every solution of the declared dimensional sine-Gordon
    equation supplies a symmetric conserved 1+1 stress without introducing a
    second action or source-normalization scale.
    """

    expression = sp.sympify(field)
    validated = _validated_coefficients(coefficients)
    speed = dimensional_sine_gordon_scales(validated).signal_speed
    time_derivative = sp.diff(expression, time)
    space_derivative = sp.diff(expression, coordinate)
    energy_density = dimensional_sine_gordon_hamiltonian_density(
        expression,
        coordinate,
        time,
        validated,
    )
    energy_flux = sp.simplify(
        -validated.gradient * time_derivative * space_derivative
    )
    longitudinal_pressure = sp.simplify(
        validated.inertia * time_derivative**2 / 2
        + validated.gradient * space_derivative**2 / 2
        - validated.onsite * (1 - sp.cos(expression))
    )
    contravariant = sp.ImmutableMatrix(
        [
            [energy_density, energy_flux / speed],
            [energy_flux / speed, longitudinal_pressure],
        ]
    ).applyfunc(sp.simplify)
    residual = dimensional_sine_gordon_residual(
        expression,
        coordinate,
        time,
        validated,
    )
    divergence = sp.ImmutableMatrix(
        [
            time_derivative * residual / speed,
            -space_derivative * residual,
        ]
    )
    direct_divergence = sp.ImmutableMatrix(
        [
            sp.diff(contravariant[0, 0], time) / speed
            + sp.diff(contravariant[1, 0], coordinate),
            sp.diff(contravariant[0, 1], time) / speed
            + sp.diff(contravariant[1, 1], coordinate),
        ]
    ).applyfunc(sp.simplify)
    if (direct_divergence - divergence).applyfunc(sp.simplify) != sp.zeros(2, 1):
        raise AssertionError("dimensional sine-Gordon stress divergence failed")
    return DimensionalSineGordonStress(
        energy_density=sp.simplify(energy_density),
        energy_flux=energy_flux,
        longitudinal_pressure=longitudinal_pressure,
        contravariant=contravariant,
        divergence=divergence,
        field_equation_residual=sp.simplify(residual),
    )


def dimensional_sine_gordon_residual(
    field: Any,
    coordinate: sp.Symbol,
    time: sp.Symbol,
    coefficients: DimensionalSineGordonCoefficients,
) -> sp.Expr:
    """Return ``lambda*u_tt-T*u_xx+mu*sin(u)``."""

    expression = sp.sympify(field)
    validated = _validated_coefficients(coefficients)
    return (
        validated.inertia * sp.diff(expression, time, 2)
        - validated.gradient * sp.diff(expression, coordinate, 2)
        + validated.onsite * sp.sin(expression)
    )


def dimensional_sine_gordon_linearized_residual(
    field: Any,
    coordinate: sp.Symbol,
    time: sp.Symbol,
    coefficients: DimensionalSineGordonCoefficients,
) -> sp.Expr:
    """Return the vacuum linearization ``lambda*u_tt-T*u_xx+mu*u``."""

    expression = sp.sympify(field)
    validated = _validated_coefficients(coefficients)
    return (
        validated.inertia * sp.diff(expression, time, 2)
        - validated.gradient * sp.diff(expression, coordinate, 2)
        + validated.onsite * expression
    )


def dimensional_sine_gordon_linear_spectrum(
    wavenumber: Any,
    coefficients: DimensionalSineGordonCoefficients,
) -> DimensionalSineGordonLinearSpectrum:
    r"""Return the positive branch of ``Omega^2=omega_0^2+c^2*k^2``.

    The signed phase velocity ``Omega/k`` is undefined at concrete ``k=0``;
    the record uses ``None`` there. For a symbolic wave number, callers are
    responsible for declaring ``k != 0`` before using the phase velocity.
    """

    k = _real_exact(wavenumber, "wavenumber")
    scales = dimensional_sine_gordon_scales(coefficients)
    squared = sp.simplify(
        scales.gap_frequency**2 + scales.signal_speed**2 * k**2
    )
    angular = sp.sqrt(squared)
    group = sp.simplify(scales.signal_speed**2 * k / angular)
    phase = None if k.is_zero is True else sp.simplify(angular / k)
    return DimensionalSineGordonLinearSpectrum(
        wavenumber=k,
        angular_frequency_squared=squared,
        angular_frequency=angular,
        group_velocity=group,
        phase_velocity=phase,
    )


def dimensional_sine_gordon_tail_coefficient(
    angular_frequency: Any,
    coefficients: DimensionalSineGordonCoefficients,
) -> sp.Expr:
    r"""Return ``(omega_0^2-Omega^2)/c^2`` in ``a''=coefficient*a``."""

    angular = _nonnegative_exact(angular_frequency, "angular_frequency")
    scales = dimensional_sine_gordon_scales(coefficients)
    return sp.simplify(
        (scales.gap_frequency**2 - angular**2) / scales.signal_speed**2
    )


def evanescent_half_line_matching_matrix(rate: Any) -> sp.ImmutableMatrix:
    r"""Return the origin matching matrix for two decaying half-line branches.

    The right branch is ``A*exp(-kappa*x)`` and the left branch is
    ``B*exp(kappa*x)``. Continuity of the field and first derivative at zero
    gives this matrix on ``(A,B)``. Its determinant is ``-2*kappa``.
    """

    kappa = _positive_exact(rate, "rate")
    return sp.ImmutableMatrix([[1, -1], [-kappa, -kappa]])


def dimensional_sine_gordon_time_harmonic_ledger(
    angular_frequency: Any,
    coefficients: DimensionalSineGordonCoefficients,
) -> DimensionalSineGordonTimeHarmonicLedger:
    """Classify exterior branches and global L2 modes on the homogeneous line.

    Sub-gap equations have one decaying solution on each exterior half-line,
    but their smooth whole-line matching matrix is full rank. Threshold and
    oscillatory constant-coefficient solutions have no nonzero half-line L2
    branch. Thus the global homogeneous linear equation has no nonzero L2
    real-frequency separated mode in any branch. A nonlinear or defect core is
    additional structure and is not inferred by this ledger.
    """

    angular = _nonnegative_exact(angular_frequency, "angular_frequency")
    coefficient = dimensional_sine_gordon_tail_coefficient(
        angular,
        coefficients,
    )
    if coefficient.is_positive is True:
        behavior: TailBehavior = "evanescent"
        rate = sp.sqrt(coefficient)
        matching = evanescent_half_line_matching_matrix(rate)
        whole_line_dimension = 2 - matching.rank()
        right_dimension = left_dimension = 1
    elif coefficient.is_zero is True:
        behavior = "threshold"
        rate = sp.Integer(0)
        whole_line_dimension = 0
        right_dimension = left_dimension = 0
    elif coefficient.is_negative is True:
        behavior = "oscillatory"
        rate = sp.sqrt(-coefficient)
        whole_line_dimension = 0
        right_dimension = left_dimension = 0
    else:
        raise ValueError(
            "angular_frequency must have a decidable relation to the gap frequency"
        )
    return DimensionalSineGordonTimeHarmonicLedger(
        angular_frequency=angular,
        spatial_coefficient=coefficient,
        behavior=behavior,
        spatial_rate=rate,
        right_half_line_l2_dimension=right_dimension,
        left_half_line_l2_dimension=left_dimension,
        whole_line_l2_dimension=whole_line_dimension,
    )


def linear_wave_traveling_field(
    profile: Any,
    coordinate: Any,
    time: Any,
    signal_speed: Any,
    *,
    direction: int = 1,
) -> sp.Expr:
    """Return ``F(x-direction*c*t)`` for a callable symbolic profile ``F``."""

    if direction not in {-1, 1}:
        raise ValueError("direction must be -1 or 1")
    speed = _positive_exact(signal_speed, "signal_speed")
    if not callable(profile):
        raise TypeError("profile must be callable")
    argument = sp.sympify(coordinate) - direction * speed * sp.sympify(time)
    return sp.sympify(profile(argument))


def linear_wave_residual(
    field: Any,
    coordinate: sp.Symbol,
    time: sp.Symbol,
    signal_speed: Any,
) -> sp.Expr:
    """Return the gapless wave residual ``u_tt-c^2*u_xx``."""

    speed = _positive_exact(signal_speed, "signal_speed")
    expression = sp.sympify(field)
    return sp.diff(expression, time, 2) - speed**2 * sp.diff(
        expression,
        coordinate,
        2,
    )


def linear_wave_energy_density(
    field: Any,
    coordinate: sp.Symbol,
    time: sp.Symbol,
    inertia: Any,
    gradient: Any,
) -> sp.Expr:
    """Return ``lambda*u_t^2/2+T*u_x^2/2`` for a gapless linear medium."""

    lam = _positive_exact(inertia, "inertia")
    tension = _positive_exact(gradient, "gradient")
    expression = sp.sympify(field)
    return (
        lam * sp.diff(expression, time) ** 2 / 2
        + tension * sp.diff(expression, coordinate) ** 2 / 2
    )


def dimensional_sine_gordon_normalized_coordinates(
    coordinate: Any,
    time: Any,
    coefficients: DimensionalSineGordonCoefficients,
) -> tuple[sp.Expr, sp.Expr]:
    """Return ``(X,tau)=(x/ell,omega_0*t)``."""

    scales = dimensional_sine_gordon_scales(coefficients)
    return (
        sp.sympify(coordinate) / scales.length,
        scales.gap_frequency * sp.sympify(time),
    )


def dimensional_sine_gordon_physical_coordinates(
    normalized_coordinate: Any,
    normalized_time: Any,
    coefficients: DimensionalSineGordonCoefficients,
) -> tuple[sp.Expr, sp.Expr]:
    """Return ``(x,t)=(ell*X,tau/omega_0)``."""

    scales = dimensional_sine_gordon_scales(coefficients)
    return (
        scales.length * sp.sympify(normalized_coordinate),
        sp.sympify(normalized_time) / scales.gap_frequency,
    )


def dimensional_breather_field(
    coordinate: Any,
    time: Any,
    frequency: Any,
    coefficients: DimensionalSineGordonCoefficients,
) -> sp.Expr:
    """Pull the accepted normalized rest breather into physical coordinates."""

    normalized_frequency = _normalized_frequency(frequency)
    normalized_x, normalized_t = dimensional_sine_gordon_normalized_coordinates(
        coordinate,
        time,
        coefficients,
    )
    return breather_field(normalized_x, normalized_t, normalized_frequency)


def dimensional_breather_observables(
    frequency: Any,
    coefficients: DimensionalSineGordonCoefficients,
) -> DimensionalBreatherObservables:
    """Return physical frequency, period, profile scales, energy, and action."""

    normalized_frequency = _normalized_frequency(frequency)
    scales = dimensional_sine_gordon_scales(coefficients)
    eta = breather_inverse_width(normalized_frequency)
    angular_frequency = normalized_frequency * scales.gap_frequency
    return DimensionalBreatherObservables(
        angular_frequency=angular_frequency,
        period=2 * sp.pi / angular_frequency,
        inverse_width=eta / scales.length,
        profile_length=scales.length / eta,
        energy=scales.energy * breather_energy(normalized_frequency),
        action=scales.action * breather_action(normalized_frequency),
    )


def dimensional_sine_gordon_log_ratio_jacobian() -> sp.ImmutableMatrix:
    r"""Return derivatives of ``log(c,omega_0,ell)`` by coefficient logs.

    Column order is ``(lambda,T,mu)``. The matrix has rank two and right
    nullspace spanned by ``(1,1,1)``, the common-multiplier direction.
    """

    return sp.ImmutableMatrix(
        [
            [-sp.Rational(1, 2), sp.Rational(1, 2), 0],
            [-sp.Rational(1, 2), 0, sp.Rational(1, 2)],
            [0, sp.Rational(1, 2), -sp.Rational(1, 2)],
        ]
    )


def dimensional_sine_gordon_coefficient_dimension_matrix() -> sp.ImmutableMatrix:
    r"""Return coefficient dimensions over rows ``(energy,length,time)``.

    Columns are ``(lambda,T,mu)`` with dimensions
    ``E*time**2/length``, ``E*length``, and ``E/length``. This matrix is
    full rank; dimensional independence and parameter identifiability remain
    different questions.
    """

    return sp.ImmutableMatrix(
        [
            [1, 1, 1],
            [-1, 1, -1],
            [2, 0, 0],
        ]
    )

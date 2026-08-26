"""Spherical stress and fixed-charge kinematics for the P240 M5 clock.

The module lifts the certified three-profile P240 ansatz into an orthonormal
spherical frame without assuming an equation of state.  It supplies the local
objects needed by an Einstein continuation; it does not claim that a supplied
profile solves either the matter or metric equations.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.chebyshev import chebder, chebval
from numpy.polynomial.legendre import legder, leggauss, legval
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import minimize

from .numerics import NumericalFailure, trapezoid_integral

FloatArray = NDArray[np.float64]

# Conventions and dimensions (c = 1, mostly plus):
#   ds^2 = -N(r)^2 dt^2 + dr^2/f(r) + r^2 dOmega^2.
#   q, tangent and split are dimensionless eigenvalue profiles; r has length.
#   The six curvature entries below are orthonormal pair energies before the
#   clock frequency is inserted.  Radial pairs acquire f, time pairs acquire
#   omega^2/N^2, and the proper radial measure is 4*pi*r^2/sqrt(f).
#   L_m = E_r + E_theta + E_phi
#         - B_rtheta - B_rphi - B_thetaphi - V.
#   T_ab = -2 dL_m/dg^(ab) + g_ab L_m fixes every pressure sign.


@dataclass(frozen=True)
class M5ClockProfiles:
    """Three radial P240 profiles and their coordinate derivatives."""

    q: FloatArray
    tangent: FloatArray
    split: FloatArray
    q_r: FloatArray
    tangent_r: FloatArray
    split_r: FloatArray


@dataclass(frozen=True)
class M5ClockKinematics:
    """Angular-averaged pair energies before metric and frequency factors."""

    clock_r: FloatArray
    clock_theta: FloatArray
    clock_phi: FloatArray
    curvature_rtheta: FloatArray
    curvature_rphi: FloatArray
    curvature_thetaphi: FloatArray
    potential: FloatArray

    @property
    def flat_inertia_density(self) -> FloatArray:
        return self.clock_r + self.clock_theta + self.clock_phi

    @property
    def flat_static_density(self) -> FloatArray:
        return (
            self.curvature_rtheta
            + self.curvature_rphi
            + self.curvature_thetaphi
            + self.potential
        )


@dataclass(frozen=True)
class M5ClockAxisymmetricKinematics:
    """Unaveraged pair energies on the radial-polar quadrature grid."""

    mu: FloatArray
    angular_weights: FloatArray
    clock_r: FloatArray
    clock_theta: FloatArray
    clock_phi: FloatArray
    curvature_rtheta: FloatArray
    curvature_rphi: FloatArray
    curvature_thetaphi: FloatArray
    potential: FloatArray
    pair_inner_products: FloatArray

    @property
    def flat_inertia_density(self) -> FloatArray:
        return self.clock_r + self.clock_theta + self.clock_phi

    @property
    def flat_static_density(self) -> FloatArray:
        return (
            self.curvature_rtheta
            + self.curvature_rphi
            + self.curvature_thetaphi
            + self.potential
        )


@dataclass(frozen=True)
class M5ClockStress:
    """Orthonormal stress components for a stationary clock profile."""

    energy_density: FloatArray
    radial_pressure: FloatArray
    polar_pressure: FloatArray
    azimuthal_pressure: FloatArray
    clock_energy: FloatArray
    curvature_energy: FloatArray
    potential_energy: FloatArray


@dataclass(frozen=True)
class M5ClockAxisymmetricStress:
    """Complete orthonormal stress on the radial-polar quadrature grid."""

    tensor: FloatArray

    @property
    def energy_density(self) -> FloatArray:
        return self.tensor[..., 0, 0]

    @property
    def momentum_density(self) -> FloatArray:
        return -self.tensor[..., 0, 1:4]


@dataclass(frozen=True)
class M5AngularRelaxation:
    """One regular even-angular fixed-J profile relaxation."""

    coefficients: FloatArray
    base_total_energy: float
    total_energy: float
    curvature_energy: float
    potential_energy: float
    inertia: float
    frequency: float
    optimizer_success: bool
    optimizer_message: str
    iterations: int
    gradient_inf_norm: float
    gradient_scale_relative: float
    refined_total_energy: float


def _finite_vector(
    value: ArrayLike, name: str, *, size: int | None = None
) -> FloatArray:
    result = np.asarray(value, dtype=np.float64)
    if result.ndim != 1 or (size is not None and result.size != size):
        expected = "a vector" if size is None else f"a vector of length {size}"
        raise ValueError(f"{name} must be {expected}")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return result


def m5_clock_profiles_from_chebyshev(
    coefficients: ArrayLike,
    radius: ArrayLike,
    domain_radius: float,
) -> M5ClockProfiles:
    """Evaluate the certified P240 ``(q,tangent,split)`` radial basis exactly."""

    points = _finite_vector(radius, "radius")
    outer = float(domain_radius)
    if not np.isfinite(outer) or outer <= 0.0:
        raise ValueError("domain_radius must be positive and finite")
    if np.any(points < 0.0) or np.any(points > outer):
        raise ValueError("radius must lie in [0, domain_radius]")
    modal_coefficients = np.asarray(coefficients, dtype=np.float64)
    if modal_coefficients.ndim == 1:
        if modal_coefficients.size % 3:
            raise ValueError(
                "flat coefficients must contain three equal profile blocks"
            )
        modal_coefficients = modal_coefficients.reshape(3, -1)
    if modal_coefficients.ndim != 2 or modal_coefficients.shape[0] != 3:
        raise ValueError("coefficients must have shape (3, order)")
    if not np.all(np.isfinite(modal_coefficients)):
        raise ValueError("coefficients must be finite")

    x = points / outer
    coordinate = 2.0 * x**2 - 1.0
    coordinate_r = 4.0 * x / outer
    modal = np.stack([chebval(coordinate, row) for row in modal_coefficients])
    modal_r = (
        np.stack([chebval(coordinate, chebder(row)) for row in modal_coefficients])
        * coordinate_r
    )

    one_minus_x2 = 1.0 - x**2
    q = x**2 * (1.0 + one_minus_x2 * modal[0])
    q_r = 2.0 * x / outer * (1.0 + one_minus_x2 * modal[0]) + x**2 * (
        -2.0 * x / outer * modal[0] + one_minus_x2 * modal_r[0]
    )
    tangent = one_minus_x2 * (1.0 / 3.0 + modal[1])
    tangent_r = -2.0 * x / outer * (1.0 / 3.0 + modal[1]) + one_minus_x2 * modal_r[1]
    split = x**4 * one_minus_x2 * modal[2]
    split_r = (4.0 * x**3 * one_minus_x2 - 2.0 * x**5) * modal[
        2
    ] / outer + x**4 * one_minus_x2 * modal_r[2]
    return M5ClockProfiles(q, tangent, split, q_r, tangent_r, split_r)


def _outer(vector: FloatArray) -> FloatArray:
    return np.einsum("...i,...j->...ij", vector, vector)


def _outer_derivative(vector: FloatArray, derivative: FloatArray) -> FloatArray:
    return np.einsum("...i,...j->...ij", derivative, vector) + np.einsum(
        "...i,...j->...ij", vector, derivative
    )


def _commutator(left: FloatArray, right: FloatArray) -> FloatArray:
    return left @ right - right @ left


def m5_clock_axisymmetric_kinematics(
    radius: ArrayLike,
    profiles: M5ClockProfiles,
    *,
    angular_quadrature_count: int = 32,
    split_angular_coefficients: ArrayLike | None = None,
) -> M5ClockAxisymmetricKinematics:
    """Return pair energies with optional regular even split deformations."""

    points = _finite_vector(radius, "radius")
    if np.any(points <= 0.0) or np.any(np.diff(points) <= 0.0):
        raise ValueError("radius must be strictly positive and increasing")
    if angular_quadrature_count < 8:
        raise ValueError("angular_quadrature_count must be at least eight")
    size = points.size
    values = tuple(
        _finite_vector(getattr(profiles, name), name, size=size)
        for name in ("q", "tangent", "split", "q_r", "tangent_r", "split_r")
    )
    q, tangent, split, q_r, tangent_r, split_r = (value[:, None] for value in values)

    mu, weights = leggauss(angular_quadrature_count)
    sine = np.sqrt(1.0 - mu**2)
    zero = np.zeros_like(mu)
    director = np.stack((sine, zero, mu), axis=-1)
    polar = np.stack((mu, zero, -sine), axis=-1)
    azimuthal = np.stack((zero, np.ones_like(mu), zero), axis=-1)
    director_mu = np.stack((-mu / sine, zero, np.ones_like(mu)), axis=-1)
    polar_mu = np.stack((np.ones_like(mu), zero, mu / sine), axis=-1)
    nn = _outer(director)[None]
    pp = _outer(polar)[None]
    aa = _outer(azimuthal)[None]
    nn_mu = _outer_derivative(director, director_mu)[None]
    pp_mu = _outer_derivative(polar, polar_mu)[None]

    angular_coefficients = (
        np.zeros(0, dtype=np.float64)
        if split_angular_coefficients is None
        else np.asarray(split_angular_coefficients, dtype=np.float64)
    )
    if angular_coefficients.ndim != 1 or not np.all(np.isfinite(angular_coefficients)):
        raise ValueError("split_angular_coefficients must be a finite vector")
    split_factor = np.ones_like(mu)
    split_factor_mu = np.zeros_like(mu)
    for index, coefficient in enumerate(angular_coefficients, start=1):
        polynomial = np.zeros(2 * index + 1)
        polynomial[-1] = 1.0
        split_factor += coefficient * legval(mu, polynomial)
        split_factor_mu += coefficient * legval(mu, legder(polynomial))
    angular_shape = sine**2 * split_factor
    angular_shape_mu = -2.0 * mu * split_factor + sine**2 * split_factor_mu
    delta = split * angular_shape[None, :]
    delta_r = split_r * angular_shape[None, :]
    delta_mu = split * angular_shape_mu[None, :]
    lambda_n = tangent + q
    lambda_n_r = tangent_r + q_r
    spatial = (
        lambda_n[..., None, None] * nn
        + (tangent + delta)[..., None, None] * pp
        + (tangent - delta)[..., None, None] * aa
    )
    derivative_r = (
        lambda_n_r[..., None, None] * nn
        + (tangent_r + delta_r)[..., None, None] * pp
        + (tangent_r - delta_r)[..., None, None] * aa
    )
    derivative_mu = (
        lambda_n[..., None, None] * nn_mu
        + (tangent + delta)[..., None, None] * pp_mu
        + delta_mu[..., None, None] * (pp - aa)
    )
    radial_grid = points[:, None, None, None]
    derivative_theta = -sine[None, :, None, None] * derivative_mu / radial_grid
    rotation_z = np.array(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        dtype=np.float64,
    )
    derivative_phi = (rotation_z @ spatial + spatial @ rotation_z.T) / (
        radial_grid * sine[None, :, None, None]
    )

    nx, ny, nz = np.moveaxis(director, -1, 0)
    clock_generator = np.stack(
        (
            np.stack((zero, -nz, ny), axis=-1),
            np.stack((nz, zero, -nx), axis=-1),
            np.stack((-ny, nx, zero), axis=-1),
        ),
        axis=-2,
    )
    clock_response = clock_generator[None] @ spatial + spatial @ np.swapaxes(
        clock_generator[None], -1, -2
    )
    pair_matrices = np.stack(
        (
            _commutator(clock_response, derivative_r),
            _commutator(clock_response, derivative_theta),
            _commutator(clock_response, derivative_phi),
            _commutator(derivative_r, derivative_theta),
            _commutator(derivative_r, derivative_phi),
            _commutator(derivative_theta, derivative_phi),
        ),
        axis=-3,
    )
    pair_inner_products = 4.0 * np.einsum(
        "...pij,...qij->...pq", pair_matrices, pair_matrices
    )
    spatial_two = spatial @ spatial
    trace_two = np.trace(spatial_two, axis1=-2, axis2=-1)
    trace_three = np.trace(spatial_two @ spatial, axis1=-2, axis2=-1)
    potential = -0.5 * trace_two - trace_three + trace_two**2 + 0.5
    result = M5ClockAxisymmetricKinematics(
        mu=mu,
        angular_weights=weights,
        clock_r=pair_inner_products[..., 0, 0],
        clock_theta=pair_inner_products[..., 1, 1],
        clock_phi=pair_inner_products[..., 2, 2],
        curvature_rtheta=pair_inner_products[..., 3, 3],
        curvature_rphi=pair_inner_products[..., 4, 4],
        curvature_thetaphi=pair_inner_products[..., 5, 5],
        potential=potential,
        pair_inner_products=pair_inner_products,
    )
    if not all(
        np.all(np.isfinite(getattr(result, field)))
        for field in result.__dataclass_fields__
    ):
        raise NumericalFailure("axisymmetric M5 clock kinematics became non-finite")
    return result


def _axisymmetric_metric_array(
    value: ArrayLike, name: str, shape: tuple[int, int]
) -> FloatArray:
    result = np.asarray(value, dtype=np.float64)
    if result.shape == (shape[0],):
        result = result[:, None]
    try:
        result = np.broadcast_to(result, shape)
    except ValueError as error:
        raise ValueError(
            f"{name} must broadcast to the radial-polar grid {shape}"
        ) from error
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be finite")
    return np.asarray(result)


def m5_clock_axisymmetric_stress(
    kinematics: M5ClockAxisymmetricKinematics,
    lapse: ArrayLike,
    radial_metric: ArrayLike,
    frequency: float,
) -> M5ClockAxisymmetricStress:
    """Return the complete diagonal-metric orthonormal M5 stress tensor."""

    shape = kinematics.potential.shape
    if len(shape) != 2:
        raise ValueError("axisymmetric kinematics must use a radial-polar grid")
    lapse_values = _axisymmetric_metric_array(lapse, "lapse", shape)
    metric_values = _axisymmetric_metric_array(radial_metric, "radial_metric", shape)
    omega = float(frequency)
    if not np.isfinite(omega):
        raise ValueError("frequency must be finite")
    if np.any(lapse_values <= 0.0) or np.any(metric_values <= 0.0):
        raise ValueError("lapse and radial_metric must be positive")

    root_f = np.sqrt(metric_values)
    clock_scale = omega / lapse_values
    scales = np.stack(
        (
            clock_scale * root_f,
            clock_scale,
            clock_scale,
            root_f,
            root_f,
            np.ones(shape),
        ),
        axis=-1,
    )
    gram = kinematics.pair_inner_products * scales[..., :, None] * scales[..., None, :]
    pair_diagonal = np.diagonal(gram, axis1=-2, axis2=-1)
    lagrangian = (
        np.sum(pair_diagonal[..., :3], axis=-1)
        - np.sum(pair_diagonal[..., 3:], axis=-1)
        - kinematics.potential
    )
    pair_indices = {
        (0, 1): 0,
        (0, 2): 1,
        (0, 3): 2,
        (1, 2): 3,
        (1, 3): 4,
        (2, 3): 5,
    }

    def oriented_pair(left: int, right: int) -> tuple[int, int]:
        if left < right:
            return pair_indices[(left, right)], 1
        return pair_indices[(right, left)], -1

    signature = np.array([-1.0, 1.0, 1.0, 1.0])
    tensor = np.zeros(shape + (4, 4), dtype=np.float64)
    for left in range(4):
        for right in range(4):
            contraction = np.zeros(shape, dtype=np.float64)
            for contracted in range(4):
                if left == contracted or right == contracted:
                    continue
                left_pair, left_sign = oriented_pair(left, contracted)
                right_pair, right_sign = oriented_pair(right, contracted)
                contraction += (
                    signature[contracted]
                    * left_sign
                    * right_sign
                    * gram[..., left_pair, right_pair]
                )
            tensor[..., left, right] = 2.0 * contraction
            if left == right:
                tensor[..., left, right] += signature[left] * lagrangian
    if not np.all(np.isfinite(tensor)):
        raise NumericalFailure("axisymmetric M5 clock stress became non-finite")
    return M5ClockAxisymmetricStress(tensor=tensor)


def m5_clock_kinematics(
    radius: ArrayLike,
    profiles: M5ClockProfiles,
    *,
    angular_quadrature_count: int = 32,
) -> M5ClockKinematics:
    """Return angular-averaged radial, angular, and clock pair energies."""

    angular = m5_clock_axisymmetric_kinematics(
        radius,
        profiles,
        angular_quadrature_count=angular_quadrature_count,
    )

    def average(value: FloatArray) -> FloatArray:
        return np.einsum("ra,a->r", value, angular.angular_weights) / 2.0

    result = M5ClockKinematics(
        clock_r=average(angular.clock_r),
        clock_theta=average(angular.clock_theta),
        clock_phi=average(angular.clock_phi),
        curvature_rtheta=average(angular.curvature_rtheta),
        curvature_rphi=average(angular.curvature_rphi),
        curvature_thetaphi=average(angular.curvature_thetaphi),
        potential=average(angular.potential),
    )
    if not all(
        np.all(np.isfinite(getattr(result, field)))
        for field in result.__dataclass_fields__
    ):
        raise NumericalFailure("M5 clock kinematics became non-finite")
    return result


def m5_clock_fixed_j_frequency(
    radius: ArrayLike,
    kinematics: M5ClockKinematics,
    lapse: ArrayLike,
    radial_metric: ArrayLike,
    angular_momentum: float,
) -> tuple[float, float]:
    """Return ``(omega, I)`` from the curved fixed-J Routhian coefficient."""

    points = _finite_vector(radius, "radius")
    lapse_values = _finite_vector(lapse, "lapse", size=points.size)
    metric_values = _finite_vector(radial_metric, "radial_metric", size=points.size)
    momentum = float(angular_momentum)
    if momentum <= 0.0 or not np.isfinite(momentum):
        raise ValueError("angular_momentum must be positive and finite")
    if np.any(lapse_values <= 0.0) or np.any(metric_values <= 0.0):
        raise ValueError("lapse and radial_metric must be positive")
    inertia_density = (
        metric_values * kinematics.clock_r
        + kinematics.clock_theta
        + kinematics.clock_phi
    )
    integrand = (
        4.0
        * np.pi
        * points**2
        * inertia_density
        / (lapse_values * np.sqrt(metric_values))
    )
    inertia = trapezoid_integral(integrand, points)
    if inertia <= 0.0 or not np.isfinite(inertia):
        raise NumericalFailure("curved M5 clock inertia must be positive and finite")
    return momentum / (2.0 * inertia), inertia


def m5_clock_stress(
    kinematics: M5ClockKinematics,
    lapse: ArrayLike,
    radial_metric: ArrayLike,
    frequency: float,
) -> M5ClockStress:
    """Metric-vary the local M5 density into orthonormal ``rho,p_r,p_theta,p_phi``."""

    size = kinematics.potential.size
    lapse_values = _finite_vector(lapse, "lapse", size=size)
    metric_values = _finite_vector(radial_metric, "radial_metric", size=size)
    omega = float(frequency)
    if not np.isfinite(omega):
        raise ValueError("frequency must be finite")
    if np.any(lapse_values <= 0.0) or np.any(metric_values <= 0.0):
        raise ValueError("lapse and radial_metric must be positive")

    redshift = omega**2 / lapse_values**2
    e_r = redshift * metric_values * kinematics.clock_r
    e_theta = redshift * kinematics.clock_theta
    e_phi = redshift * kinematics.clock_phi
    b_rtheta = metric_values * kinematics.curvature_rtheta
    b_rphi = metric_values * kinematics.curvature_rphi
    b_thetaphi = kinematics.curvature_thetaphi
    potential = kinematics.potential
    clock = e_r + e_theta + e_phi
    curvature = b_rtheta + b_rphi + b_thetaphi
    rho = clock + curvature + potential
    p_r = -e_r + e_theta + e_phi + b_rtheta + b_rphi - b_thetaphi - potential
    p_theta = e_r - e_theta + e_phi + b_rtheta - b_rphi + b_thetaphi - potential
    p_phi = e_r + e_theta - e_phi - b_rtheta + b_rphi + b_thetaphi - potential
    result = M5ClockStress(rho, p_r, p_theta, p_phi, clock, curvature, potential)
    if not all(
        np.all(np.isfinite(getattr(result, field)))
        for field in result.__dataclass_fields__
    ):
        raise NumericalFailure("M5 clock stress became non-finite")
    return result


def relax_m5_clock_split_angular_profile(
    radial_coefficients: ArrayLike,
    *,
    domain_radius: float,
    mode_count: int = 4,
    radial_quadrature_count: int = 96,
    angular_quadrature_count: int = 72,
    refinement_angular_count: int = 96,
    initial_coefficients: ArrayLike | None = None,
    maximum_iterations: int = 500,
) -> M5AngularRelaxation:
    """Relax regular even Legendre deformations of the tangent split."""

    coefficients = np.asarray(radial_coefficients, dtype=np.float64)
    if coefficients.ndim == 1:
        if coefficients.size % 3:
            raise ValueError("radial_coefficients must have three equal blocks")
        coefficients = coefficients.reshape(3, -1)
    if coefficients.ndim != 2 or coefficients.shape[0] != 3:
        raise ValueError("radial_coefficients must have shape (3, order)")
    if not np.all(np.isfinite(coefficients)):
        raise ValueError("radial_coefficients must be finite")
    outer = float(domain_radius)
    if not np.isfinite(outer) or outer <= 0.0:
        raise ValueError("domain_radius must be positive and finite")
    if mode_count <= 0:
        raise ValueError("mode_count must be positive")
    if radial_quadrature_count < 16 or angular_quadrature_count < 8:
        raise ValueError("quadrature counts are too small")
    if refinement_angular_count < angular_quadrature_count:
        raise ValueError("refinement_angular_count must not be smaller")
    if maximum_iterations <= 0:
        raise ValueError("maximum_iterations must be positive")

    nodes, raw_weights = leggauss(radial_quadrature_count)
    radius = 0.5 * outer * (nodes + 1.0)
    radial_weights = 0.5 * outer * raw_weights
    profiles = m5_clock_profiles_from_chebyshev(coefficients, radius, outer)

    def components(
        angular_coefficients: FloatArray, angular_count: int
    ) -> tuple[float, float, float, float]:
        kinematics = m5_clock_axisymmetric_kinematics(
            radius,
            profiles,
            angular_quadrature_count=angular_count,
            split_angular_coefficients=angular_coefficients,
        )

        def integrate(values: FloatArray) -> float:
            return float(
                2.0
                * np.pi
                * np.einsum(
                    "r,a,ra->",
                    radial_weights * radius**2,
                    kinematics.angular_weights,
                    values,
                )
            )

        curvature = integrate(
            kinematics.curvature_rtheta
            + kinematics.curvature_rphi
            + kinematics.curvature_thetaphi
        )
        potential = integrate(kinematics.potential)
        inertia = integrate(kinematics.flat_inertia_density)
        if inertia <= 0.0:
            raise NumericalFailure("angular relaxation reached zero inertia")
        total = curvature + potential + 1.0 / (4.0 * inertia)
        return total, curvature, potential, inertia

    zero = np.zeros(mode_count, dtype=np.float64)
    base_total = components(zero, angular_quadrature_count)[0]
    if initial_coefficients is None:
        initial = zero
    else:
        supplied = _finite_vector(
            initial_coefficients, "initial_coefficients", size=mode_count
        )
        initial = np.array(supplied, copy=True)

    def objective(values: FloatArray) -> float:
        return components(values, angular_quadrature_count)[0]

    starts = (
        initial,
        zero,
        np.linspace(0.15, -0.1, mode_count, dtype=np.float64),
    )
    best = None
    for start in starts:
        result = minimize(
            objective,
            start,
            method="L-BFGS-B",
            jac="3-point",
            bounds=[(-3.0, 3.0)] * mode_count,
            options={
                "ftol": 1.0e-14,
                "gtol": 1.0e-9,
                "maxiter": maximum_iterations,
                "maxls": 40,
                "finite_diff_rel_step": 1.0e-5,
            },
        )
        if best is None or result.fun < best.fun:
            best = result
    if best is None:
        raise NumericalFailure("no angular relaxation was attempted")
    final_coefficients = np.asarray(best.x, dtype=np.float64)
    total, curvature, potential, inertia = components(
        final_coefficients, angular_quadrature_count
    )
    refined_total = components(final_coefficients, refinement_angular_count)[0]
    gradient = np.empty(mode_count, dtype=np.float64)
    for index in range(mode_count):
        step = 1.0e-5 * (1.0 + abs(final_coefficients[index]))
        plus = np.array(final_coefficients, copy=True)
        minus = np.array(final_coefficients, copy=True)
        plus[index] += step
        minus[index] -= step
        gradient[index] = (objective(plus) - objective(minus)) / (2.0 * step)
    gradient_inf = float(np.max(np.abs(gradient)))
    return M5AngularRelaxation(
        coefficients=final_coefficients,
        base_total_energy=base_total,
        total_energy=total,
        curvature_energy=curvature,
        potential_energy=potential,
        inertia=inertia,
        frequency=1.0 / (2.0 * inertia),
        optimizer_success=bool(best.success),
        optimizer_message=str(best.message),
        iterations=int(best.nit),
        gradient_inf_norm=gradient_inf,
        gradient_scale_relative=gradient_inf / max(1.0, abs(total)),
        refined_total_energy=refined_total,
    )

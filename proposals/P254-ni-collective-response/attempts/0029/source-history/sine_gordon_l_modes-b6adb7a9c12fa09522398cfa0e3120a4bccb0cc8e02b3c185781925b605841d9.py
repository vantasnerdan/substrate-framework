"""Regular angular perturbations of a radial 3+1 sine-Gordon field.

For a radial background ``P(r,t)`` satisfying

``P_tt = P_rr + 2*P_r/r - sin(P)``,

the infinitesimal expansion ``u=P+epsilon*psi_l(r,t)*Y_lm`` gives

``psi_tt = psi_rr + 2*psi_r/r - l*(l+1)*psi/r**2 - cos(P)*psi``.

The numerical evolution uses ``v=r*psi`` so the first radial derivative is
removed.  A regular mode obeys ``psi=O(r**l)`` and hence ``v=O(r**(l+1))`` at
the origin.  Results are finite-grid, finite-time linearized evidence; this
module does not assert nonlinear stability, exact periodicity, gravitational
radiation, or an absolute physical scale.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .numerics import NumericalFailure, trapezoid_integral
from .radial_sine_gordon import (
    radial_gradient,
    radial_laplacian,
    radial_sine_gordon_energy,
)

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class LinearizedAngularModeEvolution:
    """Diagnostics from a radial background plus one linearized angular mode."""

    radius: FloatArray
    time: FloatArray
    background_center: FloatArray
    mode_norm: FloatArray
    p2_triple_stf_zz_coefficient: FloatArray
    background_energy: FloatArray
    quadratic_mode_energy: FloatArray
    final_background: FloatArray
    final_background_velocity: FloatArray
    final_mode: FloatArray
    final_mode_velocity: FloatArray
    spacing: float
    timestep: float
    outer_radius: float
    ell: int
    method: str
    completed: bool
    max_boundary_background: float
    max_boundary_mode: float


def legendre_p2(cos_theta: ArrayLike) -> FloatArray:
    """Return ``P2(mu)=(3*mu**2-1)/2`` for finite ``mu`` in ``[-1,1]``."""

    mu = np.asarray(cos_theta, dtype=np.float64)
    if not np.all(np.isfinite(mu)) or np.any(np.abs(mu) > 1.0):
        raise ValueError("cos_theta must contain finite values in [-1, 1]")
    return np.asarray(0.5 * (3.0 * np.square(mu) - 1.0), dtype=np.float64)


def multiplicative_p2_residual(
    background: ArrayLike,
    radius: ArrayLike,
    amplitude: float,
    cos_theta: ArrayLike,
) -> FloatArray:
    """Return the full-PDE residual of ``u=P*(1+a*P2)`` on ``r>0``.

    The background radial equation has already been used.  Thus the result is

    ``sin(P*(1+a*Y))-(1+a*Y)*sin(P)+6*a*P*Y/r**2``.

    Arrays follow NumPy broadcasting.  The origin is deliberately rejected:
    unless the angular part vanishes there, this multiplicative deformation
    is not a regular Cartesian scalar field.
    """

    field = np.asarray(background, dtype=np.float64)
    r = np.asarray(radius, dtype=np.float64)
    if not np.all(np.isfinite(field)):
        raise ValueError("background must contain only finite values")
    if not np.all(np.isfinite(r)) or np.any(r <= 0.0):
        raise ValueError("radius must contain finite values strictly above zero")
    if not np.isfinite(amplitude):
        raise ValueError("amplitude must be finite")
    harmonic = legendre_p2(cos_theta)
    factor = 1.0 + amplitude * harmonic
    residual = (
        np.sin(field * factor)
        - factor * np.sin(field)
        + 6.0 * amplitude * field * harmonic / np.square(r)
    )
    return np.asarray(residual, dtype=np.float64)


def multiplicative_p2_first_order_residual_coefficient(
    background: ArrayLike,
    radius: ArrayLike,
    cos_theta: ArrayLike,
) -> FloatArray:
    """Return ``d residual/da`` at ``a=0`` for the multiplicative ansatz."""

    field = np.asarray(background, dtype=np.float64)
    r = np.asarray(radius, dtype=np.float64)
    if not np.all(np.isfinite(field)):
        raise ValueError("background must contain only finite values")
    if not np.all(np.isfinite(r)) or np.any(r <= 0.0):
        raise ValueError("radius must contain finite values strictly above zero")
    harmonic = legendre_p2(cos_theta)
    return np.asarray(
        harmonic
        * (field * np.cos(field) - np.sin(field) + 6.0 * field / np.square(r)),
        dtype=np.float64,
    )


def regular_l_mode_gaussian_seed(
    radius: ArrayLike,
    *,
    ell: int,
    amplitude: float,
    width: float,
) -> FloatArray:
    """Return ``amplitude*(r/width)**ell*exp(-(r/width)**2)``.

    The seed has the required ``O(r**ell)`` origin behavior.  ``ell`` must be
    positive because the radial sector is already implemented separately.
    """

    r = _radius_vector(radius)
    mode = _positive_ell(ell)
    if not np.isfinite(amplitude):
        raise ValueError("amplitude must be finite")
    scale = _positive_finite(width, "width")
    coordinate = r / scale
    return np.asarray(
        amplitude * np.power(coordinate, mode) * np.exp(-np.square(coordinate)),
        dtype=np.float64,
    )


def transformed_l_mode_acceleration(
    background: ArrayLike,
    transformed_mode: ArrayLike,
    spacing: float,
    *,
    ell: int,
) -> FloatArray:
    """Return ``v_rr-l(l+1)v/r**2-cos(P)v`` for ``v=r*psi``.

    Both endpoints are boundary values and receive zero acceleration.  The
    origin value of a regular transformed mode is exactly zero.
    """

    field = _field_vector(background, "background")
    mode = _matching_vector(transformed_mode, field.size, "transformed_mode")
    h = _positive_finite(spacing, "spacing")
    angular_index = _positive_ell(ell)
    radius = h * np.arange(field.size, dtype=np.float64)
    acceleration = np.zeros_like(mode)
    acceleration[1:-1] = (
        (mode[2:] - 2.0 * mode[1:-1] + mode[:-2]) / h**2
        - angular_index
        * (angular_index + 1)
        * mode[1:-1]
        / np.square(radius[1:-1])
        - np.cos(field[1:-1]) * mode[1:-1]
    )
    return acceleration


def linearized_p2_energy_triple_stf(
    background: ArrayLike,
    background_velocity: ArrayLike,
    mode: ArrayLike,
    mode_velocity: ArrayLike,
    radius: ArrayLike,
) -> FloatArray:
    """Return the triple-STF energy moment at first order in P2 amplitude.

    If ``u=P+epsilon*psi*P2(cos(theta))``, the coefficient of
    ``epsilon*P2`` in ``T00`` is
    ``h=P_t*psi_t+P_r*psi_r+sin(P)*psi``.  With
    ``H=4*pi*integral r**4*h dr``, exact angular integration gives
    ``Q/epsilon=diag(-H/5,-H/5,2H/5)``.  This is moment kinematics only.
    """

    field = _field_vector(background, "background")
    velocity = _matching_vector(
        background_velocity, field.size, "background_velocity"
    )
    perturbation = _matching_vector(mode, field.size, "mode")
    perturbation_velocity = _matching_vector(
        mode_velocity, field.size, "mode_velocity"
    )
    r = _radius_vector(radius, expected_size=field.size)
    spacing = _uniform_spacing(r)
    density_coefficient = (
        velocity * perturbation_velocity
        + radial_gradient(field, spacing) * radial_gradient(perturbation, spacing)
        + np.sin(field) * perturbation
    )
    scalar = 4.0 * np.pi * trapezoid_integral(
        np.power(r, 4) * density_coefficient, r
    )
    return np.diag([-scalar / 5.0, -scalar / 5.0, 2.0 * scalar / 5.0])


def linearized_mode_quadratic_energy(
    background: ArrayLike,
    mode: ArrayLike,
    mode_velocity: ArrayLike,
    radius: ArrayLike,
    *,
    ell: int,
) -> float:
    """Return the instantaneous quadratic energy of ``psi*P_l``.

    The angular convention is the unnormalized Legendre mode, for which
    ``integral P_l(cos(theta))**2 dOmega=4*pi/(2*l+1)``.  Because ``cos(P)``
    is time dependent for a general background, this diagnostic is conserved
    only in static backgrounds such as ``P=0``.
    """

    field = _field_vector(background, "background")
    perturbation = _matching_vector(mode, field.size, "mode")
    velocity = _matching_vector(mode_velocity, field.size, "mode_velocity")
    r = _radius_vector(radius, expected_size=field.size)
    angular_index = _positive_ell(ell)
    spacing = _uniform_spacing(r)
    gradient = radial_gradient(perturbation, spacing)
    integrand = (
        np.square(r) * np.square(velocity)
        + np.square(r) * np.square(gradient)
        + angular_index * (angular_index + 1) * np.square(perturbation)
        + np.square(r) * np.cos(field) * np.square(perturbation)
    )
    return (
        2.0
        * np.pi
        / (2 * angular_index + 1)
        * trapezoid_integral(integrand, r)
    )


def evolve_radial_background_with_linearized_mode(
    background_initial: ArrayLike,
    mode_initial: ArrayLike,
    *,
    spacing: float,
    final_time: float,
    ell: int,
    background_velocity_initial: ArrayLike | None = None,
    mode_velocity_initial: ArrayLike | None = None,
    courant: float = 0.4,
    sample_interval: float = 0.2,
    boundary_monitor_fraction: float = 0.8,
) -> LinearizedAngularModeEvolution:
    """Evolve a radial background and one regular linearized angular mode.

    A velocity-Verlet discretization advances the nonlinear radial background
    and ``v=r*psi`` on the common uniform grid ``r_j=j*spacing``.  Homogeneous
    Dirichlet data are imposed at the finite outer boundary.  No sponge is
    used, so a claim must stop before reflected characteristics return to the
    region it interprets or demonstrate domain independence explicitly.
    """

    background = _field_vector(background_initial, "background_initial", copy=True)
    mode = _matching_vector(mode_initial, background.size, "mode_initial", copy=True)
    h = _positive_finite(spacing, "spacing")
    duration = _positive_finite(final_time, "final_time")
    angular_index = _positive_ell(ell)
    cfl = _positive_finite(courant, "courant")
    if cfl >= 1.0:
        raise ValueError("courant must be below one")
    interval = _positive_finite(sample_interval, "sample_interval")
    if (
        not np.isfinite(boundary_monitor_fraction)
        or boundary_monitor_fraction <= 0.0
        or boundary_monitor_fraction >= 1.0
    ):
        raise ValueError("boundary_monitor_fraction must lie strictly between zero and one")

    radius = h * np.arange(background.size, dtype=np.float64)
    if background.size < 5:
        raise ValueError("fields must contain at least five radial points")
    if abs(mode[0]) > 1.0e-14 * max(1.0, float(np.max(np.abs(mode)))):
        raise ValueError("a regular positive-l mode must vanish at the origin")
    background_velocity = _optional_velocity(
        background_velocity_initial, background.size, "background_velocity_initial"
    )
    mode_velocity = _optional_velocity(
        mode_velocity_initial, background.size, "mode_velocity_initial"
    )
    if abs(mode_velocity[0]) > 1.0e-14 * max(
        1.0, float(np.max(np.abs(mode_velocity)))
    ):
        raise ValueError("a regular positive-l mode velocity must vanish at the origin")

    background[-1] = 0.0
    background_velocity[-1] = 0.0
    mode[-1] = 0.0
    mode_velocity[-1] = 0.0
    transformed = radius * mode
    transformed_velocity = radius * mode_velocity

    requested_dt = cfl * h
    steps = int(np.ceil(duration / requested_dt))
    timestep = duration / steps
    record_stride = max(1, int(round(interval / timestep)))
    boundary_mask = radius >= boundary_monitor_fraction * radius[-1]
    max_boundary_background = 0.0
    max_boundary_mode = 0.0

    times: list[float] = []
    centers: list[float] = []
    mode_norms: list[float] = []
    qzz_coefficients: list[float] = []
    background_energies: list[float] = []
    mode_energies: list[float] = []

    def recover_mode(values: FloatArray) -> FloatArray:
        recovered = np.zeros_like(values)
        recovered[1:] = values[1:] / radius[1:]
        return recovered

    def record(time: float) -> None:
        nonlocal max_boundary_background, max_boundary_mode
        current_mode = recover_mode(transformed)
        current_mode_velocity = recover_mode(transformed_velocity)
        if not all(
            np.all(np.isfinite(values))
            for values in (
                background,
                background_velocity,
                current_mode,
                current_mode_velocity,
            )
        ):
            raise NumericalFailure("linearized angular-mode evolution became non-finite")
        tensor = linearized_p2_energy_triple_stf(
            background,
            background_velocity,
            current_mode,
            current_mode_velocity,
            radius,
        )
        norm_squared = 4.0 * np.pi / (2 * angular_index + 1) * trapezoid_integral(
            np.square(transformed), radius
        )
        times.append(float(time))
        centers.append(float(background[0]))
        mode_norms.append(float(np.sqrt(max(norm_squared, 0.0))))
        qzz_coefficients.append(float(tensor[2, 2]))
        background_energies.append(
            radial_sine_gordon_energy(background, background_velocity, radius)
        )
        mode_energies.append(
            linearized_mode_quadratic_energy(
                background,
                current_mode,
                current_mode_velocity,
                radius,
                ell=angular_index,
            )
        )
        max_boundary_background = max(
            max_boundary_background,
            float(np.max(np.abs(background[boundary_mask]))),
        )
        max_boundary_mode = max(
            max_boundary_mode,
            float(np.max(np.abs(current_mode[boundary_mask]))),
        )

    record(0.0)
    for step in range(1, steps + 1):
        background_acceleration = radial_laplacian(background, h) - np.sin(background)
        transformed_acceleration = transformed_l_mode_acceleration(
            background, transformed, h, ell=angular_index
        )
        next_background = (
            background
            + timestep * background_velocity
            + 0.5 * timestep**2 * background_acceleration
        )
        next_transformed = (
            transformed
            + timestep * transformed_velocity
            + 0.5 * timestep**2 * transformed_acceleration
        )
        next_background[-1] = 0.0
        next_transformed[0] = 0.0
        next_transformed[-1] = 0.0
        next_background_acceleration = radial_laplacian(next_background, h) - np.sin(
            next_background
        )
        next_transformed_acceleration = transformed_l_mode_acceleration(
            next_background, next_transformed, h, ell=angular_index
        )
        next_background_velocity = background_velocity + 0.5 * timestep * (
            background_acceleration + next_background_acceleration
        )
        next_transformed_velocity = transformed_velocity + 0.5 * timestep * (
            transformed_acceleration + next_transformed_acceleration
        )
        next_background_velocity[-1] = 0.0
        next_transformed_velocity[0] = 0.0
        next_transformed_velocity[-1] = 0.0
        background = next_background
        background_velocity = next_background_velocity
        transformed = next_transformed
        transformed_velocity = next_transformed_velocity
        if step % record_stride == 0 or step == steps:
            record(step * timestep)

    final_mode = recover_mode(transformed)
    final_mode_velocity = recover_mode(transformed_velocity)
    return LinearizedAngularModeEvolution(
        radius=radius,
        time=np.asarray(times, dtype=np.float64),
        background_center=np.asarray(centers, dtype=np.float64),
        mode_norm=np.asarray(mode_norms, dtype=np.float64),
        p2_triple_stf_zz_coefficient=np.asarray(qzz_coefficients, dtype=np.float64),
        background_energy=np.asarray(background_energies, dtype=np.float64),
        quadratic_mode_energy=np.asarray(mode_energies, dtype=np.float64),
        final_background=np.array(background, copy=True),
        final_background_velocity=np.array(background_velocity, copy=True),
        final_mode=final_mode,
        final_mode_velocity=final_mode_velocity,
        spacing=h,
        timestep=timestep,
        outer_radius=float(radius[-1]),
        ell=angular_index,
        method="velocity-verlet-transformed-mode",
        completed=True,
        max_boundary_background=max_boundary_background,
        max_boundary_mode=max_boundary_mode,
    )


def _field_vector(values: ArrayLike, name: str, *, copy: bool = False) -> FloatArray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size < 3:
        raise ValueError(f"{name} must be a one-dimensional array with at least three values")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return np.array(array, copy=copy)


def _matching_vector(
    values: ArrayLike,
    expected_size: int,
    name: str,
    *,
    copy: bool = False,
) -> FloatArray:
    array = _field_vector(values, name, copy=copy)
    if array.size != expected_size:
        raise ValueError(f"{name} must have the same size as the background")
    return array


def _optional_velocity(
    values: ArrayLike | None, expected_size: int, name: str
) -> FloatArray:
    if values is None:
        return np.zeros(expected_size, dtype=np.float64)
    return _matching_vector(values, expected_size, name, copy=True)


def _radius_vector(values: ArrayLike, expected_size: int | None = None) -> FloatArray:
    radius = np.asarray(values, dtype=np.float64)
    if radius.ndim != 1 or radius.size < 3:
        raise ValueError("radius must be a one-dimensional array with at least three values")
    if not np.all(np.isfinite(radius)) or radius[0] != 0.0 or np.any(np.diff(radius) <= 0.0):
        raise ValueError("radius must be finite, strictly increasing, and start at zero")
    if expected_size is not None and radius.size != expected_size:
        raise ValueError("radius must have the same size as the fields")
    return radius


def _uniform_spacing(radius: FloatArray) -> float:
    steps = np.diff(radius)
    spacing = float(np.mean(steps))
    if not np.allclose(steps, spacing, rtol=1.0e-10, atol=1.0e-13):
        raise ValueError("radius must be uniformly spaced")
    return spacing


def _positive_finite(value: float, name: str) -> float:
    converted = float(value)
    if not np.isfinite(converted) or converted <= 0.0:
        raise ValueError(f"{name} must be positive and finite")
    return converted


def _positive_ell(ell: int) -> int:
    if isinstance(ell, bool) or int(ell) != ell or ell < 1:
        raise ValueError("ell must be a positive integer")
    return int(ell)

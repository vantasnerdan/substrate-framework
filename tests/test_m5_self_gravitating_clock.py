from __future__ import annotations

import numpy as np

from substrate_framework.m5_self_gravitating_clock import (
    M5ClockKinematics,
    m5_clock_axisymmetric_kinematics,
    m5_clock_axisymmetric_stress,
    m5_clock_fixed_j_frequency,
    m5_clock_kinematics,
    m5_clock_profiles_from_chebyshev,
    m5_clock_stress,
    relax_m5_clock_split_angular_profile,
)


def test_chebyshev_profile_boundary_and_origin_regularities() -> None:
    radius = np.array([0.0, 0.25, 1.0])
    coefficients = np.array([[0.2, -0.1], [0.05, 0.02], [-0.3, 0.1]], dtype=np.float64)
    profiles = m5_clock_profiles_from_chebyshev(coefficients, radius, 1.0)

    assert profiles.q[0] == 0.0
    assert profiles.q[-1] == 1.0
    assert profiles.tangent[-1] == 0.0
    assert profiles.split[0] == 0.0
    assert profiles.split[-1] == 0.0
    assert profiles.q_r[0] == 0.0
    assert profiles.split_r[0] == 0.0


def test_split_profile_exposes_anisotropic_stress_and_recovers_flat_density() -> None:
    radius = np.linspace(0.02, 1.0, 31)
    coefficients = np.array(
        [
            [0.1, -0.03, 0.01],
            [0.02, -0.01, 0.004],
            [0.15, -0.04, 0.01],
        ],
        dtype=np.float64,
    )
    profiles = m5_clock_profiles_from_chebyshev(coefficients, radius, 1.0)
    kinematics = m5_clock_kinematics(radius, profiles, angular_quadrature_count=32)
    omega = 0.7
    stress = m5_clock_stress(
        kinematics, np.ones_like(radius), np.ones_like(radius), omega
    )

    assert np.max(np.abs(stress.polar_pressure - stress.azimuthal_pressure)) > 1.0e-4
    np.testing.assert_allclose(
        stress.energy_density,
        omega**2 * kinematics.flat_inertia_density + kinematics.flat_static_density,
        rtol=2.0e-14,
        atol=2.0e-14,
    )

    coefficients[2] = 0.0
    spherical_profiles = m5_clock_profiles_from_chebyshev(coefficients, radius, 1.0)
    spherical_kinematics = m5_clock_kinematics(
        radius, spherical_profiles, angular_quadrature_count=32
    )
    spherical_stress = m5_clock_stress(
        spherical_kinematics,
        np.ones_like(radius),
        np.ones_like(radius),
        omega,
    )
    np.testing.assert_allclose(
        spherical_stress.polar_pressure,
        spherical_stress.azimuthal_pressure,
        rtol=3.0e-13,
        atol=3.0e-13,
    )


def test_complete_axisymmetric_stress_agrees_with_diagonal_variation() -> None:
    radius = np.linspace(0.03, 1.0, 23)
    coefficients = np.array(
        [
            [0.08, -0.02, 0.01],
            [0.03, -0.01, 0.002],
            [0.12, -0.03, 0.008],
        ],
        dtype=np.float64,
    )
    profiles = m5_clock_profiles_from_chebyshev(coefficients, radius, 1.0)
    angular = m5_clock_axisymmetric_kinematics(
        radius, profiles, angular_quadrature_count=24
    )
    complete = m5_clock_axisymmetric_stress(
        angular,
        lapse=np.linspace(0.8, 1.0, radius.size),
        radial_metric=np.linspace(0.6, 1.0, radius.size),
        frequency=0.7,
    )
    averaged = m5_clock_kinematics(radius, profiles, angular_quadrature_count=24)
    diagonal = m5_clock_stress(
        averaged,
        lapse=np.linspace(0.8, 1.0, radius.size),
        radial_metric=np.linspace(0.6, 1.0, radius.size),
        frequency=0.7,
    )

    averaged_tensor = (
        np.einsum("raij,a->rij", complete.tensor, angular.angular_weights) / 2.0
    )
    np.testing.assert_allclose(
        averaged_tensor[..., 0, 0], diagonal.energy_density, rtol=2.0e-13
    )
    np.testing.assert_allclose(
        averaged_tensor[..., 1, 1], diagonal.radial_pressure, rtol=2.0e-13
    )
    np.testing.assert_allclose(
        averaged_tensor[..., 2, 2], diagonal.polar_pressure, rtol=2.0e-13
    )
    np.testing.assert_allclose(
        averaged_tensor[..., 3, 3],
        diagonal.azimuthal_pressure,
        rtol=2.0e-13,
    )
    np.testing.assert_allclose(
        complete.tensor,
        np.swapaxes(complete.tensor, -1, -2),
        rtol=2.0e-13,
        atol=2.0e-13,
    )


def test_regular_angular_relaxation_is_bounded_and_refinement_reported() -> None:
    coefficients = np.array(
        [
            [0.1, -0.03, 0.01],
            [0.02, -0.01, 0.004],
            [0.15, -0.04, 0.01],
        ],
        dtype=np.float64,
    )
    result = relax_m5_clock_split_angular_profile(
        coefficients,
        domain_radius=1.0,
        mode_count=1,
        radial_quadrature_count=16,
        angular_quadrature_count=8,
        refinement_angular_count=12,
        maximum_iterations=80,
    )

    assert result.coefficients.shape == (1,)
    assert result.optimizer_success
    assert result.total_energy <= result.base_total_energy + 1.0e-12
    assert np.isfinite(result.refined_total_energy)
    assert np.isfinite(result.gradient_scale_relative)


def test_metric_variation_has_quartic_trace_and_potential_trace() -> None:
    values = np.array([0.3, 0.7])
    kinematics = M5ClockKinematics(
        clock_r=values,
        clock_theta=2.0 * values,
        clock_phi=3.0 * values,
        curvature_rtheta=4.0 * values,
        curvature_rphi=5.0 * values,
        curvature_thetaphi=6.0 * values,
        potential=7.0 * values,
    )
    stress = m5_clock_stress(
        kinematics,
        lapse=np.array([0.8, 1.2]),
        radial_metric=np.array([0.6, 1.4]),
        frequency=0.9,
    )
    trace = (
        -stress.energy_density
        + stress.radial_pressure
        + stress.polar_pressure
        + stress.azimuthal_pressure
    )

    np.testing.assert_allclose(trace, -4.0 * kinematics.potential, rtol=1.0e-14)


def test_curved_fixed_j_frequency_uses_lapse_and_proper_radial_measure() -> None:
    radius = np.array([1.0, 2.0])
    ones = np.ones(2)
    kinematics = M5ClockKinematics(
        clock_r=ones,
        clock_theta=2.0 * ones,
        clock_phi=3.0 * ones,
        curvature_rtheta=np.zeros(2),
        curvature_rphi=np.zeros(2),
        curvature_thetaphi=np.zeros(2),
        potential=np.zeros(2),
    )
    omega, inertia = m5_clock_fixed_j_frequency(
        radius, kinematics, ones, ones, angular_momentum=12.0
    )

    assert np.isclose(inertia, 60.0 * np.pi)
    assert np.isclose(omega, 1.0 / (10.0 * np.pi))

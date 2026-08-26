import numpy as np
import pytest
import sympy as sp

from substrate_framework.nonlinear_clock_gravity import (
    homothetic_compactness_profile,
    integrate_spherical_density_cells,
    minimize_homothetic_max_compactness,
    painleve_gullstrand_mass_constraint,
    static_spherical_mass_constraint,
)


def test_exact_mass_constraint_closes_uniform_sphere_and_sign() -> None:
    r = sp.symbols("r", positive=True)
    density, gravity = sp.symbols("rho G", positive=True)
    mass = 4 * sp.pi * density * r**3 / 3

    result = static_spherical_mass_constraint(r, mass, density, gravity)

    assert result.mass_constraint_residual == 0
    assert result.compactness == 8 * sp.pi * gravity * density * r**2 / 3
    assert result.radial_metric_function == 1 - result.compactness


def test_painleve_gullstrand_data_crosses_horizon_without_metric_singularity() -> None:
    r = sp.symbols("r", positive=True)
    theta = sp.symbols("theta", real=True)
    density, gravity = sp.symbols("rho G", positive=True)
    mass = 4 * sp.pi * density * r**3 / 3

    result = painleve_gullstrand_mass_constraint(r, theta, mass, density, gravity)

    assert result.hamiltonian_constraint_residual == 0
    assert result.momentum_constraint_residual == 0
    assert (
        sp.simplify(result.infall_shift**2 - 8 * sp.pi * gravity * density * r**2 / 3)
        == 0
    )
    assert result.metric_determinant == -(r**4) * sp.sin(theta) ** 2
    horizon_gravity = 3 / (8 * sp.pi * density * r**2)
    horizon_metric = result.spacetime_metric.subs(gravity, horizon_gravity)
    assert sp.simplify(horizon_metric[0, 0]) == 0
    assert sp.simplify(horizon_metric.det() + r**4 * sp.sin(theta) ** 2) == 0


def test_piecewise_constant_uniform_sphere_has_exact_shell_mass() -> None:
    radius = 3.0
    density = 2.0
    gravity = 0.01
    edges = np.linspace(0.0, radius, 101)
    profile = integrate_spherical_density_cells(
        edges, np.full(edges.size - 1, density), gravity
    )

    expected_mass = 4.0 * np.pi * density * radius**3 / 3.0
    assert profile.total_mass == pytest.approx(expected_mass, rel=2.0e-15)
    assert profile.compactness[-1] == pytest.approx(
        2.0 * gravity * expected_mass / radius, rel=2.0e-15
    )
    assert profile.exterior_horizon_radius == pytest.approx(
        2.0 * gravity * expected_mass, rel=2.0e-15
    )


def test_trapped_surface_crossing_converges_to_uniform_sphere_root() -> None:
    density = 1.0
    gravity = 1.0
    expected = np.sqrt(3.0 / (8.0 * np.pi * gravity * density))
    errors = []
    for cells in (64, 128, 256):
        edges = np.linspace(0.0, 1.0, cells + 1)
        profile = integrate_spherical_density_cells(
            edges, np.full(cells, density), gravity
        )
        crossing = profile.first_trapped_surface_radius()
        assert crossing is not None
        errors.append(abs(crossing - expected))
    assert errors[2] < errors[1] < errors[0]
    assert errors[2] < 2.0e-5


def test_zero_density_is_horizonless_and_has_infinite_critical_g() -> None:
    edges = np.linspace(0.0, 2.0, 17)
    profile = integrate_spherical_density_cells(edges, np.zeros(16), 5.0)

    assert profile.total_mass == 0.0
    assert np.array_equal(profile.radial_metric_function, np.ones_like(edges))
    assert profile.first_trapped_surface_radius() is None
    assert profile.critical_newton_constant == float("inf")


def test_homothetic_minimum_matches_single_shell_closed_form() -> None:
    x = np.array([0.0, 1.0])
    curvature = np.array([0.0, 9.0])
    potential = np.array([0.0, 1.0])
    reference = 2.0
    gravity = 0.5

    # C(R)=2G*(M4*R0/R + M0*(R/R0)^3)/R
    #      = a/R^2+b*R^2, minimized at R=(a/b)^(1/4).
    a = 2.0 * gravity * curvature[-1] * reference
    b = 2.0 * gravity * potential[-1] / reference**3
    expected_radius = (a / b) ** 0.25
    expected_minimum = 2.0 * np.sqrt(a * b)

    result = minimize_homothetic_max_compactness(
        x,
        curvature,
        potential,
        reference_radius=reference,
        newton_constant=gravity,
        scale_bounds=(0.1, 100.0),
    )

    assert result.optimizer_success
    assert result.scale_radius == pytest.approx(expected_radius, rel=2.0e-8)
    assert result.maximum_compactness == pytest.approx(expected_minimum, rel=2.0e-12)
    assert result.critical_newton_constant == pytest.approx(
        gravity / expected_minimum, rel=2.0e-12
    )
    assert result.lower_scale_compactness > result.maximum_compactness
    assert result.upper_scale_compactness > result.maximum_compactness


def test_homothetic_reference_scale_recovers_component_sum() -> None:
    x = np.array([0.0, 0.5, 1.0])
    curvature = np.array([0.0, 1.0, 3.0])
    potential = np.array([0.0, 2.0, 4.0])
    radius = 5.0
    gravity = 0.2

    compactness = homothetic_compactness_profile(
        x,
        curvature,
        potential,
        reference_radius=radius,
        scale_radius=radius,
        newton_constant=gravity,
    )

    expected = np.zeros_like(x)
    expected[1:] = 2.0 * gravity * (curvature[1:] + potential[1:]) / (radius * x[1:])
    assert np.allclose(compactness, expected, rtol=0.0, atol=1.0e-15)


def test_wrong_sign_or_negative_density_is_rejected() -> None:
    edges = np.array([0.0, 1.0])
    with pytest.raises(ValueError, match="nonnegative"):
        integrate_spherical_density_cells(edges, np.array([-1.0]), 1.0)
    with pytest.raises(ValueError, match="positive"):
        integrate_spherical_density_cells(edges, np.array([1.0]), -1.0)

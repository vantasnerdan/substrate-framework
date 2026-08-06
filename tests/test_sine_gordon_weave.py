from __future__ import annotations

import pytest
import sympy as sp

import substrate_framework as framework
from substrate_framework.dimensional_sine_gordon import dimensional_sine_gordon_coefficients
from substrate_framework.sine_gordon_weave import (
    CollectiveWeaveMetric,
    SineGordonWeaveGravity,
    TetrahedralWeaveGeometry,
    collective_weave_metric,
    sine_gordon_weave_gravity,
    tetrahedral_weave_geometry,
)


def test_public_package_exports_sine_gordon_weave_api() -> None:
    assert framework.CollectiveWeaveMetric is CollectiveWeaveMetric
    assert framework.TetrahedralWeaveGeometry is TetrahedralWeaveGeometry
    assert framework.SineGordonWeaveGravity is SineGordonWeaveGravity
    assert framework.tetrahedral_weave_geometry is tetrahedral_weave_geometry
    assert framework.collective_weave_metric is collective_weave_metric
    assert framework.sine_gordon_weave_gravity is sine_gordon_weave_gravity


def test_tetrahedral_coordination_closes_and_has_isotropic_crossing_density() -> None:
    length = sp.symbols("ell", positive=True)
    geometry = tetrahedral_weave_geometry(length)
    assert geometry.coordination == 4
    assert geometry.direction_sum == sp.zeros(3, 1)
    assert geometry.direction_second_moment == sp.eye(3) / 3
    assert geometry.shared_line_length_per_cell == 2 * length
    assert geometry.line_length_density == 2 / length**2
    assert geometry.rotational_mean_absolute_projection == sp.Rational(1, 2)
    assert geometry.crossing_density == length**-2


def test_collective_clock_and_tetrahedral_links_reconstruct_metric_and_measure() -> None:
    time_scale, x_scale, y_scale, z_scale = sp.symbols(
        "a_t a_x a_y a_z", positive=True
    )
    metric = collective_weave_metric(
        sp.diag(time_scale, x_scale, y_scale, z_scale)
    )
    assert metric.covariant_metric == sp.diag(
        time_scale**2,
        -x_scale**2,
        -y_scale**2,
        -z_scale**2,
    )
    assert metric.link_second_moment == sp.diag(
        0,
        x_scale**2 / 3,
        y_scale**2 / 3,
        z_scale**2 / 3,
    )
    assert metric.volume_density == time_scale * x_scale * y_scale * z_scale
    assert metric.contravariant_metric * metric.covariant_metric == sp.eye(4)


def test_sine_gordon_scales_fix_entropy_einstein_coupling_and_newton_constant() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    weave = sine_gordon_weave_gravity(coefficients)
    assert weave.geometry.cell_length == sp.Rational(2, 3)
    assert weave.topological_sector_count == 2
    assert weave.entropy_per_crossing == sp.log(2)
    assert weave.entropy_area_density == 9 * sp.log(2) / 4
    assert weave.horizon_equilibrium.action_scale == 4
    assert weave.horizon_equilibrium.signal_speed == 2
    assert weave.horizon_equilibrium.einstein_stress_coupling == sp.pi / (
        9 * sp.log(2)
    )
    assert weave.horizon_equilibrium.newton_constant == 2 / (9 * sp.log(2))
    assert weave.horizon_equilibrium.cosmological_constant == 0
    assert weave.equilibrium_metric.covariant_metric == sp.diag(1, -1, -1, -1)


def test_common_microscopic_multiplier_changes_gravity_but_not_weave_geometry() -> None:
    baseline = sine_gordon_weave_gravity(
        dimensional_sine_gordon_coefficients(2, 8, 18)
    )
    scaled = sine_gordon_weave_gravity(
        dimensional_sine_gordon_coefficients(10, 40, 90)
    )
    assert scaled.scales.signal_speed == baseline.scales.signal_speed
    assert scaled.geometry.crossing_density == baseline.geometry.crossing_density
    assert scaled.scales.action == 5 * baseline.scales.action
    assert scaled.horizon_equilibrium.newton_constant == (
        baseline.horizon_equilibrium.newton_constant / 5
    )


def test_weave_has_no_independent_action_cutoff_or_newton_input() -> None:
    with pytest.raises(TypeError):
        sine_gordon_weave_gravity(
            dimensional_sine_gordon_coefficients(1, 1, 1),
            action_scale=1,
        )
    with pytest.raises(ValueError, match="cell_length"):
        tetrahedral_weave_geometry(sp.Rational(-1, 2))
    with pytest.raises(ValueError, match="4 by 4"):
        collective_weave_metric(sp.eye(3))
    with pytest.raises(ValueError, match="nonzero"):
        collective_weave_metric(sp.zeros(4))
    with pytest.raises(ValueError, match="exact"):
        collective_weave_metric(sp.eye(4) * 1.0)

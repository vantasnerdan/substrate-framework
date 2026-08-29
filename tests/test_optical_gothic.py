from __future__ import annotations

import pytest
import sympy as sp

from substrate_framework.optical_gothic import (
    determinant_slaved_optical_metric,
    newtonian_optical_gradient_ledger,
    optical_adm_metric,
    optical_continuity_compatibility,
    optical_shear_decomposition,
    recover_optical_adm,
)


def test_complete_optical_adm_block_inverse_determinant_and_jacobian() -> None:
    lapse = sp.Rational(3, 2)
    spatial = sp.diag(2, 3, 5)
    flow = sp.Matrix([1, -2, 3])
    speed = sp.Integer(7)
    result = optical_adm_metric(lapse, spatial, flow, speed)

    assert result.covariant * result.contravariant == sp.eye(4)
    assert result.determinant == -lapse**2 * spatial.det()
    assert result.volume_density == lapse * sp.sqrt(spatial.det())
    assert result.component_jacobian_determinant == 2 * lapse * spatial.det()
    assert result.gothic_contravariant == (
        result.volume_density * result.contravariant
    ).applyfunc(sp.simplify)


@pytest.mark.parametrize("orientation", [-1, 1])
def test_complete_optical_adm_round_trip_for_both_shift_signs(orientation: int) -> None:
    original = optical_adm_metric(
        2,
        [[2, 1, 0], [1, 3, 1], [0, 1, 4]],
        [3, -1, 2],
        5,
        flow_orientation=orientation,
    )
    recovered = recover_optical_adm(
        original.covariant,
        5,
        flow_orientation=orientation,
    )
    assert recovered.lapse == 2
    assert recovered.spatial_tensor == original.spatial_tensor
    assert recovered.flow_velocity == sp.Matrix([3, -1, 2])
    assert recovered.reconstructed_metric.covariant == original.covariant


def test_determinant_slaved_map_matches_paper_gothic_components() -> None:
    spatial = sp.diag(2, 4, 8)
    result = determinant_slaved_optical_metric(spatial, [3, -2, 1], 5)
    mean = sp.Integer(4)
    velocity = sp.Matrix([sp.Rational(3, 5), -sp.Rational(2, 5), sp.Rational(1, 5)])

    assert result.determinant_mean == mean
    assert result.metric.lapse == sp.Rational(1, 2)
    assert result.metric.volume_density == mean
    assert result.metric.gothic_contravariant[0, 0] == -mean**2
    assert result.metric.gothic_contravariant[1:4, 0] == mean**2 * velocity
    assert result.metric.gothic_contravariant[1:4, 1:4] == (
        mean * spatial.inv() - mean**2 * velocity * velocity.T
    )
    assert result.lapse_constraint_residual == 0
    assert result.volume_constraint_residual == 0
    assert result.field_count == 9
    assert result.general_metric_component_count == 10


def test_determinant_slaved_map_does_not_cover_an_independent_lapse() -> None:
    general = optical_adm_metric(2, sp.eye(3), [0, 0, 0], 1)
    constrained = determinant_slaved_optical_metric(sp.eye(3), [0, 0, 0], 1)
    assert general.covariant == sp.diag(-4, 1, 1, 1)
    assert constrained.metric.covariant == sp.diag(-1, 1, 1, 1)
    assert general.covariant != constrained.metric.covariant


def test_material_sign_compatibility_is_necessary_and_sufficient() -> None:
    mean, advective, divergence = sp.symbols(
        "nbar advective divergence",
        positive=True,
    )
    partial_time = -advective - mean * divergence
    result = optical_continuity_compatibility(
        mean,
        partial_time,
        advective,
        divergence,
        flow_orientation=-1,
    )
    assert result.material_residual == 0
    assert result.harmonic_residual == -mean**2 * divergence
    assert result.harmonic_from_material_residual == -mean**2 * divergence


def test_inhomogeneous_static_profile_is_material_and_harmonic() -> None:
    mean = sp.symbols("nbar", positive=True)
    result = optical_continuity_compatibility(mean, 0, 0, 0)
    assert result.material_residual == 0
    assert result.harmonic_residual == 0


def test_paper_shift_sign_changes_the_compatibility_condition() -> None:
    mean, advective, divergence = sp.symbols(
        "nbar advective divergence",
        positive=True,
    )
    result = optical_continuity_compatibility(
        mean,
        -advective - mean * divergence,
        advective,
        divergence,
        flow_orientation=1,
    )
    assert result.material_residual == 0
    assert sp.simplify(
        result.harmonic_residual
        + mean * (4 * advective + 3 * mean * divergence)
    ) == 0


def test_paper_gradient_energy_has_the_opposite_newtonian_sign() -> None:
    gradient_squared = sp.symbols("q", positive=True)
    result = newtonian_optical_gradient_ledger(3, 5, gradient_squared)
    assert result.paper_energy_density == -result.required_newtonian_energy_density
    assert result.paper_energy_minus_required == 2 * result.paper_energy_density
    assert result.maxwell_stress_prefactor_in_potential_variables == 1 / (20 * sp.pi)


def test_additive_determinant_mean_deviation_is_not_trace_free() -> None:
    result = optical_shear_decomposition([2, 1, 1])
    assert result.additive_trace != 0
    assert sp.N(result.additive_trace) > 0
    assert result.logarithmic_shear_trace == 0
    assert result.normalized_determinant == 1


def test_log_shear_vanishes_exactly_in_the_isotropic_limit() -> None:
    result = optical_shear_decomposition([3, 3, 3])
    assert result.additive_deviations == (0, 0, 0)
    assert result.logarithmic_shear == (0, 0, 0)
    assert result.normalized_principal_indices == (1, 1, 1)


@pytest.mark.parametrize(
    ("call", "message"),
    [
        (lambda: optical_adm_metric(0, sp.eye(3), [0, 0, 0], 1), "lapse"),
        (lambda: optical_adm_metric(1, sp.diag(1, 1, -1), [0, 0, 0], 1), "positive definite"),
        (lambda: optical_adm_metric(1, sp.eye(3), [0, 0], 1), "shape"),
        (lambda: optical_adm_metric(1, sp.eye(3), [0, 0, 0], 1, flow_orientation=0), "orientation"),
        (lambda: optical_shear_decomposition([1, 1]), "three"),
        (lambda: newtonian_optical_gradient_ledger(1, 1, -1), "nonnegative"),
    ],
)
def test_input_guards(call: object, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        call()

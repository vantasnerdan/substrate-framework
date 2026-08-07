import pytest
import sympy as sp

from substrate_framework.dimensional_sine_gordon import (
    DimensionalSineGordonCoefficients,
    dimensional_sine_gordon_coefficients,
    rescale_dimensional_sine_gordon_coefficients,
)
from substrate_framework.sine_gordon_wall_network import (
    icosahedral_wall_directions,
    icosahedral_wall_frame_ledger,
    local_sine_gordon_wall_network,
    sine_gordon_wall_network_scales,
    spatial_tensor_from_wall_channels,
    wall_channels_from_spatial_tensor,
)


def test_declared_wall_network_is_local_in_all_three_spatial_directions() -> None:
    time, x, y, z = sp.symbols("t x y z", real=True)
    fields = tuple(
        sp.Function(f"u_{index}")(time, x, y, z) for index in range(6)
    )
    network = local_sine_gordon_wall_network(
        fields,
        (x, y, z),
        time,
        dimensional_sine_gordon_coefficients(2, 8, 18),
    )
    assert network.spatial_gradient_hessian == -18 * sp.eye(18)
    assert network.spatial_gradient_hessian.rank() == 18
    assert network.field_equation_residuals.shape == (6, 1)
    assert all(
        network.lagrangian_density.has(sp.diff(field, coordinate) ** 2)
        for field in fields
        for coordinate in (x, y, z)
    )
    missing_transverse_channel = sp.MutableDenseMatrix(network.spatial_gradient_hessian)
    missing_transverse_channel[-1, -1] = 0
    assert missing_transverse_channel.rank() == 17


def test_icosahedral_axes_are_an_exact_rank_six_tight_frame() -> None:
    ledger = icosahedral_wall_frame_ledger()
    assert ledger.projection_rank == 6
    assert ledger.direction_second_moment == 2 * sp.eye(3)
    assert ledger.trace_quadratic_weight == sp.Rational(2, 5)
    assert ledger.traceless_quadratic_weight == sp.Rational(4, 5)
    assert ledger.tight_frame_identity_residual == 0
    assert all(sp.simplify(direction.dot(direction) - 1) == 0 for direction in ledger.directions)


def test_wall_channels_reconstruct_the_symmetric_tensor_without_a_fit() -> None:
    tensor = sp.ImmutableMatrix([[1, 2, 3], [2, 5, 7], [3, 7, 11]])
    channels = wall_channels_from_spatial_tensor(tensor)
    assert spatial_tensor_from_wall_channels(channels) == tensor


def test_repeating_one_axis_breaks_the_reconstruction_rank() -> None:
    directions = icosahedral_wall_directions()
    mutated = directions[:-1] + (directions[0],)
    projection = sp.Matrix(
        [
            [
                n[0] ** 2,
                n[1] ** 2,
                n[2] ** 2,
                2 * n[0] * n[1],
                2 * n[0] * n[2],
                2 * n[1] * n[2],
            ]
            for n in mutated
        ]
    )
    assert projection.rank() == 5


def test_kink_cell_fixes_the_collective_and_gravity_scales() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    ledger = sine_gordon_wall_network_scales(coefficients)
    assert ledger.signal_speed == 2
    assert ledger.profile_length == sp.Rational(2, 3)
    assert ledger.energy_scale == 12
    assert ledger.kink_gradient_integral == 12
    assert ledger.wall_tension == 216
    assert ledger.walls_per_length == sp.Rational(3, 2)
    assert ledger.wall_inertia_density == 36
    assert ledger.wall_modulus == 144
    assert ledger.spin2_time_inertia == sp.Rational(144, 5)
    assert ledger.spin2_spacetime_normalization == sp.Rational(576, 5)


def test_common_sine_gordon_scale_is_visible_and_load_bearing() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    baseline = sine_gordon_wall_network_scales(coefficients)
    scaled = sine_gordon_wall_network_scales(
        rescale_dimensional_sine_gordon_coefficients(coefficients, 7)
    )
    assert scaled.signal_speed == baseline.signal_speed
    assert scaled.profile_length == baseline.profile_length
    assert scaled.spin2_spacetime_normalization == 7 * baseline.spin2_spacetime_normalization


def test_wall_network_rejects_inexact_or_malformed_inputs() -> None:
    with pytest.raises(ValueError, match="exact"):
        wall_channels_from_spatial_tensor([[1.0, 0, 0], [0, 1, 0], [0, 0, 1]])
    with pytest.raises(ValueError, match="symmetric"):
        wall_channels_from_spatial_tensor([[1, 2, 0], [0, 1, 0], [0, 0, 1]])
    with pytest.raises(ValueError, match="six"):
        spatial_tensor_from_wall_channels([1, 2, 3])
    time, x, y, z = sp.symbols("t x y z", real=True)
    with pytest.raises(ValueError, match="six"):
        local_sine_gordon_wall_network(
            [sp.Function("u")(time, x, y, z)],
            (x, y, z),
            time,
            dimensional_sine_gordon_coefficients(1, 1, 1),
        )
    with pytest.raises(ValueError, match="onsite"):
        sine_gordon_wall_network_scales(DimensionalSineGordonCoefficients(1, 1, -1))

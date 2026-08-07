import pytest
import sympy as sp

from substrate_framework.dimensional_sine_gordon import dimensional_sine_gordon_coefficients
from substrate_framework.sine_gordon_spin2_predictions import (
    controlled_static_breather_prediction,
    sine_gordon_spin2_radiation_evidence,
)


def test_static_breather_prediction_has_a_quadratic_profile_error_enclosure() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    near = controlled_static_breather_prediction(sp.Rational(4, 5), coefficients, 10)
    far = controlled_static_breather_prediction(sp.Rational(4, 5), coefficients, 20)
    assert near.source_energy == far.source_energy == sp.Rational(576, 5)
    assert near.source_mass == far.source_mass == sp.Rational(144, 5)
    assert far.maximum_relative_profile_error == near.maximum_relative_profile_error / 4
    assert near.maximum_relative_profile_error < 1
    assert near.exact_potential_lower_bound < near.exact_potential_upper_bound < 0
    assert near.newton_constant == sp.Rational(5, 1152) / sp.pi


def test_static_prediction_rejects_an_uncontrolled_near_profile() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    with pytest.raises(ValueError, match="profile error"):
        controlled_static_breather_prediction(
            sp.Rational(4, 5),
            coefficients,
            sp.Rational(1, 2),
        )


def test_radiation_uses_the_same_derived_coupling_as_the_static_prediction() -> None:
    coefficients = dimensional_sine_gordon_coefficients(2, 8, 18)
    static = controlled_static_breather_prediction(sp.Rational(4, 5), coefficients, 10)
    radiation = sine_gordon_spin2_radiation_evidence(
        sp.Rational(4, 5),
        coefficients,
        spatial_extent_in_inverse_widths=10.0,
        spatial_points=257,
        time_points=33,
        direction_points=129,
        harmonic_count=2,
    )
    assert radiation.newton_constant == static.newton_constant
    assert radiation.signal_speed == 2
    assert radiation.dimensionless_coupling == pytest.approx(
        float(sp.Rational(5, 1024) / sp.pi),
        rel=2.0e-15,
    )
    assert radiation.power > 0
    assert radiation.field_cycle_energy_loss_fraction > 0

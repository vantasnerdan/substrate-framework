import pytest
import sympy as sp

from substrate_framework.euler_measurement_bridge import (
    first_event_probabilities,
    post_reset_transition_probabilities,
    projective_channel_intensities,
)


def test_projective_analyzer_energy_fractions_sum_to_one():
    theta = sp.symbols("theta", real=True)
    state = sp.Matrix([1, 0])
    analyzer = sp.Matrix(
        [[sp.cos(theta), sp.sin(theta)], [-sp.sin(theta), sp.cos(theta)]]
    )
    intensities = projective_channel_intensities(state, analyzer)
    assert intensities == (sp.cos(theta) ** 2, sp.sin(theta) ** 2)
    assert sp.trigsimp(sum(intensities)) == 1


def test_independent_exponential_race_gives_normalized_intensity_law():
    i0, i1, kappa = sp.symbols("i0 i1 kappa", positive=True)
    law = first_event_probabilities((i0, i1), kappa)
    assert law.probabilities == (i0 / (i0 + i1), i1 / (i0 + i1))
    assert law.total_rate == kappa * (i0 + i1)


def test_zero_intensity_channel_never_wins():
    intensity, kappa = sp.symbols("intensity kappa", positive=True)
    law = first_event_probabilities((intensity, 0), kappa)
    assert law.probabilities == (1, 0)


def test_assumed_reset_gives_repeatability_and_rotated_transition_law():
    identity = post_reset_transition_probabilities(sp.eye(2))
    assert identity == sp.eye(2)
    hadamard = sp.sqrt(2) * sp.Matrix([[1, 1], [1, -1]]) / 2
    assert post_reset_transition_probabilities(hadamard) == sp.ones(2) / 2


def test_domain_checks_reject_probability_insertions_hidden_in_bad_inputs():
    with pytest.raises(ValueError, match="normalized"):
        projective_channel_intensities((1, 1), sp.eye(2))
    with pytest.raises(ValueError, match="unitary"):
        projective_channel_intensities((1, 0), ((1, 1), (0, 1)))
    with pytest.raises(ValueError, match="nonnegative"):
        first_event_probabilities((1, -1), 1)
    with pytest.raises(ValueError, match="positive"):
        first_event_probabilities((0, 0), 1)

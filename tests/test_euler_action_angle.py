"""Expose the actual bracket and divergence of the proposed mode lift."""

import pytest
import sympy as sp

from substrate_framework.euler_action_angle import cohomological_kelvin_mode


@pytest.mark.parametrize("mode", [(1, 1), (2, 0)])
def test_curl_input_lifts_with_shear_and_preserves_volume_divergence(mode):
    action = sp.Symbol("action", real=True)
    carrier = sp.Symbol("N", real=True)
    jacobian = 1 + action**2
    omega = sp.Matrix([1 + action, 2 + action**2])
    m, n = mode
    phase = sp.exp(sp.I * carrier * action)
    amplitude = (1 + action) * phase
    # Exterior differentiation supplies a curl coefficient independently of
    # the inverse formula. This tests a variable volume Jacobian and shear.
    alpha1 = jacobian * (2 - action) * amplitude
    alpha2 = -jacobian * (1 + action**2) * amplitude
    w = sp.Matrix([
        (sp.I*m*alpha2 - sp.I*n*alpha1) / jacobian,
        -sp.diff(alpha2, action) / jacobian,
        sp.diff(alpha1, action) / jacobian,
    ])
    xi = cohomological_kelvin_mode(w, omega, mode, action, 3)
    frequency = m*omega[0] + n*omega[1]
    bracket = sp.I*frequency*xi - xi[0]*sp.Matrix([0, *omega.diff(action)])
    assert all(sp.simplify(value) == 0 for value in 3*bracket-w)
    divergence = (
        sp.diff(jacobian*xi[0], action)
        + jacobian*sp.I*(m*xi[1]+n*xi[2])
    ) / jacobian
    assert sp.simplify(divergence) == 0
    uncorrected = w / (3*sp.I*frequency)
    wrong_bracket = sp.I*frequency*uncorrected - uncorrected[0]*sp.Matrix(
        [0, *omega.diff(action)]
    )
    assert any(sp.simplify(value) != 0 for value in 3*wrong_bracket-w)


def test_resonance_and_nonconstant_curl_are_exposed():
    action = sp.Symbol("action", real=True)
    with pytest.raises(ValueError, match="nonzero"):
        cohomological_kelvin_mode([1, 2, 3], [1, 2], [2, -1], action, 1)
    with pytest.raises(ValueError, match="constant curl"):
        cohomological_kelvin_mode([1, 2, 3], [1, 2], [1, 1], action, 1+action)

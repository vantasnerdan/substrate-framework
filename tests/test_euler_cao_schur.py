import pytest
import sympy as sp

from substrate_framework.euler_cao_schur import cao_thin_ring_schur_jet


def test_cao_thin_ring_parameter_and_physical_schur_jet():
    kappa, radius, core_log, density = sp.symbols(
        "kappa radius core_log density", positive=True
    )
    jet = cao_thin_ring_schur_jet(kappa, radius, core_log, density)

    assert jet.chemical_potential == 3 * kappa * radius * core_log / (8 * sp.pi)
    assert jet.translation_speed == kappa * core_log / (4 * sp.pi * radius)
    assert jet.physical_axial_impulse == sp.pi * density * kappa * radius**2
    assert jet.parameter_jacobian == -3 * kappa * core_log**2 / (
        16 * sp.pi**2 * radius
    )
    assert jet.moment_jacobian == 2 * sp.pi * density * kappa * radius
    assert jet.circulation_response_at_fixed_speed == 4 * sp.pi / (
        3 * radius * core_log
    )
    assert jet.impulse_response_at_fixed_circulation == (
        -8 * sp.pi**2 * density * radius**3 / core_log
    )
    assert sp.simplify(
        jet.circulation_response_at_fixed_speed
        * jet.impulse_response_at_fixed_circulation
        - jet.physical_schur_determinant
    ) == 0
    assert jet.physical_schur_determinant == -32 * sp.pi**3 * density * radius**2 / (
        3 * core_log**2
    )


def test_cao_thin_ring_schur_rejects_nonpositive_physical_inputs():
    with pytest.raises(ValueError, match="circulation"):
        cao_thin_ring_schur_jet(0, 1, 2, 3)
    with pytest.raises(ValueError, match="radius"):
        cao_thin_ring_schur_jet(1, -1, 2, 3)
    with pytest.raises(ValueError, match="log_inverse_core"):
        cao_thin_ring_schur_jet(1, 1, 0, 3)
    with pytest.raises(ValueError, match="density"):
        cao_thin_ring_schur_jet(1, 1, 2, -3)

import pytest
import sympy as sp

from substrate_framework.euler_maxwell_radiation_flux import (
    doppler_shell_geometry,
    gaussian_curl_current_power,
    outgoing_power_prefactor,
    shell_sphere_power_weight,
    switched_continuity_residual,
)


def test_outgoing_prefactor_and_shell_coarea_constants():
    epsilon, c_em, speed, omega, n_z = sp.symbols(
        "epsilon c_em speed omega n_z", positive=True
    )
    geometry = doppler_shell_geometry(c_em, speed, omega, n_z)
    assert geometry.radius == omega / (c_em - speed * n_z)
    assert geometry.temporal_frequency == c_em * geometry.radius
    assert geometry.radial_derivative_magnitude == 2 * c_em * omega

    direct = (
        outgoing_power_prefactor(epsilon)
        * geometry.temporal_frequency
        * geometry.radius**2
        / geometry.radial_derivative_magnitude
    )
    assert sp.simplify(
        direct - shell_sphere_power_weight(epsilon, c_em, speed, omega, n_z)
    ) == 0


def test_full_shell_gradient_is_nonzero_with_subluminal_margin():
    c_em, speed, omega, n_z = sp.symbols(
        "c_em speed omega n_z", positive=True
    )
    geometry = doppler_shell_geometry(c_em, speed, omega, n_z)
    lower = 4 * c_em**2 * geometry.radius**2 * (c_em - speed) ** 2
    remainder = 8 * c_em**3 * geometry.radius**2 * speed * (1 - n_z)
    assert sp.simplify(geometry.gradient_squared - lower - remainder) == 0


def test_envelope_current_completion_cancels_continuity_defect():
    a, da, rho = sp.symbols("a da rho")
    assert switched_continuity_residual(a, da, 0, rho, -rho) == 0


def test_negative_frequency_and_speed_select_conjugate_shell_once():
    c_em, speed, omega = sp.symbols("c_em speed omega", positive=True)
    n_z = -sp.Rational(1, 3)
    positive = doppler_shell_geometry(c_em, -speed, omega, n_z)
    conjugate = doppler_shell_geometry(c_em, -speed, -omega, -n_z)
    assert sp.simplify(positive.radius - conjugate.radius) == 0
    assert sp.simplify(
        positive.temporal_frequency + conjugate.temporal_frequency
    ) == 0
    assert sp.simplify(
        shell_sphere_power_weight(1, c_em, -speed, omega, n_z)
        - shell_sphere_power_weight(1, c_em, -speed, -omega, -n_z)
    ) == 0


def test_gaussian_curl_current_power_is_exact_spherical_integral():
    epsilon, c_em, omega, width, j_sq, a_sq = sp.symbols(
        "epsilon c_em omega width j_sq a_sq", positive=True
    )
    radius = omega / c_em
    angular_transverse_square = sp.Rational(8, 3) * sp.pi * a_sq
    direct = (
        sp.pi
        * omega**2
        / (4 * epsilon * c_em**3)
        * j_sq
        * radius**2
        * sp.exp(-(width * radius) ** 2)
        * angular_transverse_square
    )
    assert sp.simplify(
        direct
        - gaussian_curl_current_power(
            epsilon, c_em, omega, width, j_sq, a_sq
        )
    ) == 0


def test_radiation_flux_helpers_reject_out_of_domain_inputs():
    with pytest.raises(ValueError, match="nonzero"):
        doppler_shell_geometry(2, 1, 0, 0)
    with pytest.raises(ValueError, match="strictly subluminal"):
        doppler_shell_geometry(1, 1, 2, 0)
    with pytest.raises(ValueError, match="strictly subluminal"):
        shell_sphere_power_weight(1, 1, 2, 3, 0)
    with pytest.raises(ValueError, match="unit interval"):
        doppler_shell_geometry(2, 1, 3, 2)
    with pytest.raises(ValueError, match="positive"):
        outgoing_power_prefactor(-1)

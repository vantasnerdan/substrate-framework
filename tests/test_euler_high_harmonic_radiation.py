import pytest
import sympy as sp

from substrate_framework.euler_high_harmonic_radiation import (
    cylindrical_vector_bessel_orders,
    debye_exponential_rate,
    debye_speed_ceiling,
    fixed_mode_bessel_ratio,
    high_index_limiting_bessel_ratio,
    physical_mode_ledger,
    transverse_shell_maximum,
)


def test_physical_clock_retains_circulation_and_two_pi():
    kappa, radius, wave_number, harmonic, sigma = sp.symbols(
        "kappa radius wave_number harmonic sigma", positive=True
    )
    ledger = physical_mode_ledger(
        kappa, radius, wave_number, harmonic, sigma
    )
    assert ledger.core_radius == radius * wave_number / harmonic
    assert ledger.core_clock == (
        kappa * harmonic**2 / (2 * sp.pi * radius**2 * wave_number**2)
    )
    assert ledger.physical_frequency == ledger.core_clock * sigma


def test_transverse_shell_maximum_is_sharper_than_radial_bound():
    c_em, speed, omega = sp.symbols("c_em speed omega", positive=True)
    shell = transverse_shell_maximum(c_em, speed, omega, 1)
    t = shell.maximizing_direction_z
    k_perp = omega * sp.sqrt(1 - t**2) / (c_em - speed * t)
    assert sp.simplify(k_perp - shell.transverse_wave_number) == 0
    assert sp.simplify(
        shell.transverse_wave_number - omega / sp.sqrt(c_em**2 - speed**2)
    ) == 0


def test_fixed_mode_ratio_has_core_support_and_n_minus_one_order():
    kappa, sigma, n, k, radius, c_em, speed, support = sp.symbols(
        "kappa sigma n k radius c_em speed support", positive=True
    )
    ratio = fixed_mode_bessel_ratio(
        kappa, sigma, n, k, radius, c_em, speed, support
    )
    expected = (
        kappa
        * sigma
        * n
        * (n + support * k)
        / (2 * sp.pi * radius * k**2 * sp.sqrt(c_em**2 - speed**2) * (n - 1))
    )
    assert sp.simplify(ratio - expected) == 0


def test_finite_j_leading_predictor_at_supplied_speed_and_ceiling():
    kappa, ell, l_phi, radius, k, c_em, speed = sp.symbols(
        "kappa ell l_phi radius k c_em speed", positive=True
    )
    eta = sp.Rational(1, 4)
    ratio = high_index_limiting_bessel_ratio(
        kappa, ell, l_phi, radius, k, c_em, speed
    )
    c_zero = kappa * ell * l_phi / (2 * sp.pi**2 * radius * k)
    assert sp.simplify(ratio - c_zero / sp.sqrt(c_em**2 - speed**2)) == 0
    ceiling = debye_speed_ceiling(c_em, c_zero, eta)
    assert sp.simplify(
        c_zero / sp.sqrt(c_em**2 - ceiling**2) - (1 - eta)
    ) == 0


def test_debye_rate_reduces_on_hyperbolic_parameterization():
    q = sp.Rational(4, 5)
    rate = debye_exponential_rate(q)
    assert sp.simplify(rate - (sp.acosh(sp.Rational(5, 4)) - sp.Rational(3, 5))) == 0


def test_cylindrical_frame_has_only_neighboring_bessel_orders():
    n = sp.symbols("n", integer=True)
    assert cylindrical_vector_bessel_orders(n) == (n - 1, n, n + 1)


def test_radiation_helpers_reject_invalid_domains():
    with pytest.raises(ValueError, match="harmonic"):
        fixed_mode_bessel_ratio(1, 1, 1, 1, 1, 2, 0)
    with pytest.raises(ValueError, match="frequency_sign"):
        transverse_shell_maximum(2, 1, 1, 0)
    with pytest.raises(ValueError, match="strictly subluminal"):
        transverse_shell_maximum(1, 1, 1)
    with pytest.raises(ValueError, match="strictly subluminal"):
        transverse_shell_maximum(1, 2, 1)
    with pytest.raises(ValueError, match="strictly subluminal"):
        high_index_limiting_bessel_ratio(1, 1, 1, 1, 1, 1, 2)
    with pytest.raises(ValueError, match="support_constant"):
        fixed_mode_bessel_ratio(1, 1, 2, 1, 1, 2, 0, -1)
    with pytest.raises(ValueError, match="strictly between"):
        debye_speed_ceiling(2, 1, 0)
    with pytest.raises(ValueError, match="impossible"):
        debye_speed_ceiling(2, sp.Rational(3, 2), sp.Rational(1, 4))
    with pytest.raises(ValueError, match="impossible"):
        debye_speed_ceiling(2, 2, sp.Rational(1, 4))
    with pytest.raises(ValueError, match="strictly between"):
        debye_exponential_rate(1)

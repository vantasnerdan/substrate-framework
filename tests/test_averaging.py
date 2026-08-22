"""Tests for exact coarse-graining machinery."""

from __future__ import annotations

import pytest
import sympy as sp

from substrate_framework.averaging import (
    commutation_residual,
    filter_direct_tophat,
    filter_polynomial,
    filtered_balance_residual,
    kernel_even_moments,
    leonard_expansion_residual,
    barotropic_balance_residual,
    microscopic_balance_residual,
)


def test_tophat_moments_match_direct_integrals() -> None:
    moments = kernel_even_moments("tophat", 2)
    assert moments[0] == 1
    assert moments[1] == sp.Rational(1, 12)
    assert moments[2] == sp.Rational(1, 80)


def test_series_equals_direct_convolution() -> None:
    x, width = sp.symbols("x Delta", positive=True)
    for polynomial in (x**2, x**3, x**4):
        series = filter_polynomial(polynomial, x, width)
        direct = filter_direct_tophat(polynomial, x, width)
        assert sp.simplify(series - direct) == 0


def test_commutator_vanishes_on_polynomial_basis() -> None:
    x, width = sp.symbols("x Delta", positive=True)
    for degree in range(1, 5):
        assert commutation_residual(x**degree, x, width) == 0
    assert commutation_residual(sp.Integer(1), x, width) == 0


def test_commutator_vanishes_on_polynomial_basis() -> None:
    x, width = sp.symbols("x Delta", positive=True)
    for degree in range(5):
        assert commutation_residual(x**degree, x, width) == 0


def test_filtered_balance_closes_on_rigid_rotation() -> None:
    x, y, t, width = sp.symbols("x y t Delta", positive=True)
    omega = sp.Symbol("omega_0", positive=True)
    rho = sp.Symbol("rho_0", positive=True)
    velocity = (omega * y, -omega * x)
    pressure = rho * omega**2 * (x**2 + y**2) / 2
    micro = microscopic_balance_residual(
        velocity, pressure, (0, 0), rho, (x, y), t
    )
    assert all(sp.simplify(component) == 0 for component in micro)
    filtered = filtered_balance_residual(
        velocity, pressure, (0, 0), rho, (x, y), t, width
    )
    assert all(sp.simplify(component) == 0 for component in filtered)


def test_subfilter_flux_quadratic_in_amplitude() -> None:
    x, y, t, width = sp.symbols("x y t Delta", positive=True)
    omega = sp.Symbol("omega_0", positive=True)
    from substrate_framework.averaging import subfilter_flux

    flux = subfilter_flux((omega * y, -omega * x), (x, y), t, width)
    for entry in (flux[0, 0], flux[1, 1]):
        assert sp.simplify(entry - omega**2 * width**2 / 12) == 0
    assert sp.simplify(flux[0, 1]) == 0


def test_leonard_expansion_exact_beyond_degree() -> None:
    x, width = sp.symbols("x Delta", positive=True)
    assert leonard_expansion_residual(x**3, x, width, 2) == 0


def test_microscopic_balance_rejects_non_euler_field() -> None:
    x, y, t = sp.symbols("x y t", positive=True)
    with pytest.raises(ValueError):
        microscopic_balance_residual(
            (x, y), sp.Integer(0), (0, 0), 1.0, (x, y), t
        )


def test_barotropic_balance_closes_on_diluting_flow() -> None:
    x, _y, t = sp.symbols("x y t", real=True)
    a, rho0, k_eos = sp.symbols("a rho0 k_eos", positive=True)
    velocity = (a * x, 0)
    density = rho0 * sp.exp(-a * t)
    pressure = k_eos * density
    body = (a**2 * x, 0)
    residual = barotropic_balance_residual(
        velocity, density, pressure, body, (x, _y), t
    )
    assert all(sp.simplify(r) == 0 for r in residual)


def test_barotropic_balance_flags_wrong_eos_term() -> None:
    x, _y, t = sp.symbols("x y t", real=True)
    a, rho0, k_eos = sp.symbols("a rho0 k_eos", positive=True)
    velocity = (a * x, 0)
    density = rho0 * sp.exp(-a * t)
    body = (sp.Integer(0), 0)
    residual = barotropic_balance_residual(
        velocity, density, k_eos * density, body, (x, _y), t
    )
    assert sp.simplify(residual[0] - a**2 * x) == 0
    assert sp.simplify(residual[2]) == 0

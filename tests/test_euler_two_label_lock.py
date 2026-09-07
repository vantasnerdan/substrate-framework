import pytest
import sympy as sp

from substrate_framework.euler_two_label_lock import (
    annular_dipole_matching_shift,
    ampere_eliminated_lock_response,
    ampere_first_speed_magnetic_response,
    cao_constant_ratio_residual,
    comoving_lorenz_scalar_source,
    constant_ratio_residual,
    fixed_profile_maxwell_obstruction,
    exponential_column_profile,
    two_label_forced_lock_residual,
    weighted_contour_zero_mean,
)


def test_two_label_forced_lock_uses_phase_gradient():
    chi, dt_lambda = sp.symbols("chi dt_lambda")
    curl_f = sp.Matrix(sp.symbols("c0:3"))
    grad_theta = sp.Matrix(sp.symbols("t0:3"))
    assert two_label_forced_lock_residual(
        chi, dt_lambda, curl_f, grad_theta
    ) == chi * dt_lambda - curl_f.dot(grad_theta)


def test_undivided_constant_lock_exposes_zero_tag_failure():
    q, lam = sp.symbols("q lam", nonzero=True)
    assert constant_ratio_residual(q, 0, lam) == q


def test_cao_constant_ratio_residual_matches_first_integral():
    h, g, rho, chi, chip, phi, H, lam = sp.symbols(
        "h g rho chi chip phi H lam"
    )
    got = cao_constant_ratio_residual(h, g, rho, chi, chip, phi, H, lam)
    assert sp.expand(got - (h + g * (chip * phi - chi * H) / rho - lam * chi)) == 0


def test_weighted_contour_projector_annihilates_constants_and_is_idempotent():
    a, b, c = sp.symbols("a b c")
    weights = (1, 2, 3)
    assert weighted_contour_zero_mean((c, c, c), weights) == (0, 0, 0)
    first = weighted_contour_zero_mean((a, b, c), weights)
    second = weighted_contour_zero_mean(first, weights)
    assert all(sp.simplify(x - y) == 0 for x, y in zip(first, second, strict=True))
    assert sp.simplify(sum(w * x for w, x in zip(weights, first, strict=True))) == 0


def test_fixed_profile_obstruction_retains_only_contour_variation():
    chi, chip, phi0, h0, a, b = sp.symbols("chi chip phi0 h0 a b")
    got = fixed_profile_maxwell_obstruction(
        chi,
        chip,
        (phi0 + a, phi0 - a),
        (h0 + b, h0 - b),
        (1, 1),
    )
    expected = chip * a - chi * b
    assert tuple(map(sp.simplify, got)) == (expected, -expected)


def test_contour_projection_rejects_bad_domains():
    with pytest.raises(ValueError):
        weighted_contour_zero_mean((), ())
    with pytest.raises(ValueError):
        weighted_contour_zero_mean((1,), (1, 2))
    with pytest.raises(ValueError):
        weighted_contour_zero_mean((1, 2), (1, -1))


def test_signed_charge_factored_lorenz_source_has_physical_prefactor():
    chi, eps, c, c_em, w_z = sp.symbols(
        "chi eps c c_em w_z", nonzero=True
    )
    got = comoving_lorenz_scalar_source(chi, eps, c, c_em, w_z)
    want = chi * (1 - c**2 / c_em**2 - c * w_z / c_em**2) / eps
    assert sp.simplify(got - want) == 0


def test_ampere_first_speed_response_uses_same_em_normalization():
    phi_r, c_em, radius = sp.symbols("phi_r c_em radius", nonzero=True)
    assert ampere_first_speed_magnetic_response(
        phi_r, c_em, radius
    ) == -phi_r / (c_em**2 * radius)


def test_maxwell_response_helpers_reject_degenerate_domains():
    with pytest.raises(ValueError):
        comoving_lorenz_scalar_source(1, 0, 1, 1, 1)
    with pytest.raises(ValueError):
        comoving_lorenz_scalar_source(1, 1, 1, 0, 1)
    with pytest.raises(ValueError):
        ampere_first_speed_magnetic_response(1, 0, 1)
    with pytest.raises(ValueError):
        ampere_first_speed_magnetic_response(1, 1, 0)


def test_ampere_eliminated_lock_response_keeps_curved_radius_row():
    chi, chip, phi, phi_r, primitive, eps, c, c_em, radius = sp.symbols(
        "chi chip phi phi_r primitive eps c c_em radius", nonzero=True
    )
    a = 1 - c**2 / c_em**2
    got = ampere_eliminated_lock_response(
        chi, chip, phi, phi_r, primitive, eps, c, c_em, radius
    )
    want = (
        chip * phi
        - chi * primitive / (eps * a * c_em**2 * radius**2)
        + chi * c * phi_r / (a * c_em**2 * radius)
    )
    assert sp.simplify(got - want) == 0


def test_ampere_elimination_rejects_characteristic_and_zero_scale_rows():
    with pytest.raises(ValueError):
        ampere_eliminated_lock_response(1, 1, 1, 1, 1, 0, 0, 1, 1)
    with pytest.raises(ValueError):
        ampere_eliminated_lock_response(1, 1, 1, 1, 1, 1, 1, 1, 1)
    with pytest.raises(ValueError):
        ampere_eliminated_lock_response(1, 1, 1, 1, 1, 1, 2, 1, 1)


def test_annular_dipole_inner_outer_moments_control_independent_rows():
    inner, outer, radius = sp.symbols("inner outer radius", nonzero=True)
    a_shift, b_shift = annular_dipole_matching_shift(inner, outer)
    response = a_shift * radius + b_shift / radius
    assert sp.diff(response, inner) == -1 / (2 * radius)
    assert sp.diff(response, outer) == -radius / 2
    assert sp.det(sp.Matrix([[0, -sp.Rational(1, 2)], [-sp.Rational(1, 2), 0]])) != 0


def test_explicit_exponential_column_solves_poisson_and_profile_law():
    s, p0, exponent, scale, radius = sp.symbols(
        "s p0 exponent scale radius", positive=True
    )
    profile, zeta = exponential_column_profile(
        s, p0, exponent, scale, radius
    )
    radial_laplacian = sp.diff(profile, s, 2) + sp.diff(profile, s) / s
    assert sp.simplify(-radial_laplacian - radius**2 * zeta) == 0
    assert sp.simplify(sp.diff(zeta, s) / sp.diff(profile, s) - exponent * zeta) == 0


def test_explicit_exponential_column_rejects_bad_scales():
    with pytest.raises(ValueError):
        exponential_column_profile(1, 1, 0, 1, 1)
    with pytest.raises(ValueError):
        exponential_column_profile(1, 1, 1, 0, 1)
    with pytest.raises(ValueError):
        exponential_column_profile(1, 1, 1, 1, 0)

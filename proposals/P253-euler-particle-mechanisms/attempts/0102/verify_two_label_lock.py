#!/usr/bin/env python3
"""Exact exposing checks for the P253/0102 two-label lock derivation."""

from __future__ import annotations

import sympy as sp

from substrate_framework.euler_two_label_lock import (
    annular_dipole_matching_shift,
    ampere_eliminated_lock_response,
    ampere_first_speed_magnetic_response,
    comoving_lorenz_scalar_source,
    constant_ratio_residual,
    fixed_profile_maxwell_obstruction,
    exponential_column_profile,
    two_label_forced_lock_residual,
    weighted_contour_zero_mean,
)


def main() -> None:
    r, z, g, rho = sp.symbols("r z g rho", nonzero=True, real=True)
    P = sp.Function("P")(r, z)
    Phi = sp.Function("Phi")(r, z)
    H = sp.Function("H")(r, z)
    chi_fun = sp.Function("chi")
    chi = chi_fun(P)

    f_r = -g * chi * (sp.diff(Phi, r) + H * sp.diff(P, r)) / rho
    f_z = -g * chi * (sp.diff(Phi, z) + H * sp.diff(P, z)) / rho
    curl_f_theta_over_r = sp.expand((sp.diff(f_r, z) - sp.diff(f_z, r)) / r)
    S = sp.diff(chi_fun(sp.Symbol("s")), sp.Symbol("s")).subs(
        sp.Symbol("s"), P
    ) * Phi - chi * H
    W_r = -sp.diff(P, z) / r
    W_z = sp.diff(P, r) / r
    rhs = g * (W_r * sp.diff(S, r) + W_z * sp.diff(S, z)) / rho
    assert sp.simplify(curl_f_theta_over_r - rhs) == 0

    tag, dt_lam = sp.symbols("tag dt_lam")
    curl_f = sp.Matrix(sp.symbols("c0:3"))
    grad_theta = sp.Matrix(sp.symbols("t0:3"))
    assert two_label_forced_lock_residual(
        tag, dt_lam, curl_f, grad_theta
    ) == tag * dt_lam - curl_f.dot(grad_theta)

    zeta, lam0 = sp.symbols("zeta lam0", nonzero=True)
    assert constant_ratio_residual(zeta, 0, lam0) == zeta

    x, y, w1, w2 = sp.symbols("x y w1 w2", nonzero=True)
    projected = weighted_contour_zero_mean((x, y), (w1, w2))
    assert sp.simplify(w1 * projected[0] + w2 * projected[1]) == 0
    assert all(
        sp.simplify(a - b) == 0
        for a, b in zip(
            weighted_contour_zero_mean(projected, (w1, w2)),
            projected,
            strict=True,
        )
    )

    chip, phi0, h0, dp, dh = sp.symbols("chip phi0 h0 dp dh")
    obstruction = fixed_profile_maxwell_obstruction(
        tag,
        chip,
        (phi0 + dp, phi0 - dp),
        (h0 + dh, h0 - dh),
        (1, 1),
    )
    assert obstruction == (chip * dp - tag * dh, -chip * dp + tag * dh)

    s, a = sp.symbols("s a", positive=True)
    radial = sp.Function("radial")(s)
    axis_difference = (1 - a) * (sp.diff(radial, s) / s - sp.diff(radial, s, 2))
    quadratic = sp.Symbol("C") * s**2 / 2 + sp.Symbol("D")
    assert sp.simplify(axis_difference.subs(radial, quadratic)) == 0

    eps_em, c, c_em, w_z = sp.symbols(
        "eps_em c c_em w_z", nonzero=True
    )
    scalar_source = comoving_lorenz_scalar_source(
        tag, eps_em, c, c_em, w_z
    )
    assert sp.simplify(
        scalar_source
        - tag * (1 - c**2 / c_em**2 - c * w_z / c_em**2) / eps_em
    ) == 0

    phi0_r, radius = sp.symbols("phi0_r radius", nonzero=True)
    assert ampere_first_speed_magnetic_response(
        phi0_r, c_em, radius
    ) == -phi0_r / (c_em**2 * radius)

    primitive = sp.symbols("primitive")
    eliminated = ampere_eliminated_lock_response(
        tag, chip, phi0, phi0_r, primitive, eps_em, c, c_em, radius
    )
    ellipticity = 1 - c**2 / c_em**2
    eliminated_expected = (
        chip * phi0
        - tag * primitive / (eps_em * ellipticity * c_em**2 * radius**2)
        + tag * c * phi0_r / (ellipticity * c_em**2 * radius)
    )
    assert sp.simplify(eliminated - eliminated_expected) == 0

    # Verify the exact m=1 radial Green formula without choosing a profile.
    # I'=s^2 f and J'=-f encode its two Volterra moments.
    f = sp.Function("f")(s)
    I = sp.Function("I")(s)
    J = sp.Function("J")(s)
    dipole = -(I / s + s * J) / 2
    dipole_operator = (
        sp.diff(dipole, s, 2)
        + sp.diff(dipole, s) / s
        - dipole / s**2
    )
    dipole_operator = dipole_operator.subs(
        {sp.diff(I, s, 2): sp.diff(s**2 * f, s), sp.diff(J, s, 2): -sp.diff(f, s)}
    ).subs({sp.diff(I, s): s**2 * f, sp.diff(J, s): -f})
    assert sp.simplify(dipole_operator - f) == 0

    # Global straight-column exponential profile: zeta_P=a0*zeta makes
    # the electric dipole a multiple of P' and cancels the Ampere term.
    zeta0, lambda0, a0, p_r = sp.symbols(
        "zeta0 lambda0 a0 p_r", nonzero=True
    )
    chi0 = zeta0 / lambda0
    chi0_p = a0 * zeta0 / lambda0
    phi0_r_exact = p_r / (lambda0 * eps_em * radius**2)
    v_exact = -p_r / (
        lambda0 * eps_em * c_em**2 * radius**3 * a0
    )
    h1_exact = ampere_first_speed_magnetic_response(
        phi0_r_exact, c_em, radius
    )
    assert sp.simplify(chi0_p * v_exact - chi0 * h1_exact) == 0

    p0, radial_scale = sp.symbols("p0 radial_scale", positive=True)
    liouville_p, liouville_zeta = exponential_column_profile(
        s, p0, a0, radial_scale, radius
    )
    liouville_laplacian = (
        sp.diff(liouville_p, s, 2) + sp.diff(liouville_p, s) / s
    )
    assert sp.simplify(-liouville_laplacian - radius**2 * liouville_zeta) == 0
    assert sp.simplify(
        sp.diff(liouville_zeta, s) / sp.diff(liouville_p, s)
        - a0 * liouville_zeta
    ) == 0

    inner_moment, outer_moment = sp.symbols(
        "inner_moment outer_moment", nonzero=True
    )
    a_shift, b_shift = annular_dipole_matching_shift(
        inner_moment, outer_moment
    )
    assert a_shift == -outer_moment / 2
    assert b_shift == -inner_moment / 2
    assert sp.det(
        sp.Matrix([[0, -sp.Rational(1, 2)], [-sp.Rational(1, 2), 0]])
    ) == -sp.Rational(1, 4)

    checks = (
        "charged Cao cylindrical forced-lock identity",
        "two-label forced quotient residual",
        "zero-tag undivided constant-lock obstruction",
        "weighted contour projection and idempotence",
        "fixed-profile Maxwell zero-mean response",
        "anisotropic radial-potential obstruction",
        "signed-charge-factored Lorenz source normalization",
        "Ampere first-speed magnetic normalization",
        "exact Ampere-eliminated lock response",
        "radial dipole Green identity",
        "global exponential-profile first-speed cancellation",
        "explicit Liouville exponential-column realization",
        "independent inner/outer annular dipole matching rows",
    )
    for index, label in enumerate(checks, 1):
        print(f"PASS {index}: {label}")
    print(f"PASS: {len(checks)}/{len(checks)} exact checks")


if __name__ == "__main__":
    main()

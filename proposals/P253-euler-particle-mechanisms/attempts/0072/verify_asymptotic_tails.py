"""Exact exposing checks for P253/0072."""

from __future__ import annotations

import sympy as sp

from substrate_framework.euler_asymptotic_tails import (
    fixed_frame_constant_norm_symbol,
    fixed_frame_angular_stress,
    gaussian_tail_cross_potential,
    l1_homogeneous_tail,
    l1_multiplicity_cross_block,
    oriented_tail_cross_kernel,
    radial_tail_fourier_coefficient,
    scalar_homogeneous_fourier_coefficient,
    steady_curl_residual,
    toroidal_tail_fourier_coefficient,
)


def main() -> None:
    assert scalar_homogeneous_fourier_coefficient(0, 1) == 4 * sp.pi
    assert toroidal_tail_fourier_coefficient(1) == -4 * sp.I * sp.pi
    assert radial_tail_fourier_coefficient(1) == sp.pi**2

    x, y, z = sp.symbols("x y z", real=True)
    ax, ay, az = sp.symbols("a_x a_y a_z", real=True)
    r2 = x**2 + y**2 + z**2
    avec = sp.Matrix([ax, ay, az])
    xvec = sp.Matrix([x, y, z])
    u = avec.cross(xvec) / (1 + r2) ** sp.Rational(3, 2)
    div_u = sum(sp.diff(u[i], (x, y, z)[i]) for i in range(3))
    assert sp.simplify(div_u) == 0

    curl_u = sp.Matrix(
        [
            sp.diff(u[2], y) - sp.diff(u[1], z),
            sp.diff(u[0], z) - sp.diff(u[2], x),
            sp.diff(u[1], x) - sp.diff(u[0], y),
        ]
    )
    expected = ((2 - r2) * avec + 3 * xvec * avec.dot(xvec)) / (
        1 + r2
    ) ** sp.Rational(5, 2)
    assert all(sp.simplify(curl_u[i] - expected[i]) == 0 for i in range(3))
    assert sp.simplify(u.dot(curl_u)) == 0

    radial_integral = sp.integrate(
        sp.Symbol("r", positive=True) ** 4
        / (1 + sp.Symbol("r", positive=True) ** 2) ** 3,
        (sp.Symbol("r", positive=True), 0, sp.oo),
    )
    assert radial_integral == 3 * sp.pi / 16
    # Angular integral of |a cross n|^2 is 8*pi*|a|^2/3.
    assert sp.simplify(radial_integral * 8 * sp.pi / 3) == sp.pi**2 / 2

    rho, d = sp.symbols("rho d", positive=True)
    parallel = oriented_tail_cross_kernel(
        (1, 0, 0), (1, 0, 0), (1, 0, 0), density=rho, separation=d
    )
    transverse = oriented_tail_cross_kernel(
        (1, 0, 0), (1, 0, 0), (0, 1, 0), density=rho, separation=d
    )
    assert parallel == 4 * sp.pi * rho / d
    assert transverse == 2 * sp.pi * rho / d
    assert parallel > transverse

    block = l1_multiplicity_cross_block(
        (1, 0, 0),
        (1, 0, 0),
        (0, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        density=rho,
        separation=d,
    )
    assert block[0, 1] == block[1, 0] == -2 * sp.pi * rho / d
    assert sp.simplify(sp.I / scalar_homogeneous_fourier_coefficient(1, 1)) == (
        -1 / (2 * sp.pi**2)
    )

    n1, n2, n3 = sp.symbols("n1 n2 n3", real=True)
    nvec = sp.ImmutableMatrix([n1, n2, n3])
    fframe = fixed_frame_constant_norm_symbol(nvec)
    sphere = {n1**2: 1 - n2**2 - n3**2}
    assert sp.simplify(nvec.dot(fframe).subs(sphere)) == 0
    norm = sum(sp.conjugate(fframe[j]) * fframe[j] for j in range(3))
    assert sp.simplify(sp.expand(norm - 1).subs(sphere)) == 0
    fminus = fixed_frame_constant_norm_symbol((-n1, -n2, -n3))
    assert all(sp.simplify(fminus[j] - sp.conjugate(fframe[j])) == 0 for j in range(3))
    rz_pi = sp.diag(-1, -1, 1)
    frot = rz_pi * fixed_frame_constant_norm_symbol(rz_pi.T * nvec)
    assert all(sp.simplify(frot[j] + fframe[j]) == 0 for j in range(3))

    sigma, q1, q2 = sp.symbols("sigma q1 q2", positive=True)
    regularized = gaussian_tail_cross_potential(q1, q2, d, sigma, density=rho)
    assert sp.limit(regularized, d, 0, dir="+") == rho * q1 * q2 / (
        4 * sp.pi ** sp.Rational(3, 2) * sigma
    )
    assert sp.limit(d * regularized, d, sp.oo) == rho * q1 * q2 / (4 * sp.pi)

    A, B, C = sp.symbols("A B C", real=True)
    l1_velocity = l1_homogeneous_tail((x, y, z), (A, 0, C), (0, 0, B))
    stationary_residual = steady_curl_residual(l1_velocity, (x, y, z))
    assert sp.factor(stationary_residual[2].subs(y, 0)) == (
        B * (5 * x**2 - 2 * z**2) * (A * x + C * z)
        / (x**2 + z**2) ** sp.Rational(9, 2)
    )
    assert sp.factor(stationary_residual[1].subs({y: 0, A: 0, C: 0})) == (
        6 * B**2 * x * z / (x**2 + z**2) ** 4
    )
    assert sp.factor(stationary_residual[1].subs({y: 0, A: 0, B: 0})) == (
        -4 * C**2 * x * z / (x**2 + z**2) ** 4
    )

    t, phi = sp.symbols("t phi", real=True)
    q22_derived = sp.pi / 4 * sp.integrate(
        (1 - t**2) ** 2 / (1 + t**2), (t, -1, 1)
    )
    q33_derived = sp.pi * sp.integrate(1 - t**2, (t, -1, 1))
    q11_derived = 4 * sp.pi - q22_derived - q33_derived
    stress = fixed_frame_angular_stress()
    q11, q22, q33 = stress[0, 0], stress[1, 1], stress[2, 2]
    assert sp.simplify(q11_derived - q11) == 0
    assert sp.simplify(q22_derived - q22) == 0
    assert sp.simplify(q33_derived - q33) == 0
    assert sp.simplify(q11 + q22 + q33 - 4 * sp.pi) == 0
    assert sp.simplify(q11 - q22) == sp.pi * (16 - 3 * sp.pi) / 3

    ndiag = sp.ImmutableMatrix([1, 1, 0]) / sp.sqrt(2)
    projected = stress * ndiag - ndiag * (ndiag.dot(stress * ndiag))
    assert sp.simplify(projected.dot(projected) - (q11 - q22) ** 2 / 4) == 0

    k = sp.symbols("k", positive=True)
    sine_transform = sp.integrate(
        sp.exp(-sigma**2 * k**2) * sp.sin(k * d) / k, (k, 0, sp.oo)
    )
    assert sine_transform == sp.pi * sp.erf(d / (2 * sigma)) / 2

    q = sp.symbols("q", nonzero=True)
    assert sp.integrate(q, (sp.Symbol("mu"), -1, 1), (sp.Symbol("phi"), 0, 2 * sp.pi)) == 4 * sp.pi * q

    print("P253/0072 exact checks: all passed")


if __name__ == "__main__":
    main()

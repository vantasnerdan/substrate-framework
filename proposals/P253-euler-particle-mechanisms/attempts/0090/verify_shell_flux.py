"""Exact normalization checks for P253/0090; no radiation numerics."""

from __future__ import annotations

import sympy as sp

from substrate_framework.euler_maxwell_radiation_flux import (
    doppler_shell_geometry,
    gaussian_curl_current_power,
    outgoing_power_prefactor,
    shell_sphere_power_weight,
    switched_continuity_residual,
)


def main() -> None:
    epsilon, mu, c, speed, omega = sp.symbols(
        "epsilon mu c speed omega", positive=True
    )
    n_z = sp.symbols("n_z", real=True)

    # With E*= -i Omega mu c^2 (D+i0 s)^-1 J*, the real delta
    # coefficient in J.E* is negative.  P_out=-<J.E> is positive.
    source_work_delta = -sp.pi * mu * c**2 * omega
    p_out_coefficient = -sp.Rational(1, 2) * source_work_delta
    assert p_out_coefficient == sp.pi * mu * c**2 * omega / 2
    assert sp.simplify(
        p_out_coefficient.subs(mu, 1 / (epsilon * c**2))
        - outgoing_power_prefactor(epsilon) * omega
    ) == 0
    print("PASS 1: outgoing Poynting sign and phasor one-half")

    geometry = doppler_shell_geometry(c, speed, omega, n_z)
    radius = geometry.radius
    Omega = omega + speed * radius * n_z
    D = sp.expand(c**2 * radius**2 - Omega**2)
    assert sp.factor(D) == 0
    assert sp.simplify(Omega - c * radius) == 0
    print("PASS 2: positive-frequency star-shaped Doppler shell")

    r_symbol = sp.Symbol("r", positive=True)
    radial = sp.diff(
        c**2 * r_symbol**2 - (omega + speed * r_symbol * n_z) ** 2,
        r_symbol,
    )
    radial_on_shell = sp.factor(radial.subs(r_symbol, radius))
    assert sp.simplify(radial_on_shell - 2 * c * omega) == 0
    assert geometry.radial_derivative_magnitude == 2 * c * omega
    print("PASS 3: exact radial coarea derivative")

    direct_weight = (
        outgoing_power_prefactor(epsilon)
        * geometry.temporal_frequency
        * radius**2
        / geometry.radial_derivative_magnitude
    )
    assert sp.simplify(
        direct_weight
        - shell_sphere_power_weight(epsilon, c, speed, omega, n_z)
    ) == 0
    print("PASS 4: exact sphere power weight")

    # The full gradient bound is exposed algebraically on the shell.
    lower = 4 * c**2 * radius**2 * (c - speed) ** 2
    remainder = 8 * c**3 * radius**2 * speed * (1 - n_z)
    assert sp.simplify(geometry.gradient_squared - lower - remainder) == 0
    print("PASS 5: full shell gradient has a subluminal lower bound")

    a, da, rho = sp.symbols("a da rho")
    assert switched_continuity_residual(a, da, 0, rho, -rho) == 0
    print("PASS 6: a-prime K switching current cancels continuity defect")

    # Independently integrate the angular square for k cross e_z.
    theta, phi = sp.symbols("theta phi", real=True)
    angular = sp.integrate(
        sp.integrate(sp.sin(theta) ** 3, (theta, 0, sp.pi)),
        (phi, 0, 2 * sp.pi),
    )
    assert angular == 8 * sp.pi / 3
    width, j_sq, pol_sq = sp.symbols("width j_sq pol_sq", positive=True)
    r0 = omega / c
    gaussian_direct = (
        sp.pi
        * omega**2
        / (4 * epsilon * c**3)
        * j_sq
        * r0**2
        * sp.exp(-(width * r0) ** 2)
        * angular
        * pol_sq
    )
    assert sp.simplify(
        gaussian_direct
        - gaussian_curl_current_power(
            epsilon, c, omega, width, j_sq, pol_sq
        )
    ) == 0
    print("PASS 7: exact Gaussian curl-current power")

    frequency, action, power, gamma = sp.symbols(
        "frequency action power gamma", positive=True
    )
    energy = frequency * action
    gamma_solution = sp.solve(sp.Eq(2 * gamma * energy, power), gamma)[0]
    assert gamma_solution == power / (2 * frequency * action)
    print("PASS 8: amplitude width versus energy-decay factor two")

    # The direct conjugate relation flips both k and omega.
    positive_negative_speed = doppler_shell_geometry(c, -speed, omega, n_z)
    conjugate = doppler_shell_geometry(c, -speed, -omega, -n_z)
    assert sp.simplify(positive_negative_speed.radius - conjugate.radius) == 0
    assert sp.simplify(
        positive_negative_speed.temporal_frequency
        + conjugate.temporal_frequency
    ) == 0
    print("PASS 9: negative-frequency conjugate shell is not double counted")


if __name__ == "__main__":
    main()

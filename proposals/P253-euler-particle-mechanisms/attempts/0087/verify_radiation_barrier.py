"""Exact algebra checks for P253/0087; no spectral or radiation numerics."""

from __future__ import annotations

import sympy as sp

from substrate_framework.euler_high_harmonic_radiation import (
    cylindrical_vector_bessel_orders,
    debye_exponential_rate,
    fixed_mode_bessel_ratio,
    high_index_limiting_bessel_ratio,
    physical_mode_ledger,
    transverse_shell_maximum,
)


def main() -> None:
    kappa, radius, k, n, sigma = sp.symbols(
        "kappa radius k n sigma", positive=True
    )
    c_em, speed = sp.symbols("c_em speed", positive=True)
    ledger = physical_mode_ledger(kappa, radius, k, n, sigma)
    assert sp.simplify(
        ledger.physical_frequency
        - kappa * sigma * n**2 / (2 * sp.pi * radius**2 * k**2)
    ) == 0
    print("PASS 1: Gamma=kappa physical core clock")

    shell = transverse_shell_maximum(c_em, speed, ledger.physical_frequency)
    t = speed / c_em
    direct = ledger.physical_frequency * sp.sqrt(1 - t**2) / (c_em - speed * t)
    assert sp.simplify(direct - shell.transverse_wave_number) == 0
    print("PASS 2: sharp transverse Doppler-shell maximum")

    assert cylindrical_vector_bessel_orders(n) == (n - 1, n, n + 1)
    print("PASS 3: cylindrical vector Bessel orders")

    ratio = fixed_mode_bessel_ratio(kappa, sigma, n, k, radius, c_em, speed)
    target = (
        kappa
        * sigma
        * n
        * (n + k)
        / (2 * sp.pi * radius * k**2 * sp.sqrt(c_em**2 - speed**2) * (n - 1))
    )
    assert sp.simplify(ratio - target) == 0
    print("PASS 4: fixed-J finite radiation ratio")

    ell, l_phi = sp.symbols("ell l_phi", positive=True)
    high = high_index_limiting_bessel_ratio(
        kappa, ell, l_phi, radius, k, c_em, speed
    )
    assert sp.simplify(
        high
        - kappa
        * ell
        * l_phi
        / (2 * sp.pi**2 * radius * k * sp.sqrt(c_em**2 - speed**2))
    ) == 0
    print("PASS 5: finite-J leading predictor at supplied carrier speed")

    q = sp.Rational(4, 5)
    assert sp.simplify(
        debye_exponential_rate(q)
        - (sp.acosh(sp.Rational(5, 4)) - sp.Rational(3, 5))
    ) == 0
    print("PASS 6: optimized Debye contour rate")


if __name__ == "__main__":
    main()

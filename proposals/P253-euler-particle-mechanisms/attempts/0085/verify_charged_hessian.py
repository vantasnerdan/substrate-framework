"""Independent exact local identities for P253/0085."""

from __future__ import annotations

import sympy as sp

from substrate_framework.euler_charged_hessian import (
    comoving_maxwell_quadratic,
    comoving_radiation_shell,
    comoving_maxwell_radiation_denominator,
    constrained_charge_schur_symbol,
    material_tag_locking_coefficients,
)


def main() -> None:
    epsilon, mu, c = sp.symbols("epsilon mu c", positive=True)
    electric = sp.symbols("E0:3", real=True)
    magnetic = sp.symbols("B0:3", real=True)
    field = comoving_maxwell_quadratic(epsilon, mu, c, electric, magnetic)
    assert sp.simplify(field.relative_density - field.completed_density) == 0
    print("PASS 1: comoving Maxwell energy-momentum square is exact")

    k_perp, k_z, rho2 = sp.symbols("k_perp k_z rho2", positive=True)
    beta2 = epsilon * mu * c**2
    schur = constrained_charge_schur_symbol(
        epsilon, mu, c, [k_perp, 0, k_z], rho2
    )
    expected = rho2 * (1 - beta2) / (
        2 * epsilon * (k_perp**2 + (1 - beta2) * k_z**2)
    )
    assert sp.simplify(schur - expected) == 0
    print("PASS 2: Gauss-fiber minimization has the anisotropic H-minus-one symbol")

    f_prime, zeta, zeta_prime, n = sp.symbols(
        "f_prime zeta zeta_prime n", nonzero=True
    )
    lock = material_tag_locking_coefficients(f_prime, zeta, zeta_prime, n)
    assert lock.axisymmetric_from_toroidal_vorticity == f_prime / zeta_prime
    assert sp.simplify(
        lock.nonaxisymmetric_from_radial_vorticity
        + f_prime / (sp.I * n * zeta)
    ) == 0
    print("PASS 3: regular-band tag variations are locked to vorticity modes")

    c_em, omega = sp.symbols("c_em omega", positive=True)
    denominator = comoving_maxwell_radiation_denominator(
        c_em, c, omega, [k_perp, 0, k_z]
    )
    assert sp.simplify(
        denominator
        - (k_perp**2 + k_z**2 - (omega + c * k_z) ** 2 / c_em**2)
    ) == 0
    print("PASS 4: translating Maxwell radiation shell has the declared sign")

    gap = sp.symbols("gap", positive=True)
    bounded = expected.subs(beta2, 1 - gap)
    assert sp.simplify(
        bounded
        - rho2 * gap / (2 * epsilon * (k_perp**2 + gap * k_z**2))
    ) == 0
    print("PASS 5: fixed subluminal gap gives a positive charge Schur form")

    direction_z = sp.symbols("direction_z", real=True)
    shell = comoving_radiation_shell(c_em, c, omega, 1, direction_z)
    radial = sp.symbols("radial", positive=True)
    shell_denominator = (
        radial**2
        - (omega + c * radial * direction_z) ** 2 / c_em**2
    )
    shell_derivative = sp.factor(
        sp.diff(shell_denominator, radial).subs(radial, shell.radius)
    )
    assert sp.simplify(
        shell_derivative - shell.radial_derivative_magnitude
    ) == 0
    assert sp.simplify(
        shell.coarea_weight
        - omega * c_em / (2 * (c_em - c * direction_z) ** 2)
    ) == 0
    print("PASS 6: radiation-shell root and coarea weight are derived")


if __name__ == "__main__":
    main()

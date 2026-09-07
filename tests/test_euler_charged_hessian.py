import pytest
import sympy as sp

from substrate_framework.euler_charged_hessian import (
    comoving_maxwell_quadratic,
    comoving_radiation_shell,
    comoving_maxwell_radiation_denominator,
    constrained_charge_schur_symbol,
    material_tag_locking_coefficients,
)


def test_comoving_maxwell_completion_is_exact():
    epsilon, mu, c = sp.symbols("epsilon mu c", positive=True)
    electric = sp.symbols("E0:3", real=True)
    magnetic = sp.symbols("B0:3", real=True)
    quadratic = comoving_maxwell_quadratic(
        epsilon, mu, c, electric, magnetic
    )
    assert sp.simplify(
        quadratic.relative_density - quadratic.completed_density
    ) == 0
    assert quadratic.speed_ratio_squared == epsilon * mu * c**2


def test_charge_schur_symbol_has_exact_parallel_and_transverse_limits():
    epsilon, mu, c, k, rho2 = sp.symbols(
        "epsilon mu c k rho2", positive=True
    )
    beta2 = epsilon * mu * c**2
    parallel = constrained_charge_schur_symbol(
        epsilon, mu, c, [0, 0, k], rho2
    )
    transverse = constrained_charge_schur_symbol(
        epsilon, mu, c, [k, 0, 0], rho2
    )
    assert sp.simplify(parallel - rho2 / (2 * epsilon * k**2)) == 0
    assert sp.simplify(
        transverse - rho2 * (1 - beta2) / (2 * epsilon * k**2)
    ) == 0


def test_charge_schur_symbol_is_the_direct_mode_minimum():
    epsilon, mu, c = sp.symbols("epsilon mu c", positive=True)
    sine, cosine, electric_normal = sp.symbols(
        "sine cosine electric_normal", real=True
    )
    electric_tangent, magnetic_tangent = sp.symbols(
        "electric_tangent magnetic_tangent", real=True
    )
    beta2 = epsilon * mu * c**2
    mode_energy = (
        epsilon * (electric_normal**2 + electric_tangent**2) / 2
        + magnetic_tangent**2 / (2 * mu)
        - c
        * epsilon
        * (cosine * electric_tangent + sine * electric_normal)
        * magnetic_tangent
    )
    e_stationary = sp.solve(
        sp.diff(mode_energy, electric_tangent), electric_tangent
    )[0]
    reduced = sp.factor(mode_energy.subs(electric_tangent, e_stationary))
    b_stationary = sp.solve(
        sp.diff(reduced, magnetic_tangent), magnetic_tangent
    )[0]
    minimum = sp.factor(reduced.subs(magnetic_tangent, b_stationary))
    expected = sp.factor(
        epsilon
        * electric_normal**2
        * (1 - beta2 * (sine**2 + cosine**2))
        / (2 * (1 - beta2 * cosine**2))
    )
    assert sp.simplify(minimum - expected) == 0
    assert sp.simplify(
        minimum.subs(sine**2, 1 - cosine**2)
        - epsilon
        * electric_normal**2
        * (1 - beta2)
        / (2 * (1 - beta2 * cosine**2))
    ) == 0


def test_material_tag_is_locked_to_vorticity_on_regular_modes():
    f_prime, zeta, zeta_prime = sp.symbols(
        "f_prime zeta zeta_prime", nonzero=True
    )
    n = sp.symbols("n", integer=True, nonzero=True)
    locking = material_tag_locking_coefficients(
        f_prime, zeta, zeta_prime, n
    )
    assert locking.axisymmetric_from_toroidal_vorticity == f_prime / zeta_prime
    assert sp.simplify(
        locking.nonaxisymmetric_from_radial_vorticity
        + f_prime / (sp.I * n * zeta)
    ) == 0


def test_comoving_radiation_shell_and_steady_subluminal_limit():
    c_em, c, omega, kz = sp.symbols(
        "c_em c omega kz", positive=True
    )
    denominator = comoving_maxwell_radiation_denominator(
        c_em, c, omega, [0, 0, kz]
    )
    assert sp.simplify(
        denominator - (kz**2 - (c * kz + omega) ** 2 / c_em**2)
    ) == 0
    steady = sp.factor(denominator.subs(omega, 0))
    assert sp.simplify(
        steady - kz**2 * (c_em - c) * (c_em + c) / c_em**2
    ) == 0


def test_comoving_radiation_shell_coarea_weight_is_derived():
    c_em, c, omega_abs, n_z = sp.symbols(
        "c_em c omega_abs n_z", positive=True
    )
    shell = comoving_radiation_shell(c_em, c, omega_abs, 1, n_z)
    assert sp.simplify(
        shell.radius - omega_abs / (c_em - c * n_z)
    ) == 0
    radial = sp.symbols("radial", positive=True)
    denominator = radial**2 - (omega_abs + c * radial * n_z) ** 2 / c_em**2
    derivative_on_shell = sp.factor(
        sp.diff(denominator, radial).subs(radial, shell.radius)
    )
    assert sp.simplify(
        derivative_on_shell - shell.radial_derivative_magnitude
    ) == 0
    assert sp.simplify(
        shell.coarea_weight
        - omega_abs * c_em / (2 * (c_em - c * n_z) ** 2)
    ) == 0


def test_quadratic_helpers_reject_invalid_domains():
    with pytest.raises(ValueError, match="permittivity"):
        comoving_maxwell_quadratic(0, 1, 0, [0, 0, 0], [0, 0, 0])
    with pytest.raises(ValueError, match="wavevector"):
        constrained_charge_schur_symbol(1, 1, 0, [0, 0, 0], 1)
    with pytest.raises(ValueError, match="frequency_sign"):
        comoving_radiation_shell(1, 0, 1, 0, 0)

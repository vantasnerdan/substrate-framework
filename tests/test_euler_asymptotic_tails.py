import sympy as sp

from substrate_framework.euler_asymptotic_tails import (
    degree_minus_two_domain_ledger,
    fixed_frame_constant_norm_symbol,
    fixed_frame_angular_stress,
    gaussian_tail_cross_potential,
    l1_multiplicity_cross_block,
    l1_homogeneous_tail,
    newton_pressure_kernel,
    oriented_tail_cross_kernel,
    pressure_stress_gradient,
    radial_tail_fourier_coefficient,
    scalar_homogeneous_fourier_coefficient,
    toroidal_tail_fourier_coefficient,
    steady_curl_residual,
    tail_evolution_order_ledger,
)


def test_homogeneous_transform_coefficients_recover_newton_and_l1_toroidal():
    assert scalar_homogeneous_fourier_coefficient(0, 1) == 4 * sp.pi
    assert toroidal_tail_fourier_coefficient(1) == -4 * sp.I * sp.pi
    assert radial_tail_fourier_coefficient(1) == sp.pi**2


def test_source_free_tail_rejects_l0_radial_and_toroidal_modes():
    for fn in (radial_tail_fourier_coefficient, toroidal_tail_fourier_coefficient):
        try:
            fn(0)
        except ValueError:
            pass
        else:
            raise AssertionError("l=0 tail must be rejected")


def test_oriented_kernel_is_repulsive_for_like_amplitude_but_anisotropic():
    rho, d = sp.symbols("rho d", positive=True)
    a = (1, 0, 0)
    parallel = oriented_tail_cross_kernel(a, a, a, density=rho, separation=d)
    transverse = oriented_tail_cross_kernel(
        a, a, (0, 1, 0), density=rho, separation=d
    )
    assert parallel == 4 * sp.pi * rho / d
    assert transverse == 2 * sp.pi * rho / d
    assert sp.simplify(parallel - transverse) == 2 * sp.pi * rho / d


def test_domain_ledger_separates_energy_from_moment_action():
    ledger = degree_minus_two_domain_ledger("toroidal_l1")
    assert ledger.finite_kinetic_energy
    assert ledger.finite_helicity_at_infinity
    assert not ledger.absolutely_integrable_velocity
    assert not ledger.absolutely_integrable_vorticity
    assert not ledger.finite_ordinary_angular_momentum
    assert not ledger.finite_absolute_vorticity_impulse
    radial = degree_minus_two_domain_ledger("radial")
    assert radial.finite_ordinary_angular_momentum
    assert not radial.finite_absolute_vorticity_impulse


def test_fixed_frame_escape_is_transverse_real_and_constant_norm():
    n1, n2, n3 = sp.symbols("n1 n2 n3", real=True)
    n = sp.ImmutableMatrix([n1, n2, n3])
    f = fixed_frame_constant_norm_symbol(n)
    sphere = {n1**2: 1 - n2**2 - n3**2}
    assert sp.simplify((n.dot(f)).subs(sphere)) == 0
    norm = sum(sp.conjugate(f[j]) * f[j] for j in range(3))
    assert sp.simplify(sp.expand(norm - 1).subs(sphere)) == 0
    f_minus = fixed_frame_constant_norm_symbol((-n1, -n2, -n3))
    assert all(sp.simplify(f_minus[j] - sp.conjugate(f[j])) == 0 for j in range(3))


def test_fixed_frame_pi_rotation_flips_the_symbol():
    n1, n2, n3 = sp.symbols("n1 n2 n3", real=True)
    rotation = sp.diag(-1, -1, 1)
    n = sp.ImmutableMatrix([n1, n2, n3])
    rotated = rotation * fixed_frame_constant_norm_symbol(rotation.T * n)
    original = fixed_frame_constant_norm_symbol(n)
    assert all(sp.simplify(rotated[j] + original[j]) == 0 for j in range(3))


def test_l1_multiplicity_block_has_even_and_chiral_entries():
    rho, d = sp.symbols("rho d", positive=True)
    block = l1_multiplicity_cross_block(
        (1, 0, 0),
        (1, 0, 0),
        (0, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        density=rho,
        separation=d,
    )
    assert block[0, 0] == 0
    assert block[1, 1] == 0
    assert block[0, 1] == -2 * sp.pi * rho / d
    assert block[1, 0] == -2 * sp.pi * rho / d


def test_gaussian_fixed_frame_cross_energy_is_finite_and_coulombic():
    q1, q2, rho, sigma, d = sp.symbols("q1 q2 rho sigma d", positive=True)
    potential = gaussian_tail_cross_potential(q1, q2, d, sigma, density=rho)
    assert sp.limit(potential, d, 0, dir="+") == rho * q1 * q2 / (
        4 * sp.pi ** sp.Rational(3, 2) * sigma
    )
    assert sp.limit(d * potential, d, sp.oo) == rho * q1 * q2 / (4 * sp.pi)


def test_complete_l1_multiplicity_has_no_nonzero_steady_tail():
    x, y, z, A, B, C = sp.symbols("x y z A B C", real=True)
    u = l1_homogeneous_tail((x, y, z), (A, 0, C), (0, 0, B))
    residual = steady_curl_residual(u, (x, y, z))
    mixed = sp.factor(residual[2].subs(y, 0))
    assert mixed == B * (5 * x**2 - 2 * z**2) * (A * x + C * z) / (
        x**2 + z**2
    ) ** sp.Rational(9, 2)
    pure_toroidal = sp.factor(residual[1].subs({y: 0, A: 0, C: 0}))
    assert pure_toroidal == 6 * B**2 * x * z / (x**2 + z**2) ** 4
    pure_radial = sp.factor(residual[1].subs({y: 0, A: 0, B: 0}))
    assert pure_radial == -4 * C**2 * x * z / (x**2 + z**2) ** 4


def test_fixed_frame_minimal_realization_violates_stationary_stress_row():
    stress = fixed_frame_angular_stress()
    assert sp.simplify(sp.trace(stress) - 4 * sp.pi) == 0
    assert sp.simplify(stress[0, 0] - stress[1, 1]) == sp.pi * (
        16 - 3 * sp.pi
    ) / 3


def test_pressure_stress_kernel_and_gradient_have_exact_far_field_orders():
    x, y, z, scale = sp.symbols("x y z scale", positive=True)
    point = (x, y, z)
    kernel = newton_pressure_kernel(point)
    scaled_kernel = newton_pressure_kernel((scale * x, scale * y, scale * z))
    assert all(
        sp.simplify(scaled_kernel[i, j] - kernel[i, j] / scale**3) == 0
        for i in range(3)
        for j in range(3)
    )
    stress = sp.diag(1, 2, 4)
    gradient = pressure_stress_gradient(point, stress)
    scaled_gradient = gradient.subs(
        {x: scale * x, y: scale * y, z: scale * z}, simultaneous=True
    )
    assert all(
        sp.simplify(scaled_gradient[i] - gradient[i] / scale**4) == 0
        for i in range(3)
    )
    assert sp.simplify(sp.trace(kernel)) == 0


def test_tail_evolution_order_ledger_retains_borderline_moment():
    ledger = tail_evolution_order_ledger()
    assert ledger.stress == 4
    assert ledger.pressure == 3
    assert ledger.pressure_gradient == 4
    assert ledger.transport == 5
    assert ledger.velocity_time_derivative == 4
    assert ledger.stress_first_moment_is_logarithmic

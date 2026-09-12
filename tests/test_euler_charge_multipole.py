import pytest
import sympy as sp

from substrate_framework.euler_charge_multipole import (
    internal_stress_source,
    pressure_quadrupole_far_field,
    radial_monopole_flux,
    transverse_static_green,
    transverse_scalar_source,
    vortex_dipole_cross_energy,
    vortex_impulse_far_velocity,
)


def test_vortex_impulse_velocity_has_dipole_sign_and_cubic_decay():
    impulse = sp.Matrix([0, 0, 2])
    assert vortex_impulse_far_velocity(impulse, [0, 0, 3]) == sp.Matrix(
        [0, 0, sp.Rational(1, 27) / sp.pi]
    )
    assert vortex_impulse_far_velocity(impulse, [3, 0, 0]) == sp.Matrix(
        [0, 0, -sp.Rational(1, 54) / sp.pi]
    )


def test_cross_energy_has_dipole_power_and_orientation():
    aligned = vortex_dipole_cross_energy(2, [0, 0, 1], [0, 0, 1], [0, 0, 3])
    transverse = vortex_dipole_cross_energy(2, [0, 0, 1], [0, 0, 1], [3, 0, 0])
    assert aligned == sp.Rational(1, 27) / sp.pi
    assert transverse == -sp.Rational(1, 54) / sp.pi
    scaled = vortex_dipole_cross_energy(2, [0, 0, 1], [0, 0, 1], [0, 0, 6])
    assert sp.simplify(scaled / aligned) == sp.Rational(1, 8)


def test_physical_pressure_quadrupole_keeps_density_and_trace_cancellation():
    moment = sp.diag(1, 1, 0)
    assert pressure_quadrupole_far_field(3, moment, [0, 0, 2]) == -sp.Rational(3, 16) / sp.pi
    assert pressure_quadrupole_far_field(3, sp.eye(3), [0, 0, 2]) == 0


def test_radial_monopole_has_nonzero_flux_and_domains_are_checked():
    q = sp.symbols("q")
    assert radial_monopole_flux(q) == 4 * sp.pi * q
    with pytest.raises(ValueError):
        vortex_impulse_far_velocity([1, 2, 3], [0, 0, 0])
    with pytest.raises(ValueError):
        vortex_dipole_cross_energy(0, [1, 0, 0], [1, 0, 0], [1, 0, 0])
    with pytest.raises(ValueError):
        pressure_quadrupole_far_field(1, [[1, 0], [0, 1]], [1, 0, 0])


def test_translation_invariant_internal_stress_has_no_force_monopole():
    stress = sp.Matrix(sp.symbols("s0:9")).reshape(3, 3)
    assert internal_stress_source([0, 0, 0], stress) == sp.zeros(3, 1)
    assert internal_stress_source([1, 0, 0], sp.eye(3)) == sp.Matrix([sp.I, 0, 0])
    green_near = transverse_static_green(2, [0, 0, 3])
    green_far = transverse_static_green(2, [0, 0, 6])
    assert green_far == green_near / 2


def test_isotropic_scalar_source_is_killed_by_transverse_projection():
    k1, k2, k3, factor = sp.symbols("k1 k2 k3 factor", nonzero=True)
    assert transverse_scalar_source([k1, k2, k3], factor) == sp.zeros(3, 1)
    with pytest.raises(ValueError):
        transverse_scalar_source([0, 0, 0])

import sympy as sp

from substrate_framework.euler_neutrino_suppliers import (
    current_parity,
    helicity_flux,
    material_frame_hamiltonian,
    packet_transition_probability,
    relativistic_phase_leading,
    two_flavor_transition_probability,
)


def test_two_flavor_probability_limits():
    theta, dnu, t = sp.symbols("theta dnu t", real=True)
    p = two_flavor_transition_probability(theta, dnu, t)
    assert sp.simplify(p.subs(t, 0)) == 0
    assert sp.simplify(p.subs(theta, 0)) == 0
    assert sp.simplify(p.subs({theta: sp.pi / 4, dnu * t: sp.pi})) == 1


def test_packet_probability_reduces_to_fixed_momentum():
    theta, phi = sp.symbols("theta phi", real=True)
    got = packet_transition_probability(theta, sp.exp(-sp.I * phi))
    want = sp.sin(2 * theta) ** 2 * sp.sin(phi / 2) ** 2
    assert sp.simplify(sp.expand_trig(got - want)) == 0


def test_relativistic_phase_has_shared_action_denominator():
    dm2, c, length, energy, action = sp.symbols(
        "dm2 c length energy action", nonzero=True
    )
    assert relativistic_phase_leading(dm2, c, length, energy, action) == (
        dm2 * c**3 * length / (2 * energy * action)
    )


def test_material_frame_connection_is_inserted_by_frame_derivative():
    a = sp.symbols("a", real=True)
    dt_u = sp.Matrix([[0, -a], [a, 0]])
    identity = sp.eye(2)
    assert material_frame_hamiltonian(dt_u, identity) == sp.I * dt_u


def test_helicity_flux_pressure_coefficient():
    u = sp.Matrix([1, 2, 0])
    omega = sp.Matrix([0, 3, 4])
    h, flux = helicity_flux(u, omega, sp.Symbol("p"), u.dot(u))
    assert h == 6
    assert flux == h * u + (sp.Symbol("p") - sp.Rational(5, 2)) * omega


def test_scalar_and_pseudoscalar_currents_have_distinct_parity():
    assert current_parity("scalar") == {
        "density": "scalar",
        "spatial_current": "polar",
    }
    assert current_parity("pseudoscalar") == {
        "density": "pseudoscalar",
        "spatial_current": "axial",
    }


def test_current_parity_rejects_untyped_tag():
    try:
        current_parity("unknown")
    except ValueError:
        pass
    else:
        raise AssertionError("untyped parity must be rejected")

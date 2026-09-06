"""Exposing current-definition and scalar-action transfer tests."""

import pytest
import sympy as sp

from substrate_framework.euler_field_maps import (
    euler_lamb_map, madelung_kinetic_terms, spinor_euler_terms,
)


def test_exact_steady_shear_has_defined_current_not_source_free_maxwell():
    x, y, z, t = sp.symbols("x y z t", real=True)
    c = sp.symbols("c", positive=True)
    result = euler_lamb_map([sp.sin(y), 0, 0], 5, [x, y, z], t, density=3, speed=c)
    assert result.euler_residual == sp.zeros(3, 1)
    assert result.faraday_residual == sp.zeros(3, 1)
    assert result.hydrodynamic_charge == -sp.cos(2*y)
    assert result.defined_current == sp.Matrix([c*c*sp.sin(y), 0, 0])
    assert result.charge_continuity_residual == 0
    assert sp.integrate(result.hydrodynamic_charge, (y, -sp.pi, sp.pi)) == 0


def test_defined_maxwell_sources_do_not_false_certify_euler():
    x, y, z, t = sp.symbols("x y z t", real=True)
    result = euler_lamb_map([t*x, 0, 0], 0, [x, y, z], t, density=1, speed=2)
    assert result.lamb_vector == sp.zeros(3, 1)
    assert result.defined_current == sp.zeros(3, 1)
    assert result.charge_continuity_residual == 0
    assert result.euler_residual == sp.Matrix([x*(1+t*t), 0, 0])
    assert result.divergence_residual == t
    with pytest.raises(ValueError, match="constants"):
        euler_lamb_map([0, 0, 0], 0, [x, y, z], t, density=1, speed=sp.exp(x))


def test_madelung_energy_against_direct_wavefunction_derivatives():
    x, y, z = sp.symbols("x y z", real=True)
    m, hbar = sp.symbols("m hbar", positive=True)
    n = sp.exp(-x*x-y*y-z*z)
    S = x*y+z*z
    result = madelung_kinetic_terms(n, S, [x, y, z], mass=m, action=hbar)
    psi = sp.sqrt(n)*sp.exp(sp.I*S/hbar)
    wave_energy = hbar*hbar*sum(sp.conjugate(psi.diff(q))*psi.diff(q) for q in [x, y, z])/(2*m)
    assert sp.simplify(wave_energy-result.classical_kinetic_density-result.density_gradient_energy) == 0
    assert sp.simplify(result.quantum_potential-hbar*hbar*(3-x*x-y*y-z*z)/(2*m)) == 0
    assert result.density_gradient_energy != 0


def test_quantum_potential_is_actual_gradient_energy_variation():
    x, y, z = sp.symbols("x y z", real=True)
    m, hbar = sp.symbols("m hbar", positive=True)
    n = sp.Function("n", positive=True)(x)
    result = madelung_kinetic_terms(n, x*y, [x, y, z], mass=m, action=hbar)
    energy = result.density_gradient_energy
    variational_derivative = sp.diff(energy, n)-sp.diff(sp.diff(energy, n.diff(x)), x)
    assert sp.simplify(variational_derivative-result.quantum_potential) == 0
    uniform = madelung_kinetic_terms(2, x*y, [x, y, z], mass=m, action=hbar)
    assert uniform.density_gradient_energy == 0
    assert uniform.quantum_potential == 0


def test_unit_spinor_velocity_and_texture_energy_are_distinct():
    x, y, z = sp.symbols("x y z", real=True)
    rho, kappa = sp.symbols("rho kappa", positive=True)
    spinor = [sp.cos(y/2), sp.exp(sp.I*x)*sp.sin(y/2)]
    result = spinor_euler_terms(spinor, [x, y, z], density=rho, circulation_scale=kappa)
    assert result.normalization_residual == 0
    assert result.velocity == sp.Matrix([kappa*sp.sin(y/2)**2, 0, 0])
    assert result.texture == sp.Matrix([sp.sin(y)*sp.cos(x), sp.sin(y)*sp.sin(x), sp.cos(y)])
    assert result.energy_identity_residual == 0
    assert result.texture_energy_density != 0
    assert sp.simplify(result.euler_kinetic_density-rho*kappa*kappa*sp.sin(y/2)**4/2) == 0
    with pytest.raises(ValueError, match="unit"):
        spinor_euler_terms([2, 0], [x, y, z], density=1, circulation_scale=1)


def test_advected_spinor_labels_execute_actual_vortical_euler():
    # Local regular chart 0<y<1; no global/finite-energy assertion.
    x, y, z, t = sp.symbols("x y z t", real=True)
    theta, beta, f = y*y*t/2, x-y*t, y
    u = sp.Matrix([sp.diff(theta, q)+f*sp.diff(beta, q) for q in [x, y, z]])
    assert u == sp.Matrix([y, 0, 0])
    spinor = sp.exp(sp.I*theta)*sp.Matrix([sp.sqrt(1-y), sp.sqrt(y)*sp.exp(sp.I*beta)])
    material = spinor.diff(t)+y*spinor.diff(x)
    assert (sp.I*material+y*y*spinor/2).applyfunc(sp.simplify) == sp.zeros(2, 1)
    physical = euler_lamb_map(u, 0, [x, y, z], t, density=1, speed=1)
    assert physical.euler_residual == sp.zeros(3, 1)
    assert physical.vorticity == sp.Matrix([0, 0, -1])

"""Direct field and full-pressure initial-action checks for P251/0167."""

import pytest
import sympy as s

from substrate_framework.euler_acoustic import (
    axial_kelvin_initial_phase,
    triangular_euler_array,
)


def test_actual_cartesian_field_and_pressure():
    x, y = s.symbols("x y", real=True)
    rho = s.Symbol("rho", positive=True)
    field = triangular_euler_array(2, 3, rho, (x, y))
    v = field.velocity
    assert s.trigsimp(s.diff(v[0], x)+s.diff(v[1], y)) == 0
    assert s.trigsimp(field.vorticity+9*field.streamfunction) == 0
    residual = rho*v.jacobian((x, y))*v+s.Matrix([
        s.diff(field.pressure, x), s.diff(field.pressure, y)])
    assert all(s.trigsimp(entry) == 0 for entry in residual)
    assert field.domain_condition is s.true


def test_actual_complete_cell_covariance_and_separatrix_dual():
    x, y, a, b = s.symbols("x y a b", real=True)
    amp, lam = s.symbols("Psi lambda", positive=True)
    field = triangular_euler_array(amp, lam, 1, (x, y))
    # Independent normalized oblique-cell integration of the actual field.
    substituted = field.velocity.subs({x: a/lam, y: (2*b+a)/(s.sqrt(3)*lam)},
                                     simultaneous=True).applyfunc(s.trigsimp)
    measured = (substituted*substituted.T).applyfunc(lambda value: s.integrate(
        s.integrate(s.expand_trig(value), (a, -s.pi, s.pi)),
        (b, -s.pi, s.pi))/(4*s.pi**2))
    assert s.simplify(measured-field.covariance) == s.zeros(2)
    assert field.covariance == 3*amp**2*lam**2*s.eye(2)/4
    assert field.separatrix_level == -amp
    assert field.translation_dual.det() == amp**2*lam**4


def test_full_kelvin_phase_pressure_factor_and_fluid_limit():
    x, y = s.symbols("x y", real=True)
    amp, lam, rho, k = s.symbols("Psi lambda rho k", positive=True)
    field = triangular_euler_array(amp, lam, rho, (x, y))
    phase = axial_kelvin_initial_phase(field.wavevectors, field.sine_velocities, k, rho)
    assert s.simplify(phase.stiffness-rho*k**2*lam**2/(lam**2+k**2)
                      *field.covariance) == s.zeros(2)
    assert phase.mass == rho*s.eye(2)
    assert phase.symplectic+phase.symplectic.T == s.zeros(4)
    assert phase.hamiltonian == s.diag(phase.stiffness, phase.mass)
    assert phase.stiffness.subs(k, 0) == s.zeros(2)
    assert s.simplify(phase.stiffness-rho*k**2*field.covariance) != s.zeros(2)


def test_general_spectrum_projects_each_actual_mode():
    k = s.Symbol("k", real=True)
    wave = s.Matrix([[1, 0], [0, 2]])
    coefficients = s.Matrix([[0, 3], [2, 0]])
    phase = axial_kelvin_initial_phase(wave, coefficients, k, 5)
    expected = s.diag(90*k**2/(4+k**2), 10*k**2/(1+k**2))
    assert s.simplify(phase.stiffness-expected) == s.zeros(2)
    assert phase.domain_condition is s.true


def test_input_domains_and_fourier_orthogonality():
    x, y = s.symbols("x y", real=True)
    for bad in (0, -1, s.oo, s.nan, s.I):
        with pytest.raises(ValueError):
            triangular_euler_array(bad, 1, 1, (x, y))
    with pytest.raises(ValueError, match="spatial constants"):
        triangular_euler_array(x, 1, 1, (x, y))
    with pytest.raises(ValueError, match="distinct"):
        triangular_euler_array(1, 1, 1, (x, x))
    with pytest.raises(ValueError, match="perpendicular"):
        axial_kelvin_initial_phase([[1], [0]], [[1], [0]], 1, 1)
    with pytest.raises(ValueError, match="duplicate or opposite"):
        axial_kelvin_initial_phase([[1, -1], [0, 0]], [[0, 0], [1, 2]], 1, 1)
    with pytest.raises(ValueError, match="harmonic"):
        axial_kelvin_initial_phase([[0], [0]], [[1], [0]], 1, 1)
    with pytest.raises(ValueError, match="finite and real"):
        axial_kelvin_initial_phase([[1], [0]], [[0], [1]], s.oo, 1)

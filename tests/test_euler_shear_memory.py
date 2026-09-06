"""Exact full Euler residual and projected-memory exposing identities."""

import pytest
import sympy as sp

from substrate_framework.euler_shear_memory import (
    sinusoidal_shear_memory,
    triangular_euler_velocity,
)


def test_finite_amplitude_triangular_solution_satisfies_all_euler_components():
    x, y, z, t = sp.symbols("x y z t", real=True)
    U = 2*sp.sin(y)+sp.cos(2*y)
    w0 = 3*sp.sin(x)*sp.cos(y)+sp.cos(2*x+y)
    u = triangular_euler_velocity(U, w0, [x, y, z], t)
    jacobian = u.jacobian([x, y, z])
    assert sp.simplify(sp.trace(jacobian)) == 0
    assert sp.simplify(u.diff(t)+jacobian*u) == sp.zeros(3, 1)
    assert u.subs(t, 0) == sp.Matrix([U, 0, w0])
    # A time-varying shear or z-dependent vertical field breaks the exact
    # invariant class and is rejected rather than silently transported.
    with pytest.raises(ValueError, match="shear"):
        triangular_euler_velocity(t*sp.sin(y), w0, [x, y, z], t)
    with pytest.raises(ValueError, match="vertical"):
        triangular_euler_velocity(U, z+sp.sin(x), [x, y, z], t)


def test_bessel_transfer_matches_exact_sine_moments_and_unresolved_initial_data():
    y, t, beta, eps = sp.symbols("y t beta eps", real=True)
    response = sinusoidal_shear_memory(beta, t, unresolved_sine=eps)
    actual_series = sp.series(response.amplitude, t, 0, 7).removeO().expand()
    # Integrate the actual transported field's Taylor coefficients; this
    # independently exposes parity, initial-noise sign and normalization.
    expected = 0
    for n in range(7):
        coefficient = sp.integrate((-sp.I*beta*sp.sin(y))**n*(1+sp.I*eps*sp.sin(y)),
                                   (y, -sp.pi, sp.pi))/(2*sp.pi*sp.factorial(n))
        expected += coefficient*t**n
    assert sp.simplify(actual_series-expected) == 0
    assert actual_series.diff(t).subs(t, 0) == beta*eps/2
    opposite = sinusoidal_shear_memory(beta, 0, unresolved_sine=-eps)
    assert opposite.amplitude == 1
    assert opposite.unresolved_forcing == -beta*eps/2


def test_volterra_identity_from_exact_resolvent_and_nonzero_initial_force():
    s, beta, eps = sp.symbols("s beta eps", positive=True)
    # Hodge block resolvent: average 1/(s+i*beta*sin(y)) is
    # 1/sqrt(s^2+beta^2); the Bessel derivative supplies LJ1 below.
    root = sp.sqrt(s*s+beta*beta)
    laplace_j0 = 1/root
    laplace_j1 = (1-s/root)/beta
    Khat = root-s
    Fhat = eps*Khat/beta
    ahat = laplace_j0+eps*laplace_j1
    assert sp.simplify(s*ahat-1-Fhat+Khat*ahat) == 0
    assert sp.simplify(sp.diff(Khat/beta, s)+laplace_j1) == 0
    assert sp.limit(s*Khat, s, sp.oo) == beta**2/2
    assert sinusoidal_shear_memory(beta, 0).kernel == beta**2/2
    assert sinusoidal_shear_memory(0, s).kernel == 0


def test_initial_memory_energy_is_retained_and_invalid_complex_inputs_rejected():
    y, eps = sp.symbols("y eps", real=True)
    profile = 1+sp.I*eps*sp.sin(y)
    energy = sp.integrate(sp.conjugate(profile)*profile, (y, -sp.pi, sp.pi))/(2*sp.pi)
    assert sp.simplify(energy-1-eps**2/2) == 0
    assert sp.simplify(energy-energy.subs(eps, -eps)) == 0
    with pytest.raises(ValueError, match="real"):
        sinusoidal_shear_memory(sp.I, 1)

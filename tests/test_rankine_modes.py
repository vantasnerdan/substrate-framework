"""Euler equations and amplitude elimination, independent of branch fitting."""
import mpmath
import pytest
import sympy as s

from substrate_framework.rankine_modes import (
    boundary_determinant,
    core_velocity,
    rankine_residual,
)


def test_velocity_matches_cartesian_advection_with_polar_basis():
    r, theta, z, t = s.symbols("r theta z t", real=True)
    Om, rho, m, k, w = s.symbols("Omega rho m k w", nonzero=True)
    pressure = s.Function("P")(r)
    er = s.Matrix([s.cos(theta), s.sin(theta), 0])
    et = s.diff(er, theta)
    ez = s.Matrix([0, 0, 1])
    phase = s.exp(s.I*(m*theta+k*z-w*t))
    vr, vt, vz = core_velocity(pressure, r, k, m, w-m*Om, Om, rho)
    field = phase*(vr*er+vt*et+vz*ez)
    background = Om*r*et
    acceleration = (s.diff(field, t)+Om*s.diff(field, theta)
                    +phase*(vr*s.diff(background, r)+vt/r*s.diff(background, theta)))
    gradient = phase*(s.diff(pressure, r)*er+s.I*m*pressure/r*et+s.I*k*pressure*ez)
    assert s.simplify(acceleration+gradient/rho) == s.zeros(3, 1)


def test_boundary_determinant_follows_from_pressure_amplitude_match():
    wt, Om, m, J, K, rho, eta, a = s.symbols("wt Om m J K rho eta a", nonzero=True)
    # J=a*Pin'/Pin; K=a*potential'/potential outside. Radial kinematics.
    pin = -wt*eta*rho*a*(4*Om**2-wt**2)/(wt*J-2*Om*m)
    pout = rho*wt**2*eta*a/K
    scaled_jump = s.factor((pin-pout)*(-K*(wt*J-2*Om*m)/(rho*wt*eta*a)))
    assert s.simplify(scaled_jump-boundary_determinant(wt, Om, m, J, K)) == 0


def test_bessel_residual_respects_precision_context_and_neutral_limit():
    before = mpmath.mp.dps
    context = mpmath.mp.clone()
    context.dps = 50
    values = [abs(rankine_residual(context, x, 1, -1)) for x in ("1e-4", "1e-6")]
    assert values[1] < values[0]/1000
    assert mpmath.mp.dps == before


@pytest.mark.parametrize("x,m,doppler", [(0, 1, -1), (1, 0, -1), (1, 1, 0), (1, 1, -2)])
def test_bessel_residual_rejects_outside_declared_branch(x, m, doppler):
    with pytest.raises(ValueError):
        rankine_residual(mpmath.mp, x, m, doppler)

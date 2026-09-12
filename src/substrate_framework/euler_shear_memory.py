"""Exact unpromoted Euler shear histories and their resolved memory (P253/0009).

The triangular, z-independent family is a full finite-amplitude Euler
solution, not an approximation to another trajectory. Periodic inputs give
finite energy per cell; isolated localization and generic 3D stability are
not supplied. The Bessel response is derived from U_s*sin(y/L_y), observed
at one x Fourier wave, including an explicit unresolved initial sine mode.
It is not a prescribed oscillator, quantum amplitude or particle dispersion.
"""

from dataclasses import dataclass

import sympy as sp


def triangular_euler_velocity(shear, initial_vertical, coordinates, time):
    """Return (U(y),0,w0(x-t*U(y),y)) with constant Euler pressure.

    Coordinates and time are distinct symbols. U and w0 are independent of
    time and z; U is also x-independent. Smoothness, reality and periodicity
    of supplied profiles are caller hypotheses where undecidable here.
    """
    if len(coordinates) != 3:
        raise ValueError("three coordinate symbols are required")
    x, y, z = coordinates
    if not all(isinstance(q, sp.Symbol) for q in (x, y, z, time)) or len({x, y, z, time}) != 4:
        raise ValueError("coordinates and time must be distinct symbols")
    U, w0 = sp.sympify(shear), sp.sympify(initial_vertical)
    if any(q in U.free_symbols for q in (x, z, time)):
        raise ValueError("shear must be independent of x, z and time")
    if any(q in w0.free_symbols for q in (z, time)):
        raise ValueError("initial vertical profile must be independent of z and time")
    return sp.ImmutableMatrix([U, 0, w0.subs(x, x-time*U)])


@dataclass(frozen=True)
class EulerShearMemory:
    """Observed amplitude and exact Volterra kernel/initial-state forcing."""

    amplitude: sp.Expr
    kernel: sp.Expr
    unresolved_forcing: sp.Expr


def sinusoidal_shear_memory(frequency, time, *, unresolved_sine=0):
    """Return exact memory for g0(y)=1+i*epsilon*sin(y).

    beta=k*U_s is real and has inverse-time units. The y-average obeys
    a'=F-integral_0^t K(t-s)*a(s)ds, a(0)=1. K(0)=beta^2/2 and
    F(0)=epsilon*beta/2 are the removable limits. Epsilon=0 is one
    explicit initial field, not a universal absence of unresolved state.
    No numerical approximation or history fitting is performed.
    """
    beta, t, epsilon = map(sp.sympify, (frequency, time, unresolved_sine))
    if any(q.is_real is False or q.is_finite is False for q in (beta, t, epsilon)):
        raise ValueError("frequency, time and unresolved amplitude must be finite real")
    if beta.is_zero is True:
        return EulerShearMemory(sp.S.One, sp.S.Zero, sp.S.Zero)
    if t.is_zero is True:
        return EulerShearMemory(sp.S.One, beta**2/2, epsilon*beta/2)
    regular_ratio = sp.Piecewise((beta/2, sp.Eq(t, 0)), (sp.besselj(1, beta*t)/t, True))
    return EulerShearMemory(sp.besselj(0, beta*t)+epsilon*sp.besselj(1, beta*t),
                            beta*regular_ratio, epsilon*regular_ratio)

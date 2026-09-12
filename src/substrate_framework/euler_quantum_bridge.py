"""Exact algebraic boundaries for Euler-orbit quantization candidates.

The functions encode the finite-dimensional calculations in P253/0043. They
do not select an action unit, a deck character, a measurement rule, or a
persistent Euler carrier.
"""

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class SpinOrbitQuantization:
    """Conditional geometric quantization data for the physical KKS sphere."""

    cern_number: sp.Expr
    hilbert_dimension: sp.Expr
    spin: sp.Expr
    rotation_phase: sp.Expr


def spin_orbit_quantization(angular_momentum, action_unit):
    """Return ``N=2*j/hbar`` data for a nondegenerate positive KKS sphere."""
    j = sp.sympify(angular_momentum)
    hbar = sp.sympify(action_unit)
    if j.is_positive is not True:
        raise ValueError("angular_momentum must be positive for the KKS sphere")
    if hbar.is_positive is not True:
        raise ValueError("action_unit must be positive")
    number = sp.simplify(2 * j / hbar)
    if number.is_integer is not True or number.is_nonnegative is not True:
        raise ValueError("2*angular_momentum/action_unit must be a nonnegative integer")
    return SpinOrbitQuantization(
        cern_number=number,
        hilbert_dimension=number + 1,
        spin=number / 2,
        rotation_phase=sp.Integer(-1) ** number,
    )


def scaled_euler_angular_momentum(angular_momentum, velocity_scale, inverse_length_scale):
    """Scale ``j`` under ``u_AB(x,t)=A*u(B*x,A*B*t)`` at fixed density."""
    j = sp.sympify(angular_momentum)
    a = sp.sympify(velocity_scale)
    b = sp.sympify(inverse_length_scale)
    if a.is_positive is not True or b.is_positive is not True:
        raise ValueError("Euler scales must be positive")
    return sp.simplify(a * j / b**4)


def exchange_characters():
    """Return both unitary characters of the two-carrier ``Z2`` exchange loop."""
    return (sp.Integer(1), sp.Integer(-1))


def hopf_fr_phase(hopf_charge, *, deck_character):
    """Return the conditional FR phase ``chi**(Q mod 2)``.

    The caller supplies the deck character; topology permits both signs and
    this function deliberately does not choose one.  The source proves this
    class for spatial rotation.  Using it for exchange additionally assumes
    that the physical exchange loop maps to the same Hopf-space class.
    """
    charge = sp.sympify(hopf_charge)
    character = sp.sympify(deck_character)
    if charge.is_integer is not True:
        raise ValueError("hopf_charge must be an integer")
    if character not in (sp.Integer(1), sp.Integer(-1)):
        raise ValueError("deck_character must be +1 or -1")
    return character ** (charge % 2)


def uniform_euler_frequency(wavevector, background_velocity):
    """Return the bare incompressible Euler convective frequency ``U dot k``."""
    k = sp.Matrix(wavevector)
    velocity = sp.Matrix(background_velocity)
    if k.cols != 1 or velocity.cols != 1 or k.rows != velocity.rows or k.rows == 0:
        raise ValueError("wavevector and background_velocity must be equal nonempty vectors")
    return sp.simplify((velocity.T * k)[0])

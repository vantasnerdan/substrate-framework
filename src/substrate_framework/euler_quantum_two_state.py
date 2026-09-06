"""Exact finite-dimensional tests for Euler quantum two-state candidates.

These utilities expose the classical KKS plane supplied by the accepted
compact-ring response construction.  They do not quantize Euler, select an
action unit, or provide a measurement rule.
"""

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class KKSPlaneDynamics:
    """Classical Hamiltonian data on one real two-dimensional KKS plane."""

    symplectic_matrix: sp.ImmutableMatrix
    hessian: sp.ImmutableMatrix
    generator: sp.ImmutableMatrix
    frequency_squared: sp.Expr


def kks_plane_dynamics(symplectic_scale, hessian):
    """Return the linear flow for ``Omega=B*dq wedge dp`` and ``i_X Omega=dH``."""
    scale = sp.sympify(symplectic_scale)
    matrix = sp.ImmutableMatrix(hessian)
    if scale.is_zero is True:
        raise ValueError("symplectic_scale must be nonzero")
    if matrix.shape != (2, 2) or matrix != matrix.T:
        raise ValueError("hessian must be a symmetric 2 by 2 matrix")
    omega = sp.ImmutableMatrix([[0, scale], [-scale, 0]])
    generator = sp.ImmutableMatrix(-omega.inv() * matrix)
    return KKSPlaneDynamics(
        symplectic_matrix=omega,
        hessian=matrix,
        generator=generator,
        frequency_squared=sp.simplify(matrix.det() / scale**2),
    )


def affine_kks_poisson_bracket(first_gradient, second_gradient, symplectic_scale):
    """Return the bracket of affine functions in ``(q,p)`` coordinates.

    The convention gives ``{q,B*p}=1`` for ``Omega=B*dq wedge dp``.
    """
    first = sp.Matrix(first_gradient)
    second = sp.Matrix(second_gradient)
    scale = sp.sympify(symplectic_scale)
    if first.shape != (2, 1) or second.shape != (2, 1):
        raise ValueError("affine gradients must have two components")
    if scale.is_zero is True:
        raise ValueError("symplectic_scale must be nonzero")
    canonical = sp.Matrix([[0, 1], [-1, 0]])
    return sp.simplify((first.T * canonical * second)[0] / scale)


def finite_ccr_trace_obstruction(dimension, action_unit):
    """Return ``tr([Q,P]-i*hbar*I)`` for any finite square matrices ``Q,P``."""
    dim = sp.sympify(dimension)
    hbar = sp.sympify(action_unit)
    if dim.is_integer is not True or dim.is_positive is not True:
        raise ValueError("dimension must be a positive integer")
    if hbar.is_positive is not True:
        raise ValueError("action_unit must be positive")
    return sp.simplify(-sp.I * hbar * dim)


def euler_circulation_scale(circulation, velocity_scale, inverse_length_scale):
    """Scale circulation under ``u_AB(x,t)=A*u(B*x,A*B*t)``."""
    gamma = sp.sympify(circulation)
    a = sp.sympify(velocity_scale)
    b = sp.sympify(inverse_length_scale)
    if a.is_positive is not True or b.is_positive is not True:
        raise ValueError("Euler scales must be positive")
    return sp.simplify(a * gamma / b)

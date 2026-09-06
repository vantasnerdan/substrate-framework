"""Exact pressure-source and bounded linear retained-memory identities.

The nonlinear Euler retained-state/history theorem is derived in P253/0042.
This module exposes two reusable algebraic pieces.  It does not solve Euler,
construct an orthogonal-dynamics semigroup, or close material moments by
themselves.
"""

from dataclasses import dataclass

import sympy as sp


def euler_pressure_poisson_source(velocity_gradient, *, density):
    """Return the source in ``-Delta p = rho*tr((Du)**2)``.

    ``velocity_gradient[i,j]`` is ``partial_j u_i``.  Incompressibility and
    same-field membership are caller hypotheses; this function retains the
    contraction with the transposed index order rather than a Frobenius norm.
    With the positive inverse convention, ``p=(-Delta)^-1*source``.
    """
    gradient = sp.Matrix(velocity_gradient)
    if gradient.rows != gradient.cols or gradient.rows == 0:
        raise ValueError("velocity_gradient must be a nonempty square matrix")
    rho = sp.sympify(density)
    if (rho.is_positive is False or rho.is_real is False
            or rho.is_finite is False
            or (rho.is_number and rho.is_positive is not True)):
        raise ValueError("density must be positive finite real")
    return sp.simplify(rho * sp.trace(gradient * gradient))


@dataclass(frozen=True)
class LinearRetainedMemory:
    """Markov, unresolved-initial, and convolution operators for ``x'=Lx``."""

    markov: sp.ImmutableMatrix
    unresolved_propagator: sp.ImmutableMatrix
    memory_kernel: sp.ImmutableMatrix


def linear_retained_memory(generator, projection, time):
    """Eliminate the complementary block of a bounded linear generator.

    With ``y=P*x``, ``z=Q*x`` and ``Q=I-P``, the returned matrices give

    ``y' = markov*y + unresolved_propagator*z0
           + integral memory_kernel(t-s)*y(s) ds``.

    The formula is exact for finite matrices.  It is an oracle for the Dyson
    signs and the bounded shear specialization, not a claim that generic
    infinite-dimensional Euler orthogonal dynamics is already generated.
    """
    operator = sp.Matrix(generator)
    resolved = sp.Matrix(projection)
    if operator.rows != operator.cols or operator.rows == 0:
        raise ValueError("generator must be a nonempty square matrix")
    if resolved.shape != operator.shape:
        raise ValueError("projection must have the generator shape")
    if resolved.applyfunc(sp.simplify) != (resolved * resolved).applyfunc(sp.simplify):
        raise ValueError("projection must be idempotent")
    t = sp.sympify(time)
    if t.is_real is False or t.is_finite is False:
        raise ValueError("time must be finite real")
    identity = sp.eye(operator.rows)
    complement = identity - resolved
    unresolved = complement * operator * complement
    propagator = (t * unresolved).exp()
    markov = resolved * operator * resolved
    noise = resolved * operator * complement * propagator * complement
    kernel = noise * operator * resolved
    return LinearRetainedMemory(*(
        sp.ImmutableMatrix(value.applyfunc(sp.simplify))
        for value in (markov, noise, kernel)
    ))

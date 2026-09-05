"""Conditional time-dependent symplectic pullbacks, P251 stronger frontier.

These exact algebraic APIs are unpromoted infrastructure. They keep moving
frame terms and physical observations explicit; they do not construct an
Euler invariant phase family or establish a Floquet energy's physical sign.
The ambient symplectic form is constant. In the convention Omega(q,s)=B,
the action one-form is -x.T*Omega*xdot/2, equivalent to B*s*qdot.
"""

from dataclasses import dataclass

import sympy as sp


def _matrix(value, name):
    result = sp.Matrix(value)
    if any(x.is_finite is False or x.has(sp.nan, sp.zoo) for x in result):
        raise ValueError(f"{name} must have finite entries")
    return result


def _zero(matrix):
    return all(sp.simplify(x) == 0 for x in matrix)


def _immutable(matrix):
    return sp.ImmutableMatrix(matrix.applyfunc(sp.simplify))


@dataclass(frozen=True)
class MovingPhasePullback:
    """Same-action forms, generator, projection and ambient residual.

    The coordinate equation is Omega_E*zdot +
    (H_eff + dotOmega_E/2)*z = 0. ``residual`` is
    Edot + E*generator - ambient_generator*E, in ambient coordinates.
    It is symplectically orthogonal to E, not necessarily zero or gauge.
    An actual observable is O*E*z plus its complementary observation.
    """

    symplectic: sp.ImmutableMatrix
    symplectic_rate: sp.ImmutableMatrix
    hamiltonian: sp.ImmutableMatrix
    generator: sp.ImmutableMatrix
    coordinates: sp.ImmutableMatrix
    projection: sp.ImmutableMatrix
    residual: sp.ImmutableMatrix


def moving_phase_pullback(symplectic, hamiltonian, embedding, embedding_rate):
    """Derive a differentiable phase restriction from its complete action.

    Inputs Omega,H,E,Edot are actual simultaneous matrices, with Omega
    invertible and skew, H symmetric, and E.T*Omega*E invertible. The
    caller supplies Edot as the actual time derivative of E and retains
    any undecidable symbolic nondegeneracy hypotheses. No positivity or
    microscopic realization is inferred. Complex bilinear representations
    may be used; a real physical encoding remains a caller obligation.
    """
    omega = _matrix(symplectic, "symplectic form")
    h = _matrix(hamiltonian, "Hamiltonian")
    e = _matrix(embedding, "embedding")
    ed = _matrix(embedding_rate, "embedding rate")
    if omega.rows != omega.cols or not _zero(omega + omega.T):
        raise ValueError("ambient symplectic form must be square and skew")
    if h.shape != omega.shape or not _zero(h - h.T):
        raise ValueError("Hamiltonian must be symmetric with the ambient shape")
    if e.rows != omega.rows or ed.shape != e.shape or not e.cols:
        raise ValueError("embedding and its rate must have compatible nonempty shapes")
    reduced_omega = sp.simplify(e.T * omega * e)
    if sp.simplify(omega.det()).is_zero is True:
        raise ValueError("ambient symplectic form must be nondegenerate")
    if sp.simplify(reduced_omega.det()).is_zero is True:
        raise ValueError("restricted symplectic form must be nondegenerate")
    connection = e.T * omega * ed
    omega_rate = connection - connection.T
    effective_h = e.T * h * e + (connection + connection.T) / 2
    coordinates = reduced_omega.inv() * e.T * omega
    ambient_generator = -omega.inv() * h
    generator = -reduced_omega.inv() * (effective_h + omega_rate / 2)
    projection = e * coordinates
    residual = ed + e * generator - ambient_generator * e
    return MovingPhasePullback(*map(_immutable, (
        reduced_omega, omega_rate, effective_h, generator,
        coordinates, projection, residual)))

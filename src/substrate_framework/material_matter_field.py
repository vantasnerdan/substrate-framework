"""Exact finite quadratic matter/field blocks for material-response reductions.

These definitions do not assign a Ni effective charge, derive a nonlinear phase,
or supply a nuclear transition. Material parameters must be earned separately.
"""

from __future__ import annotations

import sympy as sp


def minimal_coupling_hessian(mass, stiffness, field_stiffness, coupling):
    """Hessian of a single matter coordinate and a canonical field coordinate.

    Coordinate order is ``(u, p, A, Pi)`` and
    ``H=(p-coupling*A)**2/(2*mass) + stiffness*u**2/2
         + Pi**2/2 + field_stiffness*A**2/2``.
    Field normalization is canonical, not an assumed dielectric/effective charge.
    The square retains its complete field self-term.
    """
    m, k, w2, c = map(sp.sympify, (mass, stiffness, field_stiffness, coupling))
    u, p, a, pi = sp.symbols("u p A Pi", real=True)
    energy = (p - c * a) ** 2 / (2 * m) + k * u**2 / 2
    energy += (pi**2 + w2 * a**2) / 2
    return sp.ImmutableMatrix(sp.hessian(energy, (u, p, a, pi)))


def static_matched_bath_hessian(retained_stiffness, bath_stiffness, coupling):
    """Full mass-weighted stiffness with the specified relaxed static response.

    ``coupling`` has retained rows and bath columns. Completing the bath square
    adds ``C*B**-1*C.T`` to the retained block. Eliminating the bath at zero
    frequency then returns precisely ``retained_stiffness``. This definition
    applies to finite nonsingular bath blocks; zero modes need a licensed gauge.
    """
    k = sp.Matrix(retained_stiffness)
    b = sp.Matrix(bath_stiffness)
    c = sp.Matrix(coupling)
    if k.rows != k.cols or b.rows != b.cols or not k.rows or not b.rows:
        raise ValueError("stiffness blocks must be nonempty square matrices")
    if c.shape != (k.rows, b.rows):
        raise ValueError("coupling must have retained rows and bath columns")
    correction = c * b.inv() * c.T
    return sp.ImmutableMatrix((k + correction).row_join(c).col_join(
        c.T.row_join(b)))


def eliminate_bath(stiffness, retained_size, frequency):
    """Exact retained inverse response of a finite mass-weighted harmonic system.

    Returns the Schur complement of ``K-z**2*I``. Retarded response requires
    ``Im(z)>0`` before taking a boundary limit. This operation alone does not
    establish a continuum bath or physical irreversible decay.
    """
    matrix = sp.Matrix(stiffness)
    n = int(retained_size)
    if matrix.rows != matrix.cols or not 0 < n < matrix.rows:
        raise ValueError("retained_size must split a square stiffness matrix")
    z = sp.sympify(frequency)
    inverse_response = matrix - z**2 * sp.eye(matrix.rows)
    a = inverse_response[:n, :n]
    c = inverse_response[:n, n:]
    d = inverse_response[n:, :n]
    b = inverse_response[n:, n:]
    return sp.ImmutableMatrix(a - c * b.inv() * d)

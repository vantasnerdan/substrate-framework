"""Exact finite-dimensional utilities for the P253/0079 Euler analyzer test.

The functions are symbolic and contain no model-specific fitted constants.
"""

from __future__ import annotations

import sympy as sp


def standard_complex_structure(two_mode_dimension: int = 2) -> sp.Matrix:
    """Real matrix for multiplication by i on C^n, ordered (Re z, Im z)."""
    eye = sp.eye(two_mode_dimension)
    zero = sp.zeros(two_mode_dimension)
    return zero.row_join(-eye).col_join(eye.row_join(zero))


def complex_linear_split(v: sp.Matrix, j: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix]:
    """Return V_C=(V-JVJ)/2 and V_A=(V+JVJ)/2."""
    if v.rows != v.cols or j.shape != v.shape:
        raise ValueError("V and J must be square matrices of the same size")
    vc = sp.simplify((v - j * v * j) / 2)
    va = sp.simplify((v + j * v * j) / 2)
    return vc, va


def realify_complex(a: sp.Matrix) -> sp.Matrix:
    """Realification of a complex matrix acting complex linearly."""
    re = a.applyfunc(sp.re)
    im = a.applyfunc(sp.im)
    return re.row_join(-im).col_join(im.row_join(re))


def pauli_coordinates(h: sp.Matrix) -> tuple[sp.Expr, sp.Matrix]:
    """Return h0 and the Pauli vector of a 2x2 Hermitian matrix."""
    if h.shape != (2, 2):
        raise ValueError("Pauli coordinates require a 2x2 matrix")
    h0 = sp.simplify((h[0, 0] + h[1, 1]) / 2)
    hx = sp.simplify((h[0, 1] + h[1, 0]) / 2)
    hy = sp.simplify((h[1, 0] - h[0, 1]) / (2 * sp.I))
    hz = sp.simplify((h[0, 0] - h[1, 1]) / 2)
    return h0, sp.Matrix([hx, hy, hz])


def pauli_noncollinearity(h1: sp.Matrix, h2: sp.Matrix) -> sp.Matrix:
    """Cross product of the traceless Pauli vectors."""
    _, p1 = pauli_coordinates(h1)
    _, p2 = pauli_coordinates(h2)
    return sp.simplify(p1.cross(p2))


def coupled_gate_bounds(
    eps_a: sp.Expr, eps_out: sp.Expr, kappa_back: sp.Expr
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """Bogoliubov, P gain, Q leakage, and absolute P-action bounds.

    ``eps_out`` controls the one-way P-to-Q Duhamel kernel.  The distinct
    ``kappa_back`` is the norm of the Q-return Volterra kernel and must be
    strictly below one when the formula is evaluated as an estimate.
    """
    beta = sp.sinh(eps_a)
    p_gain = sp.exp(eps_a) / (1 - kappa_back)
    q_gain = eps_out * p_gain
    action_change = p_gain**2 - 1
    return beta, p_gain, q_gain, sp.simplify(action_change)


def rational_ray_crossing_shift(
    k_star: sp.Expr,
    p: sp.Expr,
    large_n: sp.Expr,
    delta_c: sp.Expr,
    f_prime: sp.Expr,
) -> tuple[sp.Expr, sp.Expr]:
    """Optional C1 refinement: leading parameter and O(N^-2) correction."""
    delta_0 = k_star / (large_n * p)
    shift = -k_star * delta_c / (large_n**2 * p**2 * f_prime)
    return delta_0, shift


def massive_bracketing_interval(
    k_star: sp.Expr, p: sp.Expr, large_n: sp.Expr, error_envelope: sp.Expr
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """C0-IVT center and endpoints with h=sqrt(error envelope)."""
    h = sp.sqrt(error_envelope)
    scale = large_n * p
    return k_star / scale, (k_star - h) / scale, (k_star + h) / scale

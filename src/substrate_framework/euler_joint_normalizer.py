"""Finite form algebra used by the proposed Euler normalizer.

These coefficient matrices require an independently constructed physical
common null frame and its duals. They do not construct an Euler field or
assert that arbitrary supplied forms occur on a given background.
"""

from __future__ import annotations

import sympy as sp


def common_null_gram_lift(energy, phase, scale=1):
    """Return coefficients in the ordered basis ``(e+, e-, p, q)``.

    The basis has energy matrix ``diag(I, -I, 0, 0)`` and phase matrix
    with only ``Omega(p,q)=I`` and ``Omega(q,p)=-I`` nonzero. The returned
    columns have exactly the Hermitian energy and skew-Hermitian phase
    matrices supplied by the caller. A positive real scale changes the
    coefficients without changing either Gram matrix.

    Symbolic assumptions not decidable by SymPy remain caller hypotheses.
    The return value is immutable; no spectral inversion is performed.
    """
    energy = sp.ImmutableMatrix(energy)
    phase = sp.ImmutableMatrix(phase)
    scale = sp.sympify(scale)
    n = energy.rows
    if n == 0 or energy.cols != n or phase.shape != (n, n):
        raise ValueError("energy and phase must be nonempty square matrices of equal size")
    if energy.equals(energy.conjugate().T) is False:
        raise ValueError("energy must be Hermitian")
    if phase.equals(-phase.conjugate().T) is False:
        raise ValueError("phase must be skew-Hermitian")
    if scale.is_real is False or scale.is_positive is False:
        raise ValueError("scale must be positive and real")
    identity = sp.eye(n)
    return sp.ImmutableMatrix.vstack(
        scale * identity + energy / (4 * scale),
        scale * identity - energy / (4 * scale),
        scale * identity,
        phase / (2 * scale),
    )

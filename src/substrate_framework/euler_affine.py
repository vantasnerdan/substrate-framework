"""Conditional, unpromoted affine Euler energy (P251 attempt 0043).

The vorticity and periodic lattice are pushed forward together by a real
volume-preserving affine map. The energy is obtained by physical curl
inversion on that lattice. This is a constrained energy calculation, not
an assertion that an affine image of a stationary field remains stationary.
It supplies no independent bulk modulus, material inertia, or cell-gluing
theorem. Fourier coefficients use the cell-average Parseval convention;
real fields include both members of each conjugate wave pair. Harmonic
velocity energy is excluded and is supplied separately when present.
"""

from __future__ import annotations

import sympy as sp


def affine_vorticity_energy(deformation, vorticity_modes, density):
    """Return the exact Euler kinetic energy density on the deformed lattice.

    ``vorticity_modes`` maps real three-dimensional nonzero wave covectors
    to complex three-component vorticity Fourier coefficients. Each mode
    is solenoidal. ``deformation`` is a real 3x3 matrix with determinant one.
    All equalities must be symbolically verifiable; positivity of symbolic
    density is a caller hypothesis unless its sign is decidable here.
    Coordinates and wave covectors are physical (not silently rescaled).
    """
    f, rho = sp.Matrix(deformation), sp.sympify(density)
    if f.shape != (3, 3) or sp.simplify(f.det()-1) != 0:
        raise ValueError("deformation must be a 3x3 determinant-one matrix")
    if any(value.is_real is False for value in f):
        raise ValueError("deformation must be real")
    if rho.is_positive is False or (rho.is_number and rho.is_positive is not True):
        raise ValueError("density must be positive")
    inverse_transpose = f.inv().T
    total = sp.S.Zero
    for wave, coefficient in vorticity_modes.items():
        k, omega = sp.Matrix(wave), sp.Matrix(coefficient)
        if k.shape != (3, 1) or omega.shape != (3, 1):
            raise ValueError("waves and coefficients must have three components")
        if any(value.is_real is False for value in k):
            raise ValueError("wave covectors must be real")
        if sp.simplify(k.dot(k)) == 0:
            raise ValueError("zero vorticity modes have no periodic curl inverse")
        if sp.simplify(k.dot(omega)) != 0:
            raise ValueError("vorticity modes must be solenoidal")
        transformed_wave = inverse_transpose*k
        transformed_omega = f*omega
        total += (sp.conjugate(transformed_omega).dot(transformed_omega)
                  / transformed_wave.dot(transformed_wave))
    return sp.simplify(rho*total/2)

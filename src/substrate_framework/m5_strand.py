"""Exact conditional tension of a fixed-slot, straight M5 transverse-pair strand.

Two *different* source functionals share the pair block
``S=s0 I+b(rho)(cos(m phi) sigma_z+sin(m phi) sigma_x)``.  The
Report-018 field-Euclidean contraction is ``2 Tr(F_xy F_xy.T)``; the
OpenWave R25 trace-power action uses ``4 Tr(F_xy F_xy.T)`` on this
positive spatial chart.  Their potentials also differ.  Neither action
is identified here with the accepted auxiliary-frame M5 action of P239.

Every tension below fixes all other eigenvalues and the pair mean, pins
``b(0)=0`` and ``b(infinity)=b0>0``, and only varies the radial amplitude.
This restricted infimum is not a lower bound on the full field theory.
"""

from __future__ import annotations

from typing import cast

import sympy as sp

Scalar = sp.Expr | int | float


def pair_block(
    amplitude: sp.Expr, angle: Scalar, winding: Scalar, mean: Scalar
) -> sp.ImmutableMatrix:
    """Return the real symmetric transverse pair block; ``winding`` is integral."""

    phase = winding * angle
    cosine, sine = sp.cos(phase), sp.sin(phase)
    return sp.ImmutableMatrix(
        [[mean + amplitude * cosine, amplitude * sine],
         [amplitude * sine, mean - amplitude * cosine]]
    )


def pair_curvature_xy(
    amplitude: sp.Expr,
    radius: sp.Symbol,
    angle: sp.Symbol,
    winding: Scalar,
    mean: Scalar,
) -> sp.ImmutableMatrix:
    """Compute the spatial commutator ``[D_x S,D_y S]`` from polar derivatives."""

    block = pair_block(amplitude, angle, winding, mean)
    radial, angular = block.diff(radius), block.diff(angle)
    commutator = cast(sp.MatrixBase, radial * angular - angular * radial)
    return sp.ImmutableMatrix(
        [[sp.trigsimp(commutator[i, j] / radius) for j in range(2)] for i in range(2)]
    )


def spectral_potential(q: sp.Expr, b0: Scalar) -> sp.Expr:
    """Report-018 potential on ``0 <= q <= b0**2``: eigenvalue distance."""

    return 2 * (b0 - sp.sqrt(q)) ** 2


def trace_stiffness(q: sp.Expr, q0: Scalar, mean: Scalar) -> sp.Expr:
    """Exact positive R25 trace-power stiffness on the fixed-slot chart."""

    return 4 + 36 * mean**2 + (12 * mean**2 + 2 * (q + q0)) ** 2


def trace_potential(
    q: sp.Expr, q0: Scalar, mean: Scalar, weight: Scalar
) -> sp.Expr:
    """OpenWave R25 sum of squared trace differences, with fixed other slots."""

    return weight * (q0 - q) ** 2 * trace_stiffness(q, q0, mean)


def trace_director_first_variation(
    q: sp.Expr, q0: Scalar, mean: Scalar, weight: Scalar
) -> sp.Expr:
    """Derivative of R25's trace potential at a free 1-eigenvalue d=1.

    Other exterior slots stay fixed.  For ``0 <= q < q0`` and positive
    mean/weight this is negative: a positive compactly supported change in
    d lowers the potential, while its diagonal gradient commutes with pair
    gradients and does not change the curvature of this spatial chart.
    """

    return 2 * weight * (q - q0) * (
        4 + 18 * mean + 4 * (12 * mean**2 + 2 * (q + q0))
    )


def radial_tension_infimum(
    curvature_coefficient: Scalar,
    winding: Scalar,
    potential: sp.Expr,
    q: sp.Symbol,
    q0: Scalar,
) -> sp.Expr:
    """Sharp radial *infimum* ``2pi sqrt(A)|m| int_0^q0 sqrt(V) dq``.

    Assume A>0, integer |m|>=1, 0<=q<=q0, V(q)>=0, V(q0)=0 and
    V(q)>0 for q<q0.  The radial density in s=rho**2/2 is
    ``A*m**2*q_s**2/4 + V(q)``; q(0)=0 and q(infinity)=q0.
    Square completion proves the bound even for nonmonotone profiles.
    If the first-order solution has the specified asymptotic endpoint,
    its energy realizes this infimum in the finite-energy class.  In the
    smooth matrix-field class winding |m|>1 generally has *no smooth
    minimizer*: the first-order core b~rho is angularly nondifferentiable;
    smooth cores b=O(rho**|m|) approach the same infimum by shrinking the
    smoothing disk.  A symbolic Integral avoids silently integrating an
    unrelated potential or promoting this bound to the full field theory.
    """

    return 2 * sp.pi * sp.sqrt(curvature_coefficient) * sp.Abs(winding) * sp.Integral(
        sp.sqrt(potential), (q, 0, q0)
    )


def spectral_strand_infimum(beta: Scalar, winding: Scalar) -> sp.Expr:
    """Report-018: ``sqrt(2)*pi*|m|*beta**3/3``, ``m=2k``.

    A=16 is the coefficient of ``m**2 (b b_r/rho)**2`` in that action.
    For beta>0 this is sharp among radial finite-energy profiles, not
    necessarily attained among smooth fields at the origin.
    """

    return sp.sqrt(2) * sp.pi * sp.Abs(winding) * beta**3 / 3


def trace_strand_infimum(
    b0: Scalar, mean: Scalar, weight: Scalar, winding: Scalar
) -> sp.Expr:
    """R25 exact variable-K radial infimum; delta=2*b0 when mean=b0.

    A=32 and weight>0.  The often-quoted constant ``K_leading`` gives
    only a lower bound on this exact integral at finite splitting.
    """

    q = sp.Symbol("q", real=True, nonnegative=True)
    q0 = b0**2
    return 2 * sp.pi * sp.sqrt(32 * weight) * sp.Abs(winding) * sp.Integral(
        (q0 - q) * sp.sqrt(trace_stiffness(q, q0, mean)), (q, 0, q0)
    )


def trace_strand_leading_bound(
    b0: Scalar, mean: Scalar, weight: Scalar, winding: Scalar
) -> sp.Expr:
    """R25 looser bound using ``K_leading=4+36*mean**2+144*mean**4``.

    The pinned exterior in R25 has ``mean=b0``.  The leading bound is
    strictly weaker than the variable-K infimum for b0>0, weight>0.
    """

    leading = 4 + 36 * mean**2 + 144 * mean**4
    return sp.pi * sp.sqrt(32 * weight * leading) * sp.Abs(winding) * b0**4

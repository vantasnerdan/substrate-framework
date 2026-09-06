"""Exact algebra for finite-energy degree-minus-two Euler velocity tails.

The helpers expose only the asymptotic harmonic and cross-energy identities.
They do not certify that a tail is a persistent particle or a stationary Euler
solution.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

import sympy as sp


@dataclass(frozen=True)
class TailDomainLedger:
    """Functional-domain consequences of a nonzero ``r**-2`` tail."""

    sector: str
    finite_kinetic_energy: bool
    absolutely_integrable_velocity: bool
    absolutely_integrable_vorticity: bool
    finite_ordinary_angular_momentum: bool
    finite_absolute_vorticity_impulse: bool
    finite_helicity_at_infinity: bool


def scalar_homogeneous_fourier_coefficient(l: int, degree: int) -> sp.Expr:
    """Coefficient in ``FT[r**(-degree) Y_l]`` for the stated convention.

    The transform is

    ``c(l,degree) |k|**(degree-3) Y_l(k_hat)``

    with ``fhat(k)=integral exp(-I*k*x) f(x) dx``.  The degree-three formula
    is used only for ``l>=1``; its ``l=0`` member needs a logarithmic/contact
    extension and is outside the source-free tail class.
    """

    if not isinstance(l, int) or l < 0:
        raise ValueError("l must be a nonnegative integer")
    if degree not in (1, 2, 3):
        raise ValueError("degree must be 1, 2, or 3")
    if degree == 3 and l == 0:
        raise ValueError("the degree-three l=0 distribution needs a separate extension")
    return sp.simplify(
        sp.I ** (-l)
        * 2 ** (3 - degree)
        * sp.pi ** sp.Rational(3, 2)
        * sp.gamma(sp.Rational(l + 3 - degree, 2))
        / sp.gamma(sp.Rational(l + degree, 2))
    )


def radial_tail_fourier_coefficient(l: int) -> sp.Expr:
    """Coefficient of ``|k|^-1 grad_S Y_l`` for ``r^-2 Y_l n``."""

    if l < 1:
        raise ValueError("the source-free radial tail has l>=1")
    return sp.simplify(sp.I * scalar_homogeneous_fourier_coefficient(l, 3))


def toroidal_tail_fourier_coefficient(l: int) -> sp.Expr:
    """Coefficient of ``|k|^-1 k_hat cross grad_S Y_l``."""

    if l < 1:
        raise ValueError("toroidal vector harmonics start at l=1")
    return scalar_homogeneous_fourier_coefficient(l, 2)


def oriented_tail_cross_kernel(
    a: Sequence[sp.Expr],
    b: Sequence[sp.Expr],
    separation_direction: Sequence[sp.Expr],
    *,
    density: sp.Expr = sp.Integer(1),
    separation: sp.Expr = sp.Integer(1),
) -> sp.Expr:
    """Leading kinetic cross energy for the smooth-core oriented tail.

    The tail is ``u_a ~ a cross x / |x|^3``.  The returned expression is
    ``2*pi*density/separation * (a.b + (a.n)(b.n))`` and assumes ``n.n=1``.
    """

    av = sp.ImmutableMatrix(a)
    bv = sp.ImmutableMatrix(b)
    nv = sp.ImmutableMatrix(separation_direction)
    if av.shape != (3, 1) or bv.shape != (3, 1) or nv.shape != (3, 1):
        raise ValueError("a, b, and separation_direction must have three components")
    if sp.simplify(separation) == 0:
        raise ValueError("separation must be nonzero")
    return sp.simplify(
        2
        * sp.pi
        * density
        / separation
        * ((av.dot(bv)) + (av.dot(nv)) * (bv.dot(nv)))
    )


def fixed_frame_constant_norm_symbol(
    direction: Sequence[sp.Expr],
) -> sp.ImmutableMatrix:
    """A single-field transverse symbol with constant Hermitian norm.

    The construction uses the fixed frame ``(e1,e2,e3)`` and therefore is not
    an SO(3)-scalar state.  For a unit real direction ``n`` it returns

    ``(P_n e1 + I*n cross e2)/sqrt(1+n3**2)``.
    """

    n = sp.ImmutableMatrix(direction)
    if n.shape != (3, 1):
        raise ValueError("direction must have three components")
    e1 = sp.ImmutableMatrix([1, 0, 0])
    e2 = sp.ImmutableMatrix([0, 1, 0])
    projector = sp.eye(3) - n * n.T
    return sp.ImmutableMatrix(
        (projector * e1 + sp.I * n.cross(e2)) / sp.sqrt(1 + n[2] ** 2)
    )


def gaussian_tail_cross_potential(
    q1: sp.Expr,
    q2: sp.Expr,
    separation: sp.Expr,
    sigma: sp.Expr,
    *,
    density: sp.Expr = sp.Integer(1),
) -> sp.Expr:
    """Exact cross energy for the unit-norm fixed-frame Gaussian symbol."""

    if sp.simplify(separation) == 0:
        return sp.simplify(density * q1 * q2 / (4 * sp.pi ** sp.Rational(3, 2) * sigma))
    return sp.simplify(
        density
        * q1
        * q2
        * sp.erf(separation / (2 * sigma))
        / (4 * sp.pi * separation)
    )


def fixed_frame_angular_stress() -> sp.ImmutableMatrix:
    """Exact real angular stress of :func:`fixed_frame_constant_norm_symbol`."""

    return sp.ImmutableMatrix.diag(
        sp.pi * (8 - sp.pi) / 2,
        sp.pi * (-8 + 3 * sp.pi) / 6,
        4 * sp.pi / 3,
    )


def degree_minus_two_domain_ledger(
    sector: Literal["radial", "toroidal_l1"] = "toroidal_l1",
) -> TailDomainLedger:
    """Return the representation-specific integrability ledger.

    A radial tail has ``x cross u=0`` in the exterior and therefore need not
    inherit the toroidal example's angular-momentum divergence.  Absolute
    vorticity impulse remains unguaranteed in both sectors.
    """

    if sector not in ("radial", "toroidal_l1"):
        raise ValueError("sector must be 'radial' or 'toroidal_l1'")
    return TailDomainLedger(
        sector=sector,
        finite_kinetic_energy=True,
        absolutely_integrable_velocity=False,
        absolutely_integrable_vorticity=False,
        finite_ordinary_angular_momentum=sector == "radial",
        finite_absolute_vorticity_impulse=False,
        finite_helicity_at_infinity=True,
    )


def l1_multiplicity_cross_block(
    radial_a: Sequence[sp.Expr],
    toroidal_a: Sequence[sp.Expr],
    radial_b: Sequence[sp.Expr],
    toroidal_b: Sequence[sp.Expr],
    separation_direction: Sequence[sp.Expr],
    *,
    density: sp.Expr = sp.Integer(1),
    separation: sp.Expr = sp.Integer(1),
) -> sp.ImmutableMatrix:
    """The radial/toroidal ``l=1`` translated cross-energy block.

    Rows refer to the first carrier's radial/toroidal copy and columns to the
    second carrier's copies. Radial coefficients are polar and toroidal
    coefficients are axial under O(3). The off-diagonal entries are
    parity-sensitive chiral kernels when the second carrier is translated as
    ``u2(x-d)``.
    """

    ar = sp.ImmutableMatrix(radial_a)
    at = sp.ImmutableMatrix(toroidal_a)
    br = sp.ImmutableMatrix(radial_b)
    bt = sp.ImmutableMatrix(toroidal_b)
    nv = sp.ImmutableMatrix(separation_direction)
    if any(v.shape != (3, 1) for v in (ar, at, br, bt, nv)):
        raise ValueError("all axes and separation_direction must have three components")
    if sp.simplify(separation) == 0:
        raise ValueError("separation must be nonzero")
    def symmetric(a: sp.ImmutableMatrix, b: sp.ImmutableMatrix) -> sp.Expr:
        return sp.simplify(a.dot(b) + a.dot(nv) * b.dot(nv))

    def chiral(a: sp.ImmutableMatrix, b: sp.ImmutableMatrix) -> sp.Expr:
        return sp.simplify(nv.dot(a.cross(b)))

    prefactor = sp.simplify(density / separation)
    return sp.ImmutableMatrix(
        [
            [
                sp.pi**3 * prefactor * symmetric(ar, br) / 8,
                -2 * sp.pi * prefactor * chiral(ar, bt),
            ],
            [
                -2 * sp.pi * prefactor * chiral(at, br),
                2 * sp.pi * prefactor * symmetric(at, bt),
            ],
        ]
    )


def l1_homogeneous_tail(
    position: Sequence[sp.Expr],
    radial_axis: Sequence[sp.Expr],
    toroidal_axis: Sequence[sp.Expr],
) -> sp.ImmutableMatrix:
    """Return the complete radial/toroidal source-free ``l=1`` tail."""

    x = sp.ImmutableMatrix(position)
    a = sp.ImmutableMatrix(radial_axis)
    b = sp.ImmutableMatrix(toroidal_axis)
    if x.shape != (3, 1) or a.shape != (3, 1) or b.shape != (3, 1):
        raise ValueError("position and axes must have three components")
    r2 = sp.simplify(x.dot(x))
    return sp.ImmutableMatrix(
        (a.dot(x)) * x / r2**2 + b.cross(x) / r2 ** sp.Rational(3, 2)
    )


def steady_curl_residual(
    velocity: Sequence[sp.Expr], coordinates: Sequence[sp.Symbol]
) -> sp.ImmutableMatrix:
    """Return ``curl((u dot grad)u)`` for an exact symbolic velocity."""

    u = sp.ImmutableMatrix(velocity)
    coords = tuple(coordinates)
    if u.shape != (3, 1) or len(coords) != 3:
        raise ValueError("velocity and coordinates must have three components")
    acceleration = sp.ImmutableMatrix(
        [sum(u[j] * sp.diff(u[i], coords[j]) for j in range(3)) for i in range(3)]
    )
    return sp.ImmutableMatrix(
        [
            sp.diff(acceleration[2], coords[1])
            - sp.diff(acceleration[1], coords[2]),
            sp.diff(acceleration[0], coords[2])
            - sp.diff(acceleration[2], coords[0]),
            sp.diff(acceleration[1], coords[0])
            - sp.diff(acceleration[0], coords[1]),
        ]
    )

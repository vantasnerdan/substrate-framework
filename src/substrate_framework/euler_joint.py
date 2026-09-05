"""Unpromoted joint-observation tools for P251 / issue 200.

These APIs derive finite material-momentum jets and solve finite passive
output moments. They do not supply a stationary Euler ensemble, analytic
continuum error bounds, or an inherited action. A continuum use needs the
actual positive material measures and regular streamline-band moments,
including ambient reaction, preparation costs and independent action matching.
"""

from dataclasses import dataclass

import sympy as sp

from substrate_framework.euler_observation import (
    material_tag_fourier_dipole,
    material_tag_moments,
)


@dataclass(frozen=True)
class MaterialMomentumJet:
    """Point-minus-centroid momentum Taylor coefficients for exp(-i epsilon k.x).

    Coefficients include their Taylor factorials. ``local_second`` retains
    absolute parcel velocities: a centroid boost multiplies the finite-size
    mass form factor and cannot be dropped as an intrinsic velocity.
    """

    first: sp.ImmutableMatrix
    local_second: sp.ImmutableMatrix
    second: sp.ImmutableMatrix


def material_momentum_second_jet(masses, positions, velocities, wavevector):
    """Derive both momentum jets from the actual simultaneous material data.

    For a continuous material tag, replace the sums by integrals. The local
    cubic remainder is bounded by |epsilon|^3 sum m |v| |k.r|^3 / 6
    for positive masses and real k,r. Differentiating a material response
    additionally requires the corresponding differentiated moment bound.
    This bound is a Taylor estimate, not a PDE regularity assertion.
    """
    moments = material_tag_moments(masses, positions, velocities)
    k = sp.Matrix(wavevector)
    if k.shape == (1, 3):
        k = k.T
    if k.shape != (3, 1):
        raise ValueError("wavevector must have three components")
    first = material_tag_fourier_dipole(k, moments.spin, moments.shape_rate)
    local_second = sp.zeros(3, 1)
    for mass, position, velocity in zip(masses, positions, velocities):
        r = sp.Matrix(position).reshape(3, 1)-moments.centroid
        v = sp.Matrix(velocity).reshape(3, 1)
        local_second -= sp.sympify(mass)*v*k.dot(r)**2/2
    second = local_second-sp.I*k.dot(moments.centroid)*first
    return MaterialMomentumJet(first, sp.ImmutableMatrix(sp.simplify(local_second)),
                               sp.ImmutableMatrix(sp.simplify(second)))


@dataclass(frozen=True)
class PassiveOutputWeights:
    """Signed coherent source weights; the band densities stay positive.

    No phase or Jacobi energy cancellation is inferred from these weights.
    Those full forms and their normalization are a separate supplier.
    """

    exponents: tuple[int, ...]
    weights: sp.ImmutableMatrix
    residuals: sp.ImmutableMatrix


def passive_output_weights(moment_matrix, coefficients, *, parity):
    """Match an even cosine or odd sine output polynomial on actual bands.

    Row j contains integral eta_l omega**(2*j+parity), where parity=0
    selects g proportional to sin(theta) and parity=1 selects cos(theta).
    The physical observation is <sin(theta) g(theta-omega*t)>.
    The coefficients multiply t**(2*j+parity). Ordered disjoint positive
    frequency bands give an invertible integrated Vandermonde. Arbitrary
    matrices exercise only this finite algebra, not Euler existence.
    """
    if parity not in (0, 1):
        raise ValueError("parity must be 0 (even) or 1 (odd)")
    coefficients = tuple(map(sp.sympify, coefficients))
    if not coefficients:
        raise ValueError("at least one polynomial coefficient is required")
    matrix = sp.Matrix(moment_matrix)
    count = len(coefficients)
    if matrix.shape != (count, count):
        raise ValueError("one band per polynomial coefficient is required")
    if sp.simplify(matrix.det()).is_zero is True:
        raise ValueError("the band moment matrix must be nonsingular")
    exponents = tuple(2*j+parity for j in range(count))
    target = sp.Matrix([(-1)**j*sp.factorial(power)*value
                        for j, (power, value) in enumerate(zip(exponents, coefficients))])
    weights = (matrix.inv()*target).applyfunc(sp.simplify)
    return PassiveOutputWeights(exponents, sp.ImmutableMatrix(weights),
                                sp.ImmutableMatrix(sp.simplify(matrix*weights-target)))

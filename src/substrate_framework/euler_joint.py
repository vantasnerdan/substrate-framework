"""Joint-observation tools for P251 / issue 200.

These APIs derive finite material-momentum jets, finite passive output
moments and the physical branch residual. They do not supply a stationary Euler ensemble, analytic
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
from substrate_framework.micropolar import (
    MicropolarCoefficients,
    micropolar_fourier_stiffness,
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


@dataclass(frozen=True)
class PreparedJointSymbol:
    """Five-position reference-axis symbol for the prepared periodic bridge.

    For physical Y=T z+e, branch residual r=z_tt+D z and constant T,
    M Y_tt+K Y = M T r + defect z + M e_tt + K e. The defect is retained
    exactly, including its cubic terms; this is not an invariant-manifold
    construction. Arrays use (Ux, Uy, Phix, Phiy, Phiz), with Kvector=(0,0,k).
    """

    coefficients: MicropolarCoefficients
    observation: sp.ImmutableMatrix
    mass: sp.ImmutableMatrix
    branch_frequency_squared: sp.ImmutableMatrix
    stiffness: sp.ImmutableMatrix
    defect: sp.ImmutableMatrix


def prepared_joint_symbol(wavenumber, density, spin_density, optical_frequency,
                          acoustic_speed_squared, transverse_dispersion,
                          longitudinal_dispersion):
    """Derive the physical residual from the canonical micropolar operator.

    Real k and positive rho,j,nu,a,cL with cT>j*nu**2/(4*rho) give the
    reviewed positive bulk sector. Signs remain available for counterexample
    probes, as in the canonical coefficient API; rho must be nonzero.
    Parameters are measured/prepared source inputs, not universal material
    constants. Proper rotation transfers this reference-axis representation.
    At k=0 it is its continuous five-position limit with the axis retained.
    """
    k, rho, j, nu, a, ct, cl = map(sp.sympify, (
        wavenumber, density, spin_density, optical_frequency,
        acoustic_speed_squared, transverse_dispersion, longitudinal_dispersion))
    if rho.is_zero is True:
        raise ValueError("density must be nonzero")
    alpha = j*nu**2/4
    coefficients = MicropolarCoefficients(
        0, rho*a, alpha, j*cl/2, 0, j*(ct-alpha/rho))
    curl = sp.I*sp.Matrix([[0, -k, 0], [k, 0, 0], [0, 0, 0]])
    retained = [0, 1, 3, 4, 5]
    observation = sp.BlockMatrix([
        [sp.eye(3), -j*curl/(2*rho)], [curl/2, sp.eye(3)]
    ]).as_explicit().extract(retained, retained)
    mass = sp.diag(rho, rho, j, j, j)
    frequency = sp.diag(a*k**2, a*k**2, nu**2+ct*k**2,
                        nu**2+ct*k**2, nu**2+cl*k**2)
    stiffness = micropolar_fourier_stiffness([0, 0, k], coefficients).extract(
        retained, retained)
    defect = sp.simplify(stiffness*observation-mass*observation*frequency)
    return PreparedJointSymbol(coefficients, *(sp.ImmutableMatrix(value) for value in
        (observation, mass, frequency, stiffness, defect)))


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

"""Conditional, unpromoted joint Euler-orbit reduction (P251/0049).

Inputs are complete Hessian/KKS integrals of an independently constructed
Euler orbit, not freely supplied physical moduli. These algebraic helpers
do not construct that orbit, enforce its Kelvin/boundary conditions, or
prove an invariant PDE truncation. Positivity of symbolic input matrices
is a caller hypothesis when it cannot be decided exactly here.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class EulerRotorReduction:
    """L=dot(B,q).M.dot(B,q)/2 + q*g.dot(B,q) - K*q²/2."""

    kinetic: sp.ImmutableMatrix
    gyro: sp.ImmutableMatrix
    stiffness: sp.Expr
    momentum_hessian: sp.ImmutableMatrix


@dataclass(frozen=True)
class AffineCageRotation:
    """Psi=beta+factor*q, with B=beta+q; all angles are dimensionless."""

    factor: sp.Expr
    spin_inertia: sp.Expr
    cage_inertia: sp.Expr
    stiffness: sp.Expr


def _positive_symmetric(matrix, size):
    value = sp.Matrix(matrix)
    if value.shape != (size, size):
        raise ValueError(f"matrix must be symmetric {size}x{size}")
    if any(entry.is_real is False or entry.is_finite is False
           or entry.has(sp.nan, sp.zoo) for entry in value):
        raise ValueError("matrix must be real and finite")
    if sp.simplify(value-value.T) != sp.zeros(size):
        raise ValueError(f"matrix must be symmetric {size}x{size}")
    # Sylvester's exact criterion; undecidable symbolic signs are hypotheses.
    for count in range(1, size+1):
        minor = sp.simplify(value[:count, :count].det())
        if minor.is_positive is False:
            raise ValueError("matrix must be positive definite")
    return value


def reduce_euler_rotor_block(compact_hessian, body_pairing, internal_pairing):
    """Eliminate both reaction momenta of one complete Euler orbit action.

    H is ordered (r,q,s); Omega has body_pairing*dB^dr and
    internal_pairing*dq^ds. Starting action is
    body_pairing*r*Bdot + internal_pairing*s*qdot - (r,q,s).H.(r,q,s)/2.
    The returned gyro includes its q*qdot total derivative so no mixed term
    disappears silently. Reversing both pairings preserves M,K and reverses
    gyro. An ensemble cancels it only when its reaction momenta are varied
    independently before averaging the reduced actions.
    """
    h = _positive_symmetric(compact_hessian, 3)
    pairings = [sp.sympify(body_pairing), sp.sympify(internal_pairing)]
    if any(value.is_zero is True or value.is_real is False or value.is_finite is False
           or value.has(sp.nan, sp.zoo) for value in pairings):
        raise ValueError("KKS pairings must be real, finite and nonzero")
    momentum = h.extract([0, 2], [0, 2])
    coupling = h.extract([0, 2], [1])
    c = sp.diag(*pairings)
    inverse = momentum.inv()
    kinetic = sp.simplify(c*inverse*c)
    gyro = sp.simplify(-c*inverse*coupling)
    stiffness = sp.simplify(h[1, 1]-(coupling.T*inverse*coupling)[0])
    return EulerRotorReduction(sp.ImmutableMatrix(kinetic), sp.ImmutableMatrix(gyro),
                               stiffness, sp.ImmutableMatrix(momentum))


def affine_cage_rotation_map(kinetic, stiffness):
    """Diagonalize the time-even kinetic term in physical affine-cage fields.

    Requires the odd gyroscopic term to be retained separately or canceled
    by the explicitly paired ensemble. Input M is ordered (B,q). The physical
    cage angle is beta=B-q. If c=M_BB+M_Bq vanishes, this particular absolute
    angle map is singular; the function does not invent a different cage.
    """
    m = _positive_symmetric(kinetic, 2)
    k = sp.sympify(stiffness)
    if k.is_real is False or k.is_positive is False or k.is_finite is False or k.has(sp.nan, sp.zoo):
        raise ValueError("stiffness must be real, finite and positive")
    c = sp.simplify(m[0, 0]+m[0, 1])
    if c.is_zero is True:
        raise ValueError("the declared affine-cage angle map is singular")
    d = sp.simplify(m[0, 0]+2*m[0, 1]+m[1, 1])
    factor = sp.simplify(d/c)
    return AffineCageRotation(factor, sp.simplify(c*c/d), sp.simplify(m.det()/d),
                              sp.simplify(k/factor**2))


@dataclass(frozen=True)
class IsotropicAxisGradient:
    """W=(A*||G||²+B*((tr G)²+tr(G²)))/2, including cell density."""

    norm_coefficient: sp.Expr
    mixed_coefficient: sp.Expr
    trace_modulus: sp.Expr
    symmetric_modulus: sp.Expr
    skew_modulus: sp.Expr


def _positive_scalar(value, name):
    value = sp.sympify(value)
    if (value.is_real is False or value.is_positive is False
            or value.is_finite is False or value.has(sp.nan, sp.zoo)):
        raise ValueError(f"{name} must be real, finite and positive")
    return value


def isotropic_axis_gradient(gradient_matrix, axis, cell_density=1):
    """Haar-average a cell's scalar-axis gradient action over simultaneous rotations.

    The cell energy is (grad(a.Phi)).C.(grad(a.Phi))/2. C must be positive
    symmetric, a real unit vector; cell_density counts cells per action volume.
    This computes an ensemble average, not a homogenization or existence proof.
    Modulus convention: W=c_tr*(tr G)²+c_s*||sym G||²+c_a*||skew G||².
    The trace coefficient may be negative; coercivity instead requires
    c_s>0, c_a>0, and 3*c_tr+c_s>0.
    """
    c = _positive_symmetric(gradient_matrix, 3)
    a = sp.Matrix(axis)
    if a.shape != (3, 1) or any(entry.is_real is False or entry.is_finite is False
                               or entry.has(sp.nan, sp.zoo) for entry in a):
        raise ValueError("axis must be a real finite unit 3-vector")
    if sp.simplify(a.dot(a)-1) != 0:
        raise ValueError("axis must be a real finite unit 3-vector")
    density = _positive_scalar(cell_density, "cell density")
    trace, longitudinal = sp.trace(c), (a.T*c*a)[0]
    norm = sp.simplify(density*(2*trace-longitudinal)/15)
    mixed = sp.simplify(density*(3*longitudinal-trace)/30)
    return IsotropicAxisGradient(norm, mixed, mixed/2,
                                 sp.simplify((norm+mixed)/2),
                                 sp.simplify((norm-mixed)/2))


@dataclass(frozen=True)
class MicropolarKineticNormalForm:
    """Physical fields=field_map*normal fields, modulo spatial order three."""

    field_map: sp.ImmutableMatrix
    transverse_curvature: sp.Expr
    residual_spin_gradient_inertia: sp.Expr


def micropolar_kinetic_normal_form(
        density, spin_inertia, locking, curvature, translation_gradient_inertia,
        spin_gradient_inertia, mixed_inertia, wave_number, helicity=1):
    """Mass-normalize an isotropic transverse action through spatial order two.

    In each curl helicity h=+/-1 the physical kinetic matrix is
    [[rho+m_u*k², b*h*k], [b*h*k, j+m_phi*k²]], and the potential matrix is
    [[A*k², -2*alpha*h*k], [-2*alpha*h*k, 4*alpha+C*k²]].
    Return the SAME field transformation for both forms and its corrected C.
    rho,j,alpha are positive; other coefficients are real finite input data.
    Equality is coefficientwise through k², not an all-k equality or an
    assertion that physical centroid displacement is unchanged. The j=0
    structure-free branch requires its own unreduced action, not this map.
    """
    rho = _positive_scalar(density, "density")
    j = _positive_scalar(spin_inertia, "spin inertia")
    alpha = _positive_scalar(locking, "locking")
    values = list(map(sp.sympify, (curvature, translation_gradient_inertia,
                                 spin_gradient_inertia, mixed_inertia, wave_number)))
    if any(value.is_real is False or value.is_finite is False
           or value.has(sp.nan, sp.zoo) for value in values):
        raise ValueError("gradient coefficients and wave number must be real and finite")
    c, mu, mp, b, k = values
    h = sp.sympify(helicity)
    if h not in (-1, 1):
        raise ValueError("helicity must be +1 or -1")
    residual = sp.simplify(mp-b*b/rho)
    transform = sp.ImmutableMatrix([[1-mu*k*k/(2*rho), -b*h*k/rho],
                                    [0, 1-residual*k*k/(2*j)]])
    corrected = sp.simplify(c+4*alpha*b/rho-4*alpha*residual/j)
    return MicropolarKineticNormalForm(transform, corrected, residual)

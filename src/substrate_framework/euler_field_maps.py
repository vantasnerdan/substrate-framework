"""Conditional, unpromoted fluid field-map audits for P253/0013.

Lamb/Maxwell-shaped identities and scalar Madelung energy terms are computed
from supplied fields. Defined currents do not supply physical electric charge,
weak charge, source-free propagation or a quantum measurement law. The actual
Euler and divergence residuals are returned separately to expose that gap.
"""

from dataclasses import dataclass

import sympy as sp


def _coordinates(coordinates):
    q = tuple(coordinates)
    if len(q) != 3 or len(set(q)) != 3 or not all(isinstance(x, sp.Symbol) for x in q):
        raise ValueError("coordinates must be three distinct symbols")
    return q


def _positive(value, name):
    x = sp.sympify(value)
    if (x.is_positive is False or x.is_real is False or x.is_finite is False
            or (x.is_number and x.is_positive is not True)):
        raise ValueError(f"{name} must be positive finite real")
    return x


def _gradient(scalar, coordinates):
    return sp.Matrix([sp.diff(scalar, q) for q in coordinates])


def _curl(vector, coordinates):
    x, y, z = coordinates
    return sp.Matrix([vector[2].diff(y)-vector[1].diff(z),
                      vector[0].diff(z)-vector[2].diff(x),
                      vector[1].diff(x)-vector[0].diff(y)])


def _immutable(matrix):
    return sp.ImmutableMatrix(matrix.applyfunc(sp.simplify))


@dataclass(frozen=True)
class EulerLambMap:
    """Physical Euler residuals and the separately defined field analogy."""

    vorticity: sp.ImmutableMatrix
    lamb_vector: sp.ImmutableMatrix
    bernoulli: sp.Expr
    hydrodynamic_charge: sp.Expr
    defined_current: sp.ImmutableMatrix
    euler_residual: sp.ImmutableMatrix
    divergence_residual: sp.Expr
    faraday_residual: sp.ImmutableMatrix
    charge_continuity_residual: sp.Expr


def euler_lamb_map(velocity, pressure, coordinates, time, *, density, speed):
    """Compute B=curl(u), E=B cross u, q=div(E), J=c^2 curl(B)-E_t.

    The supplied positive c is a constant comparison scale; changing it changes
    the defined current. It is not derived from Euler. Density is physical,
    positive and constant. Reality and smoothness of undecidable symbolic
    fields remain caller hypotheses. A zero charge-continuity residual is an
    identity and does not certify the separate Euler/divergence residuals.
    """
    q = _coordinates(coordinates)
    if not isinstance(time, sp.Symbol) or time in q:
        raise ValueError("time must be a distinct symbol")
    rho, c = _positive(density, "density"), _positive(speed, "speed")
    if set((*q, time)) & (rho.free_symbols | c.free_symbols):
        raise ValueError("density and comparison speed must be spacetime constants")
    u = sp.Matrix(velocity)
    if u.shape == (1, 3):
        u = u.T
    if u.shape != (3, 1):
        raise ValueError("velocity must have three components")
    p = sp.sympify(pressure)
    B = _curl(u, q)
    E = B.cross(u)
    charge = sp.simplify(sum(E[i].diff(q[i]) for i in range(3)))
    current = c*c*_curl(B, q)-E.diff(time)
    residual = u.diff(time)+u.jacobian(q)*u+_gradient(p, q)/rho
    continuity = sp.diff(charge, time)+sum(current[i].diff(q[i]) for i in range(3))
    return EulerLambMap(_immutable(B), _immutable(E), sp.simplify(p/rho+u.dot(u)/2),
                        charge, _immutable(current), _immutable(residual),
                        sp.simplify(sum(u[i].diff(q[i]) for i in range(3))),
                        _immutable(B.diff(time)+_curl(E, q)), sp.simplify(continuity))


@dataclass(frozen=True)
class MadelungKineticTerms:
    """Scalar-phase velocity, classical kinetic density and extra gradient energy."""

    phase_velocity: sp.ImmutableMatrix
    classical_kinetic_density: sp.Expr
    density_gradient_energy: sp.Expr
    quantum_potential: sp.Expr


def madelung_kinetic_terms(number_density, phase, coordinates, *, mass, action):
    """Split hbar^2*|grad(sqrt(n)*exp(i*S/hbar))|^2/(2m).

    ``action`` is the supplied positive hbar, not a selected action quantum;
    n is a positive number/probability density, not silently material rho.
    Smoothness and positivity that symbolic analysis cannot decide remain
    hypotheses. A single regular scalar phase has zero local curl.
    """
    q = _coordinates(coordinates)
    n = _positive(number_density, "number density")
    m, hbar = _positive(mass, "mass"), _positive(action, "action")
    if set(q) & (m.free_symbols | hbar.free_symbols):
        raise ValueError("mass and action scale must be spatial constants")
    S = sp.sympify(phase)
    dn, dS = _gradient(n, q), _gradient(S, q)
    amplitude = sp.sqrt(n)
    laplacian = sum(amplitude.diff(x, 2) for x in q)
    return MadelungKineticTerms(_immutable(dS/m), sp.simplify(n*dS.dot(dS)/(2*m)),
                                sp.simplify(hbar*hbar*dn.dot(dn)/(8*m*n)),
                                sp.simplify(-hbar*hbar*laplacian/(2*m*amplitude)))


@dataclass(frozen=True)
class SpinorEulerTerms:
    """Local unit-spinor coordinates and their physical kinetic-energy split."""

    velocity: sp.ImmutableMatrix
    texture: sp.ImmutableMatrix
    normalization_residual: sp.Expr
    euler_kinetic_density: sp.Expr
    dirichlet_energy_density: sp.Expr
    texture_energy_density: sp.Expr
    energy_identity_residual: sp.Expr


def spinor_euler_terms(spinor, coordinates, *, density, circulation_scale):
    """Compute u=-i*kappa*z^dagger grad(z) and the exact Euler energy split.

    Unit norm, regular chart coverage and real velocity are caller hypotheses
    when symbolic simplification cannot decide them. Known nonunit input is
    rejected and any undecided normalization residual is exposed. Kappa is a
    positive circulation coordinate scale, not a derived action quantum.
    The texture Dirichlet energy is SUBTRACTED from the spinor Dirichlet
    energy to recover Euler; adding a texture restoring law changes dynamics.
    """
    q = _coordinates(coordinates)
    rho = _positive(density, "density")
    kappa = _positive(circulation_scale, "circulation scale")
    if set(q) & (rho.free_symbols | kappa.free_symbols):
        raise ValueError("density and circulation scale must be spatial constants")
    z = sp.Matrix(spinor)
    if z.shape == (1, 2):
        z = z.T
    if z.shape != (2, 1):
        raise ValueError("spinor must have two components")
    norm_residual = sp.simplify((z.conjugate().T*z)[0]-1)
    if norm_residual.is_zero is False:
        raise ValueError("spinor must have unit pointwise norm")
    pauli = (sp.Matrix([[0, 1], [1, 0]]),
             sp.Matrix([[0, -sp.I], [sp.I, 0]]),
             sp.diag(1, -1))
    texture = _immutable(sp.Matrix([(z.conjugate().T*s*z)[0] for s in pauli]))
    velocity = _immutable(sp.Matrix([-sp.I*kappa*(z.conjugate().T*z.diff(x))[0] for x in q]))
    grad_norm = sum((z.diff(x).conjugate().T*z.diff(x))[0] for x in q)
    texture_norm = sum(texture.diff(x).dot(texture.diff(x)) for x in q)
    kinetic = sp.simplify(rho*velocity.dot(velocity)/2)
    dirichlet = sp.simplify(rho*kappa*kappa*grad_norm/2)
    texture_energy = sp.simplify(rho*kappa*kappa*texture_norm/8)
    return SpinorEulerTerms(velocity, texture, norm_residual, kinetic,
                           dirichlet, texture_energy,
                           sp.trigsimp(sp.simplify(dirichlet-texture_energy-kinetic), method="fu"))

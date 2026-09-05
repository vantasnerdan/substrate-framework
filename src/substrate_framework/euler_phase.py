"""Conditional time-dependent symplectic pullbacks, P251 stronger frontier.

These exact algebraic APIs are unpromoted infrastructure. They keep moving
frame terms and physical observations explicit; they do not construct an
Euler invariant phase family or establish a Floquet energy's physical sign.
The ambient symplectic form is constant. In the convention Omega(q,s)=B,
the action one-form is -x.T*Omega*xdot/2, equivalent to B*s*qdot.
"""

from dataclasses import dataclass

import sympy as sp


def _matrix(value, name):
    result = sp.Matrix(value)
    if any(x.is_finite is False or x.has(sp.nan, sp.zoo) for x in result):
        raise ValueError(f"{name} must have finite entries")
    return result


def _zero(matrix):
    return all(sp.simplify(x) == 0 for x in matrix)


def _immutable(matrix):
    return sp.ImmutableMatrix(matrix.applyfunc(sp.simplify))


@dataclass(frozen=True)
class MovingPhasePullback:
    """Same-action forms, generator, projection and ambient residual.

    The coordinate equation is Omega_E*zdot +
    (H_eff + dotOmega_E/2)*z = 0. ``residual`` is
    Edot + E*generator - ambient_generator*E, in ambient coordinates.
    It is symplectically orthogonal to E, not necessarily zero or gauge.
    An actual observable is O*E*z plus its complementary observation.
    """

    symplectic: sp.ImmutableMatrix
    symplectic_rate: sp.ImmutableMatrix
    hamiltonian: sp.ImmutableMatrix
    generator: sp.ImmutableMatrix
    coordinates: sp.ImmutableMatrix
    projection: sp.ImmutableMatrix
    residual: sp.ImmutableMatrix


def moving_phase_pullback(symplectic, hamiltonian, embedding, embedding_rate):
    """Derive a differentiable phase restriction from its complete action.

    Inputs Omega,H,E,Edot are actual simultaneous matrices, with Omega
    invertible and skew, H symmetric, and E.T*Omega*E invertible. The
    caller supplies Edot as the actual time derivative of E and retains
    any undecidable symbolic nondegeneracy hypotheses. No positivity or
    microscopic realization is inferred. Complex bilinear representations
    may be used; a real physical encoding remains a caller obligation.
    """
    omega = _matrix(symplectic, "symplectic form")
    h = _matrix(hamiltonian, "Hamiltonian")
    e = _matrix(embedding, "embedding")
    ed = _matrix(embedding_rate, "embedding rate")
    if omega.rows != omega.cols or not _zero(omega + omega.T):
        raise ValueError("ambient symplectic form must be square and skew")
    if h.shape != omega.shape or not _zero(h - h.T):
        raise ValueError("Hamiltonian must be symmetric with the ambient shape")
    if e.rows != omega.rows or ed.shape != e.shape or not e.cols:
        raise ValueError("embedding and its rate must have compatible nonempty shapes")
    reduced_omega = sp.simplify(e.T * omega * e)
    if sp.simplify(omega.det()).is_zero is True:
        raise ValueError("ambient symplectic form must be nondegenerate")
    if sp.simplify(reduced_omega.det()).is_zero is True:
        raise ValueError("restricted symplectic form must be nondegenerate")
    connection = e.T * omega * ed
    omega_rate = connection - connection.T
    effective_h = e.T * h * e + (connection + connection.T) / 2
    coordinates = reduced_omega.inv() * e.T * omega
    ambient_generator = -omega.inv() * h
    generator = -reduced_omega.inv() * (effective_h + omega_rate / 2)
    projection = e * coordinates
    residual = ed + e * generator - ambient_generator * e
    return MovingPhasePullback(*map(_immutable, (
        reduced_omega, omega_rate, effective_h, generator,
        coordinates, projection, residual)))


@dataclass(frozen=True)
class PhysicalScalarChart:
    """Actual scalar angle/rate action and separately measured spin.

    ``coordinates`` maps the original phase to (theta, theta_dot).
    The scalar action is mass*theta_dot**2/2-stiffness*theta**2/2.
    Its equation includes mass_rate*theta_dot. The measured spin is
    spin_inertia*theta_dot+spin_connection*theta; it equals canonical
    momentum only with the independently checked normalization/current.
    """

    coordinates: sp.ImmutableMatrix
    generator: sp.ImmutableMatrix
    symplectic: sp.ImmutableMatrix
    hamiltonian: sp.ImmutableMatrix
    wronskian: sp.Expr
    mass: sp.Expr
    mass_rate: sp.Expr
    stiffness: sp.Expr
    spin_inertia: sp.Expr
    spin_connection: sp.Expr
    angle_spin_bracket: sp.Expr


def physical_scalar_chart(symplectic, generator, angle, *, angle_rate,
                          angle_acceleration, generator_rate, spin):
    """Derive a two-dimensional physical observation chart, P251 frontier.

    Omega is constant, invertible and skew; B and Bdot are Hamiltonian
    for Omega. c, cdot, cddot and s are actual 1x2 angle/derivative/spin
    rows on zdot=Bz. The caller supplies their true derivatives and
    microscopic meaning; this algebra cannot license those observations.
    Undecidable symbolic nondegeneracy remains a caller hypothesis.

    A nonzero angle row need not give a chart: det[c;cdot+cB] must be
    nonzero. Positivity of mass or spin_inertia is not inferred, and their
    equality is not imposed. All moving-action and measured-current terms
    are returned, including when a winding changes the physical clock.
    """
    omega = _matrix(symplectic, "symplectic form")
    b = _matrix(generator, "generator")
    bd = _matrix(generator_rate, "generator rate")
    if omega.shape != (2, 2) or not _zero(omega+omega.T):
        raise ValueError("symplectic form must be two by two and skew")
    if sp.simplify(omega.det()).is_zero is True:
        raise ValueError("symplectic form must be nondegenerate")
    for matrix in (b, bd):
        if matrix.shape != (2, 2) or not _zero(matrix.T*omega+omega*matrix):
            raise ValueError("generator and its rate must be Hamiltonian for Omega")
    c, cd, cdd, measured = [_matrix(value, name) for value, name in (
        (angle, "angle"), (angle_rate, "angle rate"),
        (angle_acceleration, "angle acceleration"), (spin, "spin"))]
    if any(row.shape != (1, 2) for row in (c, cd, cdd, measured)):
        raise ValueError("angle, its derivatives and spin must be one by two rows")
    d = cd+c*b
    dd = cdd+cd*b+c*bd
    coordinates = c.col_join(d)
    coordinates_rate = cd.col_join(dd)
    wronskian = sp.simplify(coordinates.det())
    if wronskian.is_zero is True:
        raise ValueError("physical angle/rate Wronskian must be nonzero")
    inverse = coordinates.inv()
    phase = inverse.T*omega*inverse
    chart_generator = (coordinates_rate+coordinates*b)*inverse
    mass = sp.simplify(phase[0, 1])
    # Differentiating T^-T*Omega*T^-1 at constant Omega gives the complete
    # physical mass rate. The scalar damping is -mass_rate/mass.
    inverse_rate = -inverse*coordinates_rate*inverse
    phase_rate = inverse_rate.T*omega*inverse+inverse.T*omega*inverse_rate
    mass_rate = sp.simplify(phase_rate[0, 1])
    stiffness = sp.simplify(-mass*chart_generator[1, 0])
    hamiltonian = -phase*chart_generator-phase_rate/2
    spin_row = measured*inverse
    bracket = -(c*omega.inv()*measured.T)[0]
    return PhysicalScalarChart(
        *map(_immutable, (coordinates, chart_generator, phase, hamiltonian)),
        *map(sp.simplify, (wronskian, mass, mass_rate, stiffness,
                          spin_row[0, 1], spin_row[0, 0], bracket)))

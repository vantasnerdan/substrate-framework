"""Actual passive Euler controls and prepared acoustic action, P251/0196.

Moment inputs must be actual smooth regular streamline-band moments to make
an Euler existence claim. No numeric moment matrix or finite scale alone
establishes the fixed-cell continuum limit proved in the source construction.
"""

from dataclasses import dataclass

import sympy as sp

from substrate_framework.euler_phase import physical_scalar_chart


def _scalar(value, name, *, positive=False):
    value = sp.sympify(value)
    if (not isinstance(value, sp.Expr) or value.is_real is False
            or value.is_finite is False or value.has(sp.nan, sp.zoo, sp.oo, -sp.oo)
            or (positive and value.is_positive is False)):
        raise ValueError(f"{name} must be finite real" + (" and positive" if positive else ""))
    return value


@dataclass(frozen=True)
class WrappedStreamline:
    """Actual period/coarea on Y in (0,pi), psi=B cosY+A cosZ=c."""

    coordinate: sp.Symbol
    speed_squared: sp.Expr
    period: sp.Expr
    frequency: sp.Expr
    coarea_density: sp.Expr
    domain_condition: sp.Expr


def wrapped_streamline(amplitude, dominant, level):
    """Exact integral definitions; no near-saddle numerical quadrature is used."""
    aa = _scalar(amplitude, "A", positive=True)
    bb = _scalar(dominant, "B", positive=True)
    cc = _scalar(level, "streamline level")
    margins = [bb-aa, bb-aa-cc, bb-aa+cc]
    for margin in margins:
        _scalar(margin, "regular wrapped-band margin", positive=True)
    z = sp.Dummy("Z", real=True)
    speed = bb**2-(cc-aa*sp.cos(z))**2
    period = sp.Integral(1/sp.sqrt(speed), (z, 0, 2*sp.pi))
    condition = sp.And(*(sp.Gt(value, 0) for value in [aa, bb, *margins]))
    return WrappedStreamline(z, speed, period, 2*sp.pi/period,
                             period/(4*sp.pi**2), condition)


@dataclass(frozen=True)
class PassivePacket:
    """Actual axial Euler velocity and Eulerian-coordinate Lin displacement.

    The scalar fields multiply the fixed axial unit vector and the declared
    whole-field correlation. G and omega may depend on the first integral c,
    but are constant in angle/time. ``coordinate_rate`` is partial_t xi,
    not the distinct physical Euler velocity partial_t xi+T xi.
    """

    initial_velocity: sp.Expr
    initial_configuration: sp.Expr
    euler_velocity: sp.Expr
    lin_displacement: sp.Expr
    coordinate_rate: sp.Expr
    initial_energy_density: sp.Expr
    initial_energy_angle_average: sp.Expr


def passive_packet(amplitude, frequency, angle, time, *, configuration_fraction=sp.Rational(1, 2)):
    """Construct g=G cos(theta), h=alpha G sin(theta)/omega and their exact flow."""
    if (not isinstance(angle, sp.Symbol) or not isinstance(time, sp.Symbol)
            or angle == time):
        raise ValueError("angle and time must be distinct symbols")
    amplitude = _scalar(amplitude, "packet amplitude")
    frequency = _scalar(frequency, "frequency", positive=True)
    fraction = _scalar(configuration_fraction, "configuration fraction")
    if any(value.has(angle, time) for value in (amplitude, frequency, fraction)):
        raise ValueError("packet parameters must be constant along angle and time")
    g = amplitude*sp.cos(angle)
    h = fraction*amplitude*sp.sin(angle)/frequency
    advected = angle-frequency*time
    velocity = g.subs(angle, advected)
    displacement = (h+time*g).subs(angle, advected)
    rate = sp.diff(displacement, time)
    transport_h = frequency*sp.diff(h, angle)
    energy = sp.expand(((g-transport_h)**2-transport_h**2)/2)
    average = sp.simplify(sp.integrate(energy, (angle, 0, 2*sp.pi))/(2*sp.pi))
    return PassivePacket(g, h, velocity, displacement, rate, energy, average)


@dataclass(frozen=True)
class PhaseMatchedControls:
    """Signed preparation weights on positive actual bands, never signed probabilities."""

    exponents: tuple[int, ...]
    moment_matrix: sp.ImmutableMatrix
    weights: sp.ImmutableMatrix
    phase_residual: sp.Expr
    output_residuals: tuple[sp.Expr, ...]
    domain_condition: sp.Expr


def phase_matched_controls(moment_matrix, odd_coefficients):
    """Solve exact negative-frequency phase and finite odd output moments.

    For polynomial sum_j p_j t^(2j+1), rows of the actual moment matrix are
    integral eta_l omega^r for r=-1,1,3,... . There are len(p)+1 bands.
    Ordered disjoint positive-frequency bands guarantee nonsingularity in
    the analytic construction; unknown symbolic input domains stay explicit.
    """
    coefficients = tuple(_scalar(value, "odd polynomial coefficient")
                         for value in odd_coefficients)
    if not coefficients:
        raise ValueError("at least one odd polynomial coefficient is required")
    matrix = sp.Matrix(moment_matrix)
    count = len(coefficients)+1
    if matrix.shape != (count, count):
        raise ValueError("one extra band and phase row are required")
    for value in matrix:
        _scalar(value, "positive band moment", positive=True)
    determinant = sp.simplify(matrix.det())
    if determinant.is_zero is True:
        raise ValueError("the actual band moment matrix must be nonsingular")
    rhs = sp.Matrix([0]+[(-1)**j*sp.factorial(2*j+1)*value
                        for j, value in enumerate(coefficients)])
    weights = (matrix.inv()*rhs).applyfunc(sp.simplify)
    residual = (matrix*weights-rhs).applyfunc(sp.simplify)
    condition = sp.And(sp.Ne(determinant, 0), *(sp.Gt(value, 0) for value in matrix))
    return PhaseMatchedControls(tuple([-1]+[2*j+1 for j in range(len(coefficients))]),
                                sp.ImmutableMatrix(matrix), sp.ImmutableMatrix(weights),
                                residual[0], tuple(residual[1:]), condition)


@dataclass(frozen=True)
class AcousticActionJet:
    """Physical fixed-time second spatial asymptotic; q denotes |K| squared."""

    spatial_square: sp.Symbol
    displacement_row: sp.ImmutableMatrix
    wronskian: sp.Expr
    mass: sp.Expr
    mass_rate: sp.Expr
    stiffness: sp.Expr
    conserved_energy: sp.ImmutableMatrix


def acoustic_action_jet(restoring, density, time, error_integral):
    """Restrict the complete phase via the actual mean rows, including connection.

    H=error_integral has H(0)=H'(0)=0 and H'' equal to the actual controlled
    acceleration error. The conserved energy matrix is not silently replaced
    by the moving chart's mechanical energy. Unknown initial identities are
    mathematical input hypotheses, not numerically tested assumptions.
    """
    if not isinstance(time, sp.Symbol):
        raise ValueError("time must be a symbol")
    aa = _scalar(restoring, "restoring coefficient", positive=True)
    rho = _scalar(density, "density", positive=True)
    hh = _scalar(error_integral, "integrated actual current error")
    if aa.has(time) or rho.has(time):
        raise ValueError("density and restoring coefficient are constant")
    for value in (hh.subs(time, 0), sp.diff(hh, time).subs(time, 0)):
        if sp.simplify(value).is_zero is False:
            raise ValueError("the integrated error and its first derivative start at zero")
    q = sp.Dummy("spatial_square", real=True)
    row = sp.Matrix([[1-aa*q*time**2/2, time-aa*q*time**3/6+q*hh]])
    chart = physical_scalar_chart(
        sp.Matrix([[0, rho], [-rho, 0]]), sp.zeros(2), row,
        angle_rate=row.diff(time), angle_acceleration=row.diff(time, 2),
        generator_rate=sp.zeros(2), spin=rho*row.diff(time))

    def truncate(value):
        return sp.simplify(sp.series(value, q, 0, 2).removeO())

    inverse = chart.coordinates.inv()
    energy = inverse.T*sp.diag(rho*aa*q, rho)*inverse
    return AcousticActionJet(q, sp.ImmutableMatrix(row), truncate(chart.wronskian),
                             truncate(chart.mass), truncate(chart.mass_rate),
                             truncate(chart.stiffness),
                             sp.ImmutableMatrix(energy.applyfunc(truncate)))


def diagonal_wave_number(accuracy, control_norm, remainder_constant, previous=1):
    """Choose an actual positive scale AFTER all finite control/remainder constants.

    For accuracy tending to zero this enforces k*N ->0 and k*C ->0.
    The source proof supplies the actual C for full Euler/Lin derivatives;
    passing an arbitrary constant here does not prove that analytic bound.
    """
    epsilon, norm, constant, prior = (
        _scalar(value, name, positive=True) for value, name in (
            (accuracy, "accuracy"), (control_norm, "control norm"),
            (remainder_constant, "remainder constant"), (previous, "previous scale")))
    return sp.Min(sp.Rational(1, 4), prior/2, epsilon/(1+norm+constant))/2

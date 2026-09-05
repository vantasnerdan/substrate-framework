"""Actual finite Euler displacement preparation (P251/0188, reviewed0197).

The physical observation is X=D+integral(mean Euler velocity), not the mean
material displacement. A positive D column and initial phase do not supply
the common-V history or an autonomous acoustic/EPS continuum.
"""

from dataclasses import dataclass

import sympy as sp

from substrate_framework import euler_fourier as ef


def _clean(field):
    return ef.add({wave: sp.cancel(value) for wave, value in field.items()})


def _add(*vectors):
    return tuple(_clean(ef.add(*(vector[i] for vector in vectors))) for i in range(3))


def _scale(vector, value):
    return tuple(ef.scale(component, value) for component in vector)


def _laplacian(field):
    return _clean({wave: sum(q*q for q in wave)*value for wave, value in field.items()})


def _pair(left, right):
    return sp.cancel(ef.mul(left, right).get(ef.ZERO, 0))


def transverse_pair_average(polynomial, direction, polarization):
    """Exact Haar average for a unit orthogonal pair, quadratic in polarization.

    Inputs are two disjoint triples of independent symbols. Direction degree
    may be arbitrary finite; no sampling or fitted isotropy tensor is used.
    """
    direction, polarization = tuple(direction), tuple(polarization)
    variables = direction+polarization
    if (len(direction) != 3 or len(polarization) != 3
            or any(not isinstance(v, sp.Symbol) for v in variables)
            or len(set(variables)) != 6):
        raise ValueError("direction and polarization must be disjoint triples of symbols")
    conditional = 0
    for powers, coefficient in sp.Poly(sp.expand(polynomial), *polarization).terms():
        if coefficient == 0:
            continue
        indices = [i for i in range(3) for _ in range(powers[i])]
        if len(indices) != 2:
            raise ValueError("the observable must be homogeneous quadratic in polarization")
        i, j = indices
        conditional += coefficient*(sp.KroneckerDelta(i, j)-direction[i]*direction[j])/2
    result = 0
    for powers, coefficient in sp.Poly(sp.expand(conditional), *direction).terms():
        if any(power % 2 for power in powers):
            continue
        result += (coefficient*sp.prod(sp.factorial2(power-1) for power in powers)
                   /sp.factorial2(sum(powers)+1))
    return sp.simplify(result)


@dataclass(frozen=True)
class DisplacementCell:
    """Full finite Fourier supplier in density-one, unit-curl-length units."""

    background: tuple
    pressure: dict
    first_integral: dict
    scalar_corrector: dict
    velocity_corrector: tuple
    forcing: tuple
    range_energy: sp.Expr
    background_energy: sp.Expr
    matching_amplitude: sp.Expr
    restoring_coefficient: sp.Expr


def finite_displacement_cell():
    """Construct the actual A/B=1/100 degree-seven cell and exact matching root.

    Every Fourier product is retained. The range energy is computed from the
    constructed velocity, rather than supplied as a rational fixture. The
    full first-shell energy/current coefficients used in the matching quadratic
    are independently derived in0197 and exposed by direct-field tests.
    """
    aa = sp.Rational(1, 100)
    psi = ef.add(ef.trig(1), ef.scale(ef.trig(2), aa))
    alpha = ef.add(ef.trig(1), ef.scale(ef.trig(2), -aa))
    background = (psi, ef.scale(ef.trig(2, kind="sin"), aa),
                  ef.scale(ef.trig(1, kind="sin"), -1))
    argument = sp.Dummy("argument")

    def evaluate(expression):
        polynomial = sp.Poly(expression, argument)
        result, power = {}, {ef.ZERO: sp.S.One}
        for order in range(polynomial.degree()+1):
            result = ef.add(result, ef.scale(power, polynomial.nth(order)))
            power = ef.mul(power, psi)
        return _clean(result)

    harmonics = (3, 5, 7)
    weights = {n: sp.Rational(n*n*(n*n+1), 2*(n*n-1)**2) for n in harmonics}
    denominator = sum(sp.Rational(n*n)/weights[n] for n in harmonics)
    coefficients = {n: 2*n/(weights[n]*denominator) for n in harmonics}
    cubic = evaluate(sp.chebyshevt(3, argument))
    fixed = evaluate(sum(coefficients[n]*sp.chebyshevt(n, argument) for n in (5, 7)))
    matrix = sp.Matrix([[_pair(ef.trig(axis), basis) for basis in (psi, cubic)]
                        for axis in (1, 2)])
    rhs = sp.Matrix([-_pair(ef.trig(axis), ef.add(alpha, fixed)) for axis in (1, 2)])
    solved = matrix.inv()*rhs
    first_integral = _clean(ef.add(fixed, ef.scale(psi, solved[0]), ef.scale(cubic, solved[1])))
    forcing_scalar = _clean(ef.add(alpha, first_integral))
    if any(sum(q*q for q in wave) == 1 for wave in forcing_scalar):
        raise ArithmeticError("constructed forcing has a nonzero kernel component")
    phi = _clean({wave: value/(sum(q*q for q in wave)-1)
                  for wave, value in forcing_scalar.items()})
    corrector = (ef.add(ef.scale(alpha, -sp.Rational(1, 2)), _laplacian(phi)),
                 ef.scale(ef.derivative(phi, 2), -1), ef.derivative(phi, 1))
    wave_y = (ef.trig(1), {}, ef.scale(ef.trig(1, kind="sin"), -1))
    wave_z = (ef.scale(ef.trig(2), aa), ef.scale(ef.trig(2, kind="sin"), aa), {})
    force = ef.leray(ef.cross(wave_z, wave_y))
    remainder = tuple({wave: value for wave, value in component.items()
                       if sum(q*q for q in wave) > 1} for component in corrector)
    cost = sp.cancel(ef.inner(remainder, remainder)/5)
    energy = ef.inner(background, background)
    radicand = 100-960*cost/energy
    if not (0 < cost < 13*energy/1280):
        raise ArithmeticError("the actual finite cell lacks its positive matching window")
    amplitude = -(4+sp.sqrt(radicand))/8
    restoring = energy*(sp.sqrt(radicand)-9)/120
    pressure = ef.scale(ef.add(*(ef.mul(c, c) for c in background)), -sp.Rational(1, 2))
    return DisplacementCell(background, pressure, first_integral, phi, corrector,
                            force, cost, energy, amplitude, restoring)


def negative_helicity_shell(vector):
    """Orthogonal curl=-1 projection of the complete unit Fourier shell."""
    if len(vector) != 3:
        raise ValueError("a vector has three Fourier components")
    shell = tuple(ef.add({wave: value for wave, value in component.items()
                          if sum(q*q for q in wave) == 1}) for component in vector)
    return _scale(_add(ef.leray(shell), _scale(ef.curl(shell), -1)), sp.Rational(1, 2))


@dataclass(frozen=True)
class PreparedDisplacement:
    """Polynomial extension of actual transverse initial-field rows.

    Physical input has unit direction perpendicular to displacement. The
    polynomial extension is retained to perform exact whole-frame averaging.
    Density multiplies energy and phase, not the derived acceleration.
    """

    translation: tuple
    lift: tuple
    returned: tuple
    first_velocity: tuple
    material_rate: tuple
    current_test: tuple
    energy_coefficient: sp.Expr
    acceleration_contraction: sp.Expr


def prepared_displacement(cell, direction, displacement, *, amplitude=None):
    """Construct the actual correlated first-K data; no V trajectory is inferred."""
    if not isinstance(cell, DisplacementCell):
        raise ValueError("an actual DisplacementCell is required")
    kap, disp = sp.Matrix(direction), sp.Matrix(displacement)
    if kap.shape != (3, 1) or disp.shape != (3, 1):
        raise ValueError("direction and displacement have three components")
    amplitude = cell.matching_amplitude if amplitude is None else sp.sympify(amplitude)
    if (not isinstance(amplitude, sp.Expr) or amplitude.is_real is False
            or amplitude.is_finite is False or amplitude.has(sp.nan, sp.zoo, sp.oo, -sp.oo)):
        raise ValueError("amplitude must be finite and real")
    u = cell.background
    translation = tuple(ef.add(*(ef.scale(ef.derivative(u[i], j), -disp[j])
                                for j in range(3))) for i in range(3))
    lift = _scale(ef.cross(tuple({ef.ZERO: kap[i]} for i in range(3)), translation), -1)
    au = ef.add(*(ef.scale(u[i], kap[i]) for i in range(3)))
    ud = ef.add(*(ef.scale(u[i], disp[i]) for i in range(3)))
    d = kap[1]*disp[1]-kap[2]*disp[2]
    rate = _add(lift, _scale(cell.velocity_corrector, d),
                tuple(ef.scale(au, -disp[i]) for i in range(3)))
    current = tuple(ef.add(ef.scale(au, disp[i]), ef.scale(ud, kap[i])) for i in range(3))
    returned = _add(_scale(negative_helicity_shell(rate), -1),
                    _scale(negative_helicity_shell(current), amplitude))
    final = _add(rate, returned)
    first = _add(lift, _scale(cell.velocity_corrector, d), returned)
    return PreparedDisplacement(translation, lift, returned, first, final, current,
        sp.expand(ef.inner(final, final)-_pair(au, au)*disp.dot(disp)),
        sp.expand(_pair(au, au)*disp.dot(disp)+ef.inner(current, final)))

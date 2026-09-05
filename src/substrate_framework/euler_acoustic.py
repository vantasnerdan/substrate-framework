"""Actual planar Euler array and full axial Kelvin preparation (P251/0167).

The triangular array's acoustic-window theorem is analytic evidence in0161.
These definitions derive its field and initial phase, not an Euler trajectory.
In particular no stability theorem for its constant-curl lift is implied.
"""

from dataclasses import dataclass

import sympy as sp


def _positive(value, name):
    value = sp.sympify(value)
    if (value.is_positive is False or value.is_real is False
            or value.is_finite is False
            or value.has(sp.nan, sp.zoo, sp.oo, -sp.oo)):
        raise ValueError(f"{name} must be finite and positive")
    return value


@dataclass(frozen=True)
class TriangularEulerArray:
    """Same-cell field, physical pressure and exact averaged moments.

    Columns of ``wavevectors`` and ``sine_velocities`` specify
    v(x)=sum_j sine_velocities[:,j]*sin(wavevectors[:,j].x).
    The separatrix dual uses centered bounded coordinates on the actual
    invariant polygons, not globally single-valued Cartesian coordinates.
    """

    wavevectors: sp.ImmutableMatrix
    sine_velocities: sp.ImmutableMatrix
    streamfunction: sp.Expr
    velocity: sp.ImmutableMatrix
    vorticity: sp.Expr
    pressure: sp.Expr
    covariance: sp.ImmutableMatrix
    separatrix_level: sp.Expr
    translation_dual: sp.ImmutableMatrix
    domain_condition: sp.Expr


def triangular_euler_array(amplitude, wavenumber, density, coordinates):
    """Construct the three-mode stationary array in Cartesian coordinates.

    Amplitude is streamfunction amplitude (length²/time), wavenumber is
    positive inverse length and density is positive mass/volume. Their
    unresolved symbolic positivity remains in ``domain_condition``.
    The physical periodic cell is dual to the first two wavevectors.
    All averages use its complete area, including separatrix cells.
    """
    coords = tuple(coordinates)
    if (len(coords) != 2 or len(set(coords)) != 2
            or not all(isinstance(c, sp.Symbol) for c in coords)):
        raise ValueError("two distinct Cartesian coordinate symbols are required")
    amp, lam, rho = (_positive(value, name) for value, name in (
        (amplitude, "amplitude"), (wavenumber, "wavenumber"), (density, "density")))
    if any(value.has(*coords) for value in (amp, lam, rho)):
        raise ValueError("field parameters must be spatial constants")
    wavevectors = sp.ImmutableMatrix([
        [lam, -lam/2, -lam/2],
        [0, sp.sqrt(3)*lam/2, -sp.sqrt(3)*lam/2],
    ])
    j = sp.Matrix([[0, -1], [1, 0]])
    sine = sp.ImmutableMatrix(-amp*j*wavevectors)
    position = sp.Matrix(coords)
    psi = amp*sum(sp.cos(wavevectors[:, i].dot(position)) for i in range(3))
    gradient = sp.Matrix([sp.diff(psi, c) for c in coords])
    velocity = sp.ImmutableMatrix(j*gradient)
    zeta = sum(sp.diff(psi, c, 2) for c in coords)
    pressure = -rho*(velocity.dot(velocity)+lam**2*psi**2)/2
    covariance = sp.ImmutableMatrix(sp.simplify(sine*sine.T/2))
    level = sp.simplify(psi.subs({coords[0]: sp.pi/lam,
                                 coords[1]: sp.pi/(sp.sqrt(3)*lam)}, simultaneous=True))
    dual = sp.ImmutableMatrix(-lam**2*level*j)
    return TriangularEulerArray(
        wavevectors, sine, psi, velocity, sp.trigsimp(zeta), pressure,
        covariance, level, dual, sp.And(*(sp.Gt(v, 0) for v in (amp, lam, rho))),
    )


@dataclass(frozen=True)
class AxialKelvinInitialPhase:
    """Exact initial phase blocks; not an autonomous all-time reduction."""

    mass: sp.ImmutableMatrix
    stiffness: sp.ImmutableMatrix
    symplectic: sp.ImmutableMatrix
    hamiltonian: sp.ImmutableMatrix
    domain_condition: sp.Expr


def axial_kelvin_initial_phase(wavevectors, sine_velocities, axial_wavenumber, density):
    """Derive full-pressure initial phase from a declared planar Fourier field.

    The columns specify distinct real sine modes on the same periodic cell,
    each nonzero and divergence free. Commensurability and stationarity of
    the background, and a periodic Euler pressure, are caller hypotheses;
    this algebraic function proves none of them for an arbitrary spectrum.
    The triangular field returned above supplies them through0161.

    The common velocity has its own circulation data. The displacement
    is Kelvin prepared with BOTH horizontal and vertical pressure returns.
    The phase has physical initial coordinates (X0,V0); its Hamiltonian
    here is the initial restriction, not the later moving-chart generator.
    """
    q, coefficients = map(sp.Matrix, (wavevectors, sine_velocities))
    if q.rows != 2 or not q.cols or coefficients.shape != q.shape:
        raise ValueError("wavevectors and sine velocities must be matching 2 by N matrices")
    k, rho = sp.sympify(axial_wavenumber), _positive(density, "density")
    if any(entry.is_real is False or entry.is_finite is False
           or entry.has(sp.nan, sp.zoo, sp.oo, -sp.oo) for entry in (*q, *coefficients, k)):
        raise ValueError("wavevectors, velocities and axial wavenumber must be finite and real")
    stiffness = sp.zeros(2)
    conditions = [sp.Gt(rho, 0)]
    vertical = sp.Matrix([0, 0, 1])
    for i in range(q.cols):
        qi, ci = q[:, i], coefficients[:, i]
        norm = sp.simplify(qi.dot(qi))
        if norm.is_zero is True:
            raise ValueError("horizontal harmonic modes are not sine modes")
        if sp.simplify(qi.dot(ci)) != 0:
            raise ValueError("each sine velocity must be perpendicular to its wavevector")
        for previous in range(i):
            if any(all(sp.simplify(entry) == 0 for entry in qi+sign*q[:, previous])
                   for sign in (-1, 1)):
                raise ValueError("combine duplicate or opposite Fourier wavevectors first")
            conditions.extend(sp.Ne(sp.simplify((qi+sign*q[:, previous]).dot(
                qi+sign*q[:, previous])), 0) for sign in (-1, 1))
        conditions.append(sp.Gt(norm, 0))
        full_wave = qi.col_join(sp.Matrix([k]))
        projector = sp.eye(3)-full_wave*full_wave.T/(norm+k**2)
        returned = projector*vertical
        stiffness += rho*k**2*ci*ci.T*returned.dot(returned)/2
    mass = rho*sp.eye(2)
    stiffness = sp.simplify(stiffness)
    omega = sp.BlockMatrix([[sp.zeros(2), mass], [-mass, sp.zeros(2)]]).as_explicit()
    hamiltonian = sp.diag(stiffness, mass)
    return AxialKelvinInitialPhase(
        *map(sp.ImmutableMatrix, (mass, stiffness, omega, hamiltonian)),
        sp.And(*conditions),
    )

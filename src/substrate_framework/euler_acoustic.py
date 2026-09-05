"""Actual planar Euler array and full axial Kelvin preparation (P251/0167).

The triangular array's acoustic-window theorem is analytic evidence in0161.
These definitions derive its field and initial phase, not an Euler trajectory.
In particular no stability theorem for its constant-curl lift is implied.
"""

from dataclasses import dataclass

import sympy as sp

from substrate_framework import euler_fourier as ef


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


def _stationary_response_inputs(background, pressure, direction, *macro_vectors):
    """Validate exact periodic inputs; pressure is pressure divided by density."""
    if len(background) != 3 or ef.divergence(background):
        raise ValueError("background must have three solenoidal components")
    if any(sp.simplify(v.get(ef.ZERO, 0)) != 0 for v in background):
        raise ValueError("the prepared response requires a mean-zero background")
    for field in (*background, pressure):
        for wave, value in field.items():
            if (len(wave) != 3 or any(sp.sympify(q).is_Rational is not True for q in wave)
                    or sp.simplify(sp.sympify(value).conjugate()
                                   - field.get(tuple(-q for q in wave), 0)) != 0):
                raise ValueError("fields must be real, commensurate finite Fourier fields")
    gradient = tuple(ef.derivative(pressure, i) for i in range(3))
    advected = ef.transport(background, background)
    if any(ef.add(advected[i], gradient[i]) for i in range(3)):
        raise ValueError("the supplied pressure per density must satisfy stationary Euler")
    kappa = sp.Matrix(direction)
    vectors = [sp.Matrix(v) for v in macro_vectors]
    if kappa.shape != (3, 1) or any(v.shape != (3, 1) for v in vectors):
        raise ValueError("direction and macro vectors must have three components")
    if any(q.is_real is False or q.is_finite is False
           or q.has(sp.nan, sp.zoo, sp.oo, -sp.oo) for v in [kappa, *vectors] for q in v):
        raise ValueError("macro vectors must be finite and real")
    if sp.simplify(kappa.dot(kappa)-1) != 0:
        raise ValueError("the Bloch direction must be a unit vector")
    if any(sp.simplify(kappa.dot(v)) != 0 for v in vectors):
        raise ValueError("macro preparation must be transverse to the common direction")
    return kappa, vectors, gradient


@dataclass(frozen=True)
class PreparedAcousticCellRows:
    """First-cell data: chi(0)=0, chi_t(0), and F(t)=F0+t*F1.

    These rows enter chi_tt+2P(u.grad)chi_t+
    P[(u.grad)^2+Hess(pressure/rho)]chi=F(t).
    The pressure return is the full mean-preserving microscopic Leray.
    No PDE evolution or autonomous constitutive closure is returned.
    """

    initial_rate: tuple
    forcing_constant: tuple
    forcing_rate: tuple


def prepared_acoustic_cell_rows(background, pressure_per_density, direction, displacement, velocity):
    """Derive the actual Kelvin-D / independent-common-V first-cell rows.

    The Fourier convention is exp(i*q.x), K=k*direction and
    eta=D+t*V+i*k*chi+O(k²). The two macro vectors are transverse;
    their initial circulation classes are independently specified.
    Inputs are exact real finite Fourier fields on a common periodic cell.
    Symbolic coefficients denote fixed real physical parameters; the
    defining reality and stationary Euler identities must simplify exactly.
    """
    kappa, (d, v), grad_p = _stationary_response_inputs(
        background, pressure_per_density, direction, displacement, velocity)
    a = ef.add(*(ef.scale(background[i], kappa[i]) for i in range(3)))
    ud = ef.add(*(ef.scale(background[i], d[i]) for i in range(3)))
    kp = ef.add(*(ef.scale(grad_p[i], kappa[i]) for i in range(3)))

    def pressure_row(vector):
        vp = ef.add(*(ef.scale(grad_p[i], vector[i]) for i in range(3)))
        return tuple(ef.add(ef.scale(kp, vector[i]), ef.scale(vp, kappa[i]))
                     for i in range(3))

    initial = ef.leray(tuple(ef.add(ef.scale(a, -d[i]), ef.scale(ud, -kappa[i]))
                             for i in range(3)))
    fd = pressure_row(d)
    constant = ef.leray(tuple(ef.add(fd[i], ef.scale(a, -2*v[i])) for i in range(3)))
    return PreparedAcousticCellRows(initial, constant, ef.leray(pressure_row(v)))


@dataclass(frozen=True)
class ObservedAcousticCellRows:
    """Coefficients of m_t/k² and (m-<eta>_t)/k², respectively."""

    acceleration: sp.ImmutableMatrix
    current_correction: sp.ImmutableMatrix


def observed_acoustic_cell_rows(background, pressure_per_density, direction,
                                zeroth_displacement, cell, cell_rate):
    """Compute the physical second-jet stress and Lin current, not a closure.

    ``cell`` and ``cell_rate`` must be mean-zero solenoidal Fourier fields.
    For trajectory interpretation they must solve the actual first-cell
    problem above; this algebraic function does not assume or verify that
    evolution. The second unknown cell cancels by periodic integration by
    parts. Both pressure rows and the separate harmonic slow projector are
    retained. Current correction is to the mean MATERIAL displacement rate,
    not to X_t when X is defined by integrating the Eulerian velocity mean.
    """
    kappa, (u0,), grad_p = _stationary_response_inputs(
        background, pressure_per_density, direction, zeroth_displacement)
    for field in (cell, cell_rate):
        if len(field) != 3 or ef.divergence(field) or any(
                sp.simplify(component.get(ef.ZERO, 0)) != 0 for component in field):
            raise ValueError("first cell and its rate must be mean-zero and solenoidal")
    a = ef.add(*(ef.scale(background[i], kappa[i]) for i in range(3)))
    kchi = ef.add(*(ef.scale(cell[i], kappa[i]) for i in range(3)))
    krate = ef.add(*(ef.scale(cell_rate[i], kappa[i]) for i in range(3)))
    kp = ef.add(*(ef.scale(grad_p[i], kappa[i]) for i in range(3)))

    def mean_product(left, right):
        return ef.mul(left, right).get(ef.ZERO, 0)

    force = sp.Matrix([
        mean_product(a, a)*u0[i] + mean_product(a, cell_rate[i])
        + mean_product(background[i], krate) + mean_product(kp, cell[i])
        + mean_product(grad_p[i], kchi) for i in range(3)])
    correction = sp.Matrix([-mean_product(a, cell[i])+mean_product(background[i], kchi)
                            for i in range(3)])
    slow = sp.eye(3)-kappa*kappa.T
    return ObservedAcousticCellRows(
        sp.ImmutableMatrix(sp.simplify(slow*force)),
        sp.ImmutableMatrix(sp.simplify(slow*correction)),
    )

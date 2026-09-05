"""Conditional inherited-energy optical preparation (P251/0187, corrected0191).

The jets must come from actual positive-action Euler mode families when making
an Euler claim. This algebra does not prove existence, material-tag identities,
or a common-background EPS/acoustic embedding for arbitrary input numbers.
The response is a spatial second jet on fixed time windows, not an all-K law.
"""

from dataclasses import dataclass

import sympy as sp

from substrate_framework.micropolar import MicropolarCoefficients


def _real(value, name, *, positive=False):
    value = sp.sympify(value)
    if (not isinstance(value, sp.Expr) or value.is_real is False
            or value.is_finite is False or value.has(sp.nan, sp.zoo, sp.oo, -sp.oo)
            or (positive and value.is_positive is False)):
        raise ValueError(f"{name} must be finite, real" + (" and positive" if positive else ""))
    return value


@dataclass(frozen=True)
class OpticalModeJet:
    """Carrier derivatives at a fixed physical frequency, not derivatives of its square.

    ``phase_mass`` and its derivatives are whole-frame densities j=M_raw/3.
    The physical angle is reconstructed as 3 E[n theta], whereas the inherited
    energy is averaged without that factor. Use ``from_raw_tilt`` for raw data.
    Unresolved realness/positivity are mathematical input hypotheses; explicit
    non-real, non-finite or nonpositive frequency/mass values are rejected.
    """

    frequency: sp.Expr
    frequency_first: sp.Expr
    frequency_second: sp.Expr
    phase_mass: sp.Expr
    phase_mass_first: sp.Expr = sp.S.Zero
    phase_mass_second: sp.Expr = sp.S.Zero

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            object.__setattr__(self, name, _real(
                getattr(self, name), name, positive=name in ("frequency", "phase_mass")))

    @classmethod
    def from_raw_tilt(cls, frequency, frequency_first, frequency_second,
                      phase_mass, phase_mass_first=0, phase_mass_second=0):
        """Normalize all three raw transverse-tilt phase-density jets by three."""
        return cls(frequency, frequency_first, frequency_second,
                   *(_real(v, "raw phase mass jet")/3
                     for v in (phase_mass, phase_mass_first, phase_mass_second)))


@dataclass(frozen=True)
class CorrelatedOpticalPreparation:
    """Actual-input second-jet construction; probabilities stay positive.

    ``amplitudes`` are signed prepared initial angles/rates, and ``slopes``
    are their carrier derivatives. Second amplitude derivatives are zero.
    ``J0`` and ``J2`` are the inherited phase-density mass and its second jet,
    including all mode-mass and prepared-amplitude derivatives. J2 need not
    be positive; the truncated mass is used only near zero wavevector.
    """

    modes: tuple[OpticalModeJet, OpticalModeJet]
    weights: tuple[sp.Expr, sp.Expr]
    amplitudes: tuple[sp.Expr, sp.Expr]
    slopes: tuple[sp.Expr, sp.Expr]
    observed_curvature: sp.Expr
    J0: sp.Expr
    J2: sp.Expr
    normalization_residual: sp.Expr
    variance_residual: sp.Expr
    curvature_residual: sp.Expr
    energy_residual: sp.Expr
    slope_matrix: sp.ImmutableMatrix
    slope_rhs: sp.ImmutableMatrix
    domain_condition: sp.Expr


def correlated_optical_preparation(modes, weights, observed_curvature):
    """Solve normalization/variance and then curvature/full-energy equations.

    Two actual families have common frequency nu>0, nonzero slopes v_r with
    unequal squares, positive phase densities j_r, and positive probabilities
    summing to one. The freely declared B_star>0 specifies a spatial-gradient
    PREPARATION, not an empirically fitted force or a unique intrinsic modulus.

    For A_r(p)=a_r+d_r p, set h_r=j_r A_r². The two slope equations are
    B=sum w(a*b+2*d*v)=B_star and
    sum w[h*(v²+nu*b)+2*nu*h'*v] = J0*nu*B_star.
    Both systems are derived and solved here. Input p denotes local carrier
    increment. Symbolically undecidable domains remain in ``domain_condition``.
    """
    modes, weights = tuple(modes), tuple(weights)
    if len(modes) != 2 or any(not isinstance(m, OpticalModeJet) for m in modes):
        raise ValueError("exactly two OpticalModeJet inputs are required")
    if len(weights) != 2:
        raise ValueError("exactly two positive probabilities are required")
    weights = tuple(_real(w, "probability", positive=True) for w in weights)
    if sp.simplify(sum(weights)-1) != 0:
        raise ValueError("probabilities must sum exactly to one")
    nu = modes[0].frequency
    if sp.simplify(nu-modes[1].frequency) != 0:
        raise ValueError("mode frequencies must agree exactly")
    target = _real(observed_curvature, "observed curvature", positive=True)
    velocities = tuple(m.frequency_first for m in modes)
    difference = sp.simplify(velocities[0]**2-velocities[1]**2)
    if any(v.is_zero is True for v in velocities) or difference.is_zero is True:
        raise ValueError("frequency slopes must be nonzero with unequal squares")

    amplitude_matrix = sp.Matrix([weights, [w*v**2 for w, v in zip(weights, velocities)]])
    amplitudes = tuple(sp.simplify(v) for v in amplitude_matrix.inv()*sp.Matrix([1, 0]))
    x = sp.Dummy("carrier_increment", real=True)
    d = sp.symbols("d1 d2", real=True, cls=sp.Dummy)
    mass_polynomials = [m.phase_mass+m.phase_mass_first*x+m.phase_mass_second*x**2/2
                        for m in modes]
    h = [j*(a+di*x)**2 for j, a, di in zip(mass_polynomials, amplitudes, d)]
    h0 = [hi.subs(x, 0) for hi in h]
    hp = [sp.diff(hi, x).subs(x, 0) for hi in h]
    J0 = sp.simplify(sum(w*hi for w, hi in zip(weights, h0)))
    curvature = sum(w*(a*m.frequency_second+2*di*m.frequency_first)
                    for w, a, di, m in zip(weights, amplitudes, d, modes))
    energy = sum(w*(hi*(m.frequency_first**2+nu*m.frequency_second)
                    +2*nu*dhi*m.frequency_first)
                 for w, hi, dhi, m in zip(weights, h0, hp, modes)) - J0*nu*target
    matrix, rhs = sp.linear_eq_to_matrix([curvature-target, energy], d)
    slopes = tuple(sp.simplify(v) for v in matrix.inv()*rhs)
    solved = dict(zip(d, slopes))
    J2 = sp.simplify(sum(w*sp.diff(hi, x, 2).subs(x, 0)
                        for w, hi in zip(weights, h)).subs(solved))
    conditions = [sp.Gt(nu, 0), sp.Gt(target, 0), sp.Ne(difference, 0)]
    conditions += [sp.Gt(w, 0) for w in weights]
    conditions += [sp.Gt(m.phase_mass, 0) for m in modes]
    conditions += [sp.Ne(v, 0) for v in velocities]
    return CorrelatedOpticalPreparation(
        modes, weights, amplitudes, slopes, target, J0, J2,
        sp.simplify(sum(w*a for w, a in zip(weights, amplitudes))-1),
        sp.simplify(sum(w*a*v**2 for w, a, v in zip(weights, amplitudes, velocities))),
        sp.simplify((curvature-target).subs(solved)), sp.simplify(energy.subs(solved)),
        sp.ImmutableMatrix(matrix), sp.ImmutableMatrix(rhs), sp.And(*conditions),
    )


def transverse_tilt_tensor(wave_vector):
    """Return 3 E[n n^T (t.K)²] for Haar frames with measured tilt n perpendicular t.

    The carrier axis is t, not the tilt n. Eigenvalues are 2|K|²/5 on
    transverse vectors and |K|²/5 longitudinally. No orientation integral
    is implied to construct the missing same-cell Euler background.
    """
    k = sp.Matrix(wave_vector)
    if k.shape != (3, 1):
        raise ValueError("wave vector must have three components")
    for component in k:
        _real(component, "wave-vector component")
    return sp.ImmutableMatrix((2*k.dot(k)*sp.eye(3)-k*k.T)/5)


@dataclass(frozen=True)
class OpticalActionJet:
    """Second-order matrices in L=(Phi_t^T M Phi_t-Phi^T K Phi)/2.

    The gradient mass and its frequency-squared contribution are retained.
    Matrix products are truncated consistently through degree two in K.
    """

    mass: sp.ImmutableMatrix
    stiffness: sp.ImmutableMatrix
    frequency_squared: sp.ImmutableMatrix


def transverse_optical_action(preparation, wave_vector):
    """Construct corrected whole-frame action matrices, not an all-K dispersion."""
    if not isinstance(preparation, CorrelatedOpticalPreparation):
        raise ValueError("a CorrelatedOpticalPreparation is required")
    tensor = transverse_tilt_tensor(wave_vector)
    nu = preparation.modes[0].frequency
    mass = preparation.J0*sp.eye(3)+preparation.J2*tensor/2
    frequency = nu**2*sp.eye(3)+nu*preparation.observed_curvature*tensor
    stiffness = nu**2*mass+preparation.J0*nu*preparation.observed_curvature*tensor
    return OpticalActionJet(*map(sp.ImmutableMatrix, (mass, stiffness, frequency)))


def positive_optical_curvature(transverse, longitudinal):
    """Positive pointwise representative with specified positive bulk symbols.

    Returns canonical MicropolarCoefficients with other sectors zero. Its
    trace coefficient may be negative; the positive trace-sector coefficient
    is c_tr+c_s/3. Changing representatives also changes the explicit boundary
    null-Lagrangian flux; this function supplies no acoustic/current join.
    """
    ct = _real(transverse, "transverse curvature", positive=True)
    cl = _real(longitudinal, "longitudinal curvature", positive=True)
    denominator = 4*ct+3*cl
    symmetric = 3*ct*cl/denominator
    skew = 4*ct**2/denominator
    trace = cl/2-symmetric
    return MicropolarCoefficients(0, 0, 0, trace, symmetric, skew)

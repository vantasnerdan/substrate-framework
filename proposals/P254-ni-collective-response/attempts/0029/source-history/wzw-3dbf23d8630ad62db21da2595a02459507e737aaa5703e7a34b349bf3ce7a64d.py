"""Exact SU(3) trace forms, generator periods, and filling algebra.

This module contains mathematical differential-form statements only.  Its
normalized winding current and ``S^5`` coefficient lattice use explicit
oriented generators.  Calling a winding a baryon number, coupling the current
to an external source, or identifying a WZW coefficient with a color or
anomaly coefficient requires separate governed physics.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
from itertools import combinations, permutations
from typing import Any, Sequence

import sympy as sp

from .su3 import fundamental_generators, structure_constant


CochainIndex = tuple[int, ...]
StructureConstants = tuple[tuple[tuple[sp.Expr, ...], ...], ...]


@dataclass(frozen=True)
class SU3TraceFiveCohomology:
    """Exact finite-dimensional evidence for the invariant trace five-form."""

    algebra_dimension: int
    four_cochain_dimension: int
    five_cochain_dimension: int
    six_cochain_dimension: int
    d4_rank: int
    d5_rank: int
    five_cocycle_dimension: int
    fifth_cohomology_dimension: int
    trace_nonzero_components: int
    trace_norm_squared: sp.Expr
    augmented_d4_trace_rank: int
    dual_separating_pairing: sp.Expr
    differential_squares_to_zero: bool
    trace_is_closed: bool
    trace_is_exact: bool
    trace_annihilates_coboundaries: bool


@dataclass(frozen=True)
class SU3Pi5PeriodEvidence:
    """Exact normalization evidence for the oriented Puttmann--Rigas map.

    The unit sphere has the boundary orientation inherited from
    ``(Re z1, Im z1, Re z2, Im z2, Re z3, Im z3)``.  With that convention the
    explicit map below has column-projection degree ``+2`` and represents the
    positive generator selected by that degree.  The accepted real form is
    negative on the corresponding oriented tangent frame, so its primitive
    period is negative; reversing either orientation reverses the sign.
    """

    projection_degree: int
    positive_preimage_jacobian: sp.Expr
    negative_preimage_jacobian: sp.Expr
    raw_trace_density: sp.Expr
    real_trace_density: sp.Expr
    sphere_volume: sp.Expr
    raw_trace_period: sp.Expr
    real_trace_period: sp.Expr
    primitive_period_magnitude: sp.Expr
    coefficient_lattice_step: sp.Expr


@dataclass(frozen=True)
class SU3WindingThreeEvidence:
    """Exact cohomology and generator evidence for the trace three-form.

    The quaternion sphere has the boundary orientation inherited from
    ``(a0,a1,a2,a3)``.  Its explicit SU(2) block embedding has first-column
    degree ``+1`` and raw trace-three period ``+24*pi**2``.  The normalized
    form used by :func:`su3_winding_current` includes a minus sign, so that
    this positive generator has period ``-1`` and the decreasing hedgehog
    convention has charge ``+1``.
    """

    d2_rank: int
    d3_rank: int
    three_cocycle_dimension: int
    third_cohomology_dimension: int
    trace_nonzero_components: int
    trace_norm_squared: sp.Expr
    augmented_d2_trace_rank: int
    differential_squares_to_zero: bool
    trace_is_closed: bool
    trace_is_exact: bool
    column_projection_jacobian: sp.Expr
    column_projection_degree: int
    raw_generator_density: sp.Expr
    sphere_volume: sp.Expr
    raw_generator_period: sp.Expr
    normalized_generator_period: sp.Expr
    current_coefficient: sp.Expr


def cochain_basis(degree: int, dimension: int = 8) -> tuple[CochainIndex, ...]:
    """Return the ordered exterior-cochain basis of a declared dimension."""

    if not isinstance(degree, int) or not isinstance(dimension, int):
        raise TypeError("degree and dimension must be integers")
    if dimension < 0 or degree < 0 or degree > dimension:
        raise ValueError("require 0 <= degree <= dimension")
    return tuple(combinations(range(dimension), degree))


def _alternating_sign(indices: Sequence[int]) -> int:
    if len(set(indices)) != len(indices):
        return 0
    inversions = sum(
        indices[first] > indices[second]
        for first in range(len(indices))
        for second in range(first + 1, len(indices))
    )
    return -1 if inversions % 2 else 1


def alternating_trace(matrices: Sequence[Any]) -> sp.Expr:
    """Return the unnormalized alternating trace of square matrices.

    For matrix-valued one-forms this is the coefficient convention
    ``sum_perm sgn(perm) Tr(M_p1 ... M_pk)``.  No hidden ``1/k!`` is applied.
    """

    values = tuple(sp.Matrix(matrix) for matrix in matrices)
    if not values:
        raise ValueError("at least one matrix is required")
    shape = values[0].shape
    if shape[0] != shape[1] or any(value.shape != shape for value in values):
        raise ValueError("matrices must be square and share one shape")
    total = sp.Integer(0)
    for order in permutations(range(len(values))):
        product = sp.eye(shape[0])
        for index in order:
            product *= values[index]
        total += _alternating_sign(order) * sp.trace(product)
    return sp.simplify(total)


def su2_quaternion_embedding(point: Sequence[Any]) -> sp.ImmutableMatrix:
    r"""Embed a unit quaternion in the upper ``SU(2)`` block of ``SU(3)``.

    For ``a=(a0,a1,a2,a3)``, the upper block is
    ``a0*I + i*(a1*sigma1+a2*sigma2+a3*sigma3)`` and the last diagonal entry
    is one.  Exact algebra gives determinant ``sum(a_j**2)`` and unitarity on
    the real unit three-sphere.
    """

    values = tuple(sp.sympify(value) for value in point)
    if len(values) != 4:
        raise ValueError("a quaternion point must have four real coordinates")
    a0, a1, a2, a3 = values
    return sp.ImmutableMatrix(
        [
            [a0 + sp.I * a3, a2 + sp.I * a1, 0],
            [-a2 + sp.I * a1, a0 - sp.I * a3, 0],
            [0, 0, 1],
        ]
    )


def su2_quaternion_embedding_differential(
    tangent: Sequence[Any],
) -> sp.ImmutableMatrix:
    """Return the differential of :func:`su2_quaternion_embedding`."""

    values = tuple(sp.sympify(value) for value in tangent)
    if len(values) != 4:
        raise ValueError("a quaternion tangent must have four real coordinates")
    v0, v1, v2, v3 = values
    return sp.ImmutableMatrix(
        [
            [v0 + sp.I * v3, v2 + sp.I * v1, 0],
            [-v2 + sp.I * v1, v0 - sp.I * v3, 0],
            [0, 0, 0],
        ]
    )


@cache
def su2_quaternion_column_projection_jacobian() -> sp.Expr:
    """Derive the oriented Jacobian of the quaternion first-column map.

    Real coordinates on the first complex column are ordered as
    ``(Re z0, Im z0, Re z1, Im z1)``.  The exact determinant fixes the
    degree sign independently of the trace-three period.
    """

    a0, a1, a2, a3 = sp.symbols("a0 a1 a2 a3", real=True)
    first_column = su2_quaternion_embedding((a0, a1, a2, a3))[:, 0]
    target_coordinates = sp.ImmutableMatrix(
        [
            sp.re(first_column[0]),
            sp.im(first_column[0]),
            sp.re(first_column[1]),
            sp.im(first_column[1]),
        ]
    )
    return sp.simplify(
        target_coordinates.jacobian((a0, a1, a2, a3)).det()
    )


@cache
def su2_quaternion_trace_three_period() -> sp.Expr:
    """Derive the raw trace-three period on the oriented quaternion sphere."""

    north = (1, 0, 0, 0)
    positive_frame = (
        (0, 1, 0, 0),
        (0, 0, 1, 0),
        (0, 0, 0, 1),
    )
    group_value = su2_quaternion_embedding(north)
    maurer_cartan_values = tuple(
        sp.ImmutableMatrix(
            group_value.H * su2_quaternion_embedding_differential(tangent)
        )
        for tangent in positive_frame
    )
    density = alternating_trace(maurer_cartan_values)
    sphere_volume = 2 * sp.pi**2 / sp.gamma(2)
    return sp.simplify(density * sphere_volume)


@cache
def su3_winding_current_coefficient() -> sp.Expr:
    """Return the coefficient selecting charge ``+1`` for a decreasing hedgehog.

    The coefficient is derived as minus the reciprocal of the explicit
    positive-generator raw period rather than inserted as an expected answer.
    """

    return sp.simplify(-1 / su2_quaternion_trace_three_period())


def su3_winding_current(
    left_currents: Sequence[Any], orientation: int = 1
) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    r"""Return the normalized topological current from four left currents.

    With ``epsilon**0123=orientation`` and ``L_mu=U^dagger*d_mu U``, this is

    ``J^mu=c*epsilon^(mu nu rho sigma)*Tr(L_nu L_rho L_sigma)``,

    where ``c`` is derived by :func:`su3_winding_current_coefficient`.
    The default sign makes a hedgehog with ``F(0)=pi`` and ``F(infinity)=0``
    carry charge ``+1``.  This is a mathematical winding current, not by
    itself a Noether, gauged-WZW-response, or physical baryon current.
    """

    if orientation not in (-1, 1):
        raise ValueError("orientation must be +1 or -1")
    values = tuple(sp.Matrix(value) for value in left_currents)
    if len(values) != 4:
        raise ValueError("exactly four spacetime left currents are required")
    shape = values[0].shape
    if shape[0] != shape[1] or any(value.shape != shape for value in values):
        raise ValueError("left currents must be square matrices of one shape")
    coefficient = su3_winding_current_coefficient()
    components: list[sp.Expr] = []
    for mu in range(4):
        complement = tuple(index for index in range(4) if index != mu)
        contraction = sp.Integer(0)
        for order in permutations(complement):
            sign = orientation * sp.LeviCivita(mu, *order)
            contraction += sign * sp.trace(
                values[order[0]] * values[order[1]] * values[order[2]]
            )
        components.append(sp.simplify(coefficient * contraction))
    return tuple(components)  # type: ignore[return-value]


def hedgehog_winding_density(profile: Any, radius: sp.Symbol) -> sp.Expr:
    """Return the local static hedgehog density for positive ``radius``.

    The convention is the default one in :func:`su3_winding_current`.  The
    apparent ``radius**-2`` is a spherical-coordinate density; smooth
    hedgehog regularity makes its volume-weighted radial form finite.
    """

    if not isinstance(radius, sp.Symbol):
        raise TypeError("radius must be a SymPy symbol")
    field = sp.sympify(profile)
    return sp.simplify(
        -sp.sin(field) ** 2 * sp.diff(field, radius)
        / (2 * sp.pi**2 * radius**2)
    )


def hedgehog_winding_radial_density(profile: Any, radius: sp.Symbol) -> sp.Expr:
    """Return ``4*pi*r**2`` times :func:`hedgehog_winding_density`."""

    return sp.simplify(4 * sp.pi * radius**2 * hedgehog_winding_density(profile, radius))


def hedgehog_winding_charge(inner_value: Any, outer_value: Any) -> sp.Expr:
    """Return the exact hedgehog charge determined by its endpoint values."""

    inner = sp.sympify(inner_value)
    outer = sp.sympify(outer_value)

    def primitive(value: sp.Expr) -> sp.Expr:
        return value - sp.sin(value) * sp.cos(value)

    return sp.simplify((primitive(inner) - primitive(outer)) / sp.pi)


def _complex_cross_matrix(vector: Sequence[Any]) -> sp.ImmutableMatrix:
    """Return the conjugate-linear matrix used by the complex cross product."""

    values = tuple(sp.sympify(value) for value in vector)
    if len(values) != 3:
        raise ValueError("the complex cross-product vector must have length three")
    first, second, third = (sp.conjugate(value) for value in values)
    return sp.ImmutableMatrix(
        [
            [0, -third, second],
            [third, 0, -first],
            [-second, first, 0],
        ]
    )


def su3_pi5_generator(point: Sequence[Any]) -> sp.ImmutableMatrix:
    r"""Return the Puttmann--Rigas map ``eta:S^5 -> SU(3)``.

    For a unit vector ``z`` in ``C^3``,

    ``eta(z) = z*z^T + cross_matrix(conjugate(z))``.

    The transpose is deliberately not a Hermitian transpose.  Exact algebra
    gives ``det eta(z)=|z|^4`` and ``eta(z)^dagger eta(z)=I`` on ``|z|=1``.
    Puttmann and Rigas, Theorem 2.1, prove that this smooth embedding generates
    ``pi_5(SU(3))``; their independent criterion is degree ``+2`` of a column
    projection, reproduced exactly by :func:`su3_pi5_period_evidence`.
    """

    values = tuple(sp.sympify(value) for value in point)
    if len(values) != 3:
        raise ValueError("an S^5 point must have three complex coordinates")
    column = sp.ImmutableMatrix(3, 1, values)
    return sp.ImmutableMatrix(column * column.T + _complex_cross_matrix(values))


def su3_pi5_generator_differential(
    point: Sequence[Any], tangent: Sequence[Any]
) -> sp.ImmutableMatrix:
    """Return ``D eta_z(v)`` for the explicit ``pi_5`` generator map."""

    point_values = tuple(sp.sympify(value) for value in point)
    tangent_values = tuple(sp.sympify(value) for value in tangent)
    if len(point_values) != 3 or len(tangent_values) != 3:
        raise ValueError("point and tangent must each have three complex coordinates")
    z = sp.ImmutableMatrix(3, 1, point_values)
    v = sp.ImmutableMatrix(3, 1, tangent_values)
    return sp.ImmutableMatrix(v * z.T + z * v.T + _complex_cross_matrix(tangent_values))


def _first_column_differential(
    point: Sequence[Any], tangent: Sequence[Any]
) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(su3_pi5_generator_differential(point, tangent)[:, 0])


def _oriented_target_coordinates(vector: Sequence[Any]) -> sp.ImmutableMatrix:
    """Coordinates in the positive tangent frame at ``(1,0,0) in S^5``."""

    values = tuple(sp.sympify(value) for value in vector)
    return sp.ImmutableMatrix(
        [
            sp.im(values[0]),
            sp.re(values[1]),
            sp.im(values[1]),
            sp.re(values[2]),
            sp.im(values[2]),
        ]
    )


@cache
def su3_pi5_period_evidence() -> SU3Pi5PeriodEvidence:
    """Derive the exact primitive ``S^5`` period in C-WZW-001's convention.

    The first column has the regular value ``(1,0,0)`` at exactly the two
    points ``(+/-1,0,0)``.  Its two oriented real Jacobians are computed rather
    than inserted, so their positive signs give degree ``+2`` independently of
    the trace-five integral.  Equivariance of the explicit map makes its
    pullback of the bi-invariant trace form an invariant top form on ``S^5``;
    one exact tangent-frame value therefore fixes the integral.
    """

    target = (sp.Integer(1), sp.Integer(0), sp.Integer(0))
    jacobians: list[sp.Expr] = []
    for sign in (sp.Integer(1), sp.Integer(-1)):
        point = (sign, sp.Integer(0), sp.Integer(0))
        domain_frame = (
            (sign * sp.I, 0, 0),
            (0, 1, 0),
            (0, sp.I, 0),
            (0, 0, 1),
            (0, 0, sp.I),
        )
        jacobian = sp.Matrix.hstack(
            *[
                _oriented_target_coordinates(
                    _first_column_differential(point, tangent)
                )
                for tangent in domain_frame
            ]
        )
        jacobians.append(sp.factor(jacobian.det()))

    positive_frame = (
        (sp.I, 0, 0),
        (0, 1, 0),
        (0, sp.I, 0),
        (0, 0, 1),
        (0, 0, sp.I),
    )
    eta = su3_pi5_generator(target)
    maurer_cartan_values = tuple(
        sp.ImmutableMatrix(
            eta.H * su3_pi5_generator_differential(target, tangent)
        )
        for tangent in positive_frame
    )
    raw_density = sp.simplify(alternating_trace(maurer_cartan_values))
    real_density = sp.simplify(-sp.I * raw_density)
    sphere_volume = sp.simplify(2 * sp.pi**3 / sp.gamma(3))
    raw_period = sp.simplify(raw_density * sphere_volume)
    real_period = sp.simplify(real_density * sphere_volume)
    magnitude = sp.simplify(-real_period)
    coefficient_step = sp.simplify(2 * sp.pi / magnitude)
    degree = sum(1 if value > 0 else -1 for value in jacobians)
    return SU3Pi5PeriodEvidence(
        projection_degree=degree,
        positive_preimage_jacobian=jacobians[0],
        negative_preimage_jacobian=jacobians[1],
        raw_trace_density=raw_density,
        real_trace_density=real_density,
        sphere_volume=sphere_volume,
        raw_trace_period=raw_period,
        real_trace_period=real_period,
        primitive_period_magnitude=magnitude,
        coefficient_lattice_step=coefficient_step,
    )


def su3_sphere_trace_five_period(winding: Any = 1) -> sp.Expr:
    """Return the real trace-five period for an oriented ``pi_5`` winding.

    The result concerns maps from the oriented five-sphere.  It does not assert
    the period lattice for arbitrary closed five-manifolds or add spin/bordism
    premises that are absent from the framework.
    """

    return sp.simplify(sp.sympify(winding) * su3_pi5_period_evidence().real_trace_period)


def sphere_extension_coefficient(level: Any = 1) -> sp.Expr:
    """Return the coefficient whose ``S^5`` ambiguity is ``2*pi*level``.

    An integer ``level`` makes the phase independent of two ball fillings of a
    common oriented ``S^4`` boundary.  The name is mathematical: no physical
    identification of ``level`` is made here.
    """

    return sp.simplify(
        sp.sympify(level) * su3_pi5_period_evidence().coefficient_lattice_step
    )


def sphere_extension_phase_ratio(level: Any, winding: Any) -> sp.Expr:
    """Return the phase ratio for the declared level and sphere winding."""

    return sp.simplify(
        sp.exp(
            sp.I
            * sphere_extension_coefficient(level)
            * su3_sphere_trace_five_period(winding)
        )
    )


@cache
def antihermitian_generators() -> tuple[sp.ImmutableMatrix, ...]:
    """Return ``E_a=i*T_a`` for C-LIE-001's Hermitian generators."""

    return tuple(sp.ImmutableMatrix(sp.I * value) for value in fundamental_generators())


@cache
def antihermitian_structure_constants() -> StructureConstants:
    """Return real ``c_ab^d`` for ``[E_a,E_b]=c_ab^d E_d``.

    Since ``E_a=i*T_a`` and C-LIE-001 uses
    ``[T_a,T_b]=i*f_abc*T_c``, this basis has ``c_ab^d=-f_abd``.
    """

    return tuple(
        tuple(
            tuple(sp.simplify(-structure_constant(a, b, d)) for d in range(8))
            for b in range(8)
        )
        for a in range(8)
    )


def _validated_structure_constants(values: Any) -> StructureConstants:
    constants = tuple(
        tuple(tuple(sp.sympify(entry) for entry in row) for row in plane)
        for plane in values
    )
    dimension = len(constants)
    if dimension == 0 or any(len(plane) != dimension for plane in constants):
        raise ValueError("structure constants must be a nonempty cubic array")
    if any(len(row) != dimension for plane in constants for row in plane):
        raise ValueError("structure constants must be a nonempty cubic array")
    return constants


def _build_ce_differential(
    degree: int, constants: StructureConstants
) -> sp.ImmutableSparseMatrix:
    """Build ``d:C^degree -> C^(degree+1)`` in increasing-index bases."""

    dimension = len(constants)
    if not isinstance(degree, int):
        raise TypeError("degree must be an integer")
    if degree < 0 or degree >= dimension:
        raise ValueError("require 0 <= degree < algebra dimension")
    domain = cochain_basis(degree, dimension)
    codomain = cochain_basis(degree + 1, dimension)
    domain_position = {indices: position for position, indices in enumerate(domain)}
    entries: dict[tuple[int, int], sp.Expr] = {}

    for row, arguments in enumerate(codomain):
        for first in range(degree + 1):
            for second in range(first + 1, degree + 1):
                remaining = [
                    arguments[position]
                    for position in range(degree + 1)
                    if position not in (first, second)
                ]
                pair_sign = (-1) ** (first + second)
                for bracket_index, coefficient in enumerate(
                    constants[arguments[first]][arguments[second]]
                ):
                    if coefficient == 0:
                        continue
                    cochain_arguments = [bracket_index, *remaining]
                    reorder_sign = _alternating_sign(cochain_arguments)
                    if reorder_sign == 0:
                        continue
                    column = domain_position[tuple(sorted(cochain_arguments))]
                    key = (row, column)
                    entries[key] = sp.simplify(
                        entries.get(key, 0)
                        + pair_sign * reorder_sign * coefficient
                    )
    return sp.ImmutableSparseMatrix(len(codomain), len(domain), entries)


@cache
def _canonical_ce_differential(degree: int) -> sp.ImmutableSparseMatrix:
    return _build_ce_differential(degree, antihermitian_structure_constants())


def chevalley_eilenberg_differential(
    degree: int,
    constants: Any | None = None,
) -> sp.ImmutableSparseMatrix:
    """Return the exact invariant-cochain differential.

    The default is the anti-Hermitian SU(3) basis derived from C-LIE-001.
    An explicit cubic structure-constant array may be supplied for independent
    algebra checks and mutation-sensitive verification.
    """

    if constants is None:
        return _canonical_ce_differential(degree)
    return _build_ce_differential(degree, _validated_structure_constants(constants))


@cache
def su3_trace_power_cochain(degree: int) -> sp.ImmutableDenseMatrix:
    """Return ``Alt Tr(E_a1...E_ak)`` on increasing SU(3) basis tuples."""

    basis = cochain_basis(degree, 8)
    generators = antihermitian_generators()
    return sp.ImmutableDenseMatrix(
        [alternating_trace(tuple(generators[index] for index in item)) for item in basis]
    )


@cache
def su3_real_trace_five_cochain() -> sp.ImmutableDenseMatrix:
    """Return the real cochain ``-i Alt Tr(theta^5)``.

    The factor ``-i`` converts the imaginary trace of five anti-Hermitian
    generators into a real differential form.  Period normalization remains
    deliberately unspecified.
    """

    return sp.ImmutableDenseMatrix(
        [sp.simplify(-sp.I * value) for value in su3_trace_power_cochain(5)]
    )


@cache
def su3_winding_three_evidence() -> SU3WindingThreeEvidence:
    """Derive exact cohomology and generator normalization for the winding form."""

    d2 = chevalley_eilenberg_differential(2)
    d3 = chevalley_eilenberg_differential(3)
    trace_three = su3_trace_power_cochain(3)
    d2_rank = int(d2.rank())
    d3_rank = int(d3.rank())
    augmented_rank = int(d2.row_join(trace_three).rank())
    cocycle_dimension = len(cochain_basis(3, 8)) - d3_rank

    north = (1, 0, 0, 0)
    positive_frame = (
        (0, 1, 0, 0),
        (0, 0, 1, 0),
        (0, 0, 0, 1),
    )
    group_value = su2_quaternion_embedding(north)
    maurer_cartan_values = tuple(
        sp.ImmutableMatrix(
            group_value.H * su2_quaternion_embedding_differential(tangent)
        )
        for tangent in positive_frame
    )
    raw_density = alternating_trace(maurer_cartan_values)
    sphere_volume = sp.simplify(2 * sp.pi**2 / sp.gamma(2))
    raw_period = sp.simplify(raw_density * sphere_volume)
    coefficient = su3_winding_current_coefficient()
    column_jacobian = su2_quaternion_column_projection_jacobian()
    return SU3WindingThreeEvidence(
        d2_rank=d2_rank,
        d3_rank=d3_rank,
        three_cocycle_dimension=cocycle_dimension,
        third_cohomology_dimension=cocycle_dimension - d2_rank,
        trace_nonzero_components=sum(value != 0 for value in trace_three),
        trace_norm_squared=sp.simplify((trace_three.T * trace_three)[0]),
        augmented_d2_trace_rank=augmented_rank,
        differential_squares_to_zero=(d3 * d2 == sp.zeros(d3.rows, d2.cols)),
        trace_is_closed=(d3 * trace_three == sp.zeros(d3.rows, 1)),
        trace_is_exact=(augmented_rank == d2_rank),
        column_projection_jacobian=column_jacobian,
        column_projection_degree=int(sp.sign(column_jacobian)),
        raw_generator_density=raw_density,
        sphere_volume=sphere_volume,
        raw_generator_period=raw_period,
        normalized_generator_period=sp.simplify(coefficient * raw_period),
        current_coefficient=coefficient,
    )


def trace_power_cyclic_shift_sign(power: int) -> int:
    """Return the graded sign for moving one odd one-form past ``power-1``."""

    if not isinstance(power, int):
        raise TypeError("power must be an integer")
    if power < 1:
        raise ValueError("power must be positive")
    return (-1) ** (power - 1)


def trace_power_derivative_multiplier(power: int) -> int:
    """Return the coefficient of ``Tr(dL L^(power-1))`` in ``d Tr(L^power)``.

    The exterior derivative contributes alternating Leibniz signs.  After
    graded cyclic reordering the coefficient is ``power`` for odd powers and
    zero for even powers.  In particular, the tempting formula
    ``d Tr(L^4)=4 Tr(dL L^3)`` is false: all four terms cancel.
    """

    if not isinstance(power, int):
        raise TypeError("power must be an integer")
    if power < 1:
        raise ValueError("power must be positive")
    return sum(
        (-1) ** (position + position * (power + 1 - position))
        for position in range(power)
    )


def maurer_cartan_power_derivative_multiplier(power: int) -> int:
    """Return the coefficient in ``d(L^power)=c L^(power+1)`` under ``dL=-L^2``.

    This is a matrix-valued identity before taking a trace.  The coefficient
    is zero for even ``power`` and minus one for odd ``power``.
    """

    if not isinstance(power, int):
        raise TypeError("power must be an integer")
    if power < 1:
        raise ValueError("power must be positive")
    return -sum((-1) ** position for position in range(power))


@cache
def su3_trace_five_cohomology() -> SU3TraceFiveCohomology:
    """Derive exact invariant and global non-exactness evidence.

    The returned rank statements concern the left-invariant
    Chevalley-Eilenberg complex.  For compact SU(3), if this left-invariant
    form had a global primitive, averaging that primitive with normalized Haar
    measure would give a left-invariant primitive.  Thus ``trace_is_exact`` is
    also the global de Rham exactness verdict under that standard averaging
    theorem.
    """

    d4 = chevalley_eilenberg_differential(4)
    d5 = chevalley_eilenberg_differential(5)
    trace_five = su3_real_trace_five_cochain()
    d4_rank = int(d4.rank())
    d5_rank = int(d5.rank())
    augmented_rank = int(d4.row_join(trace_five).rank())
    trace_norm_squared = sp.simplify((trace_five.T * trace_five)[0])
    dual_pairing = trace_norm_squared
    trace_annihilates = trace_five.T * d4 == sp.zeros(1, d4.cols)
    five_cocycle_dimension = len(cochain_basis(5, 8)) - d5_rank
    return SU3TraceFiveCohomology(
        algebra_dimension=8,
        four_cochain_dimension=len(cochain_basis(4, 8)),
        five_cochain_dimension=len(cochain_basis(5, 8)),
        six_cochain_dimension=len(cochain_basis(6, 8)),
        d4_rank=d4_rank,
        d5_rank=d5_rank,
        five_cocycle_dimension=five_cocycle_dimension,
        fifth_cohomology_dimension=five_cocycle_dimension - d4_rank,
        trace_nonzero_components=sum(value != 0 for value in trace_five),
        trace_norm_squared=trace_norm_squared,
        augmented_d4_trace_rank=augmented_rank,
        dual_separating_pairing=dual_pairing,
        differential_squares_to_zero=(d5 * d4 == sp.zeros(d5.rows, d4.cols)),
        trace_is_closed=(d5 * trace_five == sp.zeros(d5.rows, 1)),
        trace_is_exact=(augmented_rank == d4_rank),
        trace_annihilates_coboundaries=trace_annihilates,
    )


def glued_filling_period(first_integral: Any, second_integral: Any) -> sp.Expr:
    """Return the oriented closed-cycle integral ``I(B)-I(B')``.

    This is a gluing identity, not a claim that the period is nonzero,
    integral, or normalized to a particular value.
    """

    return sp.simplify(sp.sympify(first_integral) - sp.sympify(second_integral))


def extension_phase_ratio(
    coefficient: Any,
    first_integral: Any,
    second_integral: Any,
) -> sp.Expr:
    """Return the phase ratio between two declared extension integrals."""

    return sp.simplify(
        sp.exp(
            sp.I
            * sp.sympify(coefficient)
            * glued_filling_period(first_integral, second_integral)
        )
    )

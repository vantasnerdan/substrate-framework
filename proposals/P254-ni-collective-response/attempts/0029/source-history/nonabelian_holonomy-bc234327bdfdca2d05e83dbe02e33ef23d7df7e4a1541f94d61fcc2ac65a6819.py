"""Exact finite-matrix parallel-transport and holonomy ledgers.

The convention is inherited from ``D = partial - i*g*W``.  Chronological
integrated connection segments ``B_1,...,B_n`` therefore transport with later
segments on the left,
``V = exp(i*B_n) ... exp(i*B_1)``.  The APIs expose ordering, reversal,
commuting collapse, endpoint gauge covariance, closed-loop conjugacy data, and
SU(2) representation dependence.  They do not construct a physical path,
gauge action, flux, matter carrier, Aharonov--Bohm experiment, weak sector, or
substrate dictionary.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp


def _zero(value: Any) -> bool:
    if isinstance(value, sp.MatrixBase):
        return all(sp.simplify(entry) == 0 for entry in value)
    return sp.simplify(value) == 0


def _exact_square_matrix(value: Any, name: str) -> sp.ImmutableMatrix:
    matrix = sp.Matrix(value)
    if matrix.rows == 0 or matrix.cols == 0 or matrix.rows != matrix.cols:
        raise ValueError(f"{name} must be a nonempty square matrix")
    if any(entry.has(sp.Float) for entry in matrix):
        raise ValueError(f"{name} must contain exact entries")
    return sp.ImmutableMatrix(matrix.applyfunc(sp.simplify))


def _exact_hermitian(value: Any, name: str) -> sp.ImmutableMatrix:
    matrix = _exact_square_matrix(value, name)
    if not _zero(matrix - matrix.H):
        raise ValueError(f"{name} must be exactly Hermitian")
    return matrix


def _exact_unitary(value: Any, name: str) -> sp.ImmutableMatrix:
    matrix = _exact_square_matrix(value, name)
    if not _zero(matrix.H * matrix - sp.eye(matrix.rows)):
        raise ValueError(f"{name} must be exactly unitary")
    return matrix


def _transport_product(
    chronological_factors: Sequence[sp.ImmutableMatrix],
) -> sp.ImmutableMatrix:
    dimension = chronological_factors[0].rows
    product = sp.ImmutableMatrix(sp.eye(dimension))
    for factor in chronological_factors:
        product = sp.ImmutableMatrix((factor * product).applyfunc(sp.simplify))
    return product


def _factor(integrated_connection: sp.ImmutableMatrix) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(
        (sp.I * integrated_connection).exp().applyfunc(sp.simplify)
    )


@dataclass(frozen=True)
class OrderedSegmentHolonomyEvidence:
    """Exact evidence for an ordered sequence of integrated connections."""

    segments: tuple[sp.ImmutableMatrix, ...]
    segment_transporters: tuple[sp.ImmutableMatrix, ...]
    transporter: sp.ImmutableMatrix
    unitarity_residual: sp.ImmutableMatrix
    determinant: sp.Expr
    determinant_formula: sp.Expr
    determinant_residual: sp.Expr
    trace: sp.Expr
    normalized_trace: sp.Expr
    pairwise_commutators: tuple[sp.ImmutableMatrix, ...]
    pairwise_commuting: bool
    exponential_of_sum: sp.ImmutableMatrix
    commuting_collapse_residual: sp.ImmutableMatrix
    reverse_segments: tuple[sp.ImmutableMatrix, ...]
    reverse_transporter: sp.ImmutableMatrix
    reverse_inverse_residual: sp.ImmutableMatrix
    composition_residuals: tuple[sp.ImmutableMatrix, ...]
    cyclic_shifted_transporter: sp.ImmutableMatrix
    cyclic_conjugation_residual: sp.ImmutableMatrix
    cyclic_trace_residual: sp.Expr

    @property
    def unitary_certified(self) -> bool:
        return _zero(self.unitarity_residual)

    @property
    def determinant_certified(self) -> bool:
        return _zero(self.determinant_residual)

    @property
    def reverse_certified(self) -> bool:
        return _zero(self.reverse_inverse_residual)

    @property
    def composition_certified(self) -> bool:
        return all(_zero(residual) for residual in self.composition_residuals)

    @property
    def cyclic_basepoint_certified(self) -> bool:
        return _zero(self.cyclic_conjugation_residual) and _zero(
            self.cyclic_trace_residual
        )


def ordered_segment_holonomy(
    segments: Sequence[Any],
) -> OrderedSegmentHolonomyEvidence:
    """Return later-left transport for exact chronological Hermitian segments.

    Each input is the already integrated matrix ``B_j`` along one declared
    oriented segment, including any coupling.  This function checks no path
    geometry and does not infer that the sequence forms a physical closed loop.
    """

    if not segments:
        raise ValueError("segments must be nonempty")
    exact_segments = tuple(
        _exact_hermitian(value, f"segment_{index}")
        for index, value in enumerate(segments)
    )
    dimension = exact_segments[0].rows
    if any(segment.shape != (dimension, dimension) for segment in exact_segments):
        raise ValueError("segments must have the same square shape")

    factors = tuple(_factor(segment) for segment in exact_segments)
    transporter = _transport_product(factors)
    identity = sp.ImmutableMatrix(sp.eye(dimension))
    unitarity_residual = sp.ImmutableMatrix(
        (transporter.H * transporter - identity).applyfunc(sp.simplify)
    )
    determinant = sp.simplify(transporter.det())
    determinant_formula = sp.simplify(
        sp.exp(sp.I * sum((sp.trace(segment) for segment in exact_segments), sp.S.Zero))
    )

    commutators = tuple(
        sp.ImmutableMatrix(
            (exact_segments[left] * exact_segments[right]
             - exact_segments[right] * exact_segments[left]).applyfunc(sp.simplify)
        )
        for left in range(len(exact_segments))
        for right in range(left + 1, len(exact_segments))
    )
    pairwise_commuting = all(_zero(commutator) for commutator in commutators)
    summed_segment = sum(exact_segments, sp.zeros(dimension))
    exponential_of_sum = _factor(sp.ImmutableMatrix(summed_segment))
    collapse_residual = sp.ImmutableMatrix(
        (transporter - exponential_of_sum).applyfunc(sp.simplify)
    )

    reverse_segments = tuple(-segment for segment in reversed(exact_segments))
    reverse_factors = tuple(_factor(segment) for segment in reverse_segments)
    reverse_transporter = _transport_product(reverse_factors)
    reverse_residual = sp.ImmutableMatrix(
        (reverse_transporter - transporter.H).applyfunc(sp.simplify)
    )

    composition_residuals = []
    for cut in range(1, len(factors)):
        prefix = _transport_product(factors[:cut])
        suffix = _transport_product(factors[cut:])
        composition_residuals.append(
            sp.ImmutableMatrix((transporter - suffix * prefix).applyfunc(sp.simplify))
        )

    shifted_factors = factors[1:] + factors[:1]
    shifted_transporter = _transport_product(shifted_factors)
    first_factor = factors[0]
    cyclic_prediction = sp.ImmutableMatrix(
        (first_factor * transporter * first_factor.H).applyfunc(sp.simplify)
    )

    return OrderedSegmentHolonomyEvidence(
        segments=exact_segments,
        segment_transporters=factors,
        transporter=transporter,
        unitarity_residual=unitarity_residual,
        determinant=determinant,
        determinant_formula=determinant_formula,
        determinant_residual=sp.simplify(determinant - determinant_formula),
        trace=sp.simplify(sp.trace(transporter)),
        normalized_trace=sp.simplify(sp.trace(transporter) / dimension),
        pairwise_commutators=commutators,
        pairwise_commuting=pairwise_commuting,
        exponential_of_sum=exponential_of_sum,
        commuting_collapse_residual=collapse_residual,
        reverse_segments=tuple(sp.ImmutableMatrix(value) for value in reverse_segments),
        reverse_transporter=reverse_transporter,
        reverse_inverse_residual=reverse_residual,
        composition_residuals=tuple(composition_residuals),
        cyclic_shifted_transporter=shifted_transporter,
        cyclic_conjugation_residual=sp.ImmutableMatrix(
            (shifted_transporter - cyclic_prediction).applyfunc(sp.simplify)
        ),
        cyclic_trace_residual=sp.simplify(
            sp.trace(shifted_transporter) - sp.trace(transporter)
        ),
    )


@dataclass(frozen=True)
class EndpointGaugeHolonomyEvidence:
    """Endpoint-factor covariance for supplied segment transporters."""

    closed_path: bool
    segment_transporters: tuple[sp.ImmutableMatrix, ...]
    node_gauges: tuple[sp.ImmutableMatrix, ...]
    transporter: sp.ImmutableMatrix
    transformed_segments: tuple[sp.ImmutableMatrix, ...]
    transformed_transporter: sp.ImmutableMatrix
    endpoint_prediction: sp.ImmutableMatrix
    endpoint_covariance_residual: sp.ImmutableMatrix
    trace_residual: sp.Expr | None
    determinant_residual: sp.Expr | None
    characteristic_polynomial_residual: sp.Expr | None

    @property
    def endpoint_covariance_certified(self) -> bool:
        return _zero(self.endpoint_covariance_residual)

    @property
    def closed_conjugacy_certified(self) -> bool:
        return (
            self.closed_path
            and self.trace_residual is not None
            and self.determinant_residual is not None
            and self.characteristic_polynomial_residual is not None
            and _zero(self.trace_residual)
            and _zero(self.determinant_residual)
            and _zero(self.characteristic_polynomial_residual)
        )


def endpoint_gauge_holonomy_evidence(
    segment_transporters: Sequence[Any],
    node_gauges: Sequence[Any],
    *,
    closed_path: bool = False,
) -> EndpointGaugeHolonomyEvidence:
    """Telescope segment endpoint gauges into the full transporter law.

    Chronological segment transporters ``E_j`` are transformed as
    ``E_j' = U_j E_j U_(j-1)^dagger``.  Their later-left product therefore
    obeys ``V' = U_n V U_0^dagger``.  For a declared closed path this function
    requires the endpoint gauge matrices to agree and also returns conjugacy-
    invariant trace, determinant, and characteristic-polynomial residuals.
    """

    if not segment_transporters:
        raise ValueError("segment_transporters must be nonempty")
    factors = tuple(
        _exact_unitary(value, f"segment_transporter_{index}")
        for index, value in enumerate(segment_transporters)
    )
    dimension = factors[0].rows
    if any(factor.shape != (dimension, dimension) for factor in factors):
        raise ValueError("segment_transporters must have the same square shape")
    if len(node_gauges) != len(factors) + 1:
        raise ValueError("node_gauges must contain one more matrix than segments")
    gauges = tuple(
        _exact_unitary(value, f"node_gauge_{index}")
        for index, value in enumerate(node_gauges)
    )
    if any(gauge.shape != (dimension, dimension) for gauge in gauges):
        raise ValueError("node_gauges must match the segment dimension")
    if closed_path and not _zero(gauges[-1] - gauges[0]):
        raise ValueError("closed_path requires equal endpoint gauge matrices")

    transporter = _transport_product(factors)
    transformed_segments = tuple(
        sp.ImmutableMatrix(
            (gauges[index + 1] * factor * gauges[index].H).applyfunc(sp.simplify)
        )
        for index, factor in enumerate(factors)
    )
    transformed = _transport_product(transformed_segments)
    prediction = sp.ImmutableMatrix(
        (gauges[-1] * transporter * gauges[0].H).applyfunc(sp.simplify)
    )
    covariance_residual = sp.ImmutableMatrix(
        (transformed - prediction).applyfunc(sp.simplify)
    )

    trace_residual: sp.Expr | None = None
    determinant_residual: sp.Expr | None = None
    characteristic_residual: sp.Expr | None = None
    if closed_path:
        spectral_parameter = sp.Symbol("lambda")
        trace_residual = sp.simplify(sp.trace(transformed) - sp.trace(transporter))
        determinant_residual = sp.simplify(transformed.det() - transporter.det())
        characteristic_residual = sp.simplify(
            transformed.charpoly(spectral_parameter).as_expr()
            - transporter.charpoly(spectral_parameter).as_expr()
        )

    return EndpointGaugeHolonomyEvidence(
        closed_path=closed_path,
        segment_transporters=factors,
        node_gauges=gauges,
        transporter=transporter,
        transformed_segments=transformed_segments,
        transformed_transporter=transformed,
        endpoint_prediction=prediction,
        endpoint_covariance_residual=covariance_residual,
        trace_residual=trace_residual,
        determinant_residual=determinant_residual,
        characteristic_polynomial_residual=characteristic_residual,
    )


@dataclass(frozen=True)
class SU2HolonomyEvidence:
    """Fundamental center, adjoint image, and noncommuting-order controls."""

    parameter: sp.Symbol
    generators: tuple[sp.ImmutableMatrix, sp.ImmutableMatrix, sp.ImmutableMatrix]
    fundamental_2pi: sp.ImmutableMatrix
    fundamental_4pi: sp.ImmutableMatrix
    adjoint_2pi: sp.ImmutableMatrix
    fundamental_2pi_trace: sp.Expr
    fundamental_2pi_normalized_trace: sp.Expr
    adjoint_2pi_trace: sp.Expr
    adjoint_2pi_normalized_trace: sp.Expr
    ordered_noncommuting: sp.ImmutableMatrix
    naive_exponential_of_sum: sp.ImmutableMatrix
    ordered_minus_naive: sp.ImmutableMatrix
    leading_quadratic_coefficient: sp.ImmutableMatrix
    expected_quadratic_coefficient: sp.ImmutableMatrix
    commuting_ordered: sp.ImmutableMatrix
    commuting_exponential_of_sum: sp.ImmutableMatrix
    commuting_residual: sp.ImmutableMatrix

    @property
    def bch_coefficient_certified(self) -> bool:
        return _zero(
            self.leading_quadratic_coefficient
            - self.expected_quadratic_coefficient
        )


def su2_holonomy_evidence(
    parameter: sp.Symbol | None = None,
) -> SU2HolonomyEvidence:
    """Return exact SU(2) controls without assigning a physical weak sector."""

    if parameter is None:
        angle = sp.Symbol("a", real=True)
    else:
        if not isinstance(parameter, sp.Symbol) or parameter.is_real is not True:
            raise ValueError("parameter must be an explicitly real Symbol")
        angle = parameter

    sigma_1 = sp.ImmutableMatrix([[0, 1], [1, 0]])
    sigma_2 = sp.ImmutableMatrix([[0, -sp.I], [sp.I, 0]])
    sigma_3 = sp.ImmutableMatrix([[1, 0], [0, -1]])
    generators = tuple(
        sp.ImmutableMatrix(matrix / 2)
        for matrix in (sigma_1, sigma_2, sigma_3)
    )
    t_1, t_2, t_3 = generators
    identity_2 = sp.ImmutableMatrix(sp.eye(2))
    fundamental_2pi = sp.ImmutableMatrix((sp.I * 2 * sp.pi * t_3).exp())
    fundamental_4pi = sp.ImmutableMatrix((sp.I * 4 * sp.pi * t_3).exp())
    adjoint_generator_3 = sp.ImmutableMatrix(sp.diag(1, 0, -1))
    adjoint_2pi = sp.ImmutableMatrix((sp.I * 2 * sp.pi * adjoint_generator_3).exp())

    factor_1 = sp.ImmutableMatrix(
        sp.cos(angle / 2) * identity_2
        + 2 * sp.I * sp.sin(angle / 2) * t_1
    )
    factor_2 = sp.ImmutableMatrix(
        sp.cos(angle / 2) * identity_2
        + 2 * sp.I * sp.sin(angle / 2) * t_2
    )
    ordered = sp.ImmutableMatrix((factor_2 * factor_1).applyfunc(sp.simplify))
    naive = sp.ImmutableMatrix(
        (
            sp.cos(angle / sp.sqrt(2)) * identity_2
            + sp.I
            * sp.sqrt(2)
            * sp.sin(angle / sp.sqrt(2))
            * (t_1 + t_2)
        ).applyfunc(sp.simplify)
    )
    difference = sp.ImmutableMatrix((ordered - naive).applyfunc(sp.simplify))
    leading = sp.ImmutableMatrix(
        difference.applyfunc(
            lambda entry: sp.simplify(sp.diff(entry, angle, 2).subs(angle, 0) / 2)
        )
    )
    expected = sp.ImmutableMatrix(sp.I * t_3 / 2)

    factor_3 = sp.ImmutableMatrix(
        sp.cos(angle / 2) * identity_2
        + 2 * sp.I * sp.sin(angle / 2) * t_3
    )
    commuting_ordered = sp.ImmutableMatrix((factor_3 * factor_3).applyfunc(sp.simplify))
    commuting_naive = sp.ImmutableMatrix(
        (sp.cos(angle) * identity_2 + 2 * sp.I * sp.sin(angle) * t_3).applyfunc(
            sp.simplify
        )
    )

    return SU2HolonomyEvidence(
        parameter=angle,
        generators=generators,
        fundamental_2pi=fundamental_2pi,
        fundamental_4pi=fundamental_4pi,
        adjoint_2pi=adjoint_2pi,
        fundamental_2pi_trace=sp.simplify(sp.trace(fundamental_2pi)),
        fundamental_2pi_normalized_trace=sp.simplify(
            sp.trace(fundamental_2pi) / 2
        ),
        adjoint_2pi_trace=sp.simplify(sp.trace(adjoint_2pi)),
        adjoint_2pi_normalized_trace=sp.simplify(sp.trace(adjoint_2pi) / 3),
        ordered_noncommuting=ordered,
        naive_exponential_of_sum=naive,
        ordered_minus_naive=difference,
        leading_quadratic_coefficient=leading,
        expected_quadratic_coefficient=expected,
        commuting_ordered=commuting_ordered,
        commuting_exponential_of_sum=commuting_naive,
        commuting_residual=sp.ImmutableMatrix(
            (commuting_ordered - commuting_naive).applyfunc(sp.simplify)
        ),
    )

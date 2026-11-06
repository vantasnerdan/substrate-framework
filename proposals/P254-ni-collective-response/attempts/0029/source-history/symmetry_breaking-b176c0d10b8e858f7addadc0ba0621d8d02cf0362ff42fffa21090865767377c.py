"""Exact finite-dimensional evidence for stationary symmetry zero modes.

This module proves a classical quadratic statement for declared scalar
coordinates.  It does not by itself construct a quantum vacuum, a Goldstone
particle, a chiral sector, or a physical-pion identification.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp


def _symbols(fields: Sequence[sp.Symbol]) -> tuple[sp.Symbol, ...]:
    result = tuple(fields)
    if not result:
        raise ValueError("fields must be non-empty")
    if any(not isinstance(field, sp.Symbol) for field in result):
        raise ValueError("fields must be SymPy symbols")
    if len(set(result)) != len(result):
        raise ValueError("fields must be unique")
    return result


def _square_matrices(
    matrices: Sequence[Any],
    dimension: int,
    *,
    name: str,
) -> tuple[sp.ImmutableMatrix, ...]:
    result = tuple(sp.ImmutableMatrix(matrix) for matrix in matrices)
    if not result:
        raise ValueError(f"{name} must be non-empty")
    if any(matrix.shape != (dimension, dimension) for matrix in result):
        raise ValueError(f"every {name} matrix must match the field dimension")
    return result


def _immutable_simplified(matrix: sp.MatrixBase) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(matrix.applyfunc(sp.simplify))


def _is_zero_matrix(matrix: sp.MatrixBase) -> bool:
    return _immutable_simplified(matrix) == sp.zeros(*matrix.shape)


@dataclass(frozen=True)
class LinearSymmetryHessianEvidence:
    """Exact objects entering the stationary symmetry-Hessian theorem.

    ``coefficient_kernel_dimension`` is the kernel dimension of the map from
    the supplied generator coefficients to vacuum tangents.  It is the
    stabilizer dimension only when the supplied generators are independent.
    """

    fields: tuple[sp.Symbol, ...]
    vacuum: sp.ImmutableMatrix
    generators: tuple[sp.ImmutableMatrix, ...]
    gradient: sp.ImmutableMatrix
    invariance_residuals: tuple[sp.Expr, ...]
    stationarity_residual: sp.ImmutableMatrix
    hessian: sp.ImmutableMatrix
    hessian_at_vacuum: sp.ImmutableMatrix
    generator_tangents: sp.ImmutableMatrix
    differentiated_invariance_at_vacuum: sp.ImmutableMatrix
    differentiated_identity_residual: sp.ImmutableMatrix
    hessian_tangent_residual: sp.ImmutableMatrix
    generator_span_rank: int
    broken_tangent_rank: int
    coefficient_kernel_dimension: int

    @property
    def invariant(self) -> bool:
        """Whether every supplied infinitesimal invariance residual is zero."""

        return all(sp.simplify(residual) == 0 for residual in self.invariance_residuals)

    @property
    def stationary(self) -> bool:
        """Whether the declared vacuum is an exact stationary point."""

        return _is_zero_matrix(self.stationarity_residual)

    @property
    def generators_independent(self) -> bool:
        """Whether the supplied matrices form an independent generator basis."""

        return self.generator_span_rank == len(self.generators)

    @property
    def theorem_hypotheses_hold(self) -> bool:
        """Whether exact invariance and exact stationarity both hold."""

        return self.invariant and self.stationary

    @property
    def tangent_kernel_certified(self) -> bool:
        """Whether the actual generator tangents lie in the Hessian kernel."""

        return _is_zero_matrix(self.hessian_tangent_residual)

    @property
    def stabilizer_dimension(self) -> int:
        """Return stabilizer dimension for an independent supplied basis."""

        if not self.generators_independent:
            raise ValueError(
                "stabilizer dimension requires an independent generator basis"
            )
        return self.coefficient_kernel_dimension


def linear_symmetry_hessian_evidence(
    potential: Any,
    fields: Sequence[sp.Symbol],
    vacuum: Sequence[Any],
    generators: Sequence[Any],
) -> LinearSymmetryHessianEvidence:
    """Evaluate the exact infinitesimal-symmetry Hessian identities.

    For a linear generator ``T`` the invariance residual is
    ``grad(V).T*T*phi``.  Its field gradient obeys
    ``grad(residual) = Hess(V)*T*phi + T.T*grad(V)``.  Thus exact invariance
    and stationarity imply that every actual vacuum tangent ``T*vacuum`` is a
    Hessian zero direction.  The tangent-matrix rank, rather than printed
    group dimensions, counts the independent directions certified this way.

    The function returns residuals even when a premise fails so mutations and
    counterexamples remain observable.
    """

    field_tuple = _symbols(fields)
    dimension = len(field_tuple)
    vacuum_tuple = tuple(sp.sympify(value) for value in vacuum)
    if len(vacuum_tuple) != dimension:
        raise ValueError("vacuum must match the field dimension")
    generator_tuple = _square_matrices(
        generators,
        dimension,
        name="generators",
    )

    expression = sp.sympify(potential)
    field_column = sp.ImmutableMatrix(field_tuple)
    vacuum_column = sp.ImmutableMatrix(vacuum_tuple)
    substitutions = dict(zip(field_tuple, vacuum_tuple, strict=True))
    gradient = _immutable_simplified(
        sp.Matrix([sp.diff(expression, field) for field in field_tuple])
    )
    hessian = _immutable_simplified(sp.hessian(expression, field_tuple))
    invariance_residuals = tuple(
        sp.simplify((gradient.T * generator * field_column)[0])
        for generator in generator_tuple
    )
    stationarity = _immutable_simplified(gradient.subs(substitutions))
    hessian_at_vacuum = _immutable_simplified(hessian.subs(substitutions))
    tangents = _immutable_simplified(
        sp.Matrix.hstack(*(generator * vacuum_column for generator in generator_tuple))
    )
    differentiated = _immutable_simplified(
        sp.Matrix.hstack(
            *(
                sp.Matrix(
                    [sp.diff(residual, field) for field in field_tuple]
                ).subs(substitutions)
                for residual in invariance_residuals
            )
        )
    )
    identity_right = _immutable_simplified(
        sp.Matrix.hstack(
            *(
                hessian_at_vacuum * tangents[:, index]
                + generator.T * stationarity
                for index, generator in enumerate(generator_tuple)
            )
        )
    )
    hessian_tangent_residual = _immutable_simplified(
        hessian_at_vacuum * tangents
    )
    flattened_generators = sp.Matrix.hstack(
        *(sp.Matrix(generator).reshape(dimension**2, 1) for generator in generator_tuple)
    )
    return LinearSymmetryHessianEvidence(
        fields=field_tuple,
        vacuum=vacuum_column,
        generators=generator_tuple,
        gradient=gradient,
        invariance_residuals=invariance_residuals,
        stationarity_residual=stationarity,
        hessian=hessian,
        hessian_at_vacuum=hessian_at_vacuum,
        generator_tangents=tangents,
        differentiated_invariance_at_vacuum=differentiated,
        differentiated_identity_residual=_immutable_simplified(
            differentiated - identity_right
        ),
        hessian_tangent_residual=hessian_tangent_residual,
        generator_span_rank=int(flattened_generators.rank()),
        broken_tangent_rank=int(tangents.rank()),
        coefficient_kernel_dimension=len(generator_tuple) - int(tangents.rank()),
    )


def orthogonal_generators(dimension: int) -> tuple[sp.ImmutableMatrix, ...]:
    """Return the standard independent antisymmetric basis of ``so(dimension)``."""

    if not isinstance(dimension, int) or isinstance(dimension, bool) or dimension < 2:
        raise ValueError("dimension must be an integer at least two")
    result: list[sp.ImmutableMatrix] = []
    for first in range(dimension):
        for second in range(first + 1, dimension):
            generator = sp.zeros(dimension)
            generator[first, second] = 1
            generator[second, first] = -1
            result.append(sp.ImmutableMatrix(generator))
    return tuple(result)


def radial_quartic_potential(
    fields: Sequence[sp.Symbol],
    coupling: Any,
    vacuum_scale: Any,
) -> sp.Expr:
    """Return ``coupling*(sum(fields**2)-vacuum_scale**2)**2`` exactly."""

    field_tuple = _symbols(fields)
    radius_squared = sum((field**2 for field in field_tuple), sp.Integer(0))
    return sp.sympify(coupling) * (radius_squared - sp.sympify(vacuum_scale) ** 2) ** 2


@dataclass(frozen=True)
class PositiveKineticMassEvidence:
    """Generalized quadratic mass operator for a proven positive kinetic metric."""

    hessian: sp.ImmutableMatrix
    kinetic_metric: sp.ImmutableMatrix
    zero_directions: sp.ImmutableMatrix
    generalized_mass_operator: sp.ImmutableMatrix
    zero_direction_residual: sp.ImmutableMatrix
    zero_direction_rank: int

    @property
    def zero_directions_certified(self) -> bool:
        """Whether all supplied directions have zero generalized mass."""

        return _is_zero_matrix(self.zero_direction_residual)


def positive_kinetic_mass_evidence(
    hessian: Any,
    kinetic_metric: Any,
    zero_directions: Any,
) -> PositiveKineticMassEvidence:
    """Convert Hessian zeros into zeros of ``K**-1*H`` under positive ``K``.

    Positivity must be decidable by SymPy; an unknown symbolic sign is rejected
    instead of being silently promoted to a physical kinetic premise.
    """

    hessian_matrix = sp.ImmutableMatrix(hessian)
    kinetic_matrix = sp.ImmutableMatrix(kinetic_metric)
    directions = sp.ImmutableMatrix(zero_directions)
    if hessian_matrix.rows == 0 or hessian_matrix.rows != hessian_matrix.cols:
        raise ValueError("hessian must be non-empty and square")
    dimension = hessian_matrix.rows
    if kinetic_matrix.shape != (dimension, dimension):
        raise ValueError("kinetic metric must match the Hessian dimension")
    if directions.rows != dimension or directions.cols == 0:
        raise ValueError("zero directions must be a non-empty matching matrix")
    if not _is_zero_matrix(hessian_matrix - hessian_matrix.T):
        raise ValueError("hessian must be symmetric")
    if not _is_zero_matrix(kinetic_matrix - kinetic_matrix.T):
        raise ValueError("kinetic metric must be symmetric")
    if kinetic_matrix.is_positive_definite is not True:
        raise ValueError("kinetic metric must be provably positive definite")
    operator = _immutable_simplified(kinetic_matrix.inv() * hessian_matrix)
    residual = _immutable_simplified(operator * directions)
    return PositiveKineticMassEvidence(
        hessian=hessian_matrix,
        kinetic_metric=kinetic_matrix,
        zero_directions=directions,
        generalized_mass_operator=operator,
        zero_direction_residual=residual,
        zero_direction_rank=int(directions.rank()),
    )


@dataclass(frozen=True)
class GroupCoordinateKineticEvidence:
    """Leading quadratic metric of a declared exponential group coordinate."""

    generators: tuple[sp.ImmutableMatrix, ...]
    coordinate_scale: sp.Expr
    action_prefactor: sp.Expr
    trace_gram: sp.ImmutableMatrix
    kinetic_metric: sp.ImmutableMatrix


def leading_exponential_kinetic_metric(
    generators: Sequence[Any],
    coordinate_scale: Any,
    action_prefactor: Any,
) -> GroupCoordinateKineticEvidence:
    """Derive the leading metric for ``U=exp(i*phi_a*T_a/scale)``.

    For ``L=prefactor*Tr(dU*dU.H)``, the leading scalar convention
    ``L=(1/2)*dphi.T*K*dphi`` gives
    ``K_ab=2*prefactor*Tr(T_a*T_b)/scale**2``.  Generators must be a nonempty
    same-size Hermitian family.  This is a coordinate-model calculation, not a
    derivation of a physical field dictionary or of the action itself.
    """

    raw_generators = tuple(sp.ImmutableMatrix(generator) for generator in generators)
    if not raw_generators:
        raise ValueError("generators must be non-empty")
    dimension = raw_generators[0].rows
    if dimension == 0 or any(
        generator.shape != (dimension, dimension) for generator in raw_generators
    ):
        raise ValueError("generators must be same-size non-empty square matrices")
    if any(not _is_zero_matrix(generator - generator.H) for generator in raw_generators):
        raise ValueError("generators must be Hermitian")
    scale = sp.sympify(coordinate_scale)
    if sp.simplify(scale) == 0:
        raise ValueError("coordinate scale must be nonzero")
    prefactor = sp.sympify(action_prefactor)
    count = len(raw_generators)
    gram = _immutable_simplified(
        sp.Matrix(
            count,
            count,
            lambda row, column: sp.trace(
                raw_generators[row] * raw_generators[column]
            ),
        )
    )
    metric = _immutable_simplified(2 * prefactor * gram / scale**2)
    return GroupCoordinateKineticEvidence(
        generators=raw_generators,
        coordinate_scale=scale,
        action_prefactor=prefactor,
        trace_gram=gram,
        kinetic_metric=metric,
    )

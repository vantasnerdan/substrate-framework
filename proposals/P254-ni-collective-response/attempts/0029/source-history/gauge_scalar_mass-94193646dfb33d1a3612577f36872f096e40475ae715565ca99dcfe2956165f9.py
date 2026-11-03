"""Exact quadratic gauge-scalar mass ledgers for declared finite carriers.

The functions in this module start from a separately declared scalar kinetic
term, Hermitian carrier generators, couplings, and a vacuum vector.  They
derive a finite-dimensional gauge-field quadratic form.  They do not derive a
scalar potential, prove that the vector is a ground state or condensate,
construct a gauge kinetic action, establish a spectral pole or Higgs particle,
identify physical photon or weak-boson fields, select the Standard Model, or
realize a substrate mass mechanism.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp

from .su2_doublets import su2_fundamental_ledger


def _immutable_simplified(matrix: sp.MatrixBase) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(sp.Matrix(matrix).applyfunc(sp.simplify))


def _zero_matrix(matrix: sp.MatrixBase) -> bool:
    return _immutable_simplified(matrix) == sp.zeros(*matrix.shape)


def _exact_matrix(value: Any, name: str) -> sp.ImmutableMatrix:
    matrix = sp.Matrix(value)
    if matrix.rows == 0 or matrix.cols == 0:
        raise ValueError(f"{name} must be nonempty")
    if any(entry.has(sp.Float) for entry in matrix):
        raise ValueError(f"{name} must contain exact entries")
    return _immutable_simplified(matrix)


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be explicitly real")
    return sp.simplify(expression)


def _realify(matrix: sp.MatrixBase) -> sp.ImmutableMatrix:
    real = sp.Matrix(matrix).applyfunc(lambda entry: sp.simplify(sp.re(entry)))
    imaginary = sp.Matrix(matrix).applyfunc(lambda entry: sp.simplify(sp.im(entry)))
    result = _immutable_simplified(sp.Matrix.vstack(real, imaginary))
    if any(entry.is_real is not True for entry in result):
        raise ValueError("realified matrix entries must be explicitly real")
    return result


def _kernel_matrix(matrix: sp.MatrixBase) -> sp.ImmutableMatrix:
    basis = matrix.nullspace()
    if not basis:
        return sp.ImmutableMatrix.zeros(matrix.cols, 0)
    return _immutable_simplified(sp.Matrix.hstack(*basis))


@dataclass(frozen=True)
class GaugeScalarMassEvidence:
    """Exact gauge-orbit Gram data from a declared scalar kinetic term.

    Gauge coefficients are real.  ``real_orbit_map`` therefore stacks real
    and imaginary carrier components before its kernel is computed; using a
    complex coefficient nullspace would answer a different stabilizer
    question.
    """

    vacuum: sp.ImmutableMatrix
    generators: tuple[sp.ImmutableMatrix, ...]
    couplings: tuple[sp.Expr, ...]
    coupled_orbit_vectors: sp.ImmutableMatrix
    real_orbit_map: sp.ImmutableMatrix
    real_generator_map: sp.ImmutableMatrix
    anticommutator_mass_matrix: sp.ImmutableMatrix
    gram_mass_matrix: sp.ImmutableMatrix
    gram_identity_residual: sp.ImmutableMatrix
    coupled_stabilizer_kernel: sp.ImmutableMatrix
    stabilizer_orbit_residual: sp.ImmutableMatrix
    mass_kernel_residual: sp.ImmutableMatrix
    generator_span_rank: int
    orbit_rank: int
    mass_rank: int
    coefficient_kernel_dimension: int
    all_couplings_nonzero: bool

    @property
    def mass_matrix(self) -> sp.ImmutableMatrix:
        """Return the matrix in ``quadratic_density=A.T*M2*A/2``."""

        return self.anticommutator_mass_matrix

    @property
    def gram_identity_certified(self) -> bool:
        """Whether the anticommutator and twice-real-Gram routes agree."""

        return _zero_matrix(self.gram_identity_residual)

    @property
    def generator_basis_independent(self) -> bool:
        """Whether the supplied generators are independent over real coefficients."""

        return self.generator_span_rank == len(self.generators)

    @property
    def stabilizer_kernel_certified(self) -> bool:
        """Whether the returned real coefficient kernel kills both maps."""

        return _zero_matrix(self.stabilizer_orbit_residual) and _zero_matrix(
            self.mass_kernel_residual
        )

    @property
    def stabilizer_dimension(self) -> int:
        """Return stabilizer dimension in an independent nonzero-coupling basis."""

        if not self.generator_basis_independent:
            raise ValueError(
                "stabilizer dimension requires an independent generator basis"
            )
        if not self.all_couplings_nonzero:
            raise ValueError("stabilizer dimension requires nonzero couplings")
        return self.coefficient_kernel_dimension

    def quadratic_density(self, gauge_fields: Sequence[Any]) -> sp.Expr:
        """Return the exact real-field quadratic density ``A.T*M2*A/2``."""

        fields = tuple(
            _exact_real(field, f"gauge_fields[{index}]")
            for index, field in enumerate(gauge_fields)
        )
        if len(fields) != len(self.generators):
            raise ValueError("gauge_fields must match the generator count")
        column = sp.ImmutableMatrix(fields)
        return sp.simplify((column.T * self.mass_matrix * column)[0] / 2)


def gauge_scalar_mass_evidence(
    generators: Sequence[Any],
    couplings: Sequence[Any],
    vacuum: Sequence[Any],
) -> GaugeScalarMassEvidence:
    """Derive the exact gauge-field mass Gram matrix at a declared vacuum.

    For ``D phi0=-i*sum_a(g_a*A_a*T_a)*phi0`` and real gauge fields, the
    scalar kinetic term gives

    ``M2_ab=g_a*g_b*phi0.H*(T_a*T_b+T_b*T_a)*phi0``

    in the convention ``L_quadratic=A.T*M2*A/2``.  Equivalently, with
    ``u_a=g_a*T_a*phi0``, ``M2`` is twice the real Gram matrix of the ``u_a``.
    This identity certifies positive semidefiniteness and makes its kernel the
    coupled vacuum stabilizer.
    """

    generator_tuple = tuple(
        _exact_matrix(generator, f"generators[{index}]")
        for index, generator in enumerate(generators)
    )
    if not generator_tuple:
        raise ValueError("generators must be nonempty")
    dimension = generator_tuple[0].rows
    if any(generator.shape != (dimension, dimension) for generator in generator_tuple):
        raise ValueError("generators must be same-size square matrices")
    if any(not _zero_matrix(generator - generator.H) for generator in generator_tuple):
        raise ValueError("generators must be Hermitian")

    coupling_tuple = tuple(
        _exact_real(coupling, f"couplings[{index}]")
        for index, coupling in enumerate(couplings)
    )
    if len(coupling_tuple) != len(generator_tuple):
        raise ValueError("couplings must match the generator count")

    vacuum_matrix = _exact_matrix(vacuum, "vacuum")
    if vacuum_matrix.shape != (dimension, 1):
        raise ValueError("vacuum must be a matching column vector")

    orbit = _immutable_simplified(
        sp.Matrix.hstack(
            *(
                coupling * generator * vacuum_matrix
                for generator, coupling in zip(
                    generator_tuple, coupling_tuple, strict=True
                )
            )
        )
    )
    real_orbit = _realify(orbit)
    real_generator_columns = tuple(
        _realify(sp.Matrix(generator).reshape(dimension**2, 1))
        for generator in generator_tuple
    )
    real_generator_map = _immutable_simplified(
        sp.Matrix.hstack(*real_generator_columns)
    )

    anticommutator = _immutable_simplified(
        sp.Matrix(
            len(generator_tuple),
            len(generator_tuple),
            lambda first, second: sp.simplify(
                coupling_tuple[first]
                * coupling_tuple[second]
                * (
                    vacuum_matrix.H
                    * (
                        generator_tuple[first] * generator_tuple[second]
                        + generator_tuple[second] * generator_tuple[first]
                    )
                    * vacuum_matrix
                )[0]
            ),
        )
    )
    if any(entry.is_real is not True for entry in anticommutator):
        raise ValueError("mass matrix entries must be explicitly real")
    gram = _immutable_simplified(2 * real_orbit.T * real_orbit)
    stabilizer_kernel = _kernel_matrix(real_orbit)

    orbit_rank = int(real_orbit.rank())
    return GaugeScalarMassEvidence(
        vacuum=vacuum_matrix,
        generators=generator_tuple,
        couplings=coupling_tuple,
        coupled_orbit_vectors=orbit,
        real_orbit_map=real_orbit,
        real_generator_map=real_generator_map,
        anticommutator_mass_matrix=anticommutator,
        gram_mass_matrix=gram,
        gram_identity_residual=_immutable_simplified(anticommutator - gram),
        coupled_stabilizer_kernel=stabilizer_kernel,
        stabilizer_orbit_residual=_immutable_simplified(
            real_orbit * stabilizer_kernel
        ),
        mass_kernel_residual=_immutable_simplified(
            anticommutator * stabilizer_kernel
        ),
        generator_span_rank=int(real_generator_map.rank()),
        orbit_rank=orbit_rank,
        mass_rank=int(anticommutator.rank()),
        coefficient_kernel_dimension=len(generator_tuple) - orbit_rank,
        all_couplings_nonzero=all(
            coupling.is_zero is False for coupling in coupling_tuple
        ),
    )


@dataclass(frozen=True)
class PositiveGaugeKineticMassEvidence:
    """Generalized mass data for a declared positive gauge kinetic metric."""

    mass_matrix: sp.ImmutableMatrix
    kinetic_metric: sp.ImmutableMatrix
    generalized_mass_operator: sp.ImmutableMatrix
    spectral_parameter: sp.Symbol
    generalized_characteristic_polynomial: sp.Expr
    mass_kernel: sp.ImmutableMatrix
    generalized_kernel_residual: sp.ImmutableMatrix

    @property
    def kernel_certified(self) -> bool:
        """Whether every algebraic mass null remains a generalized null."""

        return _zero_matrix(self.generalized_kernel_residual)


def positive_gauge_kinetic_mass_evidence(
    mass_matrix: Any,
    kinetic_metric: Any,
    *,
    spectral_parameter: sp.Symbol | None = None,
) -> PositiveGaugeKineticMassEvidence:
    """Form ``M2*x=lambda*K*x`` for an exact symmetric positive ``K``.

    The raw eigenvalues of ``M2`` are mass parameters only when the same gauge
    basis has ``K=I``.  This helper derives the generalized operator and
    characteristic polynomial; it does not turn them into physical poles.
    """

    mass = _exact_matrix(mass_matrix, "mass_matrix")
    kinetic = _exact_matrix(kinetic_metric, "kinetic_metric")
    if mass.rows != mass.cols:
        raise ValueError("mass_matrix must be square")
    if kinetic.shape != mass.shape:
        raise ValueError("kinetic_metric must match mass_matrix")
    if not _zero_matrix(mass - mass.T):
        raise ValueError("mass_matrix must be symmetric")
    if not _zero_matrix(kinetic - kinetic.T):
        raise ValueError("kinetic_metric must be symmetric")
    if any(entry.is_real is not True for entry in mass):
        raise ValueError("mass_matrix must be explicitly real")
    if any(entry.is_real is not True for entry in kinetic):
        raise ValueError("kinetic_metric must be explicitly real")
    if kinetic.is_positive_definite is not True:
        raise ValueError("kinetic_metric must be provably positive definite")

    parameter = spectral_parameter or sp.Symbol("lambda", real=True)
    if not isinstance(parameter, sp.Symbol) or parameter.is_real is not True:
        raise ValueError("spectral_parameter must be an explicitly real symbol")
    operator = _immutable_simplified(kinetic.inv() * mass)
    mass_kernel = _kernel_matrix(mass)
    return PositiveGaugeKineticMassEvidence(
        mass_matrix=mass,
        kinetic_metric=kinetic,
        generalized_mass_operator=operator,
        spectral_parameter=parameter,
        generalized_characteristic_polynomial=sp.factor(
            (mass - parameter * kinetic).det()
        ),
        mass_kernel=mass_kernel,
        generalized_kernel_residual=_immutable_simplified(operator * mass_kernel),
    )


@dataclass(frozen=True)
class GaugeQuadraticCongruence:
    """Basis-change ledger for paired gauge kinetic and mass forms."""

    mass_matrix: sp.ImmutableMatrix
    kinetic_metric: sp.ImmutableMatrix
    field_map: sp.ImmutableMatrix
    transformed_mass_matrix: sp.ImmutableMatrix
    transformed_kinetic_metric: sp.ImmutableMatrix
    spectral_parameter: sp.Symbol
    characteristic_covariance_residual: sp.Expr
    original_nullity: int
    transformed_nullity: int

    @property
    def generalized_spectrum_covariant(self) -> bool:
        """Whether the determinant transforms by the expected square factor."""

        return sp.simplify(self.characteristic_covariance_residual) == 0


def transform_gauge_quadratic_forms(
    mass_matrix: Any,
    kinetic_metric: Any,
    field_map: Any,
    *,
    spectral_parameter: sp.Symbol | None = None,
) -> GaugeQuadraticCongruence:
    """Transform ``A=S*A_new`` by congruence in both quadratic forms."""

    mass = _exact_matrix(mass_matrix, "mass_matrix")
    kinetic = _exact_matrix(kinetic_metric, "kinetic_metric")
    mapping = _exact_matrix(field_map, "field_map")
    if mass.rows != mass.cols:
        raise ValueError("mass_matrix must be square")
    if kinetic.shape != mass.shape or mapping.shape != mass.shape:
        raise ValueError("kinetic_metric and field_map must match mass_matrix")
    if not _zero_matrix(mass - mass.T) or not _zero_matrix(kinetic - kinetic.T):
        raise ValueError("quadratic-form matrices must be symmetric")
    if any(entry.is_real is not True for entry in (*mass, *kinetic, *mapping)):
        raise ValueError("quadratic-form and field-map entries must be explicitly real")
    determinant = sp.simplify(mapping.det())
    if determinant.is_zero is not False:
        raise ValueError("field_map must be provably invertible")
    parameter = spectral_parameter or sp.Symbol("lambda", real=True)
    if not isinstance(parameter, sp.Symbol) or parameter.is_real is not True:
        raise ValueError("spectral_parameter must be an explicitly real symbol")

    transformed_mass = _immutable_simplified(mapping.T * mass * mapping)
    transformed_kinetic = _immutable_simplified(mapping.T * kinetic * mapping)
    original_characteristic = sp.factor((mass - parameter * kinetic).det())
    transformed_characteristic = sp.factor(
        (transformed_mass - parameter * transformed_kinetic).det()
    )
    return GaugeQuadraticCongruence(
        mass_matrix=mass,
        kinetic_metric=kinetic,
        field_map=mapping,
        transformed_mass_matrix=transformed_mass,
        transformed_kinetic_metric=transformed_kinetic,
        spectral_parameter=parameter,
        characteristic_covariance_residual=sp.simplify(
            transformed_characteristic - determinant**2 * original_characteristic
        ),
        original_nullity=mass.cols - int(mass.rank()),
        transformed_nullity=transformed_mass.cols - int(transformed_mass.rank()),
    )


@dataclass(frozen=True)
class SU2U1LowerDoubletMassEvidence:
    """Exact canonical quadratic ledger for one declared lower doublet vacuum."""

    coupling_su2: sp.Expr
    coupling_u1: sp.Expr
    vacuum_scale: sp.Expr
    generators: tuple[sp.ImmutableMatrix, ...]
    vacuum: sp.ImmutableMatrix
    charge_operator: sp.ImmutableMatrix
    charge_vacuum_residual: sp.ImmutableMatrix
    general_evidence: GaugeScalarMassEvidence
    charged_mass_squared: sp.Expr
    neutral_mass_matrix: sp.ImmutableMatrix
    neutral_null_vector: sp.ImmutableMatrix
    neutral_massive_vector: sp.ImmutableMatrix
    neutral_mass_squared: sp.Expr
    mixing_sine: sp.Expr
    mixing_cosine: sp.Expr
    rho: sp.Expr


def su2_u1_lower_doublet_mass_evidence(
    coupling_su2: Any,
    coupling_u1: Any,
    vacuum_scale: Any,
) -> SU2U1LowerDoubletMassEvidence:
    """Specialize the Gram theorem to Pauli-half SU(2), ``Y/2=I/2``.

    The supplied couplings and vacuum scale must be exact positive quantities.
    The returned names describe neutral and charged algebraic directions only;
    they do not assign physical particles.
    """

    g = _exact_real(coupling_su2, "coupling_su2")
    gp = _exact_real(coupling_u1, "coupling_u1")
    v = _exact_real(vacuum_scale, "vacuum_scale")
    if g.is_positive is not True or gp.is_positive is not True or v.is_positive is not True:
        raise ValueError("couplings and vacuum_scale must be explicitly positive")

    su2 = su2_fundamental_ledger().generators
    abelian = sp.ImmutableMatrix(sp.eye(2) / 2)
    generators = (*su2, abelian)
    vacuum = sp.ImmutableMatrix([0, v / sp.sqrt(2)])
    general = gauge_scalar_mass_evidence(
        generators,
        (g, g, g, gp),
        vacuum,
    )
    neutral = sp.ImmutableMatrix(general.mass_matrix.extract((2, 3), (2, 3)))
    denominator = sp.sqrt(g**2 + gp**2)
    sine = sp.simplify(gp / denominator)
    cosine = sp.simplify(g / denominator)
    null_vector = sp.ImmutableMatrix([sine, cosine])
    massive_vector = sp.ImmutableMatrix([cosine, -sine])
    charged_mass = sp.simplify(g**2 * v**2 / 4)
    neutral_mass = sp.simplify((g**2 + gp**2) * v**2 / 4)
    charge = _immutable_simplified(su2[2] + abelian)
    return SU2U1LowerDoubletMassEvidence(
        coupling_su2=g,
        coupling_u1=gp,
        vacuum_scale=v,
        generators=generators,
        vacuum=vacuum,
        charge_operator=charge,
        charge_vacuum_residual=_immutable_simplified(charge * vacuum),
        general_evidence=general,
        charged_mass_squared=charged_mass,
        neutral_mass_matrix=neutral,
        neutral_null_vector=null_vector,
        neutral_massive_vector=massive_vector,
        neutral_mass_squared=neutral_mass,
        mixing_sine=sine,
        mixing_cosine=cosine,
        rho=sp.simplify(
            charged_mass / (neutral_mass * cosine**2)
        ),
    )

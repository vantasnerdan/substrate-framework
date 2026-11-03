"""Exact SU(2) current algebra and conditional leading HLS reduction.

This module separates statements that are often conflated.  The
current-wedge and Maurer--Cartan identities are exact.  Minimizing a declared
quadratic hidden-connection mass term is exact at leading derivative order.
Substituting that connection into a separately declared kinetic-curvature term
then gives the order-``p**4`` Skyrme density, while the full kinetic-vector
equation generally changes the field at order ``p**3 / M**2`` and the action at
order ``p**6 / M**2``.  A declared algebraic singlet vector can likewise be
eliminated exactly, but its local squared-current term is only the leading
term when a kinetic differential operator is restored.  Finally, ``U(2)``
invariance permits independent positive singlet and triplet quadratic
coefficients; equality follows from the additional single-trace premise, not
from the group name alone.

No API here derives an HLS sector, a physical rho meson, KSRF, a value of a
coupling, a hedgehog coefficient, an in-medium response, or a substrate map.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .effective_actions import eliminate_quadratic_field


_PAULI = (
    sp.ImmutableMatrix([[0, 1], [1, 0]]),
    sp.ImmutableMatrix([[0, -sp.I], [sp.I, 0]]),
    sp.ImmutableMatrix([[1, 0], [0, -1]]),
)


def _exact_real_matrix(value: Any, name: str) -> sp.ImmutableMatrix:
    matrix = sp.Matrix(value)
    if matrix.rows == 0 or matrix.cols != 3:
        raise ValueError(f"{name} must be a nonempty matrix with three columns")
    if any(entry.has(sp.Float) for entry in matrix):
        raise ValueError(f"{name} must contain exact entries")
    if any(entry.is_real is not True for entry in matrix):
        raise ValueError(f"{name} entries must be declared real")
    return sp.ImmutableMatrix(matrix.applyfunc(sp.simplify))


def _positive_exact(value: Any, name: str) -> sp.Expr:
    result = sp.sympify(value)
    if result.has(sp.Float):
        raise ValueError(f"{name} must be exact")
    if result.is_real is not True or result.is_positive is not True:
        raise ValueError(f"{name} must be declared positive")
    return sp.simplify(result)


def _zero_matrix(rows: int, columns: int) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(sp.zeros(rows, columns))


def _su2_current(row: sp.MatrixBase) -> sp.ImmutableMatrix:
    matrix = sp.zeros(2)
    for component, pauli in zip(row, _PAULI, strict=True):
        matrix += sp.I * component * pauli
    return sp.ImmutableMatrix(matrix.applyfunc(sp.simplify))


@dataclass(frozen=True)
class SU2CurrentQuartic:
    r"""Exact quartic data for ``L_i=i*x_i^a*sigma_a``.

    Both spatial sums are ordered sums.  Thus

    ``sum_ij |x_i cross x_j|^2 = I1-I2`` and
    ``sum_ij Tr([L_i,L_j]^2) = -8*(I1-I2)``.
    """

    components: sp.ImmutableMatrix
    gram_matrix: sp.ImmutableMatrix
    invariant_one: sp.Expr
    invariant_two: sp.Expr
    wedge_norm_squared: sp.Expr
    trace_commutator_sum: sp.Expr


@dataclass(frozen=True)
class HLSCurvaturePair:
    """One exact curvature of the half Maurer--Cartan connection."""

    first_index: int
    second_index: int
    current_commutator: sp.ImmutableMatrix
    connection_curvature: sp.ImmutableMatrix


@dataclass(frozen=True)
class HLSDerivativeOrders:
    """Formal derivative orders in the leading heavy-vector reduction."""

    current: int = 1
    connection: int = 1
    connection_curvature: int = 2
    kinetic_eom_residual: int = 3
    leading_field_correction: int = 3
    leading_quartic_energy: int = 4
    first_backreaction_energy: int = 6


@dataclass(frozen=True)
class HLSLeadingReduction:
    """Leading connection substitution for a declared HLS-type energy.

    The declared algebraic mass density is

    ``kappa*sum_(i,a) (v_i^a-x_i^a/2)^2``.

    Its exact stationary point is ``v=x/2``.  With
    ``L_i=i*x_i^a*sigma_a`` and Maurer--Cartan flatness, the connection
    ``Gamma_i=L_i/2`` has ``F_ij(Gamma)=-[L_i,L_j]/4``.  Substitution into the
    declared curvature density

    ``-sum_ij Tr(F_ij F_ij)/(2*g^2)``

    gives ``-sum_ij Tr([L_i,L_j]^2)/(32*g^2)``.  This matches the equally
    normalized Skyrme density with ``e=g``.  It is an order-``p**4`` result,
    not an exact solution of the full kinetic-vector field equation.
    """

    current_quartic: SU2CurrentQuartic
    currents: tuple[sp.ImmutableMatrix, ...]
    stationary_vector_components: sp.ImmutableMatrix
    mass_hessian: sp.ImmutableMatrix
    mass_stationarity_residual: sp.ImmutableMatrix
    connection: tuple[sp.ImmutableMatrix, ...]
    curvature_pairs: tuple[HLSCurvaturePair, ...]
    ordered_curvature_trace_sum: sp.Expr
    leading_curvature_energy: sp.Expr
    matched_skyrme_energy: sp.Expr
    mass_coefficient: sp.Expr
    gauge_coupling: sp.Expr
    matched_skyrme_coupling: sp.Expr
    derivative_orders: HLSDerivativeOrders


@dataclass(frozen=True)
class HLSKSRFMatching:
    """Conditional dimensionless matching for ``m_V^2=a*g^2*F^2``."""

    vector_mass: sp.Expr
    decay_scale: sp.Expr
    hls_parameter: sp.Expr
    gauge_coupling_squared: sp.Expr
    gauge_coupling: sp.Expr
    skyrme_coupling: sp.Expr
    inverse_skyrme_coupling_squared: sp.Expr
    relation_residual: sp.Expr


@dataclass(frozen=True)
class U2InvariantMetric:
    r"""General positive invariant metric in the Pauli-half ``u(2)`` basis.

    With ``T_0=I/2`` and ``T_a=sigma_a/2``, every real symmetric
    adjoint-invariant form is

    ``2*alpha*Tr(X*Y) + (beta-alpha)*Tr(X)*Tr(Y)``.

    Its component Gram matrix is ``diag(beta,alpha,alpha,alpha)``.  Thus
    ``alpha>0`` and ``beta>0`` give a positive metric, while the single-trace
    specialization and singlet--triplet degeneracy require ``alpha=beta``.
    """

    triplet_coefficient: sp.Expr
    singlet_coefficient: sp.Expr
    generators: tuple[sp.ImmutableMatrix, ...]
    fundamental_trace_gram: sp.ImmutableMatrix
    single_trace_coefficient: sp.Expr
    double_trace_coefficient: sp.Expr
    invariant_gram: sp.ImmutableMatrix
    singlet_triplet_degenerate: bool


@dataclass(frozen=True)
class VectorCurrentSexticMatching:
    r"""Conditional algebraic vector-current elimination and convention map.

    The declared algebraic action is
    ``m**2*w**2/2 + g*w*B``.  Exact stationary elimination gives
    ``-g**2*B**2/(2*m**2)``.  This matches ``-lambda_A**2*B**2`` with
    ``lambda_A=g/(sqrt(2)*m)``.  In C-BPS-001's convention
    ``-lambda_BPS**2*pi**4*B**2``, the same coefficient is
    ``lambda_BPS=lambda_A/pi**2``.

    The result supplies no vector field, current interpretation, mass,
    coupling, kinetic operator, or low-momentum power counting.
    """

    vector_mass: sp.Expr
    current_coupling: sp.Expr
    stationary_field_per_current: sp.Expr
    effective_current_coefficient: sp.Expr
    source_sextic_coupling: sp.Expr
    bps_sextic_coupling: sp.Expr
    convention_ratio: sp.Expr
    stationarity_residual: sp.ImmutableMatrix


def su2_current_quartic(current_components: Any) -> SU2CurrentQuartic:
    r"""Derive the exact SU(2) wedge, Gram, and commutator-square identity.

    Rows are spatial or spacetime current labels and the three columns are
    exact real Pauli components in ``L_i=i*x_i^a*sigma_a``.  The function
    attaches no field equation or physical interpretation to those rows.
    """

    components = _exact_real_matrix(current_components, "current_components")
    component_matrix = sp.Matrix(components)
    gram = component_matrix * component_matrix.T
    invariant_one = sp.expand(sp.trace(gram) ** 2)
    invariant_two = sp.expand(sp.trace(gram * gram))

    wedge_norm = sp.Integer(0)
    currents = tuple(
        _su2_current(component_matrix.row(index))
        for index in range(component_matrix.rows)
    )
    trace_sum = sp.Integer(0)
    for first in range(component_matrix.rows):
        first_vector = component_matrix.row(first).T
        for second in range(component_matrix.rows):
            second_vector = component_matrix.row(second).T
            cross = first_vector.cross(second_vector)
            wedge_norm += sp.expand(cross.dot(cross))
            commutator = sp.Matrix(currents[first] * currents[second]) - sp.Matrix(
                currents[second] * currents[first]
            )
            trace_sum += sp.trace(commutator * commutator)

    wedge_norm = sp.simplify(sp.expand(wedge_norm))
    trace_sum = sp.simplify(sp.expand(trace_sum))
    if sp.simplify(wedge_norm - (invariant_one - invariant_two)) != 0:
        raise AssertionError("current wedge and Gram invariants disagree")
    if sp.simplify(trace_sum + 8 * wedge_norm) != 0:
        raise AssertionError("Pauli commutator normalization disagrees with wedge norm")
    return SU2CurrentQuartic(
        components=components,
        gram_matrix=sp.ImmutableMatrix(gram.applyfunc(sp.simplify)),
        invariant_one=sp.simplify(invariant_one),
        invariant_two=sp.simplify(invariant_two),
        wedge_norm_squared=wedge_norm,
        trace_commutator_sum=trace_sum,
    )


def leading_hls_connection_reduction(
    current_components: Any,
    gauge_coupling: Any,
    *,
    mass_coefficient: Any = 1,
) -> HLSLeadingReduction:
    r"""Return the exact order-``p**4`` half-connection reduction.

    The vector kinetic term is retained as the source of the quartic density,
    but its contribution to the vector equation changes the connection only
    beyond leading order.  The returned order ledger makes that approximation
    boundary part of the API rather than silently calling it exact.
    """

    quartic = su2_current_quartic(current_components)
    coupling = _positive_exact(gauge_coupling, "gauge_coupling")
    coefficient = _positive_exact(mass_coefficient, "mass_coefficient")
    components = sp.Matrix(quartic.components)
    stationary_components = sp.ImmutableMatrix(
        (components / 2).applyfunc(sp.simplify)
    )
    mass_hessian = sp.ImmutableMatrix(
        2 * coefficient * sp.eye(components.rows * components.cols)
    )
    mass_residual = _zero_matrix(components.rows, components.cols)

    currents = tuple(
        _su2_current(components.row(index)) for index in range(components.rows)
    )
    connection = tuple(
        sp.ImmutableMatrix((sp.Matrix(current) / 2).applyfunc(sp.simplify))
        for current in currents
    )
    curvature_pairs: list[HLSCurvaturePair] = []
    ordered_curvature_trace_sum = sp.Integer(0)
    for first in range(components.rows):
        for second in range(first + 1, components.rows):
            commutator = sp.Matrix(currents[first] * currents[second]) - sp.Matrix(
                currents[second] * currents[first]
            )
            commutator = sp.ImmutableMatrix(commutator.applyfunc(sp.simplify))
            curvature = sp.ImmutableMatrix(
                (-sp.Matrix(commutator) / 4).applyfunc(sp.simplify)
            )
            curvature_pairs.append(
                HLSCurvaturePair(
                    first_index=first,
                    second_index=second,
                    current_commutator=commutator,
                    connection_curvature=curvature,
                )
            )
            ordered_curvature_trace_sum += 2 * sp.trace(
                sp.Matrix(curvature) * sp.Matrix(curvature)
            )

    ordered_curvature_trace_sum = sp.simplify(
        sp.expand(ordered_curvature_trace_sum)
    )
    leading_energy = sp.simplify(
        -ordered_curvature_trace_sum / (2 * coupling**2)
    )
    skyrme_energy = sp.simplify(
        -quartic.trace_commutator_sum / (32 * coupling**2)
    )
    if sp.simplify(leading_energy - skyrme_energy) != 0:
        raise AssertionError("half-connection curvature does not match Skyrme density")
    if sp.simplify(leading_energy - quartic.wedge_norm_squared / (4 * coupling**2)) != 0:
        raise AssertionError("leading HLS density has the wrong wedge coefficient")

    return HLSLeadingReduction(
        current_quartic=quartic,
        currents=currents,
        stationary_vector_components=stationary_components,
        mass_hessian=mass_hessian,
        mass_stationarity_residual=mass_residual,
        connection=connection,
        curvature_pairs=tuple(curvature_pairs),
        ordered_curvature_trace_sum=ordered_curvature_trace_sum,
        leading_curvature_energy=leading_energy,
        matched_skyrme_energy=skyrme_energy,
        mass_coefficient=coefficient,
        gauge_coupling=coupling,
        matched_skyrme_coupling=coupling,
        derivative_orders=HLSDerivativeOrders(),
    )


def conditional_hls_ksrf_matching(
    vector_mass: Any,
    decay_scale: Any,
    *,
    hls_parameter: Any = 2,
) -> HLSKSRFMatching:
    r"""Solve the separately declared relation ``m_V^2=a*g^2*F^2``.

    When ``m_V`` and ``F`` carry the same mass dimension and ``a`` is
    dimensionless, the result ``g=e=m_V/(sqrt(a)*F)`` is dimensionless.  The
    relation is a premise; this function neither derives KSRF nor selects
    ``a=2``.
    """

    mass = _positive_exact(vector_mass, "vector_mass")
    scale = _positive_exact(decay_scale, "decay_scale")
    parameter = _positive_exact(hls_parameter, "hls_parameter")
    coupling_squared = sp.simplify(mass**2 / (parameter * scale**2))
    coupling = sp.simplify(mass / (sp.sqrt(parameter) * scale))
    inverse_squared = sp.simplify(parameter * scale**2 / mass**2)
    residual = sp.simplify(
        mass**2 - parameter * coupling_squared * scale**2
    )
    return HLSKSRFMatching(
        vector_mass=mass,
        decay_scale=scale,
        hls_parameter=parameter,
        gauge_coupling_squared=coupling_squared,
        gauge_coupling=coupling,
        skyrme_coupling=coupling,
        inverse_skyrme_coupling_squared=inverse_squared,
        relation_residual=residual,
    )


def u2_invariant_metric(
    triplet_coefficient: Any,
    singlet_coefficient: Any,
) -> U2InvariantMetric:
    r"""Return the general positive adjoint-invariant ``u(2)`` metric.

    The two inputs are the component coefficients of the three traceless
    generators and the central generator.  Keeping them separate prevents a
    selected fundamental single trace from masquerading as a consequence of
    ``U(2)`` invariance.
    """

    alpha = _positive_exact(triplet_coefficient, "triplet_coefficient")
    beta = _positive_exact(singlet_coefficient, "singlet_coefficient")
    identity = sp.ImmutableMatrix(sp.eye(2))
    generators = (identity / 2,) + tuple(pauli / 2 for pauli in _PAULI)
    trace_gram = sp.ImmutableMatrix(
        4,
        4,
        lambda first, second: sp.simplify(
            sp.trace(generators[first] * generators[second])
        ),
    )
    if trace_gram != sp.ImmutableMatrix(sp.eye(4) / 2):
        raise AssertionError("Pauli-half u(2) trace normalization changed")

    single_trace = sp.simplify(2 * alpha)
    double_trace = sp.simplify(beta - alpha)
    invariant_gram = sp.ImmutableMatrix(sp.diag(beta, alpha, alpha, alpha))
    reconstructed = sp.ImmutableMatrix(
        4,
        4,
        lambda first, second: sp.simplify(
            single_trace * trace_gram[first, second]
            + double_trace
            * sp.trace(generators[first])
            * sp.trace(generators[second])
        ),
    )
    if reconstructed != invariant_gram:
        raise AssertionError("u(2) invariant-form reconstruction failed")
    return U2InvariantMetric(
        triplet_coefficient=alpha,
        singlet_coefficient=beta,
        generators=generators,
        fundamental_trace_gram=trace_gram,
        single_trace_coefficient=single_trace,
        double_trace_coefficient=double_trace,
        invariant_gram=invariant_gram,
        singlet_triplet_degenerate=sp.simplify(beta - alpha) == 0,
    )


def conditional_vector_current_sextic_matching(
    vector_mass: Any,
    current_coupling: Any,
) -> VectorCurrentSexticMatching:
    r"""Eliminate one declared algebraic vector component exactly.

    Restoring a vector kinetic operator replaces the scalar kernel ``m**2``
    by a differential kernel.  Callers must then use the inverse-expansion
    residuals in :mod:`substrate_framework.effective_actions`; this function
    does not label the local term an exact full-field result.
    """

    mass = _positive_exact(vector_mass, "vector_mass")
    coupling = _positive_exact(current_coupling, "current_coupling")
    current = sp.Symbol("B_current", real=True)
    elimination = eliminate_quadratic_field(
        sp.Matrix([[mass**2]]),
        sp.Matrix([coupling * current]),
    )
    source_coupling = sp.simplify(coupling / (sp.sqrt(2) * mass))
    bps_coupling = sp.simplify(source_coupling / sp.pi**2)
    expected_effective = sp.simplify(-source_coupling**2 * current**2)
    if sp.simplify(elimination.effective_term - expected_effective) != 0:
        raise AssertionError("vector-current elimination coefficient changed")
    return VectorCurrentSexticMatching(
        vector_mass=mass,
        current_coupling=coupling,
        stationary_field_per_current=sp.simplify(
            elimination.stationary_field[0] / current
        ),
        effective_current_coefficient=sp.simplify(
            elimination.effective_term / current**2
        ),
        source_sextic_coupling=source_coupling,
        bps_sextic_coupling=bps_coupling,
        convention_ratio=sp.simplify(source_coupling / bps_coupling),
        stationarity_residual=elimination.stationarity_residual,
    )

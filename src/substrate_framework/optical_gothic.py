"""Exact optical ADM and gothic-metric ledgers in mostly-plus signature.

The general chart uses a positive lapse ``N``, a symmetric positive-definite
spatial optical tensor ``gamma``, and dimensionless flow ``v=V/c``:

``ds^2 = -N^2 (dx^0)^2 + gamma_ij (dx^i+s*v^i*dx^0)
                                      (dx^j+s*v^j*dx^0)``.

Here ``s`` is ``-1`` for a material-flow convention and ``+1`` for the shift
sign printed in the issue-184 source.  The independent lapse is essential: a
map that fixes ``N=(det gamma)^(-1/6)`` has nine fields and covers only a
codimension-one metric slice, not all foliation-adapted Lorentzian metrics.

This module proves algebraic maps and compatibility ledgers.  It does not
derive an optical medium, a constitutive action, Newton's constant, or a
physical identification of an effective metric with gravity.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import sympy as sp

from .exact_symbolic import exact_real, positive_exact


def _exact_matrix(values: Any, rows: int, columns: int, name: str) -> sp.ImmutableMatrix:
    matrix = sp.ImmutableMatrix(values)
    if matrix.shape != (rows, columns):
        raise ValueError(f"{name} must have shape ({rows}, {columns})")
    entries = [
        exact_real(matrix[row, column], f"{name}[{row},{column}]")
        for row in range(rows)
        for column in range(columns)
    ]
    return sp.ImmutableMatrix(rows, columns, entries)


def _positive_definite(values: Any, size: int, name: str) -> sp.ImmutableMatrix:
    matrix = _exact_matrix(values, size, size, name)
    if matrix != matrix.T:
        raise ValueError(f"{name} must be symmetric")
    for order in range(1, size + 1):
        minor = sp.simplify(matrix[:order, :order].det())
        if minor.is_positive is not True:
            raise ValueError(f"{name} must be explicitly positive definite")
    return matrix


def _flow_sign(value: Any) -> sp.Integer:
    sign = sp.sympify(value)
    if sign not in (sp.Integer(-1), sp.Integer(1)):
        raise ValueError("flow_orientation must be exactly -1 or +1")
    return sign


@dataclass(frozen=True)
class OpticalADMMetric:
    """Complete exact optical ADM metric and inverse-density data."""

    lapse: sp.Expr
    spatial_tensor: sp.ImmutableMatrix
    flow_velocity: sp.ImmutableMatrix
    signal_speed: sp.Expr
    flow_orientation: sp.Integer
    dimensionless_flow: sp.ImmutableMatrix
    covariant: sp.ImmutableMatrix
    contravariant: sp.ImmutableMatrix
    determinant: sp.Expr
    volume_density: sp.Expr
    gothic_contravariant: sp.ImmutableMatrix
    component_jacobian_determinant: sp.Expr


@dataclass(frozen=True)
class RecoveredOpticalADM:
    """Optical variables reconstructed from a foliation-adapted metric."""

    lapse: sp.Expr
    spatial_tensor: sp.ImmutableMatrix
    flow_velocity: sp.ImmutableMatrix
    signal_speed: sp.Expr
    flow_orientation: sp.Integer
    reconstructed_metric: OpticalADMMetric


@dataclass(frozen=True)
class DeterminantSlavedOpticalMetric:
    """Nine-field optical map and the exact constraint defining its image."""

    metric: OpticalADMMetric
    determinant_mean: sp.Expr
    lapse_constraint_residual: sp.Expr
    volume_constraint_residual: sp.Expr
    field_count: int
    general_metric_component_count: int


@dataclass(frozen=True)
class OpticalContinuityCompatibility:
    """Material and harmonic determinant-density continuity residuals."""

    determinant_mean: sp.Expr
    time_derivative: sp.Expr
    advective_derivative: sp.Expr
    flow_divergence: sp.Expr
    flow_orientation: sp.Integer
    material_residual: sp.Expr
    harmonic_residual: sp.Expr
    harmonic_from_material_residual: sp.Expr


@dataclass(frozen=True)
class NewtonianOpticalGradientLedger:
    """Weak conformal-gradient energy and stress coefficient comparison."""

    signal_speed: sp.Expr
    newton_constant: sp.Expr
    log_index_gradient_squared: sp.Expr
    weak_potential_gradient_squared: sp.Expr
    paper_energy_density: sp.Expr
    required_newtonian_energy_density: sp.Expr
    paper_energy_minus_required: sp.Expr
    maxwell_stress_prefactor_in_potential_variables: sp.Expr


@dataclass(frozen=True)
class OpticalShearDecomposition:
    """Spectral conformal/shear decomposition for three positive indices."""

    principal_indices: tuple[sp.Expr, sp.Expr, sp.Expr]
    determinant_mean: sp.Expr
    additive_deviations: tuple[sp.Expr, sp.Expr, sp.Expr]
    additive_trace: sp.Expr
    logarithmic_mean: sp.Expr
    logarithmic_shear: tuple[sp.Expr, sp.Expr, sp.Expr]
    logarithmic_shear_trace: sp.Expr
    normalized_principal_indices: tuple[sp.Expr, sp.Expr, sp.Expr]
    normalized_determinant: sp.Expr


def optical_adm_metric(
    lapse: Any,
    spatial_tensor: Any,
    flow_velocity: Any,
    signal_speed: Any,
    *,
    flow_orientation: Any = -1,
) -> OpticalADMMetric:
    """Return the exact ten-field optical ADM metric and gothic inverse.

    Coordinates use ``x^0=c*t``.  The six spatial components, three flow
    components, and independent lapse give a nonsingular ten-to-ten point map
    with component Jacobian determinant ``-2*s*N*det(gamma)``.
    """

    lapse_value = positive_exact(lapse, "lapse")
    speed = positive_exact(signal_speed, "signal_speed")
    spatial = _positive_definite(spatial_tensor, 3, "spatial_tensor")
    flow = _exact_matrix(flow_velocity, 3, 1, "flow_velocity")
    orientation = _flow_sign(flow_orientation)
    velocity = flow.applyfunc(lambda entry: sp.simplify(entry / speed))

    covariant = sp.zeros(4)
    flow_quadratic = sp.simplify((velocity.T * spatial * velocity)[0])
    covariant[0, 0] = sp.simplify(-lapse_value**2 + flow_quadratic)
    mixed = orientation * spatial * velocity
    for index in range(3):
        covariant[0, index + 1] = mixed[index]
        covariant[index + 1, 0] = mixed[index]
        for other in range(3):
            covariant[index + 1, other + 1] = spatial[index, other]

    spatial_inverse = spatial.inv()
    contravariant = sp.zeros(4)
    contravariant[0, 0] = -1 / lapse_value**2
    for index in range(3):
        contravariant[0, index + 1] = (
            orientation * velocity[index] / lapse_value**2
        )
        contravariant[index + 1, 0] = contravariant[0, index + 1]
        for other in range(3):
            contravariant[index + 1, other + 1] = (
                spatial_inverse[index, other]
                - velocity[index] * velocity[other] / lapse_value**2
            )

    covariant_immutable = sp.ImmutableMatrix(covariant.applyfunc(sp.simplify))
    contravariant_immutable = sp.ImmutableMatrix(
        contravariant.applyfunc(sp.simplify)
    )
    if (covariant_immutable * contravariant_immutable).applyfunc(sp.simplify) != sp.eye(4):
        raise AssertionError("optical ADM block formulas are not mutual inverses")

    spatial_determinant = sp.simplify(spatial.det())
    determinant = sp.simplify(-lapse_value**2 * spatial_determinant)
    if sp.simplify(covariant_immutable.det() - determinant) != 0:
        raise AssertionError("optical ADM determinant does not match the block formula")
    volume_density = sp.simplify(lapse_value * sp.sqrt(spatial_determinant))
    gothic = sp.ImmutableMatrix(
        (volume_density * contravariant_immutable).applyfunc(sp.simplify)
    )
    jacobian = sp.simplify(
        -2 * orientation * lapse_value * spatial_determinant
    )
    return OpticalADMMetric(
        lapse=lapse_value,
        spatial_tensor=spatial,
        flow_velocity=flow,
        signal_speed=speed,
        flow_orientation=orientation,
        dimensionless_flow=sp.ImmutableMatrix(velocity),
        covariant=covariant_immutable,
        contravariant=contravariant_immutable,
        determinant=determinant,
        volume_density=volume_density,
        gothic_contravariant=gothic,
        component_jacobian_determinant=jacobian,
    )


def recover_optical_adm(
    covariant_metric: Any,
    signal_speed: Any,
    *,
    flow_orientation: Any = -1,
) -> RecoveredOpticalADM:
    """Invert :func:`optical_adm_metric` on its foliation-adapted domain."""

    metric = _exact_matrix(covariant_metric, 4, 4, "covariant_metric")
    if metric != metric.T:
        raise ValueError("covariant_metric must be symmetric")
    speed = positive_exact(signal_speed, "signal_speed")
    orientation = _flow_sign(flow_orientation)
    spatial = _positive_definite(metric[1:4, 1:4], 3, "spatial_metric")
    mixed = sp.ImmutableMatrix(metric[1:4, 0])
    dimensionless_flow = sp.ImmutableMatrix(
        (orientation * spatial.inv() * mixed).applyfunc(sp.simplify)
    )
    lapse_squared = sp.simplify(
        -(metric[0, 0] - (mixed.T * spatial.inv() * mixed)[0])
    )
    lapse = positive_exact(sp.sqrt(lapse_squared), "recovered_lapse")
    flow = sp.ImmutableMatrix(
        (speed * dimensionless_flow).applyfunc(sp.simplify)
    )
    reconstructed = optical_adm_metric(
        lapse,
        spatial,
        flow,
        speed,
        flow_orientation=orientation,
    )
    if reconstructed.covariant != metric:
        raise AssertionError("recovered optical variables do not reconstruct the metric")
    return RecoveredOpticalADM(
        lapse=lapse,
        spatial_tensor=spatial,
        flow_velocity=flow,
        signal_speed=speed,
        flow_orientation=orientation,
        reconstructed_metric=reconstructed,
    )


def determinant_slaved_optical_metric(
    spatial_tensor: Any,
    flow_velocity: Any,
    signal_speed: Any,
    *,
    flow_orientation: Any = 1,
) -> DeterminantSlavedOpticalMetric:
    """Return the paper's valid nine-field map and its image constraint."""

    spatial = _positive_definite(spatial_tensor, 3, "spatial_tensor")
    determinant = sp.simplify(spatial.det())
    determinant_mean = sp.real_root(determinant, 3)
    lapse = sp.simplify(1 / sp.sqrt(determinant_mean))
    metric = optical_adm_metric(
        lapse,
        spatial,
        flow_velocity,
        signal_speed,
        flow_orientation=flow_orientation,
    )
    lapse_residual = sp.simplify(metric.lapse**2 * determinant_mean - 1)
    volume_residual = sp.simplify(metric.volume_density - determinant_mean)
    return DeterminantSlavedOpticalMetric(
        metric=metric,
        determinant_mean=determinant_mean,
        lapse_constraint_residual=lapse_residual,
        volume_constraint_residual=volume_residual,
        field_count=9,
        general_metric_component_count=10,
    )


def optical_continuity_compatibility(
    determinant_mean: Any,
    time_derivative: Any,
    advective_derivative: Any,
    flow_divergence: Any,
    *,
    flow_orientation: Any = -1,
) -> OpticalContinuityCompatibility:
    """Compare material continuity with the gothic time-gauge equation.

    ``time_derivative`` denotes ``partial_t(nbar)`` and
    ``advective_derivative`` denotes ``V dot grad(nbar)``.  The material
    residual is ``M=partial_t nbar+div(nbar V)``.  From the exact gothic
    block, the harmonic residual is
    ``H_s=partial_t(nbar^2)-s*div(nbar^2 V)``.
    """

    mean = positive_exact(determinant_mean, "determinant_mean")
    partial_time = exact_real(time_derivative, "time_derivative")
    advective = exact_real(advective_derivative, "advective_derivative")
    divergence = exact_real(flow_divergence, "flow_divergence")
    orientation = _flow_sign(flow_orientation)
    material = sp.simplify(partial_time + advective + mean * divergence)
    harmonic = sp.simplify(
        2 * mean * partial_time
        - orientation * (2 * mean * advective + mean**2 * divergence)
    )
    harmonic_from_material = sp.simplify(
        2 * mean * material
        - 2 * mean * (1 + orientation) * advective
        - mean**2 * (2 + orientation) * divergence
    )
    if sp.simplify(harmonic - harmonic_from_material) != 0:
        raise AssertionError("continuity compatibility expansion did not close")
    return OpticalContinuityCompatibility(
        determinant_mean=mean,
        time_derivative=partial_time,
        advective_derivative=advective,
        flow_divergence=divergence,
        flow_orientation=orientation,
        material_residual=material,
        harmonic_residual=harmonic,
        harmonic_from_material_residual=harmonic_from_material,
    )


def newtonian_optical_gradient_ledger(
    signal_speed: Any,
    newton_constant: Any,
    log_index_gradient_squared: Any,
) -> NewtonianOpticalGradientLedger:
    """Expose the weak-limit sign conflict in the paper's equation (55).

    For ``Phi=(c^2/2)*(1-nbar^-1)``, weak fields give
    ``|grad Phi|^2=(c^4/4)|grad log(nbar)|^2``.  The printed positive optical
    gradient energy is therefore ``+|grad Phi|^2/(8*pi*G)``, whereas the
    paper's own equation (18) requires the negative Newtonian field energy.
    """

    speed = positive_exact(signal_speed, "signal_speed")
    newton = positive_exact(newton_constant, "newton_constant")
    gradient_squared = exact_real(
        log_index_gradient_squared,
        "log_index_gradient_squared",
    )
    if gradient_squared.is_nonnegative is not True:
        raise ValueError("log_index_gradient_squared must be explicitly nonnegative")
    potential_gradient_squared = sp.simplify(speed**4 * gradient_squared / 4)
    paper_energy = sp.simplify(
        speed**4 * gradient_squared / (32 * sp.pi * newton)
    )
    required_energy = sp.simplify(
        -potential_gradient_squared / (8 * sp.pi * newton)
    )
    stress_prefactor = sp.simplify(
        speed**4 / (16 * sp.pi * newton) * 4 / speed**4
    )
    return NewtonianOpticalGradientLedger(
        signal_speed=speed,
        newton_constant=newton,
        log_index_gradient_squared=gradient_squared,
        weak_potential_gradient_squared=potential_gradient_squared,
        paper_energy_density=paper_energy,
        required_newtonian_energy_density=required_energy,
        paper_energy_minus_required=sp.simplify(paper_energy - required_energy),
        maxwell_stress_prefactor_in_potential_variables=stress_prefactor,
    )


def optical_shear_decomposition(
    principal_indices: Iterable[Any],
) -> OpticalShearDecomposition:
    """Return a conformal/log-shear split with exactly zero shear trace.

    The additive proposal ``n_i-nbar`` is retained in the ledger so its trace
    can be tested.  The repaired coordinates subtract the arithmetic mean of
    the three logarithms.  Their sum is identically zero and their exponentials
    have unit product, so the conformal mode is not counted twice.
    """

    indices = tuple(
        positive_exact(value, f"principal_indices[{position}]")
        for position, value in enumerate(principal_indices)
    )
    if len(indices) != 3:
        raise ValueError("principal_indices must contain exactly three entries")
    determinant_mean = sp.real_root(sp.prod(indices), 3)
    additive = tuple(sp.simplify(value - determinant_mean) for value in indices)
    additive_trace = sp.simplify(sum(additive))
    logarithmic_mean = sp.simplify(sum(sp.log(value) for value in indices) / 3)
    logarithmic_shear = tuple(
        sp.simplify(sp.log(value) - logarithmic_mean) for value in indices
    )
    logarithmic_trace = sp.simplify(sum(logarithmic_shear))
    normalized = tuple(sp.exp(value) for value in logarithmic_shear)
    normalized_determinant = sp.simplify(sp.prod(normalized))
    if logarithmic_trace != 0 or normalized_determinant != 1:
        raise AssertionError("logarithmic shear did not remove the conformal mode")
    return OpticalShearDecomposition(
        principal_indices=indices,
        determinant_mean=determinant_mean,
        additive_deviations=additive,
        additive_trace=additive_trace,
        logarithmic_mean=logarithmic_mean,
        logarithmic_shear=logarithmic_shear,
        logarithmic_shear_trace=logarithmic_trace,
        normalized_principal_indices=normalized,
        normalized_determinant=normalized_determinant,
    )

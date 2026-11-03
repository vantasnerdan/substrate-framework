"""Exact conditional relations for one-loop running couplings."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Any, Literal, Sequence

import sympy as sp

from substrate_framework.linear_systems import (
    LinearSystemDiagnostics,
    diagnose_linear_system,
)


@dataclass(frozen=True)
class AffineCrossingLedger:
    """Exact crossing status for two supplied affine inverse couplings."""

    left: str
    right: str
    status: Literal["unique", "coincident", "parallel_disjoint"]
    coordinate: sp.Expr | None
    intercept_difference: sp.Expr
    coefficient_difference: sp.Expr


@dataclass(frozen=True)
class AffineUnificationDiagnostics:
    """Exact consistency ledger for ``a_i = A + B*b_i``."""

    inverse_couplings: tuple[sp.Expr, ...]
    beta_coefficients: tuple[sp.Expr, ...]
    provenance: tuple[str, ...]
    linear: LinearSystemDiagnostics
    left_nullspace: tuple[sp.ImmutableMatrix, ...]
    compatibility_residuals: tuple[sp.Expr, ...]
    pairwise_crossings: tuple[AffineCrossingLedger, ...]
    common_inverse_coupling: sp.Expr | None
    running_coordinate: sp.Expr | None


@dataclass(frozen=True)
class ElectroweakUnificationReconstruction:
    """Conditional exact inverse reconstruction from two supplied observations."""

    electromagnetic_inverse: sp.Expr
    strong_inverse: sp.Expr
    hypercharge_weight: sp.Expr
    beta_coefficients: tuple[sp.Expr, sp.Expr, sp.Expr]
    denominator: sp.Expr
    common_inverse_coupling: sp.Expr
    running_coordinate: sp.Expr
    inverse_couplings: tuple[sp.Expr, sp.Expr, sp.Expr]
    weak_angle_coordinate: sp.Expr
    boundary_weak_angle_coordinate: sp.Expr


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be provably real")
    return expression


def _exact_positive(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return expression


def _labels(values: Sequence[str], *, size: int) -> tuple[str, ...]:
    labels = tuple(values)
    if len(labels) != size:
        raise ValueError("provenance must name every inverse coupling")
    if any(not isinstance(label, str) or not label.strip() for label in labels):
        raise ValueError("provenance labels must be non-empty strings")
    return labels


def pairwise_affine_crossing(
    inverse_left: Any,
    beta_left: Any,
    inverse_right: Any,
    beta_right: Any,
    *,
    left: str = "left",
    right: str = "right",
) -> AffineCrossingLedger:
    """Classify two exact lines in the convention ``a_i=A+B*b_i``.

    Distinct coefficients give the unique coordinate
    ``B=(a_left-a_right)/(b_left-b_right)``.  Equal coefficients give
    coincident lines only when the inverse couplings also agree; otherwise the
    lines are parallel and disjoint.  The inputs and their physical meanings
    are supplied premises.
    """

    if not isinstance(left, str) or not left.strip():
        raise ValueError("left must be a non-empty label")
    if not isinstance(right, str) or not right.strip():
        raise ValueError("right must be a non-empty label")
    a_left = _exact_real(inverse_left, "inverse_left")
    a_right = _exact_real(inverse_right, "inverse_right")
    b_left = _exact_real(beta_left, "beta_left")
    b_right = _exact_real(beta_right, "beta_right")
    intercept_difference = sp.simplify(a_left - a_right)
    coefficient_difference = sp.simplify(b_left - b_right)
    if coefficient_difference == 0:
        status = "coincident" if intercept_difference == 0 else "parallel_disjoint"
        coordinate = None
    else:
        status = "unique"
        coordinate = sp.simplify(intercept_difference / coefficient_difference)
    return AffineCrossingLedger(
        left=left,
        right=right,
        status=status,
        coordinate=coordinate,
        intercept_difference=intercept_difference,
        coefficient_difference=coefficient_difference,
    )


def diagnose_affine_unification(
    inverse_couplings: Sequence[Any],
    beta_coefficients: Sequence[Any],
    *,
    provenance: Sequence[str],
) -> AffineUnificationDiagnostics:
    """Diagnose a supplied finite family ``a_i=A+B*b_i`` exactly.

    The routine composes exact rank and augmented-rank diagnostics with every
    pairwise crossing.  A returned unique solution is a conditional common
    intersection of the supplied lines; it does not derive the lines, their
    beta coefficients, a physical matching condition, or a perturbative
    domain.
    """

    inverse_values = tuple(
        _exact_real(value, f"inverse_couplings[{index}]")
        for index, value in enumerate(inverse_couplings)
    )
    beta_values = tuple(
        _exact_real(value, f"beta_coefficients[{index}]")
        for index, value in enumerate(beta_coefficients)
    )
    if len(inverse_values) < 2:
        raise ValueError("at least two inverse couplings are required")
    if len(beta_values) != len(inverse_values):
        raise ValueError("beta_coefficients must match inverse_couplings")
    labels = _labels(provenance, size=len(inverse_values))
    design = sp.Matrix([[1, coefficient] for coefficient in beta_values])
    rhs = sp.Matrix(inverse_values)
    linear = diagnose_linear_system(design, rhs)
    left_nullspace = tuple(
        sp.ImmutableMatrix(vector) for vector in design.T.nullspace()
    )
    compatibility_residuals = tuple(
        sp.simplify((vector.T * rhs)[0]) for vector in left_nullspace
    )
    pairwise = tuple(
        pairwise_affine_crossing(
            inverse_values[left_index],
            beta_values[left_index],
            inverse_values[right_index],
            beta_values[right_index],
            left=labels[left_index],
            right=labels[right_index],
        )
        for left_index, right_index in combinations(range(len(inverse_values)), 2)
    )
    common_inverse: sp.Expr | None = None
    running_coordinate: sp.Expr | None = None
    if linear.unique:
        solution = sp.linsolve((design, rhs))
        common_inverse, running_coordinate = tuple(next(iter(solution)))
        common_inverse = sp.simplify(common_inverse)
        running_coordinate = sp.simplify(running_coordinate)
    return AffineUnificationDiagnostics(
        inverse_couplings=inverse_values,
        beta_coefficients=beta_values,
        provenance=labels,
        linear=linear,
        left_nullspace=left_nullspace,
        compatibility_residuals=compatibility_residuals,
        pairwise_crossings=pairwise,
        common_inverse_coupling=common_inverse,
        running_coordinate=running_coordinate,
    )


def reconstruct_electroweak_unification(
    electromagnetic_inverse: Any,
    strong_inverse: Any,
    beta_1: Any,
    beta_2: Any,
    beta_3: Any,
    hypercharge_weight: Any,
) -> ElectroweakUnificationReconstruction:
    r"""Solve a declared two-observation exact-unification inverse problem.

    The supplied convention is

    ``a_i=A+B*b_i``, ``a_3=strong_inverse``, and
    ``electromagnetic_inverse=a_2+n*a_1``.

    Here ``n`` is the supplied hypercharge weight (``5/3`` in the legacy WM3
    coordinate).  The weak-angle coordinate is then ``a_2/E`` and its common-
    coupling boundary value is ``1/(1+n)``.  The result is an inverse
    reconstruction from two observations and an exact matching premise, not
    an ab-initio prediction or a derivation of ``n`` or the beta coefficients.
    """

    electromagnetic = _exact_positive(
        electromagnetic_inverse, "electromagnetic_inverse"
    )
    strong = _exact_positive(strong_inverse, "strong_inverse")
    weight = _exact_positive(hypercharge_weight, "hypercharge_weight")
    beta = (
        _exact_real(beta_1, "beta_1"),
        _exact_real(beta_2, "beta_2"),
        _exact_real(beta_3, "beta_3"),
    )
    denominator = sp.simplify(
        beta[1] + weight * beta[0] - (1 + weight) * beta[2]
    )
    if denominator == 0:
        raise ValueError("the reconstruction denominator must be nonzero")
    running = sp.simplify(
        (electromagnetic - (1 + weight) * strong) / denominator
    )
    common = sp.simplify(strong - running * beta[2])
    inverse_values = tuple(
        sp.simplify(common + running * coefficient) for coefficient in beta
    )
    weak_angle = sp.simplify(inverse_values[1] / electromagnetic)
    boundary = sp.simplify(1 / (1 + weight))
    return ElectroweakUnificationReconstruction(
        electromagnetic_inverse=electromagnetic,
        strong_inverse=strong,
        hypercharge_weight=weight,
        beta_coefficients=beta,
        denominator=denominator,
        common_inverse_coupling=common,
        running_coordinate=running,
        inverse_couplings=inverse_values,
        weak_angle_coordinate=weak_angle,
        boundary_weak_angle_coordinate=boundary,
    )


def shift_affine_reference(
    inverse_couplings: Sequence[Any],
    beta_coefficients: Sequence[Any],
    coordinate_shift: Any,
) -> tuple[sp.Expr, ...]:
    """Move the reference by ``delta`` in the supplied ``B`` coordinate.

    If ``a_i=A+B*b_i``, the shifted intercepts are ``a_i-delta*b_i`` and the
    same common point has coordinate ``B-delta``.  This algebraic covariance
    supplies no physical reference scale.
    """

    inverse_values = tuple(
        _exact_real(value, f"inverse_couplings[{index}]")
        for index, value in enumerate(inverse_couplings)
    )
    beta_values = tuple(
        _exact_real(value, f"beta_coefficients[{index}]")
        for index, value in enumerate(beta_coefficients)
    )
    if len(beta_values) != len(inverse_values):
        raise ValueError("beta_coefficients must match inverse_couplings")
    shift = _exact_real(coordinate_shift, "coordinate_shift")
    return tuple(
        sp.simplify(inverse - shift * coefficient)
        for inverse, coefficient in zip(inverse_values, beta_values, strict=True)
    )


def rescale_abelian_inverse_coordinate(
    inverse_coupling: Any,
    beta_coefficient: Any,
    electromagnetic_weight: Any,
    coordinate_factor: Any,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Apply a paired positive rescaling of an Abelian coupling coordinate.

    For ``alpha_1_new=q*alpha_1``, inverse coupling and its affine beta
    coefficient divide by ``q`` while the electromagnetic weight multiplying
    the inverse coupling is multiplied by ``q``.  Thus the electromagnetic
    relation is invariant, although equality of this coordinate with a
    non-Abelian coupling is not.  A physical embedding fixing ``q`` is not
    derived here.
    """

    inverse = _exact_positive(inverse_coupling, "inverse_coupling")
    beta = _exact_real(beta_coefficient, "beta_coefficient")
    weight = _exact_positive(electromagnetic_weight, "electromagnetic_weight")
    factor = _exact_positive(coordinate_factor, "coordinate_factor")
    return (
        sp.simplify(inverse / factor),
        sp.simplify(beta / factor),
        sp.simplify(weight * factor),
    )


def affine_unification_scale(reference_scale: Any, running_coordinate: Any) -> sp.Expr:
    """Return ``reference_scale*exp(2*pi*B)`` for a declared WM3 convention."""

    reference = _exact_positive(reference_scale, "reference_scale")
    coordinate = _exact_real(running_coordinate, "running_coordinate")
    return sp.simplify(reference * sp.exp(2 * sp.pi * coordinate))


def one_loop_inverse_coupling_squared(
    scale: Any,
    reference_scale: Any,
    reference_coupling: Any,
    beta_coefficient: Any,
) -> sp.Expr:
    """Return the exact solution for ``1/g(scale)**2`` of a declared flow.

    The premise is
    ``scale*dg/dscale = -beta_coefficient*g**3/(16*pi**2)``.
    This helper does not derive the beta function or its coefficient.
    """

    scale_value = _positive(scale, "scale")
    reference_value = _positive(reference_scale, "reference_scale")
    coupling_value = _positive(reference_coupling, "reference_coupling")
    coefficient_value = _positive(beta_coefficient, "beta_coefficient")
    return sp.simplify(
        1 / coupling_value**2
        + coefficient_value
        * sp.log(scale_value / reference_value)
        / (8 * sp.pi**2)
    )


def one_loop_transmutation_scale(
    reference_scale: Any,
    reference_coupling: Any,
    beta_coefficient: Any,
) -> sp.Expr:
    """Return the conditional zero of the inverse one-loop coupling.

    The returned expression is invariant under changes of reference point only
    when the reference coupling runs according to the declared one-loop flow.
    """

    reference_value = _positive(reference_scale, "reference_scale")
    coupling_value = _positive(reference_coupling, "reference_coupling")
    coefficient_value = _positive(beta_coefficient, "beta_coefficient")
    return sp.simplify(
        reference_value
        * sp.exp(-8 * sp.pi**2 / (coefficient_value * coupling_value**2))
    )


def single_scale_tension(scale: Any, dimensionless_ratio: Any) -> sp.Expr:
    """Return the conditional mass-dimension-two single-scale form.

    If an independently existing tension has mass dimension two and ``scale``
    is its only dimensionful input with mass dimension one, dimensional
    homogeneity fixes only the power: ``tension = ratio * scale**2``.  The
    dimensionless ratio remains an unconstrained, load-bearing premise.  This
    helper neither establishes that a tension exists nor supplies a
    confinement mechanism.
    """

    scale_value = _positive(scale, "scale")
    ratio_value = _positive(dimensionless_ratio, "dimensionless_ratio")
    return sp.simplify(ratio_value * scale_value**2)


def transmuted_mass_coordinate(
    coupling_squared: Any,
    beta_coefficient: Any,
    mass_energy_ratio: Any,
) -> sp.Expr:
    """Return a conditional mass coordinate tied to a transmuted scale.

    If ``mu0=S*c/a``, ``Lambda=mu0*exp(-8*pi**2/(b0*g0**2))``, and
    ``m*c**2=q*Lambda``, then ``m*c*a/S=q*exp(-8*pi**2/(b0*g0**2))``.
    The coupling squared, beta coefficient, and mass-energy ratio ``q`` remain
    independent inputs.
    """

    coupling_value = _positive(coupling_squared, "coupling_squared")
    coefficient_value = _positive(beta_coefficient, "beta_coefficient")
    ratio_value = _positive(mass_energy_ratio, "mass_energy_ratio")
    return sp.simplify(
        ratio_value
        * sp.exp(-8 * sp.pi**2 / (coefficient_value * coupling_value))
    )

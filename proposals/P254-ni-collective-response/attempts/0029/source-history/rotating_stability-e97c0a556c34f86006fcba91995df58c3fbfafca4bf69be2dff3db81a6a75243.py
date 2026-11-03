"""Exact finite rotating-frame and axisymmetric-rotor stability algebra.

The finite matrix theorem states what a co-rotating change of variables does
and keeps that reduction separate from spectral stability.  The rotor theorem
classifies a declared free symmetric top in body angular-velocity space.  It
does not construct a field action, a rotating field solution, a collective
Skyrme inertia, or a full-field Floquet operator.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _square_matrix(value: Any, name: str) -> sp.ImmutableMatrix:
    matrix = sp.ImmutableMatrix(value)
    if matrix.rows == 0 or matrix.rows != matrix.cols:
        raise ValueError(f"{name} must be a non-empty square matrix")
    return matrix


def _positive_quantity(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return expression


def _immutable_simplified(matrix: sp.MatrixBase) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(matrix.applyfunc(sp.simplify))


@dataclass(frozen=True)
class FiniteMatrixPowerEvidence:
    """Exact eigenvalue and Jordan diagnostics for powers of one matrix."""

    matrix: sp.ImmutableMatrix
    eigenvalues: tuple[sp.Expr, ...]
    algebraic_multiplicities: tuple[int, ...]
    geometric_multiplicities: tuple[int, ...]
    modulus_squared: tuple[sp.Expr, ...]
    all_inside_closed_unit_disk: bool | None
    unit_circle_eigenvalues_semisimple: bool | None
    powers_bounded: bool | None


def finite_matrix_power_evidence(matrix: Any) -> FiniteMatrixPowerEvidence:
    r"""Classify the exact finite-matrix power-boundedness criterion.

    Powers are bounded exactly when every eigenvalue lies in the closed unit
    disk and every eigenvalue on its boundary is semisimple.  ``None`` is
    returned when SymPy cannot decide a required sign or equality.
    """

    value = _square_matrix(matrix, "matrix")
    multiplicities = value.eigenvals()
    ordered = tuple(sorted(multiplicities, key=sp.sstr))
    algebraic = tuple(int(multiplicities[item]) for item in ordered)
    geometric = tuple(
        value.rows - int((value - item * sp.eye(value.rows)).rank())
        for item in ordered
    )
    modulus_squared = tuple(
        sp.simplify(sp.conjugate(item) * item) for item in ordered
    )
    offsets = tuple(sp.simplify(item - 1) for item in modulus_squared)

    inside_flags: list[bool | None] = []
    boundary_flags: list[bool | None] = []
    for modulus, offset in zip(modulus_squared, offsets, strict=True):
        if offset == 0:
            inside_flags.append(True)
            boundary_flags.append(True)
        elif modulus.func == sp.exp and modulus.args[0].is_positive is True:
            inside_flags.append(False)
            boundary_flags.append(False)
        elif modulus.func == sp.exp and modulus.args[0].is_negative is True:
            inside_flags.append(True)
            boundary_flags.append(False)
        else:
            inside_flags.append(offset.is_negative)
            boundary_flags.append(False if offset.is_zero is False else None)

    if any(flag is False for flag in inside_flags):
        all_inside: bool | None = False
    elif all(flag is True for flag in inside_flags):
        all_inside = True
    else:
        all_inside = None

    boundary_semisimple_flags: list[bool | None] = []
    for is_boundary, alg, geom in zip(
        boundary_flags, algebraic, geometric, strict=True
    ):
        if is_boundary is True:
            boundary_semisimple_flags.append(alg == geom)
        elif is_boundary is False:
            boundary_semisimple_flags.append(True)
        else:
            boundary_semisimple_flags.append(None)
    if any(flag is False for flag in boundary_semisimple_flags):
        boundary_semisimple: bool | None = False
    elif all(flag is True for flag in boundary_semisimple_flags):
        boundary_semisimple = True
    else:
        boundary_semisimple = None

    if all_inside is False or boundary_semisimple is False:
        powers_bounded: bool | None = False
    elif all_inside is True and boundary_semisimple is True:
        powers_bounded = True
    else:
        powers_bounded = None

    return FiniteMatrixPowerEvidence(
        matrix=value,
        eigenvalues=ordered,
        algebraic_multiplicities=algebraic,
        geometric_multiplicities=geometric,
        modulus_squared=modulus_squared,
        all_inside_closed_unit_disk=all_inside,
        unit_circle_eigenvalues_semisimple=boundary_semisimple,
        powers_bounded=powers_bounded,
    )


@dataclass(frozen=True)
class CoRotatingLinearSystemEvidence:
    """Exact finite-dimensional rotating-frame reduction data."""

    body_generator: sp.ImmutableMatrix
    frame_generator: sp.ImmutableMatrix
    transformed_generator: sp.ImmutableMatrix
    period: sp.Expr
    frame_at_period: sp.ImmutableMatrix
    frame_periodic: bool
    transformed_monodromy: sp.ImmutableMatrix
    laboratory_monodromy: sp.ImmutableMatrix
    generator_identity_residual: sp.ImmutableMatrix
    transformed_power_evidence: FiniteMatrixPowerEvidence
    laboratory_power_evidence: FiniteMatrixPowerEvidence


def co_rotating_linear_system_evidence(
    body_generator: Any,
    frame_generator: Any,
    period: Any,
) -> CoRotatingLinearSystemEvidence:
    r"""Return the exact finite rotating-frame fundamental matrices.

    Let ``Q(t)=exp(t*K)`` and let the laboratory equation be
    ``x_dot=Q(t)*B*Q(t)^(-1)*x``.  With ``x=Q(t)*y`` the transformed equation is
    ``y_dot=(B-K)*y``.  Therefore the laboratory fundamental matrix is
    ``Q(t)*exp(t*(B-K))``.  At a period ``T`` with ``Q(T)=I``, its monodromy is
    exactly ``exp(T*(B-K))``.  The returned Jordan diagnostics decide finite-
    matrix power boundedness separately; time independence alone supplies no
    stability verdict.
    """

    body = _square_matrix(body_generator, "body_generator")
    frame = _square_matrix(frame_generator, "frame_generator")
    if frame.shape != body.shape:
        raise ValueError("body_generator and frame_generator must have equal shape")
    cycle = _positive_quantity(period, "period")
    transformed = _immutable_simplified(body - frame)
    frame_at_period = _immutable_simplified((cycle * frame).exp())
    transformed_monodromy = _immutable_simplified((cycle * transformed).exp())
    laboratory_monodromy = _immutable_simplified(
        frame_at_period * transformed_monodromy
    )
    identity = sp.ImmutableMatrix.eye(body.rows)
    periodic = frame_at_period == identity
    residual = _immutable_simplified(body - frame - transformed)
    return CoRotatingLinearSystemEvidence(
        body_generator=body,
        frame_generator=frame,
        transformed_generator=transformed,
        period=cycle,
        frame_at_period=frame_at_period,
        frame_periodic=periodic,
        transformed_monodromy=transformed_monodromy,
        laboratory_monodromy=laboratory_monodromy,
        generator_identity_residual=residual,
        transformed_power_evidence=finite_matrix_power_evidence(
            transformed_monodromy
        ),
        laboratory_power_evidence=finite_matrix_power_evidence(
            laboratory_monodromy
        ),
    )


@dataclass(frozen=True)
class AxisymmetricTransverseRotorEvidence:
    """Exact body-angular-velocity theorem for an oblate free symmetric top."""

    transverse_inertia: sp.Expr
    symmetry_axis_inertia: sp.Expr
    angular_speed: sp.Expr
    time: sp.Symbol
    angular_velocity_symbols: tuple[sp.Symbol, sp.Symbol, sp.Symbol]
    euler_rhs: sp.ImmutableMatrix
    base_equilibrium: sp.ImmutableMatrix
    linearized_generator: sp.ImmutableMatrix
    fundamental_matrix: sp.ImmutableMatrix
    monodromy: sp.ImmutableMatrix
    monodromy_power_evidence: FiniteMatrixPowerEvidence
    angular_momentum_norm_squared: sp.Expr
    twice_rotational_energy: sp.Expr
    invariant_derivatives: tuple[sp.Expr, sp.Expr]
    exact_perturbed_solution: sp.ImmutableMatrix
    exact_solution_residual: sp.ImmutableMatrix
    equilibrium_circle_distance_squared: sp.Expr
    fixed_equilibrium_initial_distance_squared: sp.Expr
    fixed_equilibrium_witness_time: sp.Expr
    fixed_equilibrium_witness_distance_squared: sp.Expr


def axisymmetric_transverse_rotor_evidence(
    transverse_inertia: Any,
    symmetry_axis_inertia: Any,
    angular_speed: Any,
    time: sp.Symbol,
) -> AxisymmetricTransverseRotorEvidence:
    r"""Classify transverse-axis rotation of a declared oblate free rotor.

    The principal inertias are ``diag(A,A,C)`` with ``C>A>0``.  The circle
    ``(Omega*cos(phi),Omega*sin(phi),0)`` consists of Euler equilibria.  At one
    member the nonzero nilpotent linearization has an unbounded ``I+t*J``
    fundamental matrix even though every monodromy eigenvalue is one.  Exact
    nearby solutions prove that a fixed member is not Lyapunov stable.  The
    distance to the entire equilibrium circle is constant, so that set is
    stable in body angular-velocity space.  No orientation-space or field-
    theory stability is inferred.
    """

    if not isinstance(time, sp.Symbol):
        raise ValueError("time must be a SymPy Symbol")
    A = _positive_quantity(transverse_inertia, "transverse_inertia")
    C = _positive_quantity(symmetry_axis_inertia, "symmetry_axis_inertia")
    difference = sp.simplify(C - A)
    if difference.is_positive is not True:
        raise ValueError("symmetry_axis_inertia must be provably greater")
    Omega = _positive_quantity(angular_speed, "angular_speed")
    w1, w2, w3 = sp.symbols("w1 w2 w3", real=True)
    vector = sp.ImmutableMatrix([w1, w2, w3])
    rhs = sp.ImmutableMatrix(
        [
            (A - C) * w2 * w3 / A,
            (C - A) * w3 * w1 / A,
            sp.S.Zero,
        ]
    )
    base = sp.ImmutableMatrix([Omega, 0, 0])
    substitutions = dict(zip((w1, w2, w3), base, strict=True))
    linearized = _immutable_simplified(rhs.jacobian(vector).subs(substitutions))
    fundamental = _immutable_simplified(sp.eye(3) + time * linearized)
    period = sp.simplify(2 * sp.pi / Omega)
    monodromy = _immutable_simplified(fundamental.subs(time, period))

    angular_momentum = sp.simplify(
        A**2 * (w1**2 + w2**2) + C**2 * w3**2
    )
    twice_energy = sp.simplify(A * (w1**2 + w2**2) + C * w3**2)
    invariant_derivatives = tuple(
        sp.simplify(
            sum(sp.diff(invariant, item) * rhs[index] for index, item in enumerate(vector))
        )
        for invariant in (angular_momentum, twice_energy)
    )

    radius, phase, epsilon = sp.symbols(
        "transverse_radius phase axial_perturbation", positive=True
    )
    precession_rate = sp.simplify(difference * epsilon / A)
    solution = sp.ImmutableMatrix(
        [
            radius * sp.cos(phase + precession_rate * time),
            radius * sp.sin(phase + precession_rate * time),
            epsilon,
        ]
    )
    solution_rhs = rhs.subs(dict(zip((w1, w2, w3), solution, strict=True)))
    solution_residual = _immutable_simplified(
        sp.diff(solution, time) - solution_rhs
    )
    set_distance_squared = sp.simplify((radius - Omega) ** 2 + epsilon**2)
    initial_distance_squared = epsilon**2
    witness_time = sp.simplify(sp.pi * A / (2 * difference * epsilon))
    witness = solution.subs({radius: Omega, phase: 0, time: witness_time})
    witness_distance_squared = sp.simplify((witness - base).dot(witness - base))

    return AxisymmetricTransverseRotorEvidence(
        transverse_inertia=A,
        symmetry_axis_inertia=C,
        angular_speed=Omega,
        time=time,
        angular_velocity_symbols=(w1, w2, w3),
        euler_rhs=rhs,
        base_equilibrium=base,
        linearized_generator=linearized,
        fundamental_matrix=fundamental,
        monodromy=monodromy,
        monodromy_power_evidence=finite_matrix_power_evidence(monodromy),
        angular_momentum_norm_squared=angular_momentum,
        twice_rotational_energy=twice_energy,
        invariant_derivatives=invariant_derivatives,
        exact_perturbed_solution=solution,
        exact_solution_residual=solution_residual,
        equilibrium_circle_distance_squared=set_distance_squared,
        fixed_equilibrium_initial_distance_squared=initial_distance_squared,
        fixed_equilibrium_witness_time=witness_time,
        fixed_equilibrium_witness_distance_squared=witness_distance_squared,
    )


@dataclass(frozen=True)
class AxisymmetricDensityInertiaRelation:
    """Ordinary density inertia versus its normalized STF second moment."""

    radial_second_moment: sp.Expr
    axial_second_moment: sp.Expr
    transverse_axis_inertia: sp.Expr
    symmetry_axis_inertia: sp.Expr
    normalized_stf_zz: sp.Expr
    relation_residual: sp.Expr


def axisymmetric_density_inertia_relation(
    radial_second_moment: Any,
    axial_second_moment: Any,
) -> AxisymmetricDensityInertiaRelation:
    r"""Return ``I_zz-I_xx=-(3/2)*(M_STF)_zz`` for a declared density.

    ``radial_second_moment`` is ``integral rho*(x^2+y^2+z^2)`` and
    ``axial_second_moment`` is ``integral rho*z^2``.  This ordinary second-
    moment inertia need not equal a field theory's collective rotational metric.
    """

    radial = sp.sympify(radial_second_moment)
    axial = sp.sympify(axial_second_moment)
    symmetry_inertia = sp.simplify(radial - axial)
    transverse_inertia = sp.simplify((radial + axial) / 2)
    stf_zz = sp.simplify(axial - radial / 3)
    residual = sp.simplify(
        symmetry_inertia - transverse_inertia + sp.Rational(3, 2) * stf_zz
    )
    return AxisymmetricDensityInertiaRelation(
        radial_second_moment=radial,
        axial_second_moment=axial,
        transverse_axis_inertia=transverse_inertia,
        symmetry_axis_inertia=symmetry_inertia,
        normalized_stf_zz=stf_zz,
        relation_residual=residual,
    )

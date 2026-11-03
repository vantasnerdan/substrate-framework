"""Exact pure-spin-1 invariants and conditional mean-field selection.

The functions in this module concern a three-component *pure* spinor in the
ordered ``m=(+1, 0, -1)`` basis.  The phase-selection ledger assumes the
displayed fixed-density spin functional; it does not derive a material
coupling, a spatial condensate solution, or an experimentally realized phase.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, Sequence

import sympy as sp


Spin1ProjectiveOrbit = Literal["polar", "ferromagnetic", "intermediate"]


_SQRT_TWO = sp.sqrt(2)
_SPIN1_MATRICES = (
    sp.ImmutableMatrix([[0, 1, 0], [1, 0, 1], [0, 1, 0]]) / _SQRT_TWO,
    sp.ImmutableMatrix([[0, -sp.I, 0], [sp.I, 0, -sp.I], [0, sp.I, 0]])
    / _SQRT_TWO,
    sp.ImmutableMatrix([[1, 0, 0], [0, 0, 0], [0, 0, -1]]),
)


def _column_three(values: Sequence[Any] | sp.MatrixBase, name: str) -> sp.ImmutableMatrix:
    if isinstance(values, sp.MatrixBase):
        entries = tuple(values)
    else:
        try:
            entries = tuple(values)
        except TypeError as error:
            raise ValueError(f"{name} must contain exactly three components") from error
    if len(entries) != 3:
        raise ValueError(f"{name} must contain exactly three components")
    return sp.ImmutableMatrix(3, 1, [sp.sympify(value) for value in entries])


def _real_expression(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_real is not True:
        raise ValueError(f"{name} must be explicitly real")
    return expression


def _positive_expression(value: Any, name: str) -> sp.Expr:
    expression = _real_expression(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be explicitly positive")
    return expression


def _zero_decision(value: Any, name: str) -> bool:
    expression = sp.simplify(value)
    if expression == 0 or expression.is_zero is True:
        return True
    if expression.is_zero is False:
        return False
    raise ValueError(f"{name} could not be decided exactly")


def spin1_matrices() -> tuple[sp.ImmutableMatrix, sp.ImmutableMatrix, sp.ImmutableMatrix]:
    """Return ``(F_x,F_y,F_z)`` for the standard spin-1 representation."""

    return _SPIN1_MATRICES


def spin1_norm(spinor: Sequence[Any] | sp.MatrixBase) -> sp.Expr:
    """Return the Hermitian norm ``Psi^dagger Psi``."""

    state = _column_three(spinor, "spinor")
    return sp.simplify((state.conjugate().T * state)[0])


def spin1_expectation(
    spinor: Sequence[Any] | sp.MatrixBase,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return the three exact spin expectations ``Psi^dagger F_a Psi``."""

    state = _column_three(spinor, "spinor")
    return tuple(
        sp.simplify((state.conjugate().T * matrix * state)[0])
        for matrix in _SPIN1_MATRICES
    )  # type: ignore[return-value]


def spin1_magnitude_squared(spinor: Sequence[Any] | sp.MatrixBase) -> sp.Expr:
    """Return ``|Psi^dagger F Psi|^2`` without imposing unit normalization."""

    return sp.simplify(sum(value**2 for value in spin1_expectation(spinor)))


def spin1_singlet_amplitude(spinor: Sequence[Any] | sp.MatrixBase) -> sp.Expr:
    """Return ``A=Psi_0**2-2*Psi_+*Psi_-`` in the declared basis."""

    plus, zero, minus = _column_three(spinor, "spinor")
    return sp.simplify(zero**2 - 2 * plus * minus)


def spin1_singlet_magnitude_squared(
    spinor: Sequence[Any] | sp.MatrixBase,
) -> sp.Expr:
    """Return the exact squared magnitude ``conjugate(A)*A``."""

    amplitude = spin1_singlet_amplitude(spinor)
    return sp.simplify(sp.conjugate(amplitude) * amplitude)


def spin1_invariant_residual(spinor: Sequence[Any] | sp.MatrixBase) -> sp.Expr:
    """Return the residual in ``|f|^2+|A|^2=n^2``.

    The identity makes both density scaling and the sharp upper endpoint
    explicit.  For a nonzero pure spinor, ``A=0`` is the ferromagnetic orbit;
    ``f=0`` is the polar orbit.
    """

    norm = spin1_norm(spinor)
    return sp.simplify(
        spin1_magnitude_squared(spinor)
        + spin1_singlet_magnitude_squared(spinor)
        - norm**2
    )


def spin1_to_cartesian(spinor: Sequence[Any] | sp.MatrixBase) -> sp.ImmutableMatrix:
    """Map spherical spin components to the complex Cartesian vector ``d``.

    The convention is ``Psi_+=-(d_x-i*d_y)/sqrt(2)``, ``Psi_0=d_z``, and
    ``Psi_-=(d_x+i*d_y)/sqrt(2)``.
    """

    plus, zero, minus = _column_three(spinor, "spinor")
    return sp.ImmutableMatrix(
        [
            (minus - plus) / _SQRT_TWO,
            -sp.I * (plus + minus) / _SQRT_TWO,
            zero,
        ]
    )


def cartesian_to_spin1(vector: Sequence[Any] | sp.MatrixBase) -> sp.ImmutableMatrix:
    """Apply the inverse complex-Cartesian to spherical-spinor map."""

    x_value, y_value, z_value = _column_three(vector, "vector")
    return sp.ImmutableMatrix(
        [
            -(x_value - sp.I * y_value) / _SQRT_TWO,
            z_value,
            (x_value + sp.I * y_value) / _SQRT_TWO,
        ]
    )


def spin1_cartesian_spin(
    vector: Sequence[Any] | sp.MatrixBase,
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return ``2 Re(d) cross Im(d)`` for a complex Cartesian vector."""

    cartesian = _column_three(vector, "vector")
    real_part = sp.ImmutableMatrix([sp.re(value) for value in cartesian])
    imaginary_part = sp.ImmutableMatrix([sp.im(value) for value in cartesian])
    result = 2 * real_part.cross(imaginary_part)
    return tuple(sp.simplify(value) for value in result)  # type: ignore[return-value]


@dataclass(frozen=True)
class Spin1OrbitLedger:
    """Exact invariant ledger for one nonzero pure spin-1 state."""

    norm: sp.Expr
    spin_squared: sp.Expr
    singlet_squared: sp.Expr
    invariant_residual: sp.Expr
    projective_orbit: Spin1ProjectiveOrbit
    saturates_lower_endpoint: bool
    saturates_upper_endpoint: bool


def spin1_orbit_ledger(
    spinor: Sequence[Any] | sp.MatrixBase,
) -> Spin1OrbitLedger:
    """Classify an exactly decidable nonzero pure spinor by its sharp endpoint.

    Projective means that global phase is quotiented.  The polar equality orbit
    is ``RP^2`` and the ferromagnetic equality orbit is ``S^2`` under spatial
    ``SO(3)``.  This function deliberately accepts a vector, not a density
    matrix; mixed zero-spin states are a different set.
    """

    norm = spin1_norm(spinor)
    if norm.is_positive is not True:
        raise ValueError("spinor norm must be explicitly positive")
    spin_squared = spin1_magnitude_squared(spinor)
    singlet_squared = spin1_singlet_magnitude_squared(spinor)
    residual = spin1_invariant_residual(spinor)
    if not _zero_decision(residual, "spin-1 invariant residual"):
        raise ValueError("spin-1 invariant identity failed")
    lower = _zero_decision(spin_squared, "spin magnitude squared")
    upper = _zero_decision(singlet_squared, "singlet magnitude squared")
    if lower:
        orbit: Spin1ProjectiveOrbit = "polar"
    elif upper:
        orbit = "ferromagnetic"
    else:
        orbit = "intermediate"
    return Spin1OrbitLedger(
        norm=norm,
        spin_squared=spin_squared,
        singlet_squared=singlet_squared,
        invariant_residual=residual,
        projective_orbit=orbit,
        saturates_lower_endpoint=lower,
        saturates_upper_endpoint=upper,
    )


def spin1_mean_field_energy(
    spinor: Sequence[Any] | sp.MatrixBase,
    coupling: Any,
) -> sp.Expr:
    """Return the supplied spin energy ``(c2/2)|Psi^dagger F Psi|^2``."""

    coefficient = _real_expression(coupling, "coupling")
    return sp.simplify(coefficient * spin1_magnitude_squared(spinor) / 2)


@dataclass(frozen=True)
class Spin1MeanFieldSelection:
    """Complete endpoint selection for the supplied fixed-density functional."""

    density: sp.Expr
    coupling: sp.Expr
    attainable_spin_squared: tuple[sp.Expr, sp.Expr]
    minimizing_projective_orbits: tuple[str, ...]
    maximizing_projective_orbits: tuple[str, ...]
    minimum_energy: sp.Expr
    maximum_energy: sp.Expr
    polar_minus_ferromagnetic_energy: sp.Expr


def fixed_density_spin1_selection(
    density: Any,
    coupling: Any,
) -> Spin1MeanFieldSelection:
    """Minimize the supplied spin energy over pure spinors of fixed norm.

    For ``c2>0`` the projective polar orbit minimizes; for ``c2<0`` the
    projective ferromagnetic orbit minimizes; for ``c2=0`` every pure spin-1
    ray is degenerate.  The endpoint energy gap scales as ``density**2``.
    """

    norm = _positive_expression(density, "density")
    coefficient = _real_expression(coupling, "coupling")
    endpoint = sp.simplify(norm**2)
    ferro_energy = sp.simplify(coefficient * endpoint / 2)
    gap = sp.simplify(-ferro_energy)
    if coefficient.is_positive is True:
        minimizers = ("polar",)
        maximizers = ("ferromagnetic",)
        minimum_energy = sp.Integer(0)
        maximum_energy = ferro_energy
    elif coefficient.is_negative is True:
        minimizers = ("ferromagnetic",)
        maximizers = ("polar",)
        minimum_energy = ferro_energy
        maximum_energy = sp.Integer(0)
    elif coefficient.is_zero is True or coefficient == 0:
        minimizers = ("all_pure_spin1_rays",)
        maximizers = ("all_pure_spin1_rays",)
        minimum_energy = sp.Integer(0)
        maximum_energy = sp.Integer(0)
    else:
        raise ValueError("coupling sign must be exactly decidable")
    return Spin1MeanFieldSelection(
        density=norm,
        coupling=coefficient,
        attainable_spin_squared=(sp.Integer(0), endpoint),
        minimizing_projective_orbits=minimizers,
        maximizing_projective_orbits=maximizers,
        minimum_energy=sp.simplify(minimum_energy),
        maximum_energy=sp.simplify(maximum_energy),
        polar_minus_ferromagnetic_energy=gap,
    )

"""Exact conditional energy and topology ledgers for angular defects.

These APIs separate one-field annular self energy, a declared matched-shell
split model, and the projective/full polar topology.  They do not solve a
finite-domain multi-core PDE or supply material stiffnesses or core energies.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp


def _real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_real is not True:
        raise ValueError(f"{name} must be explicitly real")
    return expression


def _positive(value: Any, name: str) -> sp.Expr:
    expression = _real(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be explicitly positive")
    return expression


def _nonnegative(value: Any, name: str) -> sp.Expr:
    expression = _real(value, name)
    if expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be explicitly nonnegative")
    return expression


def _integer(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_integer is not True:
        raise ValueError(f"{name} must be an integer")
    return expression


def _positive_concrete_integer(value: Any, name: str) -> int:
    expression = _integer(value, name)
    if expression.is_number is not True or expression.is_positive is not True:
        raise ValueError(f"{name} must be a positive concrete integer")
    return int(expression)


def _radii(inner_radius: Any, outer_radius: Any) -> tuple[sp.Expr, sp.Expr]:
    inner = _positive(inner_radius, "inner_radius")
    outer = _positive(outer_radius, "outer_radius")
    if sp.simplify(outer - inner).is_positive is not True:
        raise ValueError("outer_radius must be explicitly greater than inner_radius")
    return inner, outer


def annular_angular_energy(
    stiffness: Any,
    charge: Any,
    inner_radius: Any,
    outer_radius: Any,
) -> sp.Expr:
    """Return the exact uniform-winding Dirichlet energy on an annulus.

    For ``E=(K/2) integral |grad(theta)|^2 d^2x`` and
    ``theta=q*phi``, the result is ``pi*K*q^2*log(R/xi)``.  The same value is
    the sharp lower bound at fixed degree ``q`` on every concentric circle.
    """

    coefficient = _positive(stiffness, "stiffness")
    winding = _real(charge, "charge")
    inner, outer = _radii(inner_radius, outer_radius)
    return sp.simplify(sp.pi * coefficient * winding**2 * sp.log(outer / inner))


@dataclass(frozen=True)
class EqualSplitShellLedger:
    """Near/far matched-shell energy for equal pieces of fixed total charge."""

    pieces: int
    total_charge: sp.Expr
    near_energy: sp.Expr
    far_energy: sp.Expr
    split_core_energy: sp.Expr
    split_total_energy: sp.Expr
    unsplit_field_energy: sp.Expr
    unsplit_total_energy: sp.Expr
    split_minus_unsplit: sp.Expr
    field_energy_ratio: sp.Expr
    independent_copy_ratio: sp.Expr


def equal_split_shell_ledger(
    stiffness: Any,
    total_charge: Any,
    pieces: Any,
    core_radius: Any,
    separation_radius: Any,
    outer_radius: Any,
    *,
    piece_core_energy: Any = 0,
    unsplit_core_energy: Any = 0,
) -> EqualSplitShellLedger:
    """Return an exact declared matched-shell split-energy ledger.

    Each of ``n`` equal charges ``Q/n`` occupies a near annulus from ``xi``
    to ``d``.  The common far annulus from ``d`` to ``R`` carries the fixed
    total charge ``Q``.  This is a scale-matched shell model, not an exact
    finite-domain multi-core solution.  It makes the far-field term omitted by
    independent-copy comparisons explicit.
    """

    coefficient = _positive(stiffness, "stiffness")
    charge = _real(total_charge, "total_charge")
    count = _positive_concrete_integer(pieces, "pieces")
    inner, outer = _radii(core_radius, outer_radius)
    separation = _positive(separation_radius, "separation_radius")
    if sp.simplify(separation - inner).is_nonnegative is not True:
        raise ValueError("separation_radius must be at least core_radius")
    if sp.simplify(outer - separation).is_nonnegative is not True:
        raise ValueError("separation_radius must not exceed outer_radius")
    piece_core = _nonnegative(piece_core_energy, "piece_core_energy")
    unsplit_core = _nonnegative(unsplit_core_energy, "unsplit_core_energy")
    near = sp.simplify(
        sp.pi
        * coefficient
        * charge**2
        / count
        * sp.log(separation / inner)
    )
    far = sp.simplify(
        sp.pi * coefficient * charge**2 * sp.log(outer / separation)
    )
    split_cores = sp.simplify(count * piece_core)
    split_total = sp.simplify(near + far + split_cores)
    unsplit_field = sp.simplify(
        sp.pi * coefficient * charge**2 * sp.log(outer / inner)
    )
    unsplit_total = sp.simplify(unsplit_field + unsplit_core)
    return EqualSplitShellLedger(
        pieces=count,
        total_charge=charge,
        near_energy=near,
        far_energy=far,
        split_core_energy=split_cores,
        split_total_energy=split_total,
        unsplit_field_energy=unsplit_field,
        unsplit_total_energy=unsplit_total,
        split_minus_unsplit=sp.simplify(split_total - unsplit_total),
        field_energy_ratio=sp.simplify((near + far) / unsplit_field),
        independent_copy_ratio=sp.Rational(1, count),
    )


@dataclass(frozen=True)
class PolarTopologyLedger:
    """Fundamental-group distinction for projective and full polar order."""

    projective_manifold: str
    projective_fundamental_group: str
    projective_generator_order: int
    full_polar_manifold: str
    full_polar_fundamental_group: str
    full_generator_order: None
    full_generator_square: str


def polar_topology_ledger() -> PolarTopologyLedger:
    """Return the exact loop-group distinction relevant to half defects."""

    return PolarTopologyLedger(
        projective_manifold="RP2",
        projective_fundamental_group="Z2",
        projective_generator_order=2,
        full_polar_manifold="(S2 x U1)/Z2",
        full_polar_fundamental_group="Z",
        full_generator_order=None,
        full_generator_square="nontrivial_integer_phase_vortex",
    )


@dataclass(frozen=True)
class FullPolarDeckTransformation:
    """One deck transformation of the universal cover ``S2 x R``."""

    index: sp.Expr
    director_sign: sp.Expr
    phase_shift: sp.Expr


def full_polar_deck_transformation(index: Any) -> FullPolarDeckTransformation:
    """Return ``g^k:(d,t)->((-1)^k d,t+k*pi)`` for integer ``k``."""

    step = _integer(index, "index")
    return FullPolarDeckTransformation(
        index=step,
        director_sign=sp.simplify((-1) ** step),
        phase_shift=sp.simplify(sp.pi * step),
    )


def projective_rp2_loop_class(generator_steps: Any) -> sp.Expr:
    """Return the projective loop class in ``Z2``."""

    steps = _integer(generator_steps, "generator_steps")
    return sp.Mod(steps, 2)


def full_polar_loop_class(generator_steps: Any) -> sp.Expr:
    """Return the full polar loop class in the deck group ``Z``."""

    return _integer(generator_steps, "generator_steps")


@dataclass(frozen=True)
class HalfQuantumPairLedger:
    """Isolated-self-energy comparison with explicit stiffnesses and cores."""

    one_half_field_energy: sp.Expr
    pair_field_energy: sp.Expr
    pair_total_energy: sp.Expr
    integer_field_energy: sp.Expr
    integer_total_energy: sp.Expr
    pair_minus_integer: sp.Expr
    zero_core_field_ratio: sp.Expr


def half_quantum_pair_ledger(
    phase_stiffness: Any,
    director_stiffness: Any,
    inner_radius: Any,
    outer_radius: Any,
    *,
    half_core_energy: Any = 0,
    integer_core_energy: Any = 0,
) -> HalfQuantumPairLedger:
    """Compare two isolated half textures with one integer phase vortex.

    The declared functional is
    ``(1/2) integral (K_phase|grad theta|^2+K_dir|grad d|^2)``.
    A half texture has phase and great-circle director charges ``1/2``.
    The result is still an isolated-copy comparison; use
    :func:`equal_split_shell_ledger` to retain a common far field.
    """

    phase = _positive(phase_stiffness, "phase_stiffness")
    director = _positive(director_stiffness, "director_stiffness")
    inner, outer = _radii(inner_radius, outer_radius)
    half_core = _nonnegative(half_core_energy, "half_core_energy")
    integer_core = _nonnegative(integer_core_energy, "integer_core_energy")
    logarithm = sp.log(outer / inner)
    one_half = sp.simplify(sp.pi * (phase + director) * logarithm / 4)
    pair_field = sp.simplify(2 * one_half)
    pair_total = sp.simplify(pair_field + 2 * half_core)
    integer_field = sp.simplify(sp.pi * phase * logarithm)
    integer_total = sp.simplify(integer_field + integer_core)
    return HalfQuantumPairLedger(
        one_half_field_energy=one_half,
        pair_field_energy=pair_field,
        pair_total_energy=pair_total,
        integer_field_energy=integer_field,
        integer_total_energy=integer_total,
        pair_minus_integer=sp.simplify(pair_total - integer_total),
        zero_core_field_ratio=sp.simplify(pair_field / integer_field),
    )

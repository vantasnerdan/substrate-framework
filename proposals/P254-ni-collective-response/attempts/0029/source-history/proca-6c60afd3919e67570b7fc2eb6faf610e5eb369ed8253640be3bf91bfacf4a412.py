"""Exact ledgers for a separately declared source-free Proca model.

The APIs in this module use the mostly-plus Minkowski metric
``diag(-1,+1,+1,+1)`` and the action density
``-F_{mu nu} F^{mu nu}/4 - m^2 A_mu A^mu/2``. They expose the full vector
momentum kernel, the massive divergence constraint, and one constraint-
compatible half-line boundary-value problem. They do not derive this action
from a scalar sector, supply a London current or material response, identify a
physical weak boson, select the Standard Model, or realize a substrate mass
mechanism.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp


def _exact_expression(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact")
    return sp.simplify(expression)


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = _exact_expression(value, name)
    if expression.is_real is not True:
        raise ValueError(f"{name} must be explicitly real")
    return expression


def _exact_positive(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be explicitly positive")
    return expression


def _exact_real_column(
    values: Sequence[Any],
    name: str,
    *,
    length: int,
) -> sp.ImmutableMatrix:
    entries = tuple(
        _exact_real(value, f"{name}[{index}]")
        for index, value in enumerate(values)
    )
    if len(entries) != length:
        raise ValueError(f"{name} must contain exactly {length} entries")
    return sp.ImmutableMatrix(entries)


@dataclass(frozen=True)
class MostlyPlusProcaMomentumEvidence:
    """Full Fourier-space Proca kernel in the frozen mostly-plus convention."""

    mass: sp.Expr
    frequency: sp.Expr
    spatial_momentum: tuple[sp.Expr, sp.Expr, sp.Expr]
    metric: sp.ImmutableMatrix
    momentum_covector: sp.ImmutableMatrix
    momentum_vector: sp.ImmutableMatrix
    momentum_norm_squared: sp.Expr
    spatial_momentum_norm_squared: sp.Expr
    euler_kernel: sp.ImmutableMatrix
    divergence_contraction: sp.ImmutableMatrix
    divergence_constraint_residual: sp.ImmutableMatrix
    transverse_kernel_factor: sp.Expr
    dispersion_frequency_squared: sp.Expr
    dispersion_residual: sp.Expr

    @property
    def divergence_constraint_certified(self) -> bool:
        """Whether left contraction gives ``-m^2*k_covector`` exactly."""

        return self.divergence_constraint_residual == sp.zeros(1, 4)

    def euler_residual(self, amplitude: Sequence[Any]) -> sp.ImmutableMatrix:
        """Apply the exact mixed-index Euler kernel to a real amplitude."""

        vector = _exact_real_column(amplitude, "amplitude", length=4)
        return sp.ImmutableMatrix(
            (self.euler_kernel * vector).applyfunc(sp.simplify)
        )

    def transversality_residual(self, amplitude: Sequence[Any]) -> sp.Expr:
        """Return ``k_mu*A^mu`` for a supplied real plane-wave amplitude."""

        vector = _exact_real_column(amplitude, "amplitude", length=4)
        return sp.simplify((self.momentum_covector.T * vector)[0])


def mostly_plus_proca_momentum_evidence(
    mass: Any,
    frequency: Any,
    spatial_momentum: Sequence[Any],
) -> MostlyPlusProcaMomentumEvidence:
    """Derive the full source-free Proca plane-wave kernel exactly.

    With phase ``exp(i*k_mu*x^mu)`` and covector
    ``k_mu=(-omega,k_x,k_y,k_z)``, variation of the declared action gives

    ``E^nu=[-(k^2+m^2) delta^nu_sigma+k^nu k_sigma] A^sigma``.

    Left contraction by ``k_nu`` is ``-m^2*k_sigma``. Thus the equation of
    motion derives ``k.A=0`` when ``m>0``; on that transverse subspace its
    remaining factor gives ``omega^2=|k|^2+m^2``.
    """

    mass_expression = _exact_positive(mass, "mass")
    frequency_expression = _exact_real(frequency, "frequency")
    spatial = _exact_real_column(
        spatial_momentum,
        "spatial_momentum",
        length=3,
    )
    metric = sp.ImmutableMatrix(sp.diag(-1, 1, 1, 1))
    covector = sp.ImmutableMatrix(
        [-frequency_expression, spatial[0], spatial[1], spatial[2]]
    )
    vector = sp.ImmutableMatrix(metric * covector)
    momentum_norm = sp.simplify((covector.T * vector)[0])
    spatial_norm = sp.simplify((spatial.T * spatial)[0])
    mass_squared = sp.simplify(mass_expression**2)
    kernel = sp.ImmutableMatrix(
        (
            -(momentum_norm + mass_squared) * sp.eye(4)
            + vector * covector.T
        ).applyfunc(sp.simplify)
    )
    contraction = sp.ImmutableMatrix(
        (covector.T * kernel).applyfunc(sp.simplify)
    )
    expected_contraction = sp.ImmutableMatrix(
        (-mass_squared * covector.T).applyfunc(sp.simplify)
    )
    transverse_factor = sp.simplify(-(momentum_norm + mass_squared))
    frequency_squared = sp.simplify(spatial_norm + mass_squared)
    return MostlyPlusProcaMomentumEvidence(
        mass=mass_expression,
        frequency=frequency_expression,
        spatial_momentum=(spatial[0], spatial[1], spatial[2]),
        metric=metric,
        momentum_covector=covector,
        momentum_vector=vector,
        momentum_norm_squared=momentum_norm,
        spatial_momentum_norm_squared=spatial_norm,
        euler_kernel=kernel,
        divergence_contraction=contraction,
        divergence_constraint_residual=sp.ImmutableMatrix(
            (contraction - expected_contraction).applyfunc(sp.simplify)
        ),
        transverse_kernel_factor=transverse_factor,
        dispersion_frequency_squared=frequency_squared,
        dispersion_residual=sp.simplify(
            transverse_factor
            - (frequency_expression**2 - spatial_norm - mass_squared)
        ),
    )


@dataclass(frozen=True)
class TransverseHalfLineProcaEvidence:
    """Exact static half-line BVP data for a tangential vector component."""

    mass: sp.Expr
    boundary_amplitude: sp.Expr
    coordinate: sp.Symbol
    characteristic_variable: sp.Symbol
    characteristic_polynomial: sp.Expr
    characteristic_roots: tuple[sp.Expr, sp.Expr]
    decay_coefficient: sp.Symbol
    growing_coefficient: sp.Symbol
    general_solution: sp.Expr
    decaying_profile: sp.Expr
    growing_profile: sp.Expr
    equation_residual: sp.Expr
    boundary_residual: sp.Expr
    decay_limit: sp.Expr
    decaying_basis_limit: sp.Expr
    growing_basis_absolute_limit: sp.Expr
    selected_general_solution_residual: sp.Expr
    tangential_divergence_residual: sp.Expr
    longitudinal_divergence: sp.Expr
    inverse_length: sp.Expr
    penetration_length: sp.Expr

    @property
    def bvp_certified(self) -> bool:
        """Whether the returned profile solves the equation and both BVP data."""

        return (
            self.equation_residual == 0
            and self.boundary_residual == 0
            and self.decay_limit == 0
            and self.decaying_basis_limit == 0
            and self.growing_basis_absolute_limit == sp.oo
            and self.selected_general_solution_residual == 0
            and self.tangential_divergence_residual == 0
        )


def transverse_half_line_proca_evidence(
    mass: Any,
    boundary_amplitude: Any,
    *,
    coordinate: sp.Symbol | None = None,
) -> TransverseHalfLineProcaEvidence:
    """Solve a constraint-compatible static Proca half-line problem.

    The profile is the tangential component ``A_y(x)`` on ``x>=0``. Because
    it depends on ``x`` rather than ``y``, its divergence is zero. The general
    static solution contains ``C_decay*exp(-m*x)+C_grow*exp(m*x)``. Decay at
    infinity forces ``C_grow=0`` and ``A_y(0)=A0`` then forces
    ``C_decay=A0``, yielding the unique profile ``A0*exp(-m*x)``.

    By contrast, assigning the same nonconstant profile to the longitudinal
    component ``A_x(x)`` gives the returned nonzero derivative and generally
    violates the massive divergence constraint.
    """

    mass_expression = _exact_positive(mass, "mass")
    amplitude = _exact_real(boundary_amplitude, "boundary_amplitude")
    if coordinate is None:
        x = sp.Symbol("x", nonnegative=True)
    else:
        if not isinstance(coordinate, sp.Symbol) or coordinate.is_real is not True:
            raise ValueError("coordinate must be an explicitly real symbol")
        x = coordinate

    rate = sp.Symbol("r", real=True)
    decay_coefficient = sp.Symbol("C_decay", real=True)
    growing_coefficient = sp.Symbol("C_grow", real=True)
    characteristic = sp.factor(rate**2 - mass_expression**2)
    decaying = sp.simplify(amplitude * sp.exp(-mass_expression * x))
    growing = sp.simplify(amplitude * sp.exp(mass_expression * x))
    general = sp.simplify(
        decay_coefficient * sp.exp(-mass_expression * x)
        + growing_coefficient * sp.exp(mass_expression * x)
    )
    return TransverseHalfLineProcaEvidence(
        mass=mass_expression,
        boundary_amplitude=amplitude,
        coordinate=x,
        characteristic_variable=rate,
        characteristic_polynomial=characteristic,
        characteristic_roots=(-mass_expression, mass_expression),
        decay_coefficient=decay_coefficient,
        growing_coefficient=growing_coefficient,
        general_solution=general,
        decaying_profile=decaying,
        growing_profile=growing,
        equation_residual=sp.simplify(
            sp.diff(decaying, x, 2) - mass_expression**2 * decaying
        ),
        boundary_residual=sp.simplify(decaying.subs(x, 0) - amplitude),
        decay_limit=sp.limit(decaying, x, sp.oo),
        decaying_basis_limit=sp.limit(sp.exp(-mass_expression * x), x, sp.oo),
        growing_basis_absolute_limit=sp.limit(
            sp.Abs(sp.exp(mass_expression * x)),
            x,
            sp.oo,
        ),
        selected_general_solution_residual=sp.simplify(
            general.subs(
                {
                    decay_coefficient: amplitude,
                    growing_coefficient: 0,
                }
            )
            - decaying
        ),
        tangential_divergence_residual=sp.Integer(0),
        longitudinal_divergence=sp.simplify(sp.diff(decaying, x)),
        inverse_length=mass_expression,
        penetration_length=sp.simplify(1 / mass_expression),
    )


@dataclass(frozen=True)
class NormalizedProcaModeEvidence:
    """One-mode normalization data for positive kinetic and mass coefficients."""

    kinetic_coefficient: sp.Expr
    quadratic_coefficient: sp.Expr
    mass_squared: sp.Expr
    mass: sp.Expr
    inverse_length: sp.Expr
    penetration_length: sp.Expr


def normalized_proca_mode_evidence(
    kinetic_coefficient: Any,
    quadratic_coefficient: Any,
) -> NormalizedProcaModeEvidence:
    """Normalize one declared free vector mode with positive coefficients.

    If the kinetic term has coefficient ``kappa>0`` and the quadratic form has
    coefficient ``q>0`` in the same field coordinate, canonical rescaling gives
    ``m^2=q/kappa``. This algebra alone supplies no action origin, material
    response, particle identity, or observed penetration depth.
    """

    kinetic = _exact_positive(kinetic_coefficient, "kinetic_coefficient")
    quadratic = _exact_positive(quadratic_coefficient, "quadratic_coefficient")
    mass_squared = sp.simplify(quadratic / kinetic)
    mass = sp.simplify(sp.sqrt(mass_squared))
    return NormalizedProcaModeEvidence(
        kinetic_coefficient=kinetic,
        quadratic_coefficient=quadratic,
        mass_squared=mass_squared,
        mass=mass,
        inverse_length=mass,
        penetration_length=sp.simplify(1 / mass),
    )

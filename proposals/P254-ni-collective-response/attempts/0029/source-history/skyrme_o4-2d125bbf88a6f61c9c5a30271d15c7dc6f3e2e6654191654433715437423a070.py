"""Exact conditional pointwise identities for an O(4) Skyrme field.

This module starts from declared spatial gradients in ``R**4``.  It proves
algebraic positivity identities for the usual quadratic/quartic static
density and for the associated time-derivative mass operator.  It does not
derive those declarations from a physical action, establish a stationary
field, or imply local, rotating, or dynamical stability.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp


@dataclass(frozen=True)
class O4SkyrmePointwiseEvidence:
    """Exact identities for three declared O(4) spatial gradients."""

    gradients: sp.ImmutableDenseMatrix
    tangent: sp.ImmutableDenseMatrix
    strain: sp.ImmutableDenseMatrix
    quadratic_density: sp.Expr
    quartic_trace_form: sp.Expr
    quartic_minor_sos: sp.Expr
    static_density: sp.Expr
    mass_operator: sp.ImmutableDenseMatrix
    tangent_norm_squared: sp.Expr
    mass_quadratic_form: sp.Expr
    mass_lower_bound_gap: sp.Expr
    mass_minor_sos: sp.Expr
    quartic_identity_residual: sp.Expr
    mass_identity_residual: sp.Expr

    @property
    def quartic_has_sos_certificate(self) -> bool:
        """Return whether the quartic trace form equals its minor-square sum."""

        return self.quartic_identity_residual == 0

    @property
    def mass_has_sharp_lower_bound_certificate(self) -> bool:
        """Return whether ``w.T*M*w - 2*w.T*w`` equals a square sum."""

        return self.mass_identity_residual == 0


def _exact_real_matrix(
    values: Sequence[Sequence[Any]] | Sequence[Any],
    shape: tuple[int, int],
    name: str,
) -> sp.ImmutableDenseMatrix:
    matrix = sp.ImmutableDenseMatrix(values)
    if matrix.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    for entry in matrix:
        if entry.has(sp.Float):
            raise ValueError(f"{name} must use exact entries, not SymPy Float values")
        if entry.is_real is not True:
            raise ValueError(f"{name} entries must be explicitly real")
    return matrix


def o4_skyrme_pointwise_evidence(
    spatial_gradients: Sequence[Sequence[Any]],
    tangent: Sequence[Any],
) -> O4SkyrmePointwiseEvidence:
    """Return exact pointwise square certificates for a declared O(4) model.

    ``spatial_gradients`` has three rows ``g_i`` and four O(4) components.
    For ``D = G*G.T``, the declared densities and mass operator are

    ``e2 = tr(D)``,
    ``e4 = ((tr(D))**2 - tr(D**2))/2``, and
    ``M = 2*((1 + tr(D))*I - G.T*G)``.

    The result derives ``e4`` as a sum of squared 2-by-2 minors and derives
    ``w.T*M*w - 2*w.T*w`` as twice a sum of squared minors.  Consequently
    ``e4 >= 0`` and ``M >= 2*I`` over the reals, conditional on the declared
    formulas.  Equality of the mass bound is possible, so the coefficient
    ``2`` is sharp.
    """

    gradients = _exact_real_matrix(spatial_gradients, (3, 4), "spatial_gradients")
    tangent_column = _exact_real_matrix(tangent, (4, 1), "tangent")
    strain = sp.ImmutableDenseMatrix(gradients * gradients.T)
    quadratic = sp.expand(sp.trace(strain))
    quartic_trace = sp.expand((sp.trace(strain) ** 2 - sp.trace(strain * strain)) / 2)
    quartic_sos = sp.expand(
        sum(
            (
                gradients[i, a] * gradients[j, b]
                - gradients[i, b] * gradients[j, a]
            )
            ** 2
            for i in range(3)
            for j in range(i + 1, 3)
            for a in range(4)
            for b in range(a + 1, 4)
        )
    )
    mass = sp.ImmutableDenseMatrix(
        2 * ((1 + quadratic) * sp.eye(4) - gradients.T * gradients)
    )
    tangent_norm = sp.expand((tangent_column.T * tangent_column)[0])
    mass_form = sp.expand((tangent_column.T * mass * tangent_column)[0])
    mass_gap = sp.expand(mass_form - 2 * tangent_norm)
    mass_sos = sp.expand(
        2
        * sum(
            (
                tangent_column[a, 0] * gradients[i, b]
                - tangent_column[b, 0] * gradients[i, a]
            )
            ** 2
            for i in range(3)
            for a in range(4)
            for b in range(a + 1, 4)
        )
    )
    return O4SkyrmePointwiseEvidence(
        gradients=gradients,
        tangent=tangent_column,
        strain=strain,
        quadratic_density=quadratic,
        quartic_trace_form=quartic_trace,
        quartic_minor_sos=quartic_sos,
        static_density=sp.expand(quadratic + quartic_trace),
        mass_operator=mass,
        tangent_norm_squared=tangent_norm,
        mass_quadratic_form=mass_form,
        mass_lower_bound_gap=mass_gap,
        mass_minor_sos=mass_sos,
        quartic_identity_residual=sp.expand(quartic_trace - quartic_sos),
        mass_identity_residual=sp.expand(mass_gap - mass_sos),
    )

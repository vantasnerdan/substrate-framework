"""Leading scalar one-loop contribution to the inverse Newton coupling.

This module isolates one reusable heat-kernel result for a declared
four-dimensional real scalar. Its conventions are

    D_E = -nabla_E**2 + xi*R_E + m**2,
    Gamma_E = (1/2)*ln(det(D_E)),
    S_E,EH = -integral(sqrt(g_E)*R_E)/(16*pi*G).

In the boundaryless Gilkey--Seeley--DeWitt convention of Vassilevich
arXiv:hep-th/0306138 equations (2.2) and (4.27), the integrated coefficient
usually named a_2 is

    (4*pi)**-2 * integral(sqrt(g_E)*((1/6-xi)*R_E-m**2)).

With a sharp lower proper-time cutoff tau=Lambda**-2, its leading
UV-asymptotic quadratic contribution matches

    Delta(1/G)_leading = N*(1/6-xi)*Lambda**2/(2*pi).

The returned dimensionless coefficient is therefore the s that, when nonzero,
can be fed into induced_gravity.induced_inverse_newton_ledger for this
specific leading term. The regulator must be named explicitly. Visser
arXiv:gr-qc/0204062 equations (8)--(21) records the same cutoff-squared,
mass-logarithmic, finite, and additive-baseline structure under its stated
curvature and action conventions.

This module does not derive a field spectrum, a sine-Gordon covariant action,
a cutoff identification, a vanishing bare term, the finite-mass or logarithmic
terms, curvature-squared terms, the cosmological term, a total Newton constant,
or a stability verdict. In particular, a positive leading shift does not by
itself establish attractive gravity because the accepted ledger retains an
independent additive baseline and the omitted one-loop terms.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .exact_symbolic import exact_real as _exact_real
from .exact_symbolic import positive_exact as _positive_exact


FOUR_DIMENSIONAL_CONFORMAL_COUPLING = sp.Rational(1, 6)
SHARP_PROPER_TIME_REGULATOR = "proper_time_sharp"

_REAL_SCALAR_DETERMINANT_WEIGHT = sp.Rational(1, 2)
_FOUR_DIMENSIONAL_HEAT_KERNEL_PREFACTOR = (4 * sp.pi) ** -2
_EINSTEIN_HILBERT_INVERSE_COUPLING_FACTOR = 16 * sp.pi


@dataclass(frozen=True)
class ScalarHeatKernelA2:
    """Integrated local a_2 data before the universal (4*pi)^-2 factor."""

    non_minimal_coupling: sp.Expr
    mass_squared: sp.Expr
    conformal_coupling: sp.Expr
    curvature_weight: sp.Expr
    mass_weight: sp.Expr
    conformal_zero: bool | None


@dataclass(frozen=True)
class LeadingScalarNewtonShiftCoefficient:
    """Dimensionless coefficient of one isolated leading inverse-Newton shift."""

    field_count: sp.Expr
    non_minimal_coupling: sp.Expr
    regulator: str
    determinant_weight: sp.Expr
    heat_kernel_prefactor: sp.Expr
    inverse_newton_matching_factor: sp.Expr
    scheme_factor: sp.Expr
    curvature_weight: sp.Expr
    coefficient_per_field: sp.Expr
    coefficient: sp.Expr
    sign: int
    positive_leading_shift: bool


def _nonnegative_exact(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be provably nonnegative")
    return expression


def _positive_integer(value: Any, name: str) -> sp.Expr:
    expression = _positive_exact(value, name)
    if expression.is_integer is not True:
        raise ValueError(f"{name} must be a positive integer field count")
    return expression


def scalar_heat_kernel_a2(
    non_minimal_coupling: Any,
    mass_squared: Any = 0,
) -> ScalarHeatKernelA2:
    r"""Return the local scalar a_2 weights in the declared convention.

    Vassilevich writes a Laplace-type operator as D=-(nabla**2+E) and gives
    a_2 proportional to E+R/6. For D_E=-nabla_E**2+xi*R_E+m**2,
    substitution of E=-xi*R_E-m**2 gives the returned curvature and mass
    weights.

    The mass term does not alter the leading UV-asymptotic quadratic curvature
    weight. It does generate mass-logarithmic and finite contributions at
    later order, which are intentionally outside this function. Retaining only
    the returned leading weight therefore supplies no controlled finite-mass
    verdict when the cutoff is comparable to the mass.
    """

    xi = _exact_real(non_minimal_coupling, "non_minimal_coupling")
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    curvature_weight = sp.simplify(FOUR_DIMENSIONAL_CONFORMAL_COUPLING - xi)
    if curvature_weight.is_zero is True:
        conformal_zero: bool | None = True
    elif curvature_weight.is_zero is False:
        conformal_zero = False
    else:
        conformal_zero = None
    return ScalarHeatKernelA2(
        non_minimal_coupling=xi,
        mass_squared=mass,
        conformal_coupling=FOUR_DIMENSIONAL_CONFORMAL_COUPLING,
        curvature_weight=curvature_weight,
        mass_weight=sp.simplify(-mass),
        conformal_zero=conformal_zero,
    )


def leading_scalar_newton_shift_coefficient(
    field_count: Any,
    non_minimal_coupling: Any,
    *,
    regulator: str,
) -> LeadingScalarNewtonShiftCoefficient:
    r"""Return s=N*(1-6*xi)/(12*pi) for the named leading cutoff term.

    Only the exact sharp proper-time prescription declared by
    SHARP_PROPER_TIME_REGULATOR is implemented. The scheme factor is derived
    from the one-real-scalar determinant weight, the four-dimensional
    heat-kernel prefactor, and conversion from the Euclidean Einstein-Hilbert
    coefficient to 1/G:

        16*pi * (1/2) * (4*pi)**-2 = 1/(2*pi).

    A coupling with an undecidable relation to 1/6 is rejected because the
    returned sign is part of the API. That sign belongs only to this leading
    additive shift, not to the full inverse coupling or a stability theorem.
    """

    if regulator != SHARP_PROPER_TIME_REGULATOR:
        raise ValueError(
            f"unknown regulator {regulator!r}; pass "
            f"{SHARP_PROPER_TIME_REGULATOR!r} explicitly"
        )
    count = _positive_integer(field_count, "field_count")
    heat_kernel = scalar_heat_kernel_a2(non_minimal_coupling)
    curvature_weight = heat_kernel.curvature_weight
    scheme_factor = sp.simplify(
        _EINSTEIN_HILBERT_INVERSE_COUPLING_FACTOR
        * _REAL_SCALAR_DETERMINANT_WEIGHT
        * _FOUR_DIMENSIONAL_HEAT_KERNEL_PREFACTOR
    )
    coefficient_per_field = sp.simplify(scheme_factor * curvature_weight)
    coefficient = sp.simplify(count * coefficient_per_field)

    if curvature_weight.is_positive is True:
        sign = 1
    elif curvature_weight.is_zero is True:
        sign = 0
    elif curvature_weight.is_negative is True:
        sign = -1
    else:
        raise ValueError(
            "non_minimal_coupling must have a decidable relation to the "
            "four-dimensional conformal value 1/6"
        )

    return LeadingScalarNewtonShiftCoefficient(
        field_count=count,
        non_minimal_coupling=heat_kernel.non_minimal_coupling,
        regulator=regulator,
        determinant_weight=_REAL_SCALAR_DETERMINANT_WEIGHT,
        heat_kernel_prefactor=_FOUR_DIMENSIONAL_HEAT_KERNEL_PREFACTOR,
        inverse_newton_matching_factor=_EINSTEIN_HILBERT_INVERSE_COUPLING_FACTOR,
        scheme_factor=scheme_factor,
        curvature_weight=curvature_weight,
        coefficient_per_field=coefficient_per_field,
        coefficient=coefficient,
        sign=sign,
        positive_leading_shift=sign == 1,
    )

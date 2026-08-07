"""Induced (Sakharov) Newton coefficient from the sine-Gordon fluctuation heat-kernel.

This module discharges the single free coordinate that :mod:`induced_gravity`
leaves open: the dimensionless ``induced_coefficient`` ``s`` in the accepted
relation ``Delta(1/G) = s*hbar/(a**2*c**3)``.  It derives ``s`` by composing
the standard Gilkey--Seeley--DeWitt heat-kernel coefficient of the sine-Gordon
fluctuation operator with a declared proper-time cutoff, then feeds the result
into :func:`induced_gravity.induced_inverse_newton_ledger` together with the
accepted sine-Gordon length from :mod:`dimensional_sine_gordon`.

Route provenance: this is a Route-1 (induced / Sakharov) rung of the campaign
in issue #12, chosen after the Route-2 collective-tensor ladder (#9) closed on
a proven no-go.  In Route 1 there is no fundamental propagating collective
tensor and no microscopic graviton particle; the metric is a background field
and the Einstein--Hilbert term is *induced* by integrating out the accepted
sine-Gordon fluctuations.  The #9 spectrum gate therefore does not apply here.

Imported, cited, NOT re-derived here:

- the second Gilkey--Seeley--DeWitt coefficient ``a_2 = (1/6 - xi)*R - m**2``
  for a Laplace-type operator ``-box + xi*R + m**2`` (P. B. Gilkey,
  *Invariance Theory, the Heat Equation, and the Atiyah--Singer Index
  Theorem*; Birrell & Davies, *Quantum Fields in Curved Space*, sec. 6.3;
  Parker & Toms, *Quantum Field Theory in Curved Spacetime*);
- the Sakharov induced-gravity mechanism (A. D. Sakharov, 1967; S. L. Adler,
  Rev. Mod. Phys. 54 (1982) 729; M. Visser, Mod. Phys. Lett. A17 (2002) 977).

Derived here: the composition of that coefficient with a declared cutoff into
the dimensionless ``s``; its exact sign gate at conformal coupling; and the
induced coupling and induced Planck length written in upstream sine-Gordon
coefficients.

Regulator dependence is a DECLARED PREMISE, not a hidden fit.  The overall
O(1) factor of ``s`` depends on the cutoff scheme; the ``(1/6 - xi)``
structure, the sign flip at ``xi = 1/6``, the linearity in the field count
``N``, and the ``Lambda**2`` (i.e. ``1/ell**2``) leading scaling are
scheme-independent.  No fitted ``G``, additive baseline, or observed scale
enters this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .dimensional_sine_gordon import (
    DimensionalSineGordonCoefficients,
    dimensional_sine_gordon_coefficients,
    dimensional_sine_gordon_scales,
)
from .exact_symbolic import exact_real as _exact_real
from .exact_symbolic import positive_exact as _positive_exact
from .induced_gravity import induced_inverse_newton_ledger


CONFORMAL_COUPLING = sp.Rational(1, 6)

# Scheme-dependent overall factor applied to the scheme-independent (1/6 - xi).
# proper_time_sharp: real scalar, Gamma = (1/2)*ln det, sharp proper-time lower
# cutoff at s = 1/Lambda**2.  The s^1 heat-kernel term gives
#   1/(16 pi G) = (1/2)*(1/(4 pi)**2)*(1/6 - xi)*Lambda**2
# hence 1/G = (1/6 - xi)/(2 pi) * Lambda**2, i.e. a scheme factor of 1/(2 pi).
_REGULATOR_SCHEME_FACTOR: dict[str, sp.Expr] = {
    "proper_time_sharp": sp.Integer(1) / (2 * sp.pi),
}


@dataclass(frozen=True)
class SeeleyDeWittCurvatureCoefficient:
    """The ``a_2`` heat-kernel data of ``-box + xi*R + m**2``.

    ``curvature_coefficient`` is the coefficient of ``R`` in
    ``a_2 = (1/6 - xi)*R - m**2``; it is the scheme-independent object that
    fixes the sign of any induced Einstein--Hilbert term.
    """

    non_minimal_coupling: sp.Expr
    mass_squared: sp.Expr
    conformal_coupling: sp.Expr
    curvature_coefficient: sp.Expr
    mass_coefficient: sp.Expr
    vanishes_at_conformal_coupling: bool


@dataclass(frozen=True)
class InducedNewtonCoefficient:
    """The dimensionless induced coefficient ``s`` and its sign gate."""

    field_count: sp.Expr
    non_minimal_coupling: sp.Expr
    regulator: str
    scheme_factor: sp.Expr
    curvature_coefficient: sp.Expr
    coefficient_per_field: sp.Expr
    coefficient: sp.Expr
    conformal_coupling: sp.Expr
    sign: int
    attractive: bool


@dataclass(frozen=True)
class SineGordonInducedNewton:
    """Compose ``s`` with the accepted sine-Gordon length into a coupling.

    ``induced_inverse_newton`` is ``Delta(1/G) = s*hbar/(ell**2*c**3)`` with the
    substrate cutoff identified with the sine-Gordon length ``ell`` (a declared
    premise).  ``induced_inverse_newton_in_coefficients`` is the same quantity
    written in the model coefficients as ``s*mu*lambda**2/T**2``.
    ``induced_planck_length`` is ``ell/sqrt(s)`` when the induced coupling is
    positive, and ``None`` otherwise (no real induced Planck length exists when
    the induced Einstein--Hilbert term vanishes or has the wrong sign).
    """

    coefficient: InducedNewtonCoefficient
    sine_gordon_length: sp.Expr
    signal_speed: sp.Expr
    action_scale: sp.Expr
    induced_inverse_newton: sp.Expr
    induced_inverse_newton_in_coefficients: sp.Expr
    induced_planck_length: sp.Expr | None


def _positive_integer(value: Any, name: str) -> sp.Expr:
    expression = _positive_exact(value, name)
    if expression.is_integer is not True:
        raise ValueError(f"{name} must be a positive integer field count")
    return expression


def seeley_dewitt_curvature_coefficient(
    non_minimal_coupling: Any,
    mass_squared: Any = 0,
) -> SeeleyDeWittCurvatureCoefficient:
    r"""Return the ``a_2`` coefficient of ``-box + xi*R + m**2``.

    With ``a_2 = (1/6 - xi)*R - m**2`` (Gilkey), the coefficient of ``R`` is
    ``1/6 - xi``.  It is regulator-independent and vanishes exactly at the
    conformal value ``xi = 1/6``.  ``mass_squared`` is carried for provenance
    (for the sine-Gordon vacuum ``m**2 = U''(phi_bar) = mu*cos(phi_bar)``); it
    does not enter the leading induced ``R`` coefficient.
    """

    xi = _exact_real(non_minimal_coupling, "non_minimal_coupling")
    mass = _exact_real(mass_squared, "mass_squared")
    curvature_coefficient = sp.simplify(CONFORMAL_COUPLING - xi)
    return SeeleyDeWittCurvatureCoefficient(
        non_minimal_coupling=xi,
        mass_squared=mass,
        conformal_coupling=CONFORMAL_COUPLING,
        curvature_coefficient=curvature_coefficient,
        mass_coefficient=sp.simplify(-mass),
        vanishes_at_conformal_coupling=curvature_coefficient.is_zero is True,
    )


def induced_newton_coefficient(
    field_count: Any,
    non_minimal_coupling: Any,
    *,
    regulator: str = "proper_time_sharp",
) -> InducedNewtonCoefficient:
    r"""Return the dimensionless induced coefficient ``s`` and its sign gate.

    ``s = N * scheme_factor * (1/6 - xi)``.  For the ``proper_time_sharp``
    regulator, ``scheme_factor = 1/(2*pi)`` so ``s = N*(1 - 6*xi)/(12*pi)``.
    The sign is decided from ``(1/6 - xi)`` with ``N > 0``:

    - ``xi < 1/6``  ->  ``s > 0``: attractive induced gravity (correct sign);
    - ``xi = 1/6``  ->  ``s = 0``: the induced Einstein--Hilbert term vanishes;
    - ``xi > 1/6``  ->  ``s < 0``: wrong-sign induced gravity (a stability
      failure that redirects the field content, not GR).

    The coupling must have a decidable relation to ``1/6``; a purely symbolic
    coupling raises, matching the framework's decidability discipline.
    """

    if regulator not in _REGULATOR_SCHEME_FACTOR:
        raise ValueError(
            f"unknown regulator {regulator!r}; the scheme factor is a declared "
            f"premise and only {sorted(_REGULATOR_SCHEME_FACTOR)} is implemented"
        )
    count = _positive_integer(field_count, "field_count")
    coefficient_data = seeley_dewitt_curvature_coefficient(non_minimal_coupling)
    curvature_coefficient = coefficient_data.curvature_coefficient
    scheme_factor = _REGULATOR_SCHEME_FACTOR[regulator]
    per_field = sp.simplify(scheme_factor * curvature_coefficient)
    coefficient = sp.simplify(count * per_field)
    if curvature_coefficient.is_positive is True:
        sign = 1
    elif curvature_coefficient.is_zero is True:
        sign = 0
    elif curvature_coefficient.is_negative is True:
        sign = -1
    else:
        raise ValueError(
            "non_minimal_coupling must have a decidable relation to the "
            "conformal value 1/6"
        )
    return InducedNewtonCoefficient(
        field_count=count,
        non_minimal_coupling=coefficient_data.non_minimal_coupling,
        regulator=regulator,
        scheme_factor=scheme_factor,
        curvature_coefficient=curvature_coefficient,
        coefficient_per_field=per_field,
        coefficient=coefficient,
        conformal_coupling=CONFORMAL_COUPLING,
        sign=sign,
        attractive=sign == 1,
    )


def sine_gordon_induced_newton(
    coefficients: DimensionalSineGordonCoefficients,
    field_count: Any,
    non_minimal_coupling: Any,
    *,
    regulator: str = "proper_time_sharp",
) -> SineGordonInducedNewton:
    r"""Compose ``s`` with the accepted sine-Gordon length into a coupling.

    The substrate ultraviolet scale is identified with the sine-Gordon length
    ``ell = sqrt(T/mu)`` (a declared premise), so the cutoff is ``1/ell``.
    Feeding ``s`` and ``ell`` into the accepted induced-Newton ledger gives

    ``Delta(1/G) = s*hbar/(ell**2*c**3) = s*mu*lambda**2/T**2``

    in the model coefficients ``(lambda, T, mu) = (inertia, gradient,
    onsite)``.  When ``s > 0`` the induced Planck length is ``ell/sqrt(s)``:
    the substrate length *is* the induced gravitational length up to the pure
    number ``1/sqrt(s)``.  No observed scale selects any of these; ``G`` is an
    output in coefficients, never an input.
    """

    if not isinstance(coefficients, DimensionalSineGordonCoefficients):
        raise TypeError(
            "coefficients must be a DimensionalSineGordonCoefficients record"
        )
    validated = dimensional_sine_gordon_coefficients(
        coefficients.inertia,
        coefficients.gradient,
        coefficients.onsite,
    )
    scales = dimensional_sine_gordon_scales(validated)
    coefficient = induced_newton_coefficient(
        field_count,
        non_minimal_coupling,
        regulator=regulator,
    )
    s = coefficient.coefficient
    ell = scales.length
    speed = scales.signal_speed
    action = scales.action

    induced_direct = sp.simplify(s * action / (ell**2 * speed**3))
    induced_in_coefficients = sp.simplify(
        s * validated.onsite * validated.inertia**2 / validated.gradient**2
    )
    if sp.simplify(induced_direct - induced_in_coefficients) != 0:
        raise AssertionError(
            "induced inverse Newton mismatch between cutoff form and coefficients"
        )
    # Cross-check against the accepted ledger whenever s is nonzero (the ledger
    # requires a nonzero coefficient by construction).
    if coefficient.sign != 0:
        ledger = induced_inverse_newton_ledger(
            cutoff_length=ell,
            induced_coefficient=s,
            signal_speed=speed,
            action_scale=action,
        )
        if sp.simplify(ledger.induced_inverse_newton - induced_direct) != 0:
            raise AssertionError(
                "induced inverse Newton disagrees with induced_gravity ledger"
            )

    if coefficient.sign == 1:
        induced_planck_length: sp.Expr | None = sp.simplify(sp.sqrt(ell**2 / s))
    else:
        induced_planck_length = None

    return SineGordonInducedNewton(
        coefficient=coefficient,
        sine_gordon_length=ell,
        signal_speed=speed,
        action_scale=action,
        induced_inverse_newton=induced_direct,
        induced_inverse_newton_in_coefficients=induced_in_coefficients,
        induced_planck_length=induced_planck_length,
    )

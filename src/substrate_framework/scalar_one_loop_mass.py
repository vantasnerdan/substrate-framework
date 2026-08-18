"""Exact-mass one-loop proper-time coefficients for preregistered regulators.

This is the Route 1 (induced / Sakharov gravity) rung that follows
``covariant_sine_gordon_action``.  The landed action rung hands the
fluctuation operator

    D_E = -nabla_E**2 + xi*R_E + m**2,      m**2 = V''(phi_bg),

to the determinant.  This module evaluates the two leading proper-time
integral classes of ``Gamma_E = (1/2)*ln(det(D_E))`` with the mass retained
exactly --- no expansion of ``exp(-tau*m**2)`` inside the regulated integral
--- for three preregistered regulators:

``SHARP_PROPER_TIME_REGULATOR`` (reused from ``scalar_induced_newton``)
    A sharp lower proper-time cutoff ``tau_0 = Lambda**-2``.  Writing
    ``z = m**2/Lambda**2``, the curvature class integral is

        I_2 = integral_{tau_0}^{infty} tau**-2 * exp(-m**2*tau) dtau
            = Lambda**2 * (exp(-z) - z*E1(z)),

    with ``E1 = expint(1, .)``.  The factor ``exp(-z) - z*E1(z)`` is the
    load-bearing finite-mass correction named by the PR #13 review
    (approximately 0.1485 at z = 1): a leading UV-asymptotic ``Lambda**2``
    term alone does not establish the finite-cutoff result.

``SMOOTH_PROPER_TIME_REGULATOR``
    A smooth essential-singularity proper-time weight ``exp(-1/(Lambda**2
    * tau))``.  All power-divergent classes become convergent single
    integrals; the closed forms are modified Bessel functions:

        I_2 = 2*Lambda**2*sqrt(z)*BesselK_1(2*sqrt(z)),
        I_3 = 2*Lambda**4*z*BesselK_2(2*sqrt(z)).

``ZETA_POWER_SUBTRACTED_REGULATOR``
    The Mellin finite-part scheme in which the power divergences are
    subtracted; it requires a declared renormalization scale ``mu`` (the
    cutoff argument is rejected) and yields the pure logarithmic running

        I_2 = m**2*(ln(m**2/mu**2) + EulerGamma - 1),
        I_3 = -(m**4/2)*(ln(m**2/mu**2) + EulerGamma - 3/2).

Authority note: ``scalar_induced_newton`` and
``covariant_sine_gordon_action`` are landed conditional (unpromoted)
prior-work APIs -- the PR #14 and PR #25 harvests promoted no claims -- and
the accepted ``C-GRV-001`` supplies only the conditional dimensional and
additive-baseline ledger, leaving the coefficient, field content, and
regulator as premises.  Every public symbol in this module is likewise
conditional, unpromoted infrastructure linked to open goal #76.

Composition with that landed conditional scheme factor (one real scalar,
determinant weight 1/2, heat-kernel prefactor ``(4*pi)**-2``,
Einstein-Hilbert matching factor ``16*pi``) gives the exact-mass induced
inverse-Newton shift

    Delta(1/G) = N * coefficient_per_field(xi) * I_2(regulator),

where ``coefficient_per_field`` is taken from the landed conditional
``scalar_induced_newton.leading_scalar_newton_shift_coefficient`` API, so
the massless sharp limit reproduces that conditional API's
``s*Lambda**2 = N*(1-6*xi)*Lambda**2/(12*pi)`` exactly.  The vacuum sector
of the same mass-resummed expansion is

    Delta(rho_Gamma) = -(N/2)*(4*pi)**-2 * I_3(m**2)

per scheme.  The ``-m**2`` entry of the landed conditional
``scalar_heat_kernel_a2`` weights belongs to the *unresummed* organization,
in which the exponential is expanded and the mass survives only as that
coefficient; applying it on top of the resummed ``I_3`` double-counts the
mass.  The two organizations agree to first order through the exact
derivative identity ``d I_3/d m**2 = -I_2`` (which holds for all three
schemes and is tested), i.e. ``I_3(m**2) = I_3(0) - m**2*I_2(0) + O(m**4)``.

Declared derivation oracles (all exact, SymPy):

* sharp ``I_2``: ``d/dtau_0 I_2(tau_0) = -tau_0**-2*exp(-m**2*tau_0)`` plus
  the tail decay, so the closed form is the tail integral by uniqueness;
* sharp ``I_3``: the integration-by-parts identity
  ``I_3 = exp(-z)*Lambda**4/2 - (m**2/2)*I_2``;
* smooth closed forms: the exact differential recurrence
  ``d/dz J_p = -J_{p-1}`` of the integrals
  ``J_p(z) = integral_0^{infty} t**-p exp(-z*t - 1/t) dt`` together with the
  boundary limits ``J_2(0+) = 1`` and ``J_3(0+) = 1`` (verified against
  SymPy Bessel derivatives), plus independent high-precision quadrature;
* zeta finite parts: the cutoff-subtraction limits
  ``lim_{tau_0 -> 0} [I(tau_0) - power terms - log terms]`` evaluated
  symbolically.

This module does NOT identify the cutoff with any substrate scale, choose a
renormalization condition among the three schemes, derive a total Newton
constant (the accepted C-GRV-001 additive baseline remains independent),
evaluate the tau**-1 class (curvature-squared and m**2*R structures), or
confront any empirical comparator.  Those are later rungs of issue #76.
Scheme differences are the deliverable: at z = 1 the smooth regulator
induces about 1.88 times the sharp curvature coefficient, and the
power-subtracted scheme changes both the parametric structure and the sign
at small log; any later usable-normalization claim must cite this spread.

``ExactMassNewtonShift.curvature_weight_sign`` describes only the weight
``1/6 - xi`` and carries no total-shift verdict; ``value_sign`` is the
decidable sign of the full returned value (``None`` when symbolic inputs
leave it undecidable).  For the two cutoff schemes ``I_2 > 0`` always, so
the two agree there; for the power-subtracted scheme ``I_2`` itself changes
sign with ``m**2/mu**2`` and the two fields genuinely differ.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .exact_symbolic import exact_real as _exact_real
from .exact_symbolic import positive_exact as _positive_exact
from .scalar_induced_newton import SHARP_PROPER_TIME_REGULATOR
from .scalar_induced_newton import leading_scalar_newton_shift_coefficient  # conditional landed API

SMOOTH_PROPER_TIME_REGULATOR = "proper_time_smooth_essential"
ZETA_POWER_SUBTRACTED_REGULATOR = "zeta_power_subtracted"

KNOWN_ONE_LOOP_REGULATORS = (
    SHARP_PROPER_TIME_REGULATOR,
    SMOOTH_PROPER_TIME_REGULATOR,
    ZETA_POWER_SUBTRACTED_REGULATOR,
)

_FOUR_DIMENSIONAL_HEAT_KERNEL_PREFACTOR = (4 * sp.pi) ** -2


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


def _resolve_regulator(
    regulator: Any,
    cutoff: Any,
    renormalization_scale: Any,
) -> tuple[sp.Expr | None, sp.Expr | None]:
    """Validate the scheme contract and return (cutoff, scale) exact values."""

    if regulator not in KNOWN_ONE_LOOP_REGULATORS:
        raise ValueError(
            f"unknown regulator {regulator!r}; pass one of "
            f"{KNOWN_ONE_LOOP_REGULATORS}"
        )
    if regulator == ZETA_POWER_SUBTRACTED_REGULATOR:
        if cutoff is not None:
            raise ValueError(
                "the power-subtracted scheme has no cutoff; pass "
                "renormalization_scale only"
            )
        if renormalization_scale is None:
            raise ValueError(
                "the power-subtracted scheme requires an explicit "
                "renormalization_scale; it must not default"
            )
        return None, _positive_exact(renormalization_scale, "renormalization_scale")
    if renormalization_scale is not None:
        raise ValueError(
            "renormalization_scale belongs to the power-subtracted scheme "
            "only; a cutoff scheme must not declare one"
        )
    if cutoff is None:
        raise ValueError("this regulator requires an explicit positive cutoff")
    return _positive_exact(cutoff, "cutoff"), None


def _sharp_curvature_integral(cutoff: sp.Expr, mass_squared: sp.Expr) -> sp.Expr:
    if mass_squared.is_zero is True:
        # Massless limit: the tail integral of tau^-2 is exactly Lambda^2.
        return cutoff**2
    z = mass_squared / cutoff**2
    return sp.simplify(cutoff**2 * (sp.exp(-z) - z * sp.expint(1, z)))


def _smooth_curvature_integral(cutoff: sp.Expr, mass_squared: sp.Expr) -> sp.Expr:
    if mass_squared.is_zero is True:
        # Massless limit: 2*sqrt(z)*K_1(2*sqrt(z)) -> 1 as z -> 0+.
        return cutoff**2
    z = mass_squared / cutoff**2
    return sp.simplify(
        2 * cutoff**2 * sp.sqrt(z) * sp.besselk(1, 2 * sp.sqrt(z))
    )


def _zeta_curvature_integral(
    renormalization_scale: sp.Expr, mass_squared: sp.Expr
) -> sp.Expr:
    if mass_squared.is_zero is True:
        # Power divergences subtracted: no induced curvature term at m = 0.
        return sp.Integer(0)
    return sp.simplify(
        mass_squared
        * (sp.log(mass_squared / renormalization_scale**2) + sp.EulerGamma - 1)
    )


def _sharp_vacuum_integral(cutoff: sp.Expr, mass_squared: sp.Expr) -> sp.Expr:
    if mass_squared.is_zero is True:
        # Massless limit: the tail integral of tau^-3 is Lambda^4/2.
        return cutoff**4 / 2
    z = mass_squared / cutoff**2
    return sp.simplify(
        sp.exp(-z) * cutoff**4 / 2
        - (mass_squared / 2) * _sharp_curvature_integral(cutoff, mass_squared)
    )


def _smooth_vacuum_integral(cutoff: sp.Expr, mass_squared: sp.Expr) -> sp.Expr:
    if mass_squared.is_zero is True:
        # Massless limit: integral_0^infty t^-3 exp(-1/t) dt = 1, so Lambda^4.
        return cutoff**4
    z = mass_squared / cutoff**2
    return sp.simplify(2 * cutoff**4 * z * sp.besselk(2, 2 * sp.sqrt(z)))


def _zeta_vacuum_integral(
    renormalization_scale: sp.Expr, mass_squared: sp.Expr
) -> sp.Expr:
    if mass_squared.is_zero is True:
        return sp.Integer(0)
    return sp.simplify(
        -(mass_squared**2 / 2)
        * (
            sp.log(mass_squared / renormalization_scale**2)
            + sp.EulerGamma
            - sp.Rational(3, 2)
        )
    )


def curvature_proper_time_integral(
    regulator: Any,
    *,
    cutoff: Any = None,
    mass_squared: Any = 0,
    renormalization_scale: Any = None,
) -> sp.Expr:
    """Return the exact tau**-2-class proper-time integral for one scheme.

    This is the integral that multiplies the declared heat-kernel weights
    ``((1/6 - xi)*R_E - m**2)`` in ``Gamma_E``; with the landed conditional
    scheme factor it also multiplies the induced inverse-Newton shift.
    """

    resolved_cutoff, scale = _resolve_regulator(regulator, cutoff, renormalization_scale)
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    if regulator == SHARP_PROPER_TIME_REGULATOR:
        return _sharp_curvature_integral(resolved_cutoff, mass)
    if regulator == SMOOTH_PROPER_TIME_REGULATOR:
        return _smooth_curvature_integral(resolved_cutoff, mass)
    return _zeta_curvature_integral(scale, mass)


def vacuum_proper_time_integral(
    regulator: Any,
    *,
    cutoff: Any = None,
    mass_squared: Any = 0,
    renormalization_scale: Any = None,
) -> sp.Expr:
    """Return the exact tau**-3-class proper-time integral for one scheme.

    This is the integral that multiplies the heat-kernel weight ``1`` (the
    cosmological/vacuum sector of the one-loop action).
    """

    resolved_cutoff, scale = _resolve_regulator(regulator, cutoff, renormalization_scale)
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    if regulator == SHARP_PROPER_TIME_REGULATOR:
        return _sharp_vacuum_integral(resolved_cutoff, mass)
    if regulator == SMOOTH_PROPER_TIME_REGULATOR:
        return _smooth_vacuum_integral(resolved_cutoff, mass)
    return _zeta_vacuum_integral(scale, mass)


@dataclass(frozen=True)
class ExactMassNewtonShift:
    """Exact-mass induced inverse-Newton shift data for one scheme."""

    regulator: str
    field_count: sp.Expr
    non_minimal_coupling: sp.Expr
    mass_squared: sp.Expr
    cutoff: sp.Expr | None
    renormalization_scale: sp.Expr | None
    proper_time_value: sp.Expr
    coefficient_per_field: sp.Expr
    value: sp.Expr
    massless_leading_value: sp.Expr
    finite_mass_factor: sp.Expr | None
    curvature_weight_sign: int
    value_sign: int | None


def exact_mass_inverse_newton_shift(
    field_count: Any,
    non_minimal_coupling: Any,
    *,
    regulator: Any,
    cutoff: Any = None,
    mass_squared: Any = 0,
    renormalization_scale: Any = None,
) -> ExactMassNewtonShift:
    """Return ``Delta(1/G) = N * coefficient_per_field * I_2`` exactly.

    ``coefficient_per_field`` is the per-field scheme factor
    ``scheme_factor*(1/6 - xi)`` read from the landed conditional
    ``scalar_induced_newton.leading_scalar_newton_shift_coefficient`` API
    (its massless sharp regulator tag is used only to read that factor,
    which is regulator-independent).  The massless sharp limit reproduces
    that module's ``N*(1-6*xi)*Lambda**2/(12*pi)`` exactly.

    ``curvature_weight_sign`` is the decidable sign of ``1/6 - xi`` only.
    ``value_sign`` is the decidable sign of the full returned value, or
    ``None`` when symbolic inputs leave it undecidable; the two coincide
    for the cutoff schemes (``I_2 > 0``) but differ for the
    power-subtracted scheme whenever ``I_2 < 0``.
    """

    resolved_cutoff, scale = _resolve_regulator(regulator, cutoff, renormalization_scale)
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    count = _positive_integer(field_count, "field_count")
    xi = _exact_real(non_minimal_coupling, "non_minimal_coupling")

    landed = leading_scalar_newton_shift_coefficient(
        1, xi, regulator=SHARP_PROPER_TIME_REGULATOR
    )
    coefficient_per_field = landed.coefficient_per_field

    if regulator == SHARP_PROPER_TIME_REGULATOR:
        proper_time_value = _sharp_curvature_integral(resolved_cutoff, mass)
        massless_leading_value = sp.simplify(
            count * coefficient_per_field * resolved_cutoff**2
        )
        finite_mass_factor = sp.simplify(proper_time_value / resolved_cutoff**2)
    elif regulator == SMOOTH_PROPER_TIME_REGULATOR:
        proper_time_value = _smooth_curvature_integral(resolved_cutoff, mass)
        massless_leading_value = sp.simplify(
            count * coefficient_per_field * resolved_cutoff**2
        )
        finite_mass_factor = sp.simplify(proper_time_value / resolved_cutoff**2)
    else:
        proper_time_value = _zeta_curvature_integral(scale, mass)
        massless_leading_value = sp.Integer(0)
        finite_mass_factor = None

    value = sp.simplify(count * coefficient_per_field * proper_time_value)
    curvature_weight = sp.simplify(sp.Rational(1, 6) - xi)
    if curvature_weight.is_positive is True:
        curvature_weight_sign = 1
    elif curvature_weight.is_zero is True:
        curvature_weight_sign = 0
    elif curvature_weight.is_negative is True:
        curvature_weight_sign = -1
    else:
        raise ValueError(
            "non_minimal_coupling must have a decidable relation to the "
            "four-dimensional conformal value 1/6"
        )
    if value.is_positive is True:
        value_sign: int | None = 1
    elif value.is_zero is True:
        value_sign = 0
    elif value.is_negative is True:
        value_sign = -1
    else:
        value_sign = None

    return ExactMassNewtonShift(
        regulator=regulator,
        field_count=count,
        non_minimal_coupling=xi,
        mass_squared=mass,
        cutoff=resolved_cutoff,
        renormalization_scale=scale,
        proper_time_value=proper_time_value,
        coefficient_per_field=coefficient_per_field,
        value=value,
        massless_leading_value=massless_leading_value,
        finite_mass_factor=finite_mass_factor,
        curvature_weight_sign=curvature_weight_sign,
        value_sign=value_sign,
    )


@dataclass(frozen=True)
class ExactMassVacuumShift:
    """Exact-mass one-loop vacuum-sector data for one scheme."""

    regulator: str
    field_count: sp.Expr
    mass_squared: sp.Expr
    cutoff: sp.Expr | None
    renormalization_scale: sp.Expr | None
    tau_minus_two_value: sp.Expr
    tau_minus_three_value: sp.Expr
    value: sp.Expr
    # tau_minus_two_value is exposed for reference and for the derivative
    # identity d I_3 / d m**2 = -I_2; the returned value composes only the
    # mass-resummed tau**-3 class.


def exact_mass_vacuum_density_shift(
    field_count: Any,
    *,
    regulator: Any,
    cutoff: Any = None,
    mass_squared: Any = 0,
    renormalization_scale: Any = None,
) -> ExactMassVacuumShift:
    """Return ``Delta(rho) = -(N/2)*(4*pi)**-2*I_3(m**2)`` exactly.

    This is the vacuum (cosmological-sector) coefficient of the one-loop
    action density in the mass-resummed organization: the trace integrand is
    ``exp(-tau*m**2)*[tau**-2 + tau**-1*(1/6-xi)*R_E + ...]``, so the vacuum
    sector is the tau**-3 class with the exponential retained.  The
    ``-m**2`` heat-kernel weight of the landed conditional
    ``scalar_heat_kernel_a2`` is the first-order remnant of that same
    exponential in the unresummed organization and must not be added again;
    the exact bridge is the derivative identity ``d I_3/d m**2 = -I_2``.
    The sector is exhibited rather than omitted whenever the curvature-sector
    shift is quoted.
    """

    resolved_cutoff, scale = _resolve_regulator(regulator, cutoff, renormalization_scale)
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    count = _positive_integer(field_count, "field_count")

    tau_minus_two = curvature_proper_time_integral(
        regulator,
        cutoff=resolved_cutoff,
        mass_squared=mass,
        renormalization_scale=scale,
    )
    tau_minus_three = vacuum_proper_time_integral(
        regulator,
        cutoff=resolved_cutoff,
        mass_squared=mass,
        renormalization_scale=scale,
    )
    value = sp.simplify(
        -count
        * sp.Rational(1, 2)
        * _FOUR_DIMENSIONAL_HEAT_KERNEL_PREFACTOR
        * tau_minus_three
    )
    return ExactMassVacuumShift(
        regulator=regulator,
        field_count=count,
        mass_squared=mass,
        cutoff=resolved_cutoff,
        renormalization_scale=scale,
        tau_minus_two_value=tau_minus_two,
        tau_minus_three_value=tau_minus_three,
        value=value,
    )


@dataclass(frozen=True)
class RegulatorSchemeLedger:
    """Exact scheme spread of the curvature-class integral at common scales."""

    cutoff: sp.Expr
    mass_squared: sp.Expr
    renormalization_scale: sp.Expr
    sharp_value: sp.Expr
    smooth_value: sp.Expr
    zeta_value: sp.Expr
    sharp_over_smooth: sp.Expr
    sharp_minus_zeta: sp.Expr


def regulator_scheme_ledger(
    cutoff: Any,
    mass_squared: Any,
    renormalization_scale: Any,
) -> RegulatorSchemeLedger:
    """Return the three curvature-class values and their exact contrasts.

    The sharp and smooth schemes are evaluated at the declared cutoff; the
    power-subtracted scheme at the declared renormalization scale.  The
    contrasts are exact expressions, not numeric fits, so a later
    renormalization-condition rung can cite them as reviewable provenance.
    """

    resolved_cutoff = _positive_exact(cutoff, "cutoff")
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    scale = _positive_exact(renormalization_scale, "renormalization_scale")

    sharp_value = _sharp_curvature_integral(resolved_cutoff, mass)
    smooth_value = _smooth_curvature_integral(resolved_cutoff, mass)
    zeta_value = _zeta_curvature_integral(scale, mass)
    return RegulatorSchemeLedger(
        cutoff=resolved_cutoff,
        mass_squared=mass,
        renormalization_scale=scale,
        sharp_value=sharp_value,
        smooth_value=smooth_value,
        zeta_value=zeta_value,
        sharp_over_smooth=sp.simplify(sharp_value / smooth_value),
        sharp_minus_zeta=sp.simplify(sharp_value - zeta_value),
    )

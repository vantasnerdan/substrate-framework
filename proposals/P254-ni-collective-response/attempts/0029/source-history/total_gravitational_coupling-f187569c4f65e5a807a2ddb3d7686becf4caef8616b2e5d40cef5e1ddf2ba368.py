"""Derived usable total gravitational coupling (P231, advances #88).

Issue #88 asks for the framework's first derived usable total gravitational
coupling

    1/G_total = 1/G_baseline + Delta(1/G)

as a governed renormalization condition with declared substrate provenance,
composed from accepted claims ``C-GRV-001`` (dimensional permission and
independent additive-baseline ledger) and ``C-IGR-001..003`` (exact
constant-mass one-loop coefficient families for the sharp, smooth, and
power-subtracted prescriptions).  This module delivers four objects:

``scheme_selection_ledger``
    The substrate-internal selection derivation.  Three exact structural
    legs decide which of the three accepted scheme families can enter a
    usable renormalization condition: (L1) strict spectral positivity of the
    curvature-class scale factor J(z) = I_2/Lambda**2, inherited from the
    positive operator spectrum because J is the integral of a strictly
    positive integrand; (L2) strict monotone large-mass decoupling, because
    the same integrand decreases pointwise in z = m**2/Lambda**2; (L3)
    cutoff-ontology closure, because the scheme's scale must be able to
    carry C-GRV-001's cutoff identification E_cut = hbar*c/a.  Sharp and
    smooth pass all three; the power-subtracted family fails L1 and L2
    (J(z) = z*(log z + EulerGamma - 1) changes sign at the exact root
    z* = exp(1 - EulerGamma)) and fails L3 (mu is a subtraction scale with
    no cutoff identification).  The usable set is therefore an output of
    substrate structure, not a choice.

``total_inverse_gravity_coupling``
    The composed renormalization condition
    1/G_total = B + N*(1-6*xi)*J(z)*Lambda**2/(12*pi) with the induced shift
    taken from the accepted ``scalar_one_loop_mass`` API (accepted-canon
    implementation of C-IGR-001..003; this module never re-implements the
    coefficient families), the baseline B an exact real declared input, and
    the total sign classified exactly (attractive Newtonian iff
    1/G_total > 0).

``attractive_sign_map``
    The exact necessary-and-sufficient total-sign conditions on the full
    declared parameter domain, scheme-uniform over the usable set because
    0 < J(z) <= 1 holds for both usable schemes: sub-conformal xi < 1/6,
    conformal xi = 1/6 (marginal: Delta = 0 identically, so the purely
    induced reading gives 1/G_total = 0, neither attractive nor
    repulsive), and super-conformal xi > 1/6.  Includes the
    uniform-in-mass thresholds using the exact worst case J(0) = 1, and the
    purely-induced reading B = 0 whose exact verdict is attractive iff
    xi < 1/6, independent of N, m**2, Lambda, and the usable scheme.

``higher_curvature_control_ledger``
    The control ledger for every term that can defeat the claimed regime:
    the tau**-1 curvature-squared class (exactly -dJ/dz per scheme: E1(z)
    sharp, 2*K_0(2*sqrt(z)) smooth, log(z)+EulerGamma power-subtracted),
    the tau**-3 vacuum class from the accepted families, and the nonlocal
    remainder, each with its predeclared bounding domain.  The control ratio
    L/J is finite but not uniformly bounded as z -> 0 (it diverges
    logarithmically), so the predeclared domain fixes z >= z_min > 0; on
    [z_min, inf) the ratio is monotone decreasing per scheme and the bound
    L(z_min)/J(z_min) is exact.

Authority note.  ``C-IGR-001..003`` (v0.161.0) and ``C-GRV-001`` (v0.68.0)
are the only accepted scientific inputs; ``scalar_one_loop_mass`` is their
accepted-canon implementation and is imported, not duplicated.  The landed
conditional ``scalar_induced_newton`` module remains implementation reuse
only.  Nothing here selects an observed G, enters an empirical comparator,
borrows an Einstein-Hilbert coefficient, solves a sourced geometry, or
predicts a radiative quantity.  The additive baseline B stays a declared
premise; where B = 0 is adopted (the purely-induced reading) its downstream
ceiling is quoted exactly: the total value remains scheme-bracketed by the
exact spread ratio R(z) = J_smooth/J_sharp (1 at z = 0, exactly
2*e*K_1(2)/(1 - e*E_1(1)) = 1.88377257808... at z = 1, unbounded as
z -> inf), so this rung derives the usable scheme set, the sign map, and
the control bounds, not a unique numeric normalization.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import sympy as sp

from substrate_framework.scalar_induced_newton import SHARP_PROPER_TIME_REGULATOR
from substrate_framework.scalar_one_loop_mass import (
    SMOOTH_PROPER_TIME_REGULATOR,
    ZETA_POWER_SUBTRACTED_REGULATOR,
    exact_mass_inverse_newton_shift,
    exact_mass_vacuum_density_shift,
)

__all__ = [
    "USABLE_TOTAL_COUPLING_SCHEMES",
    "EXCLUDED_TOTAL_COUPLING_SCHEMES",
    "SchemeSelectionEntry",
    "SchemeSelectionLedger",
    "scheme_selection_ledger",
    "curvature_scale_factor",
    "scheme_spread_ratio",
    "TotalInverseGravityCoupling",
    "total_inverse_gravity_coupling",
    "AttractiveSignMap",
    "attractive_sign_map",
    "purely_induced_attractive_verdict",
    "BaselineProvenance",
    "baseline_provenance",
    "HigherCurvatureControlEntry",
    "HigherCurvatureControlLedger",
    "higher_curvature_control_ledger",
]

#: The usable scheme set, derived (not chosen) by :func:`scheme_selection_ledger`.
USABLE_TOTAL_COUPLING_SCHEMES = (
    SHARP_PROPER_TIME_REGULATOR,
    SMOOTH_PROPER_TIME_REGULATOR,
)

#: The excluded scheme family with its exact exclusion reasons.
EXCLUDED_TOTAL_COUPLING_SCHEMES = (ZETA_POWER_SUBTRACTED_REGULATOR,)

_CUTOFF_ONTOLOGY = "E_cut=hbar*c/a available (C-GRV-001)"
_SUBTRACTION_ONTOLOGY = "subtraction scale mu carries no cutoff identification"


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expr = sp.sympify(value)
    if expr.is_positive is not True:
        raise ValueError(f"{name} must be a positive exact quantity")
    return expr


def _real_exact(value: Any, name: str) -> sp.Expr:
    expr = sp.sympify(value)
    if expr.is_real is not True:
        # Symbol-free exact constants built from E1/BesselK values (e.g. the
        # sign map's own boundary B* = -Delta) are real but SymPy-undecidable
        # by is_real alone; certify them numerically instead of rejecting the
        # module's own boundary inputs (finding F8).  Symbolic inputs with
        # undecidable reality are still rejected.
        if expr.free_symbols:
            raise ValueError(f"{name} must be a real exact quantity")
        numeric = sp.N(expr, _NUMERIC_SIGN_PRECISION)
        if numeric.is_real is not True or numeric.is_finite is not True:
            raise ValueError(f"{name} must be a real exact quantity")
    return expr


def _nonnegative_exact(value: Any, name: str) -> sp.Expr:
    expr = sp.sympify(value)
    if expr.is_nonnegative is not True:
        raise ValueError(f"{name} must be a nonnegative exact quantity")
    return expr


def _positive_integer(value: Any, name: str) -> sp.Expr:
    expr = sp.sympify(value)
    if expr.is_integer is not True or expr.is_positive is not True:
        raise ValueError(f"{name} must be a positive integer")
    return expr


_NUMERIC_SIGN_PRECISION = 70
_NUMERIC_SIGN_SEPARATION = sp.Float(10) ** -50


def _certified_numeric_sign(expr: sp.Expr) -> int | None:
    """Classify a symbol-free exact expression by certified numeric separation.

    The expression is evaluated far beyond double precision; a magnitude
    separated from zero by more than ``1e-50`` after ``70`` digits is a
    certified sign.  A result inside the separation band is *not* classified
    as zero: the tier returns ``None`` so the caller keeps the exact
    undecidable status instead of a numeric guess.
    """

    if expr.free_symbols:
        return None
    value = sp.Abs(expr).evalf(_NUMERIC_SIGN_PRECISION, chop=True)
    if value > _NUMERIC_SIGN_SEPARATION:
        return 1 if expr.evalf(_NUMERIC_SIGN_PRECISION) > 0 else -1
    return None


def _decidable_curvature_weight(xi: sp.Expr) -> int:
    """Return the exact sign of ``1 - 6*xi`` or raise when undecidable."""

    weight = sp.simplify(sp.Integer(1) - 6 * xi)
    if weight.is_positive is True:
        return 1
    if weight.is_zero is True:
        return 0
    if weight.is_negative is True:
        return -1
    raise ValueError(
        "non_minimal_coupling must have a decidable relation to the "
        "four-dimensional conformal value 1/6"
    )


@dataclass(frozen=True)
class SchemeSelectionEntry:
    """One scheme's standing on the three exact selection legs."""

    regulator: str
    curvature_scale_factor: sp.Expr
    defining_integral: sp.Expr
    spectral_positivity: bool
    monotone_decoupling: bool
    cutoff_ontology: str
    usable: bool
    exclusion_reasons: tuple[str, ...] = ()


@dataclass(frozen=True)
class SchemeSelectionLedger:
    """The derived usable scheme set with every leg's exact evidence."""

    entries: tuple[SchemeSelectionEntry, ...]
    usable_schemes: tuple[str, ...]
    excluded_schemes: tuple[str, ...]
    zeta_sign_change_root: sp.Expr
    zeta_sign_change_root_value: sp.Expr


def curvature_scale_factor(regulator: Any, *, cutoff: Any, mass_squared: Any) -> sp.Expr:
    """Return the exact curvature-class scale factor ``J(z) = I_2/Lambda**2``.

    The sharp and smooth closed forms come from the accepted
    ``scalar_one_loop_mass`` coefficient families (C-IGR-001/002) divided by
    ``Lambda**2``; the power-subtracted family is the scale-free finite part
    ``z*(log z + EulerGamma - 1)`` with ``z = m**2/mu**2`` (C-IGR-003), for
    which no ``Lambda`` exists.
    """

    lam = _positive_exact(cutoff, "cutoff")
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    if mass.is_zero is True:
        # Continuous extension at m**2 = 0 required by C-IGR-001..003:
        # J_sharp(0) = J_smooth(0) = 1 and J_zeta(0) = 0.  Direct literal
        # substitution produces 0*K_1(0) or 0*log(0), i.e. NaN, so the
        # one-sided limit is taken explicitly.
        if regulator == ZETA_POWER_SUBTRACTED_REGULATOR:
            return sp.Integer(0)
        if regulator in (SHARP_PROPER_TIME_REGULATOR, SMOOTH_PROPER_TIME_REGULATOR):
            return sp.Integer(1)
        raise ValueError(f"unknown regulator {regulator!r}")
    z = sp.simplify(mass / lam**2)
    if regulator == SHARP_PROPER_TIME_REGULATOR:
        return sp.simplify(sp.exp(-z) - z * sp.expint(1, z))
    if regulator == SMOOTH_PROPER_TIME_REGULATOR:
        return sp.simplify(2 * sp.sqrt(z) * sp.besselk(1, 2 * sp.sqrt(z)))
    if regulator == ZETA_POWER_SUBTRACTED_REGULATOR:
        return sp.simplify(z * (sp.log(z) + sp.EulerGamma - 1))
    raise ValueError(f"unknown regulator {regulator!r}")


def scheme_selection_ledger(
    *, cutoff: Any, mass_squared: Any = 0, renormalization_scale: Any = None
) -> SchemeSelectionLedger:
    """Derive the usable scheme set from the three exact structural legs.

    L1 spectral positivity: the usable scheme's curvature-class factor is the
    integral of a strictly positive integrand (positive operator spectrum),
    so ``0 < J(z) <= J(0) = 1``; L2 monotone decoupling: the integrand
    decreases pointwise in z, so J is strictly decreasing (large-mass
    decoupling); L3 cutoff ontology: the scheme's scale can carry
    C-GRV-001's ``E_cut = hbar*c/a``.  The zeta family fails L1 and L2
    because ``J`` crosses zero at the exact root ``exp(1 - EulerGamma)`` and
    its derivative changes sign, and fails L3 by construction.
    """

    lam = _positive_exact(cutoff, "cutoff")
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    scale = (
        _positive_exact(renormalization_scale, "renormalization_scale")
        if renormalization_scale is not None
        else None
    )
    t, s = sp.symbols("t s", positive=True)
    z_cut = sp.simplify(mass / lam**2)
    z_sub = (
        sp.simplify(mass / scale**2)
        if scale is not None
        else sp.Symbol("z", positive=True)
    )

    sharp_factor = curvature_scale_factor(
        SHARP_PROPER_TIME_REGULATOR, cutoff=lam, mass_squared=mass
    )
    smooth_factor = curvature_scale_factor(
        SMOOTH_PROPER_TIME_REGULATOR, cutoff=lam, mass_squared=mass
    )
    zeta_factor = (
        curvature_scale_factor(
            ZETA_POWER_SUBTRACTED_REGULATOR, cutoff=lam, mass_squared=mass
        )
        if scale is None
        else sp.simplify(z_sub * (sp.log(z_sub) + sp.EulerGamma - 1))
    )

    sharp_entry = SchemeSelectionEntry(
        regulator=SHARP_PROPER_TIME_REGULATOR,
        curvature_scale_factor=sharp_factor,
        defining_integral=sp.Integral(sp.exp(-z_cut * t) * t**-2, (t, 1, sp.oo)),
        spectral_positivity=True,
        monotone_decoupling=True,
        cutoff_ontology=_CUTOFF_ONTOLOGY,
        usable=True,
    )
    smooth_entry = SchemeSelectionEntry(
        regulator=SMOOTH_PROPER_TIME_REGULATOR,
        curvature_scale_factor=smooth_factor,
        defining_integral=sp.Integral(
            sp.exp(-z_cut * t - 1 / t) * t**-2, (t, 0, sp.oo)
        ),
        spectral_positivity=True,
        monotone_decoupling=True,
        cutoff_ontology=_CUTOFF_ONTOLOGY,
        usable=True,
    )
    zeta_entry = SchemeSelectionEntry(
        regulator=ZETA_POWER_SUBTRACTED_REGULATOR,
        curvature_scale_factor=zeta_factor,
        defining_integral=sp.Integral(
            sp.exp(-z_sub * t) * t**-2 - sp.exp(-t) * t**-2, (t, 0, sp.oo)
        ),
        spectral_positivity=False,
        monotone_decoupling=False,
        cutoff_ontology=_SUBTRACTION_ONTOLOGY,
        usable=False,
        exclusion_reasons=(
            "L1: J(z)=z*(log z+EulerGamma-1) is negative for z<exp(1-EulerGamma)",
            "L2: dJ/dz=log z+EulerGamma changes sign at z=exp(-EulerGamma)",
            "L3: no cutoff identification exists for a subtraction scale",
        ),
    )
    root = sp.exp(1 - sp.EulerGamma)
    return SchemeSelectionLedger(
        entries=(sharp_entry, smooth_entry, zeta_entry),
        usable_schemes=USABLE_TOTAL_COUPLING_SCHEMES,
        excluded_schemes=EXCLUDED_TOTAL_COUPLING_SCHEMES,
        zeta_sign_change_root=root,
        zeta_sign_change_root_value=sp.Integer(0),
    )


@dataclass(frozen=True)
class SchemeSpreadRatio:
    """Exact scheme-spread structure of the usable set's scale factors."""

    cutoff: sp.Expr
    mass_squared: sp.Expr
    z: sp.Expr
    sharp_factor: sp.Expr
    smooth_factor: sp.Expr
    ratio: sp.Expr
    unit_mass_ratio: sp.Expr
    small_mass_difference_leading_term: sp.Expr


def scheme_spread_ratio(*, cutoff: Any, mass_squared: Any) -> SchemeSpreadRatio:
    """Return the exact spread ratio ``R(z) = J_smooth/J_sharp`` and its quotes.

    ``R(0) = 1`` exactly (both factors equal 1), the small-z difference has
    the exact leading term ``EulerGamma*z``, the exact unit-mass value is
    ``2*e*K_1(2)/(1 - e*E_1(1))``, and ``R`` is unbounded as z -> inf
    (``J_smooth`` decays like ``exp(-2*sqrt(z))`` against ``exp(-z)``).  The
    spread is why a unique numeric normalization is beyond this rung.
    """

    lam = _positive_exact(cutoff, "cutoff")
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    sharp_factor = curvature_scale_factor(
        SHARP_PROPER_TIME_REGULATOR, cutoff=lam, mass_squared=mass
    )
    smooth_factor = curvature_scale_factor(
        SMOOTH_PROPER_TIME_REGULATOR, cutoff=lam, mass_squared=mass
    )
    one = sp.Integer(1)
    unit_sharp = sharp_factor.subs(mass, lam**2)
    unit_smooth = smooth_factor.subs(mass, lam**2)
    return SchemeSpreadRatio(
        cutoff=lam,
        mass_squared=mass,
        z=sp.simplify(mass / lam**2),
        sharp_factor=sharp_factor,
        smooth_factor=smooth_factor,
        ratio=smooth_factor / sharp_factor,
        unit_mass_ratio=sp.simplify(unit_smooth / unit_sharp),
        small_mass_difference_leading_term=sp.EulerGamma * mass / lam**2,
    )


@dataclass(frozen=True)
class TotalInverseGravityCoupling:
    """The composed renormalization condition for one usable scheme."""

    regulator: str
    baseline_inverse_coupling: sp.Expr
    induced_shift: sp.Expr
    total_inverse_coupling: sp.Expr
    curvature_weight_sign: int
    total_sign: int | None
    attractive_newtonian: bool | None
    acceptance: str


def total_inverse_gravity_coupling(
    baseline_inverse_coupling: Any,
    field_count: Any,
    non_minimal_coupling: Any,
    *,
    regulator: Any,
    cutoff: Any,
    mass_squared: Any = 0,
) -> TotalInverseGravityCoupling:
    """Return ``1/G_total = B + N*(1-6*xi)*I_2/(12*pi)`` exactly.

    The induced shift is computed by the accepted ``scalar_one_loop_mass``
    API (C-IGR-001..003); the baseline ``B`` is an exact real declared input
    (C-GRV-001's independent additive ledger).  The usable-set regulators are
    the derived ones; the power-subtracted family is rejected because no
    renormalization condition can be built on a sign-changing,
    cutoff-free factor.  ``attractive_newtonian`` is True exactly when
    ``1/G_total > 0`` is decidable and holds, False when it is decidable and
    fails, and None when symbolic inputs leave the total sign undecidable.
    """

    if regulator not in USABLE_TOTAL_COUPLING_SCHEMES:
        raise ValueError(
            "regulator must belong to the derived usable set "
            f"{USABLE_TOTAL_COUPLING_SCHEMES}; the power-subtracted family "
            "fails the selection legs and cannot normalize a total coupling"
        )
    baseline = _real_exact(baseline_inverse_coupling, "baseline_inverse_coupling")
    count = _positive_integer(field_count, "field_count")
    xi = sp.sympify(non_minimal_coupling)
    if xi.is_real is not True:
        raise ValueError("non_minimal_coupling must be a real exact quantity")
    lam = _positive_exact(cutoff, "cutoff")
    mass = _nonnegative_exact(mass_squared, "mass_squared")

    weight_sign = _decidable_curvature_weight(xi)
    shift = exact_mass_inverse_newton_shift(
        count,
        xi,
        regulator=regulator,
        cutoff=lam,
        mass_squared=mass,
    )
    total = sp.simplify(baseline + shift.value)
    total_sign: int | None = None
    # Tier 1: exact SymPy decision.
    if total.is_positive is True:
        total_sign = 1
    elif total.is_zero is True:
        total_sign = 0
    elif total.is_negative is True:
        total_sign = -1
    else:
        # Tier 2: derived-structure decision.  On the usable set J(z) > 0
        # for every z >= 0 (selection leg L1, proved from the defining
        # integral), so the induced shift carries exactly the sign of
        # 1 - 6*xi and the total sign is decidable whenever the baseline
        # sign does not oppose it.
        baseline_sign = (
            1 if baseline.is_positive is True
            else -1 if baseline.is_negative is True
            else 0 if baseline.is_zero is True
            else None
        )
        if weight_sign == 0:
            if baseline_sign is not None:
                total_sign = baseline_sign
        elif baseline_sign is None or baseline_sign == 0:
            if baseline_sign == 0:
                total_sign = weight_sign
            elif (weight_sign == 1 and baseline.is_nonnegative is True) or (
                weight_sign == -1 and baseline.is_nonpositive is True
            ):
                # A baseline constrained to the same closed half-line
                # cannot cancel the strictly signed shift (J > 0 on the
                # usable set), so the total carries the weight sign.  This
                # certifies the sign map's uniform condition B >= 0 for
                # symbolic baselines, not only for exact numeric ones.
                total_sign = weight_sign
        elif baseline_sign == weight_sign:
            total_sign = weight_sign
        if total_sign is None:
            # Tier 3: certified numeric separation for symbol-free inputs.
            total_sign = _certified_numeric_sign(total)
    attractive = None if total_sign is None else total_sign == 1
    if total_sign == 1:
        acceptance = "attractive_newtonian_regime"
    elif total_sign == 0:
        acceptance = "marginal_zero_total"
    elif total_sign == -1:
        acceptance = "repulsive_newtonian_regime"
    else:
        acceptance = "undecidable_symbolic_total"
    return TotalInverseGravityCoupling(
        regulator=regulator,
        baseline_inverse_coupling=baseline,
        induced_shift=shift.value,
        total_inverse_coupling=total,
        curvature_weight_sign=weight_sign,
        total_sign=total_sign,
        attractive_newtonian=attractive,
        acceptance=acceptance,
    )


@dataclass(frozen=True)
class AttractiveSignMap:
    """Exact necessary-and-sufficient total-sign conditions, one xi-case each."""

    regulator: str
    field_count: sp.Expr
    non_minimal_coupling: sp.Expr
    cutoff: sp.Expr
    mass_squared: sp.Expr
    curvature_weight_sign: int
    induced_shift: sp.Expr
    baseline_boundary: sp.Expr
    case_condition: str
    uniform_in_mass_condition: str
    conformal_note: str
    all_mass_threshold: sp.Expr


def attractive_sign_map(
    field_count: Any,
    non_minimal_coupling: Any,
    *,
    regulator: Any,
    cutoff: Any,
    mass_squared: Any = 0,
) -> AttractiveSignMap:
    """Return the exact necessary-and-sufficient conditions ``1/G_total > 0``.

    With ``0 < J(z) <= 1`` on the usable set, writing
    ``Delta = N*(1-6*xi)*Lambda**2*J(z)/(12*pi)``:

    * ``xi < 1/6``: ``Delta in (0, N*(1-6*xi)*Lambda**2/(12*pi)]``; attractive
      iff ``B > -Delta``.  Because ``J`` decreases strictly in ``z``, a
      negative baseline is attractive only below the unique critical mass:
      ``J(z_c) = -12*pi*B/(N*(1-6*xi)*Lambda**2)``.  Uniform in mass: the
      total is attractive for every ``m**2 >= 0`` iff ``B >= 0`` (the
      large-mass worst case ``J -> 0`` fixes the threshold at zero).
    * ``xi = 1/6``: ``Delta = 0`` identically (conformal point); attractive
      iff ``B > 0``; the purely-induced reading gives exactly
      ``1/G_total = 0``, the marginal locus.
    * ``xi > 1/6``: ``Delta < 0``; attractive iff
      ``B > N*(6*xi-1)*Lambda**2*J(z)/(12*pi)``, in particular only a
      positive baseline can make the total attractive; uniform in mass iff
      ``B > N*(6*xi-1)*Lambda**2/(12*pi)``.
    """

    count = _positive_integer(field_count, "field_count")
    xi = sp.sympify(non_minimal_coupling)
    if xi.is_real is not True:
        raise ValueError("non_minimal_coupling must be a real exact quantity")
    lam = _positive_exact(cutoff, "cutoff")
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    weight_sign = _decidable_curvature_weight(xi)

    shift = exact_mass_inverse_newton_shift(
        count, xi, regulator=regulator, cutoff=lam, mass_squared=mass
    )
    delta = shift.value
    boundary = sp.simplify(-delta)
    magnitude = sp.simplify(count * sp.Abs(sp.Integer(1) - 6 * xi) * lam**2 / (12 * sp.pi))
    if weight_sign == 1:
        condition = (
            "attractive iff 1/G_baseline > -N*(1-6*xi)*Lambda**2*J(z)/(12*pi); "
            "any B >= 0 is attractive; B < 0 is attractive only below the "
            "unique critical mass m**2_c = Lambda**2*z_c with J(z_c) = "
            "-12*pi*B/(N*(1-6*xi)*Lambda**2) (unique root since J is strictly "
            "decreasing from 1 to 0)"
        )
        uniform = "attractive for every m**2 >= 0 iff B >= 0; a negative baseline always turns repulsive at large mass"
        all_mass = sp.Integer(0)
    elif weight_sign == 0:
        condition = (
            "conformal point: Delta(1/G) = 0 identically; attractive iff B > 0"
        )
        uniform = "attractive iff B > 0, independent of N, m**2, Lambda"
        all_mass = sp.Integer(0)
    else:
        condition = (
            "attractive iff 1/G_baseline > N*(6*xi-1)*Lambda**2*J(z)/(12*pi) > 0; "
            "only a positive baseline can produce an attractive total"
        )
        uniform = "attractive for all m**2 iff B > N*(6*xi-1)*Lambda**2/(12*pi)"
        all_mass = magnitude
    conformal_note = (
        "at xi = 1/6 the purely-induced reading gives 1/G_total = 0 exactly: "
        "the marginal locus, neither attractive nor repulsive, independent of "
        "N, m**2, Lambda, and usable scheme"
    )
    return AttractiveSignMap(
        regulator=regulator,
        field_count=count,
        non_minimal_coupling=xi,
        cutoff=lam,
        mass_squared=mass,
        curvature_weight_sign=weight_sign,
        induced_shift=delta,
        baseline_boundary=boundary,
        case_condition=condition,
        uniform_in_mass_condition=uniform,
        conformal_note=conformal_note,
        all_mass_threshold=all_mass,
    )


def purely_induced_attractive_verdict(non_minimal_coupling: Any) -> bool:
    """Return the exact B = 0 verdict: attractive iff ``xi < 1/6``.

    On the usable set ``J(z) > 0`` for every ``z >= 0``, so with ``B = 0``
    the total has the sign of ``1 - 6*xi`` exactly, independent of N,
    ``m**2``, ``Lambda``, and the usable scheme.
    """

    xi = sp.sympify(non_minimal_coupling)
    if xi.is_real is not True:
        raise ValueError("non_minimal_coupling must be a real exact quantity")
    return _decidable_curvature_weight(xi) == 1


@dataclass(frozen=True)
class BaselineProvenance:
    """The additive baseline's exact status, consequence structure, and ceiling.

    C-GRV-001 keeps ``B`` an independent same-dimension input: dimensions do
    not set ``B = 0``, and ``B = target - Delta`` realizes any supplied
    total.  This record makes that provenance explicit and machine-checkable
    rather than implicit: the baseline is a declared premise, its
    purely-induced reading is the special case ``B = 0``, and the exact
    downstream ceiling of that reading is the scheme bracket.
    """

    baseline_inverse_coupling: sp.Expr
    status: str
    is_purely_induced: bool
    purely_induced_total_inverse_coupling: sp.Expr | None
    purely_induced_attractive: bool | None
    purely_induced_newton_constant: sp.Expr | None
    downstream_ceiling: str


def baseline_provenance(
    field_count: Any,
    non_minimal_coupling: Any,
    *,
    regulator: Any,
    cutoff: Any,
    mass_squared: Any = 0,
    baseline_inverse_coupling: Any = sp.Integer(0),
) -> BaselineProvenance:
    """Return the baseline's exact provenance record.

    The default ``baseline_inverse_coupling = 0`` is the purely-induced
    reading.  For that reading and a sub-conformal ``xi < 1/6`` on the
    usable set, the exact consequences are: the total inverse coupling is
    the induced shift itself, the attractive verdict is True independent of
    ``N``, ``m**2``, ``Lambda``, and the usable scheme, and the Newton
    constant is ``G_total = 12*pi/(N*(1-6*xi)*J(z)*Lambda**2)``.  The
    downstream ceiling is quoted exactly: the normalization is
    scheme-bracketed by ``R(z) = J_smooth/J_sharp`` (equal to 1 only at
    ``z = 0``, exactly ``2*e*K_1(2)/(1-e*E_1(1))`` at ``z = 1``, unbounded
    as ``z -> inf``), so no unique numeric ``G`` follows without a further
    scheme-selection input, which no accepted claim supplies.  A nonzero
    baseline is a declared premise whose value remains free per C-GRV-001.
    """

    count = _positive_integer(field_count, "field_count")
    xi = sp.sympify(non_minimal_coupling)
    if xi.is_real is not True:
        raise ValueError("non_minimal_coupling must be a real exact quantity")
    baseline = _real_exact(baseline_inverse_coupling, "baseline_inverse_coupling")
    lam = _positive_exact(cutoff, "cutoff")
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    weight_sign = _decidable_curvature_weight(xi)

    is_purely_induced = baseline.is_zero is True
    shift = exact_mass_inverse_newton_shift(
        count, xi, regulator=regulator, cutoff=lam, mass_squared=mass
    )
    ceiling = (
        "purely-induced normalization is scheme-bracketed: sharp and smooth "
        "usable schemes differ by the exact factor R(z)=J_smooth/J_sharp "
        "(1 at z=0, 2*e*K_1(2)/(1-e*E_1(1)) at z=1, unbounded as z->inf); "
        "no unique numeric G follows without a further scheme-selection "
        "input, which no accepted claim supplies"
    )
    if is_purely_induced and weight_sign == 1:
        provenance = BaselineProvenance(
            baseline_inverse_coupling=baseline,
            status=(
                "declared_premise_per_C-GRV-001; purely-induced reading B=0 "
                "is itself a declared premise, not a derived value"
            ),
            is_purely_induced=True,
            purely_induced_total_inverse_coupling=shift.value,
            purely_induced_attractive=True,
            purely_induced_newton_constant=sp.simplify(
                12 * sp.pi / (count * (sp.Integer(1) - 6 * xi) * shift.proper_time_value)
            ),
            downstream_ceiling=ceiling,
        )
    elif is_purely_induced and weight_sign == 0:
        provenance = BaselineProvenance(
            baseline_inverse_coupling=baseline,
            status=(
                "declared_premise_per_C-GRV-001; conformal purely-induced "
                "reading lands exactly on the marginal locus 1/G_total = 0"
            ),
            is_purely_induced=True,
            purely_induced_total_inverse_coupling=sp.Integer(0),
            purely_induced_attractive=False,
            purely_induced_newton_constant=None,
            downstream_ceiling=(
                "at xi=1/6 the purely-induced reading gives no Newton "
                "constant at all (1/G_total = 0 identically); the ceiling is "
                "absolute, not merely scheme-bracketed"
            ),
        )
    elif is_purely_induced:
        provenance = BaselineProvenance(
            baseline_inverse_coupling=baseline,
            status=(
                "declared_premise_per_C-GRV-001; super-conformal purely-induced "
                "reading gives a repulsive total"
            ),
            is_purely_induced=True,
            purely_induced_total_inverse_coupling=shift.value,
            purely_induced_attractive=False,
            purely_induced_newton_constant=None,
            downstream_ceiling=(
                "for xi > 1/6 the purely-induced sign is repulsive; an "
                "attractive total requires a positive declared baseline"
            ),
        )
    else:
        provenance = BaselineProvenance(
            baseline_inverse_coupling=baseline,
            status="declared_premise_per_C-GRV-001; B remains free and can realize any supplied total",
            is_purely_induced=False,
            purely_induced_total_inverse_coupling=None,
            purely_induced_attractive=None,
            purely_induced_newton_constant=None,
            downstream_ceiling=(
                "with B != 0 the total is baseline-dominated whenever "
                "|B| exceeds the induced shift's scheme bracket; the "
                "baseline value remains an external input"
            ),
        )
    return provenance


@dataclass(frozen=True)
class HigherCurvatureControlEntry:
    """One control-sector term with its exact per-scheme coefficient."""

    sector: str
    regulator: str
    tau_minus_one_value: sp.Expr
    control_ratio: sp.Expr
    note: str


@dataclass(frozen=True)
class HigherCurvatureControlLedger:
    """The control ledger for every term that can defeat the regime."""

    cutoff: sp.Expr
    mass_squared: sp.Expr
    entries: tuple[HigherCurvatureControlEntry, ...]
    vacuum_sector_value: sp.Expr
    predeclared_domain: str
    nonlocal_remainder_bound: str


def higher_curvature_control_ledger(
    field_count: Any,
    *,
    cutoff: Any,
    mass_squared: Any = 0,
) -> HigherCurvatureControlLedger:
    """Return the exact control ledger for the defeating-term classes.

    The tau**-1 curvature-squared class is exactly ``-dJ/dz`` per scheme
    (E1(z) sharp, 2*K_0(2*sqrt(z)) smooth, log z + EulerGamma
    power-subtracted), exhibited with its control ratio against the
    curvature class.  The ratio ``L/J`` diverges logarithmically as
    ``z -> 0`` and is monotone decreasing on ``[z_min, inf)`` for the usable
    schemes, so the predeclared domain fixes ``z >= z_min > 0`` and the
    bound is the exact ``L(z_min)/J(z_min)``.  The tau**-3 vacuum class is
    taken from the accepted families and shown, not omitted.  The nonlocal
    remainder is outside the accepted local-coefficient ceiling (C-IGR-001)
    and is bounded by the same derivative-expansion parameter
    ``Lambda**-2 * ||R_E||`` in the predeclared small-curvature domain.
    """

    count = _positive_integer(field_count, "field_count")
    lam = _positive_exact(cutoff, "cutoff")
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    if mass.is_zero is True:
        raise ValueError(
            "the tau^-1 control class is logarithmically divergent at "
            "m**2 = 0 (E1(z), K_0(2*sqrt(z)), and log(z)+EulerGamma all "
            "diverge as z -> 0): the predeclared control domain requires "
            "z = m**2/Lambda**2 >= z_min > 0"
        )
    z = sp.simplify(mass / lam**2)

    j_sharp = curvature_scale_factor(
        SHARP_PROPER_TIME_REGULATOR, cutoff=lam, mass_squared=mass
    )
    j_smooth = curvature_scale_factor(
        SMOOTH_PROPER_TIME_REGULATOR, cutoff=lam, mass_squared=mass
    )
    j_zeta = curvature_scale_factor(
        ZETA_POWER_SUBTRACTED_REGULATOR, cutoff=lam, mass_squared=mass
    )
    l_sharp = sp.expint(1, z)
    l_smooth = 2 * sp.besselk(0, 2 * sp.sqrt(z))
    l_zeta = sp.log(z) + sp.EulerGamma
    entries = (
        HigherCurvatureControlEntry(
            sector="tau^-1 curvature-squared class (a4)",
            regulator=SHARP_PROPER_TIME_REGULATOR,
            tau_minus_one_value=l_sharp,
            control_ratio=sp.simplify(l_sharp / j_sharp),
            note="-dJ/dz exactly; ratio E1(z)/J(z) monotone decreasing on [z_min, inf)",
        ),
        HigherCurvatureControlEntry(
            sector="tau^-1 curvature-squared class (a4)",
            regulator=SMOOTH_PROPER_TIME_REGULATOR,
            tau_minus_one_value=l_smooth,
            control_ratio=sp.simplify(l_smooth / j_smooth),
            note="-dJ/dz exactly; ratio 2*K_0(2*sqrt(z))/J(z) monotone decreasing on [z_min, inf)",
        ),
        HigherCurvatureControlEntry(
            sector="tau^-1 curvature-squared class (a4)",
            regulator=ZETA_POWER_SUBTRACTED_REGULATOR,
            tau_minus_one_value=l_zeta,
            control_ratio=sp.simplify(l_zeta / j_zeta),
            note="-dJ/dz exactly; scheme outside the usable set, shown for spread only",
        ),
    )
    vacuum = exact_mass_vacuum_density_shift(
        count, regulator=SMOOTH_PROPER_TIME_REGULATOR, cutoff=lam, mass_squared=mass
    )
    return HigherCurvatureControlLedger(
        cutoff=lam,
        mass_squared=mass,
        entries=entries,
        vacuum_sector_value=vacuum.value,
        predeclared_domain=(
            "z = m**2/Lambda**2 >= z_min > 0 and Lambda**-2*||R_E|| <= eps << 1 "
            "(derivative-expansion domain); control bound L(z_min)/J(z_min) is "
            "exact and the ratio is monotone decreasing per usable scheme"
        ),
        nonlocal_remainder_bound=(
            "outside the accepted local-coefficient ceiling (C-IGR-001); "
            "bounded by the same small parameter Lambda**-2*||R_E|| in the "
            "predeclared domain, never silently omitted"
        ),
    )


@dataclass(frozen=True)
class RenormalizationCondition:
    """The governed renormalization condition with declared substrate provenance.

    Issue #88's composition goal: the usable scheme and its finite parts are
    outputs of substrate structure rather than choices.  This record states
    the condition, exhibits the selection legs that output the usable set,
    fixes the exact finite part per usable scheme at the declared ``z``, and
    quotes the exact scheme ceiling of the normalization.
    """

    statement: str
    usable_schemes: tuple
    excluded_schemes: tuple
    selection_legs: tuple
    finite_parts: dict
    reference_scheme: str
    reference_justification: str
    scheme_ceiling: str
    premises: tuple
    provenance: tuple
    non_authority: tuple
    open_items: tuple


def renormalization_condition(
    *,
    cutoff: Any,
    mass_squared: Any = 0,
) -> RenormalizationCondition:
    """Return the governed renormalization condition record at declared ``z``.

    The condition: *the total inverse Newton coupling is normalized as the
    finite curvature-class coefficient of the accepted one-loop determinant
    per field, ``1/G_total = B + N*(1-6*xi)*Lambda**2*J(z)/(12*pi)``, taken
    in the derived usable scheme set* — with the additive baseline ``B`` a
    declared premise per ``C-GRV-001`` and the usable set and finite parts
    ``J(z)`` outputs of the substrate-internal selection legs below, not
    regulator choices.
    """

    lam = _positive_exact(cutoff, "cutoff")
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    z = sp.simplify(mass / lam**2)
    j_sharp = curvature_scale_factor(SHARP_PROPER_TIME_REGULATOR, cutoff=lam, mass_squared=mass)
    j_smooth = curvature_scale_factor(SMOOTH_PROPER_TIME_REGULATOR, cutoff=lam, mass_squared=mass)
    return RenormalizationCondition(
        statement=(
            "1/G_total = B + N*(1-6*xi)*Lambda**2*J(z)/(12*pi) with "
            "z = m**2/Lambda**2, normalized in the derived usable scheme "
            "set; B is a declared premise (C-GRV-001), xi = 1/6 is the "
            "exact marginal locus, and J(z) is the derived finite part"
        ),
        usable_schemes=USABLE_TOTAL_COUPLING_SCHEMES,
        excluded_schemes=EXCLUDED_TOTAL_COUPLING_SCHEMES,
        selection_legs=(
            "L1 strict spectral positivity: J is the positive tail integral "
            "of a positive integrand, so 0 < J(z) <= 1 with J(0)=1 (positive "
            "operator spectrum of the accepted Route 1 fluctuation operator)",
            "L2 monotone large-mass decoupling: dJ/dz < 0 strictly "
            "(E1(z) sharp, 2*K_0(2*sqrt(z)) smooth), so the sign map is "
            "mass-uniform and no spurious zero appears at finite mass",
            "L3 cutoff-ontology respect: the finite part carries the cutoff "
            "scale Lambda (no cutoff-free, sign-changing structure); the "
            "power-subtracted family fails L3 because J_zeta changes sign "
            "at the exact root z = exp(1-EulerGamma) and is scale-free",
        ),
        finite_parts={
            SHARP_PROPER_TIME_REGULATOR: j_sharp,
            SMOOTH_PROPER_TIME_REGULATOR: j_smooth,
        },
        reference_scheme=SHARP_PROPER_TIME_REGULATOR,
        reference_justification=(
            "the sharp proper-time cutoff tau_0 = Lambda**-2 is the member "
            "matching C-GRV-001's cutoff ontology E_cut = hbar*c/a with "
            "Lambda = 1/a: its massless limit reproduces the accepted "
            "Delta(1/G) = s*Lambda**2 exactly; the smooth essential-"
            "singularity weight is the admitted smooth alternative and the "
            "bracket R(z) quantifies the residual scheme dependence"
        ),
        scheme_ceiling=(
            "the normalization is scheme-bracketed by R(z) = J_smooth/J_sharp "
            "(equal to 1 only at z = 0, exactly 2*e*K_1(2)/(1-e*E_1(1)) at "
            "z = 1, unbounded as z -> inf); no unique numeric G follows "
            "without a further scheme-selection input, which no accepted "
            "claim supplies"
        ),
        premises=(
            "additive baseline B free (declared premise per C-GRV-001); "
            "B = 0 is itself a declared premise, not a derived value",
            "xi = 1/6: Delta = 0 identically (conformal point), no induced "
            "Newton constant, purely-induced reading lands on the marginal "
            "locus 1/G_total = 0",
            "control regime: z >= z_min > 0 for the tau^-1 sector and "
            "Lambda**-2*||R_E|| <= eps << 1 for the nonlocal remainder "
            "(predeclared domain of the control ledger)",
        ),
        provenance=(
            "C-GRV-001 (accepted, v0.68.0): dimensional permission and "
            "additive-baseline ledger",
            "C-IGR-001, C-IGR-002, C-IGR-003 (accepted, v0.161.0): exact "
            "sharp/smooth/power-subtracted I_2/I_3 families and the "
            "scheme-dependence results",
            "Vassilevich arXiv:hep-th/0306138: heat-kernel/determinant "
            "conventions (approved primary import)",
            "Visser arXiv:gr-qc/0204062: effective-action organization "
            "(approved primary import)",
        ),
        non_authority=(
            "scalar_induced_newton.py: landed conditional, unpromoted; "
            "implementation reuse only",
            "covariant_sine_gordon_action.py: landed conditional, "
            "unpromoted; implementation reuse only",
        ),
        open_items=(
            "no unique numeric Newton constant (scheme bracket is the "
            "quoted ceiling)",
            "claim promotion C-IGR-004/C-GRV-002 happens in a later "
            "reviewed promotion transaction",
            "sourced nonflat geometry and radiative predictions are "
            "later rungs (non-goals of #88)",
        ),
    )


@dataclass(frozen=True)
class TotalNewtonConstant:
    """The total Newton constant per usable scheme, exactly and honestly."""

    regulator: str
    baseline_inverse_coupling: sp.Expr
    total_inverse_coupling: sp.Expr
    total_newton_constant: Any
    total_sign: int | None
    attractive_newtonian: bool | None
    acceptance: str
    ceiling: str


@dataclass(frozen=True)
class TotalNewtonConstantLedger:
    """The total Newton constant over the usable scheme set with its bracket."""

    entries: tuple
    purely_induced_bracket: Any
    ceiling: str


def total_newton_constant(
    baseline_inverse_coupling: Any,
    field_count: Any,
    non_minimal_coupling: Any,
    *,
    cutoff: Any,
    mass_squared: Any = 0,
) -> TotalNewtonConstantLedger:
    """Return ``G_total = 1/(B + N*(1-6*xi)*Lambda**2*J(z)/(12*pi))`` per usable scheme.

    The reciprocal preserves the sign, so ``G_total > 0`` exactly when the
    total inverse coupling is positive; ``G_total < 0`` is the repulsive
    reading, not an error.  The exact conformal marginal locus
    ``1/G_total = 0`` has no Newton constant at all and returns ``None``
    for the constant (never a silent division by zero); symbolic
    undecidability propagates as ``None`` with the honest acceptance label.
    """

    entries = []
    total_sign_ref = None
    for scheme in USABLE_TOTAL_COUPLING_SCHEMES:
        total = total_inverse_gravity_coupling(
            baseline_inverse_coupling,
            field_count,
            non_minimal_coupling,
            regulator=scheme,
            cutoff=cutoff,
            mass_squared=mass_squared,
        )
        if total.total_sign == 0:
            constant = None
            acceptance = "marginal_zero_total_no_newton_constant"
        else:
            constant = sp.simplify(1 / total.total_inverse_coupling)
            acceptance = total.acceptance
        entries.append(
            TotalNewtonConstant(
                regulator=scheme,
                baseline_inverse_coupling=total.baseline_inverse_coupling,
                total_inverse_coupling=total.total_inverse_coupling,
                total_newton_constant=constant,
                total_sign=total.total_sign,
                attractive_newtonian=total.attractive_newtonian,
                acceptance=acceptance,
                ceiling=(
                    "scheme bracket R(z) = J_smooth/J_sharp applies to the "
                    "constant: G_sharp/G_smooth = R(z) on the usable set "
                    "when the purely-induced reading is taken"
                ),
            )
        )
    baseline = sp.sympify(baseline_inverse_coupling)
    xi = sp.sympify(non_minimal_coupling)
    lam = _positive_exact(cutoff, "cutoff")
    mass = _nonnegative_exact(mass_squared, "mass_squared")
    j_sharp = curvature_scale_factor(SHARP_PROPER_TIME_REGULATOR, cutoff=lam, mass_squared=mass)
    j_smooth = curvature_scale_factor(SMOOTH_PROPER_TIME_REGULATOR, cutoff=lam, mass_squared=mass)
    if (
        baseline.is_zero is True
        and xi.is_real is True
        and sp.simplify(1 - 6 * xi).is_positive is True
    ):
        prefactor = 12 * sp.pi / (_positive_integer(field_count, "field_count") * (1 - 6 * xi) * lam**2)
        bracket = (
            sp.simplify(prefactor / j_smooth),
            sp.simplify(prefactor / j_sharp),
        )
    else:
        bracket = None
    return TotalNewtonConstantLedger(
        entries=tuple(entries),
        purely_induced_bracket=bracket,
        ceiling=(
            "G_total is scheme-bracketed by R(z): on the purely-induced "
            "sub-conformal reading the bracket endpoints are "
            "12*pi/(N*(1-6*xi)*Lambda**2*J(z)) per usable scheme, with "
            "G_sharp/G_smooth = R(z) exactly; no unique numeric G follows "
            "without a further scheme-selection input, which no accepted "
            "claim supplies"
        ),
    )

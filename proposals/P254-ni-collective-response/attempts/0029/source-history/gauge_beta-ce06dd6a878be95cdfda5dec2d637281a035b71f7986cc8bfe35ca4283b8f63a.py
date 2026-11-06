"""Exact gauge-only one- and two-loop coefficient ledgers.

The representation invariants and the perturbative loop weights are explicit
inputs to this module's scope.  The module does not infer a physical field
content from labels or charges, and its two-loop matrix excludes Yukawa terms.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Literal

import sympy as sp


MultipletKind = Literal["weyl_fermion", "complex_scalar"]


@dataclass(frozen=True)
class GaugeFactor:
    """One declared gauge factor in a fixed generator normalization."""

    label: str
    adjoint_casimir: Any
    is_abelian: bool = False


@dataclass(frozen=True)
class ProductMultiplet:
    """Supplied product-representation invariants for one multiplet.

    ``dynkin_indices[a]`` includes spectator-factor degeneracies, while
    ``quadratic_casimirs[a]`` is the Casimir of the representation under only
    factor ``a``. ``multiplicity`` counts identical copies not already present
    in those indices.
    """

    label: str
    kind: MultipletKind
    multiplicity: int
    dynkin_indices: tuple[Any, ...]
    quadratic_casimirs: tuple[Any, ...]


@dataclass(frozen=True)
class GaugeCoefficientLedger:
    """Exact coefficients and separately auditable contribution ledgers."""

    factors: tuple[GaugeFactor, ...]
    multiplets: tuple[ProductMultiplet, ...]
    one_loop_gauge: tuple[sp.Expr, ...]
    one_loop_weyl_fermions: tuple[sp.Expr, ...]
    one_loop_complex_scalars: tuple[sp.Expr, ...]
    one_loop: tuple[sp.Expr, ...]
    two_loop_gauge: tuple[tuple[sp.Expr, ...], ...]
    two_loop_weyl_fermions: tuple[tuple[sp.Expr, ...], ...]
    two_loop_complex_scalars: tuple[tuple[sp.Expr, ...], ...]
    two_loop_gauge_matrix: tuple[tuple[sp.Expr, ...], ...]
    beta_convention: str
    omitted_terms: tuple[str, ...]


@dataclass(frozen=True)
class AbelianGaugeRescalingLedger:
    """Covariance evidence for ``T'_a=rho_a*T_a`` and ``g'_a=g_a/rho_a``."""

    base: GaugeCoefficientLedger
    generator_rescalings: tuple[sp.Expr, ...]
    rescaled: GaugeCoefficientLedger
    expected_one_loop: tuple[sp.Expr, ...]
    expected_two_loop: tuple[tuple[sp.Expr, ...], ...]
    one_loop_residuals: tuple[sp.Expr, ...]
    two_loop_residuals: tuple[tuple[sp.Expr, ...], ...]


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be provably real")
    return expression


def _nonnegative_exact(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be nonnegative")
    return expression


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    return expression


def _normalize_factors(factors: Iterable[GaugeFactor]) -> tuple[GaugeFactor, ...]:
    normalized: list[GaugeFactor] = []
    labels: set[str] = set()
    for index, factor in enumerate(factors):
        if not isinstance(factor, GaugeFactor):
            raise TypeError("factors must contain GaugeFactor records")
        if not isinstance(factor.label, str) or not factor.label.strip():
            raise ValueError("factor labels must be non-empty strings")
        if factor.label in labels:
            raise ValueError("factor labels must be unique provenance keys")
        if not isinstance(factor.is_abelian, bool):
            raise TypeError("is_abelian must be boolean")
        labels.add(factor.label)
        casimir = _nonnegative_exact(
            factor.adjoint_casimir,
            f"factors[{index}].adjoint_casimir",
        )
        if factor.is_abelian and casimir != 0:
            raise ValueError("an Abelian factor must have zero adjoint Casimir")
        normalized.append(GaugeFactor(factor.label, casimir, factor.is_abelian))
    if not normalized:
        raise ValueError("at least one gauge factor is required")
    if sum(factor.is_abelian for factor in normalized) > 1:
        raise ValueError(
            "multiple Abelian factors require a kinetic-mixing treatment"
        )
    return tuple(normalized)


def _normalize_multiplets(
    multiplets: Iterable[ProductMultiplet],
    factor_count: int,
) -> tuple[ProductMultiplet, ...]:
    normalized: list[ProductMultiplet] = []
    labels: set[str] = set()
    for index, multiplet in enumerate(multiplets):
        if not isinstance(multiplet, ProductMultiplet):
            raise TypeError("multiplets must contain ProductMultiplet records")
        if not isinstance(multiplet.label, str) or not multiplet.label.strip():
            raise ValueError("multiplet labels must be non-empty strings")
        if multiplet.label in labels:
            raise ValueError("multiplet labels must be unique provenance keys")
        if multiplet.kind not in ("weyl_fermion", "complex_scalar"):
            raise ValueError("multiplet kind must be Weyl fermion or complex scalar")
        multiplicity = sp.sympify(multiplet.multiplicity)
        if (
            multiplicity.is_number is not True
            or multiplicity.is_integer is not True
            or multiplicity.is_positive is not True
        ):
            raise ValueError("multiplet multiplicities must be positive integers")
        if len(multiplet.dynkin_indices) != factor_count:
            raise ValueError("each multiplet needs one Dynkin index per factor")
        if len(multiplet.quadratic_casimirs) != factor_count:
            raise ValueError("each multiplet needs one quadratic Casimir per factor")
        labels.add(multiplet.label)
        normalized.append(
            ProductMultiplet(
                label=multiplet.label,
                kind=multiplet.kind,
                multiplicity=int(multiplicity),
                dynkin_indices=tuple(
                    _nonnegative_exact(value, f"multiplets[{index}].dynkin_indices[{a}]")
                    for a, value in enumerate(multiplet.dynkin_indices)
                ),
                quadratic_casimirs=tuple(
                    _nonnegative_exact(
                        value,
                        f"multiplets[{index}].quadratic_casimirs[{a}]",
                    )
                    for a, value in enumerate(multiplet.quadratic_casimirs)
                ),
            )
        )
    return tuple(normalized)


def _zero_matrix(size: int) -> list[list[sp.Expr]]:
    return [[sp.Integer(0) for _ in range(size)] for _ in range(size)]


def _freeze_matrix(matrix: list[list[sp.Expr]]) -> tuple[tuple[sp.Expr, ...], ...]:
    return tuple(tuple(sp.simplify(value) for value in row) for row in matrix)


def product_gauge_coefficients(
    factors: Iterable[GaugeFactor],
    multiplets: Iterable[ProductMultiplet] = (),
) -> GaugeCoefficientLedger:
    r"""Build exact gauge-only one- and two-loop coefficients.

    The fixed convention is

    ``mu dg_a/dmu = b_a g_a^3/L + g_a^3/L^2 sum_b B_ab g_b^2``,

    where ``L=16*pi^2``. Fermions are two-component Weyl fields, scalars are
    complex, and row ``a`` of ``B`` belongs to ``beta(g_a)``. The general
    perturbative weights are imported premises; this function performs the
    representation sums supplied by its inputs. Yukawa-dependent two-loop
    terms are deliberately absent.
    """

    factor_table = _normalize_factors(factors)
    size = len(factor_table)
    matter = _normalize_multiplets(multiplets, size)

    one_gauge = tuple(
        sp.simplify(-sp.Rational(11, 3) * factor.adjoint_casimir)
        for factor in factor_table
    )
    one_fermion = [sp.Integer(0) for _ in range(size)]
    one_scalar = [sp.Integer(0) for _ in range(size)]
    two_gauge = _zero_matrix(size)
    two_fermion = _zero_matrix(size)
    two_scalar = _zero_matrix(size)

    for a, factor in enumerate(factor_table):
        two_gauge[a][a] = sp.simplify(
            -sp.Rational(34, 3) * factor.adjoint_casimir**2
        )

    for multiplet in matter:
        multiplicity = sp.Integer(multiplet.multiplicity)
        for a, factor in enumerate(factor_table):
            s2_a = sp.sympify(multiplet.dynkin_indices[a])
            if multiplet.kind == "weyl_fermion":
                one_fermion[a] += multiplicity * sp.Rational(2, 3) * s2_a
                for b in range(size):
                    diagonal = (
                        sp.Rational(10, 3) * factor.adjoint_casimir
                        if a == b
                        else sp.Integer(0)
                    )
                    two_fermion[a][b] += multiplicity * s2_a * (
                        diagonal + 2 * multiplet.quadratic_casimirs[b]
                    )
            else:
                one_scalar[a] += multiplicity * sp.Rational(1, 3) * s2_a
                for b in range(size):
                    diagonal = (
                        sp.Rational(2, 3) * factor.adjoint_casimir
                        if a == b
                        else sp.Integer(0)
                    )
                    two_scalar[a][b] += multiplicity * s2_a * (
                        diagonal + 4 * multiplet.quadratic_casimirs[b]
                    )

    one_fermion_frozen = tuple(sp.simplify(value) for value in one_fermion)
    one_scalar_frozen = tuple(sp.simplify(value) for value in one_scalar)
    two_gauge_frozen = _freeze_matrix(two_gauge)
    two_fermion_frozen = _freeze_matrix(two_fermion)
    two_scalar_frozen = _freeze_matrix(two_scalar)
    one_loop = tuple(
        sp.simplify(one_gauge[a] + one_fermion_frozen[a] + one_scalar_frozen[a])
        for a in range(size)
    )
    two_loop = tuple(
        tuple(
            sp.simplify(
                two_gauge_frozen[a][b]
                + two_fermion_frozen[a][b]
                + two_scalar_frozen[a][b]
            )
            for b in range(size)
        )
        for a in range(size)
    )
    return GaugeCoefficientLedger(
        factors=factor_table,
        multiplets=matter,
        one_loop_gauge=one_gauge,
        one_loop_weyl_fermions=one_fermion_frozen,
        one_loop_complex_scalars=one_scalar_frozen,
        one_loop=one_loop,
        two_loop_gauge=two_gauge_frozen,
        two_loop_weyl_fermions=two_fermion_frozen,
        two_loop_complex_scalars=two_scalar_frozen,
        two_loop_gauge_matrix=two_loop,
        beta_convention=(
            "mu*dg_a/dmu = b_a*g_a^3/(16*pi^2) + "
            "g_a^3/(16*pi^2)^2*sum_b(B_ab*g_b^2); "
            "Weyl fermions; complex scalars; row a is beta(g_a)"
        ),
        omitted_terms=(
            "two-loop Yukawa contribution",
            "multiple-Abelian kinetic mixing",
            "threshold and matching corrections",
            "boundary conditions and physical field-content derivation",
        ),
    )


def gauge_only_beta(
    ledger: GaugeCoefficientLedger,
    couplings: Iterable[Any],
) -> tuple[sp.Expr, ...]:
    """Evaluate the declared gauge-only beta polynomial through two loops."""

    if not isinstance(ledger, GaugeCoefficientLedger):
        raise TypeError("ledger must be a GaugeCoefficientLedger")
    coupling_table = tuple(
        _exact_real(value, f"couplings[{index}]")
        for index, value in enumerate(couplings)
    )
    size = len(ledger.factors)
    if len(coupling_table) != size:
        raise ValueError("one coupling is required per gauge factor")
    loop = 16 * sp.pi**2
    return tuple(
        sp.simplify(
            ledger.one_loop[a] * coupling_table[a] ** 3 / loop
            + coupling_table[a] ** 3
            * sum(
                ledger.two_loop_gauge_matrix[a][b] * coupling_table[b] ** 2
                for b in range(size)
            )
            / loop**2
        )
        for a in range(size)
    )


def abelian_gauge_rescaling_ledger(
    ledger: GaugeCoefficientLedger,
    generator_rescalings: Iterable[Any],
) -> AbelianGaugeRescalingLedger:
    """Rescale Abelian generators and return exact coefficient residuals.

    A non-Abelian factor must have scale one. For an Abelian factor,
    ``S2`` and ``C2`` scale by ``rho**2``. Consequently
    ``b'_a=rho_a**2*b_a`` and ``B'_ab=rho_a**2*rho_b**2*B_ab``. Together with
    ``g'_a=g_a/rho_a``, these laws make every gauge-only beta component scale
    as ``beta'_a=beta_a/rho_a``.
    """

    if not isinstance(ledger, GaugeCoefficientLedger):
        raise TypeError("ledger must be a GaugeCoefficientLedger")
    scales = tuple(
        _positive_exact(value, f"generator_rescalings[{index}]")
        for index, value in enumerate(generator_rescalings)
    )
    size = len(ledger.factors)
    if len(scales) != size:
        raise ValueError("one generator rescaling is required per factor")
    for factor, scale in zip(ledger.factors, scales, strict=True):
        if not factor.is_abelian and scale != 1:
            raise ValueError("non-Abelian generator rescaling is outside this ledger")

    rescaled_multiplets = tuple(
        ProductMultiplet(
            label=multiplet.label,
            kind=multiplet.kind,
            multiplicity=multiplet.multiplicity,
            dynkin_indices=tuple(
                sp.simplify(scales[a] ** 2 * multiplet.dynkin_indices[a])
                for a in range(size)
            ),
            quadratic_casimirs=tuple(
                sp.simplify(scales[a] ** 2 * multiplet.quadratic_casimirs[a])
                for a in range(size)
            ),
        )
        for multiplet in ledger.multiplets
    )
    rescaled = product_gauge_coefficients(ledger.factors, rescaled_multiplets)
    expected_one = tuple(
        sp.simplify(scales[a] ** 2 * ledger.one_loop[a]) for a in range(size)
    )
    expected_two = tuple(
        tuple(
            sp.simplify(
                scales[a] ** 2
                * scales[b] ** 2
                * ledger.two_loop_gauge_matrix[a][b]
            )
            for b in range(size)
        )
        for a in range(size)
    )
    one_residuals = tuple(
        sp.simplify(rescaled.one_loop[a] - expected_one[a])
        for a in range(size)
    )
    two_residuals = tuple(
        tuple(
            sp.simplify(rescaled.two_loop_gauge_matrix[a][b] - expected_two[a][b])
            for b in range(size)
        )
        for a in range(size)
    )
    return AbelianGaugeRescalingLedger(
        base=ledger,
        generator_rescalings=scales,
        rescaled=rescaled,
        expected_one_loop=expected_one,
        expected_two_loop=expected_two,
        one_loop_residuals=one_residuals,
        two_loop_residuals=two_residuals,
    )

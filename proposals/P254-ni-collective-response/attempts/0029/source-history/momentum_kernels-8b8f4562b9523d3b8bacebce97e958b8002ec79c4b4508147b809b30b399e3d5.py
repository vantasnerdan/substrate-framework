"""Exact conditional low-momentum and static Green-kernel ledgers.

The APIs in this module derive consequences of explicitly supplied momentum
kernels.  They do not derive a charged loop species, a gauge kinetic action,
a physical mass or coupling, a spatial dimension, or an electromagnetic
dictionary.  In particular, analyticity of a correction does not establish
that its first coefficient is nonzero or that no lower fractional bare term
is present.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable

import sympy as sp


def _symbol(value: Any, name: str) -> sp.Symbol:
    if not isinstance(value, sp.Symbol):
        raise ValueError(f"{name} must be a SymPy symbol")
    return value


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number and expression.is_positive is not True:
        raise ValueError(f"{name} must be positive")
    if expression.is_number is False and expression.is_positive is False:
        raise ValueError(f"{name} must not be known nonpositive")
    return expression


def _nonzero(value: Any, name: str) -> sp.Expr:
    expression = sp.simplify(sp.sympify(value))
    if expression.is_zero is not False:
        qualifier = "provably nonzero" if expression.is_zero is None else "nonzero"
        raise ValueError(f"{name} must be {qualifier}")
    return expression


def _nonnegative_integer(value: Any, name: str) -> int:
    if not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


@dataclass(frozen=True)
class MassiveParameterKernel:
    r"""Exact ledger for a declared massive Feynman-parameter kernel.

    The dimensionless kernel is

    ``C*integral_0^1 du a(u) Q2/(m2+a(u) Q2)``,
    ``a(u)=u(1-u)``.

    ``m2`` denotes a positive mass squared.  The coefficient of ``Q2**n``
    follows from the beta integral and the Taylor series converges for
    ``|Q2| < 4*m2``.  The branch point at ``Q2=-4*m2`` is a property of this
    declared integrand, not a consequence of a mass gap alone.
    """

    transfer_variable: sp.Symbol
    parameter: sp.Symbol
    mass_squared: sp.Expr
    overall_coefficient: sp.Expr
    integrand: sp.Expr
    closed_form: sp.Expr
    coefficient_index: sp.Symbol
    coefficient_formula: sp.Expr
    convergence_radius: sp.Expr
    first_coefficient: sp.Expr
    second_coefficient: sp.Expr
    zero_transfer_limit: sp.Expr
    massless_at_fixed_positive_transfer: sp.Expr

    def coefficient(self, order: int) -> sp.Expr:
        """Return the exact coefficient of ``Q2**order`` for ``order >= 1``."""

        if not isinstance(order, int):
            raise TypeError("order must be an integer")
        if order < 1:
            raise ValueError("order must be positive")
        return sp.simplify(self.coefficient_formula.subs(self.coefficient_index, order))

    def series_polynomial(self, max_order: int) -> sp.Expr:
        """Return the exact Taylor polynomial through ``Q2**max_order``."""

        order = _nonnegative_integer(max_order, "max_order")
        return sp.simplify(
            sum(
                (self.coefficient(index) * self.transfer_variable**index for index in range(1, order + 1)),
                sp.S.Zero,
            )
        )

    def pointwise_remainder_integrand(self, max_order: int) -> sp.Expr:
        """Return the exact geometric-series remainder before integration."""

        order = _nonnegative_integer(max_order, "max_order")
        a = self.parameter * (1 - self.parameter)
        ratio = a * self.transfer_variable / self.mass_squared
        return sp.simplify(
            self.overall_coefficient
            * (-1) ** order
            * ratio ** (order + 1)
            / (1 + ratio)
        )


def massive_parameter_kernel(
    transfer_variable: sp.Symbol,
    mass_squared: Any,
    overall_coefficient: Any = 1,
    *,
    parameter: sp.Symbol | None = None,
) -> MassiveParameterKernel:
    r"""Derive the exact ledger for ``C int a Q2/(m2+a Q2) du``.

    Euclidean ``Q2`` and positive ``m2`` are used.  The returned massless
    limit is the fixed-``Q2>0`` limit.  It equals ``C``, whereas setting
    ``Q2=0`` first gives zero; callers must not interchange these limits.
    """

    q2 = _symbol(transfer_variable, "transfer variable")
    m2 = _positive(mass_squared, "mass squared")
    coefficient = sp.sympify(overall_coefficient)
    u = sp.Symbol("u", real=True) if parameter is None else _symbol(parameter, "parameter")
    a = u * (1 - u)
    integrand = sp.simplify(coefficient * a * q2 / (m2 + a * q2))
    closed = sp.simplify(
        coefficient
        * (
            1
            - 4
            * m2
            / sp.sqrt(q2 * (4 * m2 + q2))
            * sp.atanh(sp.sqrt(q2 / (4 * m2 + q2)))
        )
    )
    index = sp.Symbol("n", integer=True, positive=True)
    coefficient_formula = sp.simplify(
        coefficient
        * (-1) ** (index - 1)
        * sp.factorial(index) ** 2
        / (sp.factorial(2 * index + 1) * m2**index)
    )
    first = sp.simplify(coefficient_formula.subs(index, 1))
    second = sp.simplify(coefficient_formula.subs(index, 2))
    return MassiveParameterKernel(
        transfer_variable=q2,
        parameter=u,
        mass_squared=m2,
        overall_coefficient=coefficient,
        integrand=integrand,
        closed_form=closed,
        coefficient_index=index,
        coefficient_formula=coefficient_formula,
        convergence_radius=sp.simplify(4 * m2),
        first_coefficient=first,
        second_coefficient=second,
        zero_transfer_limit=sp.simplify(sp.limit(closed, q2, 0, dir="+")),
        massless_at_fixed_positive_transfer=sp.simplify(
            sp.limit(closed, m2, 0, dir="+")
        ),
    )


@dataclass(frozen=True)
class SpectralMomentExpansion:
    r"""Finite exact expansion of a declared subtracted Stieltjes kernel.

    For a separately supplied density ``rho(t)`` supported above ``gap``,
    this records

    ``K(Q2)=integral rho(t) Q2/(t+Q2) dt``.

    The caller must establish that the displayed improper integrals exist.
    A positive lower support bound controls the infrared denominator but does
    not by itself establish ultraviolet convergence or a nonzero moment.
    """

    transfer_variable: sp.Symbol
    spectral_variable: sp.Symbol
    spectral_density: sp.Expr
    gap: sp.Expr
    max_order: int
    exact_kernel: sp.Integral
    inverse_moments: tuple[sp.Integral, ...]
    series_polynomial: sp.Expr
    exact_remainder: sp.Integral
    pointwise_identity_residual: sp.Expr


def spectral_moment_expansion(
    transfer_variable: sp.Symbol,
    spectral_variable: sp.Symbol,
    spectral_density: Any,
    gap: Any,
    max_order: int,
) -> SpectralMomentExpansion:
    r"""Construct the exact finite inverse-moment identity through an order.

    The identity is algebraic before integration.  No convergence verdict is
    inferred from the gap.  Promoting the resulting coefficients therefore
    requires separate evidence that each displayed moment and remainder is
    finite under the caller's regulator/subtraction convention.
    """

    q2 = _symbol(transfer_variable, "transfer variable")
    t = _symbol(spectral_variable, "spectral variable")
    lower = _positive(gap, "gap")
    order = _nonnegative_integer(max_order, "max_order")
    density = sp.sympify(spectral_density)
    exact_integrand = sp.simplify(density * q2 / (t + q2))
    moments = tuple(
        sp.Integral(density / t**index, (t, lower, sp.oo))
        for index in range(1, order + 1)
    )
    polynomial = sp.expand(
        sum(
            (
                (-1) ** (index - 1)
                * q2**index
                * moments[index - 1]
                for index in range(1, order + 1)
            ),
            sp.S.Zero,
        )
    )
    remainder_integrand = sp.simplify(
        density * (-1) ** order * q2 ** (order + 1) / (t**order * (t + q2))
    )
    pointwise_polynomial = sp.simplify(
        sum(
            (
                density * (-1) ** (index - 1) * q2**index / t**index
                for index in range(1, order + 1)
            ),
            sp.S.Zero,
        )
    )
    return SpectralMomentExpansion(
        transfer_variable=q2,
        spectral_variable=t,
        spectral_density=density,
        gap=lower,
        max_order=order,
        exact_kernel=sp.Integral(exact_integrand, (t, lower, sp.oo)),
        inverse_moments=moments,
        series_polynomial=polynomial,
        exact_remainder=sp.Integral(remainder_integrand, (t, lower, sp.oo)),
        pointwise_identity_residual=sp.simplify(
            exact_integrand - pointwise_polynomial - remainder_integrand
        ),
    )


@dataclass(frozen=True)
class LeadingPowerLedger:
    """Smallest provably nonzero exponent of a supplied inverse kernel."""

    combined_terms: tuple[tuple[sp.Rational, sp.Expr], ...]
    leading_exponent: sp.Rational
    leading_coefficient: sp.Expr
    inverse_kernel: sp.Expr
    propagator_momentum_exponent: sp.Expr


def leading_power_ledger(
    momentum_squared: sp.Symbol,
    terms: Iterable[tuple[Any, Any]],
) -> LeadingPowerLedger:
    r"""Classify ``sum coefficient*(k^2)**exponent`` exactly.

    Equal exponents are combined before classification.  Every surviving
    coefficient must be provably nonzero after simplification; an undecidable
    symbolic coefficient raises ``ValueError`` instead of being silently
    treated as nonzero.  This makes cancellation premises explicit.
    """

    k2 = _symbol(momentum_squared, "momentum-squared variable")
    combined: dict[sp.Rational, sp.Expr] = {}
    supplied = list(terms)
    if not supplied:
        raise ValueError("at least one power term is required")
    for exponent_value, coefficient_value in supplied:
        exponent = sp.sympify(exponent_value)
        if exponent.is_number is not True or exponent.is_real is not True:
            raise ValueError("every exponent must be an exact real number")
        exponent = sp.Rational(exponent)
        if exponent <= 0:
            raise ValueError("every exponent must be positive")
        coefficient = sp.sympify(coefficient_value)
        combined[exponent] = sp.simplify(combined.get(exponent, sp.S.Zero) + coefficient)
    nonzero_terms: list[tuple[sp.Rational, sp.Expr]] = []
    for exponent in sorted(combined):
        coefficient = sp.simplify(combined[exponent])
        if coefficient.is_zero is True:
            continue
        _nonzero(coefficient, f"coefficient of exponent {exponent}")
        nonzero_terms.append((exponent, coefficient))
    if not nonzero_terms:
        raise ValueError("the supplied inverse kernel vanishes identically")
    leading_exponent, leading_coefficient = nonzero_terms[0]
    inverse_kernel = sp.simplify(
        sum(
            (coefficient * k2**exponent for exponent, coefficient in nonzero_terms),
            sp.S.Zero,
        )
    )
    return LeadingPowerLedger(
        combined_terms=tuple(nonzero_terms),
        leading_exponent=leading_exponent,
        leading_coefficient=leading_coefficient,
        inverse_kernel=inverse_kernel,
        propagator_momentum_exponent=sp.simplify(2 * leading_exponent),
    )


@dataclass(frozen=True)
class RieszGreenKernel:
    r"""Static inverse Fourier kernel away from the source point."""

    spatial_dimension: sp.Expr
    laplacian_power: sp.Expr
    radius: sp.Expr
    inverse_kernel_coefficient: sp.Expr
    fourier_convention: str
    convergence_condition: sp.Rel
    normalization: sp.Expr
    green_kernel: sp.Expr
    radial_power: sp.Expr
    radial_derivative: sp.Expr


def riesz_green_kernel(
    spatial_dimension: Any,
    laplacian_power: Any,
    radius: Any,
    inverse_kernel_coefficient: Any = 1,
    *,
    fourier_convention: str = "inverse_angular",
) -> RieszGreenKernel:
    r"""Return the conditional Green kernel of ``A*(k^2)**s``.

    The fixed convention is
    ``G(x)=integral d^d k/(2*pi)^d exp(i k.x)/(A*(k^2)^s)``.
    For ``r>0`` and ``0<s<d/2`` this equals
    ``Gamma(d/2-s)/(A*4^s*pi^(d/2)*Gamma(s))*r^(2s-d)``.
    Distributional contact terms and analytic continuation outside that
    domain require a separate prescription.
    """

    if fourier_convention != "inverse_angular":
        raise ValueError("unsupported Fourier convention")
    dimension = _positive(spatial_dimension, "spatial dimension")
    power = _positive(laplacian_power, "Laplacian power")
    radial = _symbol(radius, "radius")
    if radial.is_positive is False:
        raise ValueError("radius must not be known nonpositive")
    coefficient = _nonzero(inverse_kernel_coefficient, "inverse-kernel coefficient")
    if dimension.is_number and power.is_number and not bool(power < dimension / 2):
        raise ValueError("Riesz kernel requires 0 < power < dimension/2")
    normalization = sp.simplify(
        sp.gamma(dimension / 2 - power)
        / (coefficient * 4**power * sp.pi ** (dimension / 2) * sp.gamma(power))
    )
    green = sp.simplify(normalization * radial ** (2 * power - dimension))
    return RieszGreenKernel(
        spatial_dimension=dimension,
        laplacian_power=power,
        radius=radial,
        inverse_kernel_coefficient=coefficient,
        fourier_convention=fourier_convention,
        convergence_condition=sp.Lt(power, dimension / 2),
        normalization=normalization,
        green_kernel=green,
        radial_power=sp.simplify(2 * power - dimension),
        radial_derivative=sp.simplify(sp.diff(green, radial)),
    )


@dataclass(frozen=True)
class CriticalRieszLogKernel:
    r"""Reference-subtracted logarithmic limit at ``s=d/2``.

    The un-subtracted subcritical Riesz kernel has a divergent constant as
    ``s`` approaches ``d/2``. Subtracting its value at a positive reference
    radius before taking the limit leaves a finite logarithmic difference.
    For noninteger ``d`` this is an analytic continuation of the normalization,
    not a construction of a fractal metric-measure space.
    """

    spatial_dimension: sp.Expr
    critical_laplacian_power: sp.Expr
    radius: sp.Symbol
    reference_radius: sp.Expr
    inverse_kernel_coefficient: sp.Expr
    approach_parameter: sp.Symbol
    reference_subtracted_subcritical_kernel: sp.Expr
    raw_critical_limit: sp.Expr
    logarithmic_normalization: sp.Expr
    logarithmic_kernel: sp.Expr
    limit_reconstruction_residual: sp.Expr
    reference_residual: sp.Expr
    radial_derivative: sp.Expr


def critical_riesz_log_kernel(
    spatial_dimension: Any,
    radius: sp.Symbol,
    reference_radius: Any,
    inverse_kernel_coefficient: Any = 1,
) -> CriticalRieszLogKernel:
    r"""Derive the reference-subtracted ``s -> d/2`` Riesz limit.

    Starting from the same inverse-Fourier convention as
    :func:`riesz_green_kernel`, this returns

    ``2*log(r0/r)/(A*4**(d/2)*pi**(d/2)*Gamma(d/2))``.

    The subtraction fixes only a reference value; another reference changes
    the kernel by a constant. A real continuation parameter ``d`` is not by
    itself a geometric, Hausdorff, spectral, or walk dimension.
    """

    dimension = _positive(spatial_dimension, "spatial dimension")
    radial = _symbol(radius, "radius")
    if radial.is_positive is False:
        raise ValueError("radius must not be known nonpositive")
    reference = _positive(reference_radius, "reference radius")
    coefficient = _nonzero(
        inverse_kernel_coefficient, "inverse-kernel coefficient"
    )
    epsilon = sp.Symbol("epsilon", positive=True)
    subcritical_power = dimension / 2 - epsilon
    subcritical_normalization = (
        sp.gamma(epsilon)
        / (
            coefficient
            * 4**subcritical_power
            * sp.pi ** (dimension / 2)
            * sp.gamma(subcritical_power)
        )
    )
    subtracted = sp.simplify(
        subcritical_normalization
        * (radial ** (-2 * epsilon) - reference ** (-2 * epsilon))
    )
    raw_limit = sp.simplify(sp.limit(subtracted, epsilon, 0, dir="+"))
    reference_log = sp.log(reference / radial)
    normalization = sp.simplify(raw_limit / reference_log)
    logarithmic = sp.simplify(normalization * reference_log)
    return CriticalRieszLogKernel(
        spatial_dimension=dimension,
        critical_laplacian_power=sp.simplify(dimension / 2),
        radius=radial,
        reference_radius=reference,
        inverse_kernel_coefficient=coefficient,
        approach_parameter=epsilon,
        reference_subtracted_subcritical_kernel=subtracted,
        raw_critical_limit=raw_limit,
        logarithmic_normalization=normalization,
        logarithmic_kernel=logarithmic,
        limit_reconstruction_residual=sp.simplify(raw_limit - logarithmic),
        reference_residual=sp.simplify(logarithmic.subs(radial, reference)),
        radial_derivative=sp.simplify(sp.diff(logarithmic, radial)),
    )


@dataclass(frozen=True)
class RieszRadialForceLedger:
    r"""Conditional source-probe radial derivative of a subcritical kernel."""

    kernel: RieszGreenKernel
    source_strength: sp.Expr
    probe_strength: sp.Expr
    potential: sp.Expr
    potential_energy: sp.Expr
    radial_force: sp.Expr
    force_radial_power: sp.Expr
    inverse_square_dimension_family: sp.Expr
    inverse_square_residual: sp.Expr


def riesz_radial_force_law(
    spatial_dimension: Any,
    laplacian_power: Any,
    radius: sp.Symbol,
    source_strength: Any,
    probe_strength: Any,
    inverse_kernel_coefficient: Any = 1,
) -> RieszRadialForceLedger:
    r"""Return ``U=probe*source*G`` and ``F_r=-dU/dr`` conditionally.

    The subcritical kernel domain remains ``0<s<d/2``. Its radial-force power
    is ``2*s-d-1``; inverse-square behavior therefore fixes the family
    ``d=2*s+1``, not either parameter separately. Treating this derivative as
    a physical force requires the supplied source/probe dictionary and sign
    convention.
    """

    kernel = riesz_green_kernel(
        spatial_dimension,
        laplacian_power,
        radius,
        inverse_kernel_coefficient,
    )
    source = sp.sympify(source_strength)
    probe = sp.sympify(probe_strength)
    potential = sp.simplify(source * kernel.green_kernel)
    energy = sp.simplify(probe * potential)
    force = sp.simplify(-sp.diff(energy, kernel.radius))
    force_power = sp.simplify(
        2 * kernel.laplacian_power - kernel.spatial_dimension - 1
    )
    inverse_square_dimension = sp.simplify(2 * kernel.laplacian_power + 1)
    return RieszRadialForceLedger(
        kernel=kernel,
        source_strength=source,
        probe_strength=probe,
        potential=potential,
        potential_energy=energy,
        radial_force=force,
        force_radial_power=force_power,
        inverse_square_dimension_family=inverse_square_dimension,
        inverse_square_residual=sp.simplify(force_power + 2),
    )

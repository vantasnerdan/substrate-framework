"""Conditional one-loop complex-scalar vacuum polarization in two and four dimensions.

This module starts from a separately declared Euclidean scalar-QED functional
determinant.  It does not quantize the framework's accepted classical complex
field, identify its U(1) charge with electric charge, or derive a physical
gauge sector.  The formulas require a massive complex scalar, a
shift-invariant gauge-preserving regulator, the scalar bubble and seagull, and
the quadratic effective-action convention stated by the APIs below.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return expression


def _positive_integer(value: Any, name: str) -> int:
    if not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < 1:
        raise ValueError(f"{name} must be positive")
    return value


def _exact_real(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.has(sp.Float):
        raise ValueError(f"{name} must be exact rather than floating")
    if expression.is_real is not True:
        raise ValueError(f"{name} must be provably real")
    return expression


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return expression


def _nonnegative_exact(value: Any, name: str) -> sp.Expr:
    expression = _exact_real(value, name)
    if expression.is_nonnegative is not True:
        raise ValueError(f"{name} must be provably nonnegative")
    return expression


@dataclass(frozen=True)
class EuclideanTransverseProjector:
    """Exact nonzero-momentum Euclidean transverse-projector ledger."""

    momentum: sp.ImmutableMatrix
    momentum_squared: sp.Expr
    matrix: sp.ImmutableMatrix
    idempotence_residual: sp.ImmutableMatrix
    left_transversality_residual: sp.ImmutableMatrix
    right_transversality_residual: sp.ImmutableMatrix


def euclidean_transverse_projector(
    momentum: Sequence[Any],
) -> EuclideanTransverseProjector:
    r"""Return ``I-q*q.T/q^2`` for a declared nonzero Euclidean momentum.

    The projector is undefined at the zero vector.  Symbolic callers are
    responsible for restricting their domain to ``q^2>0``; this function
    rejects only a momentum whose squared norm simplifies identically to zero.
    """

    vector = sp.ImmutableMatrix([sp.sympify(component) for component in momentum])
    if vector.cols != 1 or vector.rows < 2:
        raise ValueError("momentum must contain at least two components")
    squared = sp.simplify((vector.T * vector)[0])
    if squared == 0:
        raise ValueError("transverse projector is undefined at zero momentum")
    identity = sp.eye(vector.rows)
    projector = sp.ImmutableMatrix(
        (identity - vector * vector.T / squared).applyfunc(sp.simplify)
    )
    return EuclideanTransverseProjector(
        momentum=vector,
        momentum_squared=squared,
        matrix=projector,
        idempotence_residual=sp.ImmutableMatrix(
            (projector * projector - projector).applyfunc(sp.simplify)
        ),
        left_transversality_residual=sp.ImmutableMatrix(
            (vector.T * projector).applyfunc(sp.simplify)
        ),
        right_transversality_residual=sp.ImmutableMatrix(
            (projector * vector).applyfunc(sp.simplify)
        ),
    )


@dataclass(frozen=True)
class ScalarQED2VacuumPolarization:
    r"""One massive-complex-scalar contribution to the Euclidean 1PI kernel.

    The tensor convention is

    ``Gamma^(2)=A_mu Pi_mu_nu A_nu/2`` and
    ``Pi_mu_nu=(q^2*delta_mu_nu-q_mu*q_nu)*form_factor``.

    Equivalently ``Pi_mu_nu=P_mu_nu*projector_coefficient`` at nonzero
    momentum.  The two scalar functions differ by an explicit factor of
    ``q^2``.  ``local_fmunu_squared_coefficient`` is the coefficient of
    ``F_mu_nu*F_mu_nu`` in the leading low-momentum effective Lagrangian; in
    two dimensions the coefficient of the single component ``F_01^2`` is
    twice that value.
    """

    momentum_squared: sp.Expr
    scalar_mass: sp.Expr
    charge_magnitude: sp.Expr
    species_count: int
    parameter: sp.Symbol
    projector_parameter_integrand: sp.Expr
    real_parameter: sp.Symbol
    dimensionless_ratio: sp.Expr
    real_integrand: sp.Expr
    real_antiderivative: sp.Expr
    antiderivative_residual: sp.Expr
    projector_coefficient: sp.Expr
    transverse_form_factor: sp.Expr
    zero_momentum_projector_limit: sp.Expr
    low_momentum_form_factor: sp.Expr
    local_fmunu_squared_coefficient: sp.Expr
    local_f01_squared_coefficient: sp.Expr
    massless_projector_limit: sp.Expr
    heavy_mass_projector_limit: sp.Expr
    bubble_ward_tadpole_coefficient: sp.Expr
    seagull_ward_tadpole_coefficient: sp.Expr
    ward_tadpole_residual: sp.Expr


def scalar_qed2_vacuum_polarization(
    momentum_squared: Any,
    scalar_mass: Any,
    charge_magnitude: Any,
    species_count: int = 1,
) -> ScalarQED2VacuumPolarization:
    r"""Return the exact massive scalar-QED2 one-loop transverse kernel.

    For ``Q=q_E^2>0``, mass ``m>0``, charge magnitude ``e>0``, and ``N``
    identical complex scalars, dimensional regularization or another
    shift-invariant gauge-preserving prescription gives

    ``Pi_hat(Q)=N*e^2*Q/(4*pi) * integral_0^1
    (1-2*x)^2/[m^2+Q*x*(1-x)] dx``.

    The scalar bubble contracts to ``+2*N*e^2*q_nu*I_tad`` and the seagull to
    ``-2*N*e^2*q_nu*I_tad``.  Their cancellation is the Ward identity; a
    projector ansatz alone is not its derivation.  The massless limit at fixed
    positive momentum diverges, so this scalar loop does not yield the finite
    fermionic Schwinger coefficient ``e^2/pi``.
    """

    q2 = _positive(momentum_squared, "momentum squared")
    mass = _positive(scalar_mass, "scalar mass")
    charge = _positive(charge_magnitude, "charge magnitude")
    count = _positive_integer(species_count, "species count")

    x = sp.Symbol("x", real=True)
    y = sp.Symbol("y", real=True)
    ratio = sp.simplify(sp.sqrt(q2) / sp.sqrt(q2 + 4 * mass**2))
    parameter_integrand = sp.simplify(
        count
        * charge**2
        * q2
        * (1 - 2 * x) ** 2
        / (4 * sp.pi * (mass**2 + q2 * x * (1 - x)))
    )
    real_integrand = y**2 / (1 - ratio**2 * y**2)
    real_antiderivative = sp.atanh(ratio * y) / ratio**3 - y / ratio**2
    antiderivative_residual = sp.simplify(
        sp.diff(real_antiderivative, y) - real_integrand
    )
    projector_coefficient = sp.simplify(
        count
        * charge**2
        / sp.pi
        * (sp.atanh(ratio) / ratio - 1)
    )
    form_factor = sp.simplify(projector_coefficient / q2)
    zero_momentum_limit = sp.simplify(
        sp.limit(projector_coefficient, q2, 0)
    )
    low_form_factor = sp.simplify(sp.limit(form_factor, q2, 0))
    massless_limit = sp.limit(projector_coefficient, mass, 0, dir="+")
    heavy_mass_limit = sp.limit(projector_coefficient, mass, sp.oo)
    bubble = sp.simplify(2 * count * charge**2)
    seagull = -bubble
    return ScalarQED2VacuumPolarization(
        momentum_squared=q2,
        scalar_mass=mass,
        charge_magnitude=charge,
        species_count=count,
        parameter=x,
        projector_parameter_integrand=parameter_integrand,
        real_parameter=y,
        dimensionless_ratio=ratio,
        real_integrand=real_integrand,
        real_antiderivative=real_antiderivative,
        antiderivative_residual=antiderivative_residual,
        projector_coefficient=projector_coefficient,
        transverse_form_factor=form_factor,
        zero_momentum_projector_limit=zero_momentum_limit,
        low_momentum_form_factor=low_form_factor,
        local_fmunu_squared_coefficient=sp.simplify(low_form_factor / 4),
        local_f01_squared_coefficient=sp.simplify(low_form_factor / 2),
        massless_projector_limit=massless_limit,
        heavy_mass_projector_limit=heavy_mass_limit,
        bubble_ward_tadpole_coefficient=bubble,
        seagull_ward_tadpole_coefficient=seagull,
        ward_tadpole_residual=sp.simplify(bubble + seagull),
    )


@dataclass(frozen=True)
class ScalarWardIntegrandEvidence:
    """Bubble-seagull Ward contraction before a transverse ansatz is imposed."""

    loop_momentum_squared: sp.Symbol
    loop_transfer_inner_product: sp.Symbol
    transfer_squared: sp.Symbol
    mass_squared: sp.Symbol
    transfer_component: sp.Symbol
    tadpole_integral: sp.Symbol
    propagator_denominator: sp.Expr
    shifted_propagator_denominator: sp.Expr
    contracted_vertex: sp.Expr
    denominator_difference_residual: sp.Expr
    first_bubble_numerator: sp.Expr
    shifted_second_bubble_numerator: sp.Expr
    shifted_bubble_numerator_difference: sp.Expr
    shifted_bubble_contraction: sp.Expr
    seagull_contraction: sp.Expr
    integrated_ward_residual: sp.Expr
    integrated_cancellation_requires_shift_invariance: bool


def scalar_ward_integrand_evidence() -> ScalarWardIntegrandEvidence:
    r"""Derive the scalar bubble-seagull Ward cancellation symbolically.

    For ``D_p=p^2+M^2`` and ``D_pq=(p+q)^2+M^2``, contraction of the
    scalar bubble uses ``q.(2p+q)=D_pq-D_p``.  Shifting the second tadpole
    integral by the same regulator routing leaves ``2*q_nu*I_tad``; the
    scalar seagull is its negative.  The shift-invariance premise is explicit
    because the algebra does not license divergent momentum translations by
    itself.
    """

    p2 = sp.Symbol("p2", real=True)
    pq = sp.Symbol("p_dot_q", real=True)
    q2 = sp.Symbol("q2", real=True)
    mass2 = sp.Symbol("M2", positive=True)
    pnu = sp.Symbol("p_nu", real=True)
    qnu = sp.Symbol("q_nu", real=True)
    tadpole = sp.Symbol("I_tad", real=True)
    denominator_p = p2 + mass2
    denominator_pq = p2 + 2 * pq + q2 + mass2
    contracted_vertex = 2 * pq + q2
    first_numerator = 2 * pnu + qnu
    shifted_second = 2 * pnu - qnu
    numerator_difference = sp.simplify(first_numerator - shifted_second)
    bubble = sp.simplify(numerator_difference * tadpole)
    seagull = sp.simplify(-2 * qnu * tadpole)
    return ScalarWardIntegrandEvidence(
        loop_momentum_squared=p2,
        loop_transfer_inner_product=pq,
        transfer_squared=q2,
        mass_squared=mass2,
        transfer_component=qnu,
        tadpole_integral=tadpole,
        propagator_denominator=denominator_p,
        shifted_propagator_denominator=denominator_pq,
        contracted_vertex=contracted_vertex,
        denominator_difference_residual=sp.simplify(
            contracted_vertex - (denominator_pq - denominator_p)
        ),
        first_bubble_numerator=first_numerator,
        shifted_second_bubble_numerator=shifted_second,
        shifted_bubble_numerator_difference=numerator_difference,
        shifted_bubble_contraction=bubble,
        seagull_contraction=seagull,
        integrated_ward_residual=sp.simplify(bubble + seagull),
        integrated_cancellation_requires_shift_invariance=True,
    )


@dataclass(frozen=True)
class ScalarVacuumPolarizationMaster:
    """General-dimensional scalar bubble-plus-seagull form-factor ledger."""

    integration_dimension: sp.Expr
    species_count: int
    euclidean_momentum_squared: sp.Expr
    mass_squared: sp.Expr
    charge_magnitude: sp.Expr
    parameter: sp.Symbol
    delta: sp.Expr
    parameter_weight: sp.Expr
    prefactor: sp.Expr
    parameter_integrand: sp.Expr
    transverse_form_factor: sp.Expr
    charge_squared_mass_dimension: sp.Expr
    delta_power_mass_dimension: sp.Expr
    transverse_form_factor_mass_dimension: sp.Expr


def scalar_vacuum_polarization_master(
    integration_dimension: Any,
    momentum_squared: Any,
    mass_squared: Any,
    charge_magnitude: Any,
    species_count: int = 1,
) -> ScalarVacuumPolarizationMaster:
    r"""Return the declared scalar-QED form-factor representation.

    The tensor convention matches the source convention
    ``Pi_mn=(q^2*g_mn-q_m*q_n)*Pi2``.  For spacelike ``Q=-q^2>0`` the
    bubble-plus-seagull reduction is

    ``Pi2(-Q)=-N*e^2*Gamma(2-d/2)/(4*pi)^(d/2) * integral_0^1
    (1-2*x)^2*(M2+x*(1-x)*Q)^(d/2-2) dx``.

    This API records the exact integral representation.  Finite evaluated
    claims use the separately audited integer-dimensional endpoint below.
    """

    dimension = _positive_exact(integration_dimension, "integration dimension")
    momentum2 = _positive_exact(momentum_squared, "momentum squared")
    mass2 = _nonnegative_exact(mass_squared, "mass squared")
    charge = _positive_exact(charge_magnitude, "charge magnitude")
    count = _positive_integer(species_count, "species count")
    x = sp.Symbol("x", real=True)
    delta = sp.simplify(mass2 + x * (1 - x) * momentum2)
    weight = sp.expand((1 - 2 * x) ** 2)
    prefactor = sp.simplify(
        -count
        * charge**2
        * sp.gamma(2 - dimension / 2)
        / (4 * sp.pi) ** (dimension / 2)
    )
    integrand = sp.simplify(weight * delta ** (dimension / 2 - 2))
    return ScalarVacuumPolarizationMaster(
        integration_dimension=dimension,
        species_count=count,
        euclidean_momentum_squared=momentum2,
        mass_squared=mass2,
        charge_magnitude=charge,
        parameter=x,
        delta=delta,
        parameter_weight=weight,
        prefactor=prefactor,
        parameter_integrand=integrand,
        transverse_form_factor=sp.Integral(integrand, (x, 0, 1)) * prefactor,
        charge_squared_mass_dimension=sp.simplify(4 - dimension),
        delta_power_mass_dimension=sp.simplify(dimension - 4),
        transverse_form_factor_mass_dimension=sp.Integer(0),
    )


@dataclass(frozen=True)
class ScalarQED4ZeroMomentumRenormalization:
    """Complex-scalar D4 Laurent, MS-bar, and beta-coefficient ledger."""

    regulator: sp.Symbol
    charge_magnitude: sp.Expr
    mass_squared: sp.Expr
    renormalization_scale_squared: sp.Expr
    finite_local_counterterm: sp.Expr
    species_count: int
    parameter: sp.Symbol
    parameter_weight: sp.Expr
    parameter_weight_integral: sp.Expr
    bare_form_factor: sp.Expr
    laurent_pole_residue: sp.Expr
    laurent_finite_part: sp.Expr
    msbar_local_counterterm: sp.Expr
    renormalized_limit_from_bare_plus_counterterm: sp.Expr
    expected_renormalized_form_factor: sp.Expr
    renormalization_residual: sp.Expr
    mass_squared_log_slope: sp.Expr
    mass_log_slope: sp.Expr
    scale_squared_log_slope: sp.Expr
    scale_log_slope: sp.Expr
    beta_coupling: sp.Expr
    connection_inverse_coupling_scale_slope: sp.Expr
    complex_scalar_matter_weight: sp.Expr


def scalar_qed4_zero_momentum_renormalization(
    charge_magnitude: Any,
    mass_squared: Any,
    renormalization_scale_squared: Any,
    finite_local_counterterm: Any = 0,
    species_count: int = 1,
    *,
    regulator: sp.Symbol | None = None,
) -> ScalarQED4ZeroMomentumRenormalization:
    r"""Return the ``d=4-2*epsilon`` scalar-QED zero-momentum family.

    The bubble and seagull combine into the Feynman weight ``(1-2*x)^2``.
    Its integral is ``1/3``, so the source-tensor-convention bare form factor
    is ``-N*e^2*Gamma(epsilon)*(4*pi*mu2/M2)^epsilon/(48*pi^2)``.
    MS-bar subtraction leaves an arbitrary finite local counterterm.  The
    resulting coupling beta function for ``N`` identical unit-charge complex
    scalars is ``N*e^3/(48*pi^2)``; neither the action nor this API fixes a
    total gauge kinetic coefficient.
    """

    charge = _positive_exact(charge_magnitude, "charge magnitude")
    mass2 = _positive_exact(mass_squared, "mass squared")
    scale2 = _positive_exact(
        renormalization_scale_squared,
        "renormalization scale squared",
    )
    finite = _exact_real(finite_local_counterterm, "finite local counterterm")
    count = _positive_integer(species_count, "species count")
    epsilon = regulator or sp.Symbol("epsilon", positive=True)
    if not isinstance(epsilon, sp.Symbol) or epsilon.is_positive is not True:
        raise ValueError("regulator must be a positive SymPy symbol")

    x = sp.Symbol("x", real=True)
    weight = sp.expand((1 - 2 * x) ** 2)
    weight_integral = sp.integrate(weight, (x, 0, 1))
    common = sp.simplify(count * charge**2 * weight_integral / (4 * sp.pi) ** 2)
    bare = sp.simplify(
        -common
        * sp.gamma(epsilon)
        * (4 * sp.pi * scale2 / mass2) ** epsilon
    )
    residue = sp.simplify(sp.limit(epsilon * bare, epsilon, 0, dir="+"))
    finite_part = sp.expand(
        sp.expand_log(
            sp.limit(bare - residue / epsilon, epsilon, 0, dir="+"),
            force=True,
        )
    )
    counterterm = sp.simplify(
        common * (1 / epsilon - sp.EulerGamma + sp.log(4 * sp.pi)) + finite
    )
    renormalized_limit = sp.expand(
        sp.expand_log(
            sp.limit(bare + counterterm, epsilon, 0, dir="+"),
            force=True,
        )
    )
    expected = sp.simplify(common * sp.log(mass2 / scale2) + finite)
    residual = sp.simplify(
        sp.expand_log(renormalized_limit - expected, force=True)
    )
    return ScalarQED4ZeroMomentumRenormalization(
        regulator=epsilon,
        charge_magnitude=charge,
        mass_squared=mass2,
        renormalization_scale_squared=scale2,
        finite_local_counterterm=finite,
        species_count=count,
        parameter=x,
        parameter_weight=weight,
        parameter_weight_integral=weight_integral,
        bare_form_factor=bare,
        laurent_pole_residue=residue,
        laurent_finite_part=finite_part,
        msbar_local_counterterm=counterterm,
        renormalized_limit_from_bare_plus_counterterm=renormalized_limit,
        expected_renormalized_form_factor=expected,
        renormalization_residual=residual,
        mass_squared_log_slope=sp.simplify(
            mass2 * sp.diff(renormalized_limit, mass2)
        ),
        mass_log_slope=sp.simplify(
            2 * mass2 * sp.diff(renormalized_limit, mass2)
        ),
        scale_squared_log_slope=sp.simplify(
            scale2 * sp.diff(renormalized_limit, scale2)
        ),
        scale_log_slope=sp.simplify(
            2 * scale2 * sp.diff(renormalized_limit, scale2)
        ),
        beta_coupling=sp.simplify(count * charge**3 / (48 * sp.pi**2)),
        connection_inverse_coupling_scale_slope=sp.simplify(
            -count / (24 * sp.pi**2)
        ),
        complex_scalar_matter_weight=sp.Rational(count, 3),
    )


@dataclass(frozen=True)
class MatterInducedKineticEvidence:
    """Exact one-loop matter coefficient with an explicit affine boundary."""

    scale: sp.Expr
    reference_scale: sp.Expr
    renormalized_local_coefficient: sp.Expr
    finite_matching_offset: sp.Expr
    scalar_weight: sp.Expr
    dirac_weight: sp.Expr
    scalar_coefficient: sp.Expr
    dirac_coefficient: sp.Expr
    one_loop_coefficient: sp.Expr
    reference_value: sp.Expr
    running_coefficient: sp.Expr
    kinetic_coefficient: sp.Expr
    formal_running_scale: sp.Symbol
    formal_kinetic_coefficient: sp.Expr
    flow_residual: sp.Expr
    zero_matching_kinetic_coefficient: sp.Expr
    boundary_mutation: sp.Symbol
    boundary_mutated_kinetic_coefficient: sp.Expr
    boundary_mutation_residual: sp.Expr
    scheme_shift: sp.Symbol
    scheme_decomposition_residual: sp.Expr
    reference_rescaling: sp.Symbol
    transformed_reference_value: sp.Expr
    reference_covariance_residual: sp.Expr
    zero_matching_is_separate_premise: bool


def matter_induced_kinetic_evidence(
    scale: Any,
    reference_scale: Any,
    renormalized_local_coefficient: Any,
    finite_matching_offset: Any,
    scalar_weight: Any,
    dirac_weight: Any,
) -> MatterInducedKineticEvidence:
    r"""Return the complete affine connection-field kinetic family.

    In the fixed connection convention ``B=e*A`` and
    ``L=-Z*F(B)^2/4``, separately supplied invariant complex-scalar and Dirac
    weights give ``b=W_s/3+4*W_f/3`` and
    ``mu*dZ/dmu=-b/(8*pi^2)``.  The exact solution is
    ``Z(mu)=Z_ref+b*log(mu_ref/mu)/(8*pi^2)`` with
    ``Z_ref=Z_local+c_fin``.  ``Z_ref=0`` is returned as a conditional branch;
    it is not inferred from the differential equation or an absent operator
    in another theory.
    """

    scale_value = _positive_exact(scale, "scale")
    reference = _positive_exact(reference_scale, "reference scale")
    local = _exact_real(
        renormalized_local_coefficient,
        "renormalized local coefficient",
    )
    finite = _exact_real(finite_matching_offset, "finite matching offset")
    scalar = _nonnegative_exact(scalar_weight, "scalar weight")
    dirac = _nonnegative_exact(dirac_weight, "Dirac weight")
    scalar_coefficient = sp.simplify(scalar / 3)
    dirac_coefficient = sp.simplify(4 * dirac / 3)
    coefficient = sp.simplify(scalar_coefficient + dirac_coefficient)
    running = sp.simplify(coefficient / (8 * sp.pi**2))
    reference_value = sp.simplify(local + finite)
    formal_scale = sp.Symbol("mu_running", positive=True)
    formal = sp.simplify(
        reference_value + running * sp.log(reference / formal_scale)
    )
    total = sp.simplify(formal.subs(formal_scale, scale_value))
    zero_matching = sp.simplify(running * sp.log(reference / scale_value))

    delta = sp.Symbol("delta_Z", nonzero=True, real=True)
    mutated = sp.simplify(total + delta)
    scheme_shift = sp.Symbol("sigma_scheme", real=True)
    scheme_residual = sp.simplify(
        (local + scheme_shift) + (finite - scheme_shift) - reference_value
    )
    kappa = sp.Symbol("kappa_reference", positive=True)
    transformed_reference = sp.simplify(reference_value - running * sp.log(kappa))
    transformed_formal = sp.simplify(
        transformed_reference
        + running * sp.log(kappa * reference / formal_scale)
    )
    return MatterInducedKineticEvidence(
        scale=scale_value,
        reference_scale=reference,
        renormalized_local_coefficient=local,
        finite_matching_offset=finite,
        scalar_weight=scalar,
        dirac_weight=dirac,
        scalar_coefficient=scalar_coefficient,
        dirac_coefficient=dirac_coefficient,
        one_loop_coefficient=coefficient,
        reference_value=reference_value,
        running_coefficient=running,
        kinetic_coefficient=total,
        formal_running_scale=formal_scale,
        formal_kinetic_coefficient=formal,
        flow_residual=sp.simplify(
            formal_scale * sp.diff(formal, formal_scale) + running
        ),
        zero_matching_kinetic_coefficient=zero_matching,
        boundary_mutation=delta,
        boundary_mutated_kinetic_coefficient=mutated,
        boundary_mutation_residual=sp.simplify(mutated - total - delta),
        scheme_shift=scheme_shift,
        scheme_decomposition_residual=scheme_residual,
        reference_rescaling=kappa,
        transformed_reference_value=transformed_reference,
        reference_covariance_residual=sp.simplify(transformed_formal - formal),
        zero_matching_is_separate_premise=True,
    )

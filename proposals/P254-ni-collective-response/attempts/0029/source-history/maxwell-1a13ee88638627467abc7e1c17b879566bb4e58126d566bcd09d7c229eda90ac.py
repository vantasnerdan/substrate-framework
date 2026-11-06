"""Conditional Maxwell action and static point-source consequences.

The APIs in this module require an independently supplied gauge kinetic
coefficient, current, spatial dimension, source normalization, boundary data,
and test-charge force dictionary.  They do not derive a photon, physical
electric charge, a preferred dimension, or a substrate electromagnetic sector.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import sympy as sp

from .momentum_kernels import riesz_green_kernel


def _positive(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_positive is not True:
        raise ValueError(f"{name} must be provably positive")
    return expression


def _spatial_dimension(value: Any) -> int:
    if not isinstance(value, int):
        raise TypeError("spatial dimension must be an integer")
    if value < 1:
        raise ValueError("spatial dimension must be positive")
    return value


def _coordinate(value: Any, name: str) -> sp.Symbol:
    if not isinstance(value, sp.Symbol):
        raise TypeError(f"{name} must be a SymPy symbol")
    return value


@dataclass(frozen=True)
class MaxwellEulerLagrange:
    """Exact component ledger for a declared flat-space Maxwell action."""

    coordinates: tuple[sp.Symbol, ...]
    covariant_potential: tuple[sp.Expr, ...]
    contravariant_current: tuple[sp.Expr, ...]
    inverse_metric: sp.ImmutableMatrix
    kinetic_coefficient: sp.Expr
    field_strength: sp.ImmutableMatrix
    raised_field_strength: sp.ImmutableMatrix
    lagrangian_density: sp.Expr
    euler_lagrange_residuals: tuple[sp.Expr, ...]
    expected_field_equation_residuals: tuple[sp.Expr, ...]
    derivation_residuals: tuple[sp.Expr, ...]
    source_only_euler_residuals: tuple[sp.Expr, ...]
    bianchi_residuals: tuple[sp.Expr, ...]
    continuity_identity: sp.Expr


def maxwell_euler_lagrange(
    covariant_potential: Sequence[Any],
    contravariant_current: Sequence[Any],
    coordinates: Sequence[sp.Symbol],
    inverse_metric: Any,
    kinetic_coefficient: Any = 1,
) -> MaxwellEulerLagrange:
    r"""Vary ``-kappa*F_mu_nu*F^mu_nu/4-j^mu*A_mu`` exactly.

    ``inverse_metric`` must be a constant, symmetric, invertible matrix in the
    supplied coordinates.  The returned Euler residual uses the convention
    ``partial L/partial A_nu - partial_mu[partial L/partial(partial_mu A_nu)]``
    and therefore equals ``kappa*partial_mu F^mu_nu-j^nu``.  Gauge invariance
    of the source action additionally requires the supplied current to be
    conserved after boundary terms are removed.
    """

    xs = tuple(_coordinate(value, "coordinate") for value in coordinates)
    potentials = tuple(sp.sympify(value) for value in covariant_potential)
    currents = tuple(sp.sympify(value) for value in contravariant_current)
    dimension = len(xs)
    if dimension < 2:
        raise ValueError("at least two spacetime coordinates are required")
    if len(potentials) != dimension or len(currents) != dimension:
        raise ValueError("potential, current, and coordinate lengths must agree")

    metric = sp.ImmutableMatrix(inverse_metric)
    if metric.shape != (dimension, dimension):
        raise ValueError("inverse metric shape must match the spacetime dimension")
    if metric != metric.T:
        raise ValueError("inverse metric must be symmetric")
    if any(entry.has(*xs) for entry in metric):
        raise ValueError("inverse metric must be constant in the supplied coordinates")
    determinant = sp.simplify(metric.det())
    if determinant.is_zero is not False:
        raise ValueError("inverse metric must be provably invertible")

    kappa = _positive(kinetic_coefficient, "kinetic coefficient")
    lower = sp.ImmutableMatrix(
        dimension,
        dimension,
        lambda mu, nu: sp.simplify(
            sp.diff(potentials[nu], xs[mu])
            - sp.diff(potentials[mu], xs[nu])
        ),
    )
    raised = sp.ImmutableMatrix(metric * lower * metric.T)
    kinetic = -kappa * sp.Add(
        *(lower[mu, nu] * raised[mu, nu]
          for mu in range(dimension)
          for nu in range(dimension))
    ) / 4
    source = -sp.Add(
        *(potentials[nu] * currents[nu] for nu in range(dimension))
    )
    lagrangian = sp.simplify(kinetic + source)

    euler: list[sp.Expr] = []
    expected: list[sp.Expr] = []
    for nu in range(dimension):
        direct = sp.diff(lagrangian, potentials[nu])
        divergence_terms: list[sp.Expr] = []
        for mu in range(dimension):
            derivative_variable = sp.diff(potentials[nu], xs[mu])
            if derivative_variable == 0:
                divergence_terms.append(sp.S.Zero)
            else:
                divergence_terms.append(
                    sp.diff(
                        sp.diff(lagrangian, derivative_variable),
                        xs[mu],
                    )
                )
        divergence = sp.Add(*divergence_terms)
        euler.append(sp.simplify(direct - divergence))
        expected.append(sp.simplify(
            kappa
            * sp.Add(*(sp.diff(raised[mu, nu], xs[mu]) for mu in range(dimension)))
            - currents[nu]
        ))

    bianchi = tuple(
        sp.simplify(
            sp.diff(lower[nu, rho], xs[mu])
            + sp.diff(lower[rho, mu], xs[nu])
            + sp.diff(lower[mu, nu], xs[rho])
        )
        for mu in range(dimension)
        for nu in range(dimension)
        for rho in range(dimension)
    )
    double_divergence = sp.simplify(
        sp.Add(
            *(sp.diff(
                sp.Add(*(sp.diff(raised[mu, nu], xs[mu]) for mu in range(dimension))),
                xs[nu],
            ) for nu in range(dimension))
        )
    )
    return MaxwellEulerLagrange(
        coordinates=xs,
        covariant_potential=potentials,
        contravariant_current=currents,
        inverse_metric=metric,
        kinetic_coefficient=kappa,
        field_strength=lower,
        raised_field_strength=raised,
        lagrangian_density=lagrangian,
        euler_lagrange_residuals=tuple(euler),
        expected_field_equation_residuals=tuple(expected),
        derivation_residuals=tuple(
            sp.simplify(actual - target)
            for actual, target in zip(euler, expected, strict=True)
        ),
        source_only_euler_residuals=tuple(-current for current in currents),
        bianchi_residuals=bianchi,
        continuity_identity=double_divergence,
    )


@dataclass(frozen=True)
class StaticMaxwellPointSource:
    """Exact radial point-source and test-charge ledger in a supplied dimension."""

    spatial_dimension: int
    radius: sp.Symbol
    source_charge: sp.Expr
    test_charge: sp.Expr
    kinetic_coefficient: sp.Expr
    unit_sphere_area: sp.Expr
    reference_radius: sp.Expr | None
    reference_potential: sp.Expr
    boundary_condition: str
    potential: sp.Expr
    radial_electric_field: sp.Expr
    potential_energy: sp.Expr
    radial_force: sp.Expr
    radial_harmonic_residual: sp.Expr
    electric_flux: sp.Expr
    normalized_source_flux: sp.Expr
    potential_radial_power: sp.Expr | None
    field_radial_power: sp.Expr
    decays_at_infinity: bool
    inverse_square_force: bool


def static_maxwell_point_source(
    spatial_dimension: int,
    radius: sp.Symbol,
    source_charge: Any,
    test_charge: Any,
    kinetic_coefficient: Any = 1,
    *,
    reference_radius: Any | None = None,
    reference_potential: Any = 0,
) -> StaticMaxwellPointSource:
    r"""Return a source-normalized radial solution of ``-kappa*Delta phi=rho``.

    For integer ``d>2`` the boundary condition is ``phi -> 0`` at infinity and
    the implementation reuses the accepted Riesz kernel.  In ``d=2`` a
    positive reference radius is required for the logarithmic potential.  In
    ``d=1`` the returned even full-line branch has ``phi(0)=reference_potential``
    and grows linearly in magnitude.  Homogeneous additions require separate
    boundary data.  The test-charge relations ``U=q*phi`` and ``F=q*E`` are
    declared model structure, not consequences of the field equation alone.
    """

    dimension = _spatial_dimension(spatial_dimension)
    radial = _coordinate(radius, "radius")
    if radial.is_positive is not True:
        raise ValueError("radius must be provably positive")
    charge = sp.sympify(source_charge)
    probe = sp.sympify(test_charge)
    kappa = _positive(kinetic_coefficient, "kinetic coefficient")
    reference_value = sp.sympify(reference_potential)
    sphere_area = sp.simplify(
        2 * sp.pi ** sp.Rational(dimension, 2)
        / sp.gamma(sp.Rational(dimension, 2))
    )

    reference: sp.Expr | None = None
    if dimension > 2:
        if reference_radius is not None or reference_value != 0:
            raise ValueError(
                "d>2 uses the zero-at-infinity branch; reference data are not accepted"
            )
        potential = sp.simplify(
            charge
            * riesz_green_kernel(
                dimension,
                1,
                radial,
                inverse_kernel_coefficient=kappa,
            ).green_kernel
        )
        boundary = "zero_at_infinity"
        potential_power: sp.Expr | None = sp.Integer(2 - dimension)
    elif dimension == 2:
        if reference_radius is None:
            raise ValueError("d=2 requires a positive reference radius")
        reference = _positive(reference_radius, "reference radius")
        potential = (
            reference_value
            - charge * sp.log(radial / reference) / (kappa * sphere_area)
        )
        boundary = "potential_fixed_at_reference_radius"
        potential_power = None
    else:
        if reference_radius is not None:
            raise ValueError("d=1 uses the origin as its potential reference")
        potential = sp.simplify(
            reference_value - charge * radial / (kappa * sphere_area)
        )
        boundary = "even_full_line_branch_with_potential_fixed_at_origin"
        potential_power = sp.Integer(1)

    electric_field = sp.simplify(-sp.diff(potential, radial))
    energy = sp.simplify(probe * potential)
    force = sp.simplify(probe * electric_field)
    harmonic_residual = sp.simplify(
        sp.diff(radial ** (dimension - 1) * sp.diff(potential, radial), radial)
        / radial ** (dimension - 1)
    )
    flux = sp.simplify(
        sphere_area * radial ** (dimension - 1) * electric_field
    )
    return StaticMaxwellPointSource(
        spatial_dimension=dimension,
        radius=radial,
        source_charge=charge,
        test_charge=probe,
        kinetic_coefficient=kappa,
        unit_sphere_area=sphere_area,
        reference_radius=reference,
        reference_potential=reference_value,
        boundary_condition=boundary,
        potential=potential,
        radial_electric_field=electric_field,
        potential_energy=energy,
        radial_force=force,
        radial_harmonic_residual=harmonic_residual,
        electric_flux=flux,
        normalized_source_flux=sp.simplify(kappa * flux),
        potential_radial_power=potential_power,
        field_radial_power=sp.Integer(1 - dimension),
        decays_at_infinity=dimension > 2,
        inverse_square_force=dimension == 3,
    )

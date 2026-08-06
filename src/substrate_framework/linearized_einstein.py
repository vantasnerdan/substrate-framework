"""Conditional linearized-Einstein source-response ledgers.

The declared spacetime convention is mostly plus,
``eta=diag(-1, 1, 1, 1)``, with
``box=-partial_t**2/c**2+nabla**2``.  In harmonic gauge the separately
supplied Einstein coupling gives ``box(hbar_ab)=-2*kappa*T_ab`` with
``kappa=8*pi*G/c**4``.

The APIs below derive static-monopole and no-incoming retarded line-source
consequences from that one conditional normalization.  They do not derive
``G`` or ``c``, identify a microscopic gravity mechanism, or couple a
sine-Gordon source automatically.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .exact_symbolic import exact_real as _exact_real
from .exact_symbolic import positive_exact as _positive_exact
from .tt_angular import frobenius_norm_squared, tt_project_symmetric


def _direction_cosine(value: Any) -> sp.Expr:
    cosine = _exact_real(value, "direction_cosine")
    if cosine.is_number and not bool(abs(cosine) <= 1):
        raise ValueError(
            f"direction_cosine must lie in [-1, 1]; got {value!r}"
        )
    return cosine


@dataclass(frozen=True)
class HarmonicGaugeEinsteinLedger:
    """Mostly-plus coefficients of supplied-coupling linearized Einstein theory."""

    minkowski_metric_covariant: sp.ImmutableMatrix
    newton_constant: sp.Expr
    signal_speed: sp.Expr
    d_alembertian_time_coefficient: sp.Expr
    d_alembertian_spatial_coefficient: sp.Expr
    einstein_stress_coupling: sp.Expr
    trace_reversed_wave_source_coefficient: sp.Expr
    retarded_far_zone_coefficient: sp.Expr
    isaacson_flux_coefficient: sp.Expr


@dataclass(frozen=True)
class WeakFieldMonopole:
    """Leading mostly-plus static field of a positive isolated rest mass.

    ``trace_reversed_poisson_source_coefficient`` multiplies
    ``delta**3(x)`` in the static equation for ``hbar_00``.
    ``spatial_metric_factor`` multiplies the Cartesian ``delta_ij``.
    """

    mass: sp.Expr
    radius: sp.Expr
    newtonian_potential: sp.Expr
    potential_over_speed_squared: sp.Expr
    trace_reversed_poisson_source_coefficient: sp.Expr
    trace_reversed_metric_00: sp.Expr
    metric_perturbation_00: sp.Expr
    metric_00: sp.Expr
    spatial_metric_factor: sp.Expr


@dataclass(frozen=True)
class RetardedLineSourceRadiationLedger:
    """No-incoming far-zone TT readout of a longitudinal line-source term.

    If ``H_zeta(t)=integral dx T_xx(t+zeta*x/c,x)``, then the far field is
    ``h_TT=(4*G/(c**4*R))*TT[e_x e_x]*H_zeta``.  The azimuth-integrated power
    is the integral over ``zeta`` of
    ``azimuth_power_coefficient*<dot(H_zeta)**2>``.  The source history and
    no-incoming boundary condition are premises; this ledger does not solve
    them.
    """

    direction_cosine: sp.Expr
    projected_longitudinal_tensor: sp.ImmutableMatrix
    projected_norm_squared: sp.Expr
    retarded_far_zone_coefficient: sp.Expr
    differential_solid_angle_power_coefficient: sp.Expr
    azimuth_power_coefficient: sp.Expr


def harmonic_gauge_einstein_ledger(
    newton_constant: Any,
    signal_speed: Any,
) -> HarmonicGaugeEinsteinLedger:
    r"""Return the shared coefficients in ``box(hbar_ab)=-2*kappa*T_ab``.

    The no-incoming Green function gives the leading far-zone coefficient
    ``4*G/c**4=-source_coefficient/(4*pi)``.  The declared TT flux functional
    has coefficient ``c**3/(32*pi*G)``.
    """

    newton = _positive_exact(newton_constant, "newton_constant")
    speed = _positive_exact(signal_speed, "signal_speed")
    coupling = sp.simplify(8 * sp.pi * newton / speed**4)
    source = sp.simplify(-2 * coupling)
    retarded = sp.simplify(4 * newton / speed**4)
    flux = sp.simplify(speed**3 / (32 * sp.pi * newton))
    if sp.simplify(retarded + source / (4 * sp.pi)) != 0:
        raise AssertionError(
            "retarded coefficient does not match the harmonic-gauge source sign"
        )
    return HarmonicGaugeEinsteinLedger(
        minkowski_metric_covariant=sp.ImmutableMatrix(sp.diag(-1, 1, 1, 1)),
        newton_constant=newton,
        signal_speed=speed,
        d_alembertian_time_coefficient=sp.simplify(-1 / speed**2),
        d_alembertian_spatial_coefficient=sp.Integer(1),
        einstein_stress_coupling=coupling,
        trace_reversed_wave_source_coefficient=source,
        retarded_far_zone_coefficient=retarded,
        isaacson_flux_coefficient=flux,
    )


def weak_field_monopole(
    newton_constant: Any,
    signal_speed: Any,
    mass: Any,
    radius: Any,
) -> WeakFieldMonopole:
    r"""Return the leading mostly-plus harmonic-gauge static monopole.

    For ``T_00=M*c**2*delta**3(x)``, the static equation is
    ``nabla**2(hbar_00)=-16*pi*G*M*delta**3(x)/c**2``.  Since
    ``nabla**2(1/r)=-4*pi*delta**3(x)``,
    ``hbar_00=4*G*M/(c**2*r)=-4*Phi/c**2`` for ``Phi=-G*M/r``.  Trace
    reversal then gives ``g_00=-1-2*Phi/c**2`` and
    ``g_ij=(1-2*Phi/c**2)*delta_ij``.  No Schwarzschild metric is supplied as
    an input.
    """

    ledger = harmonic_gauge_einstein_ledger(newton_constant, signal_speed)
    source_mass = _positive_exact(mass, "mass")
    distance = _positive_exact(radius, "radius")
    potential = sp.simplify(-ledger.newton_constant * source_mass / distance)
    reduced = sp.simplify(potential / ledger.signal_speed**2)
    poisson_source = sp.simplify(
        ledger.trace_reversed_wave_source_coefficient
        * source_mass
        * ledger.signal_speed**2
    )
    trace_reversed = sp.simplify(-4 * reduced)
    metric_perturbation = sp.simplify(-2 * reduced)
    if sp.simplify(-4 * sp.pi * distance * trace_reversed - poisson_source) != 0:
        raise AssertionError(
            "static Green function does not reproduce the mostly-plus source sign"
        )
    return WeakFieldMonopole(
        mass=source_mass,
        radius=distance,
        newtonian_potential=potential,
        potential_over_speed_squared=reduced,
        trace_reversed_poisson_source_coefficient=poisson_source,
        trace_reversed_metric_00=trace_reversed,
        metric_perturbation_00=metric_perturbation,
        metric_00=sp.simplify(-1 + metric_perturbation),
        spatial_metric_factor=sp.simplify(1 - 2 * reduced),
    )


def retarded_line_source_radiation_ledger(
    newton_constant: Any,
    signal_speed: Any,
    direction_cosine: Any,
) -> RetardedLineSourceRadiationLedger:
    r"""Return the exact finite-wavelength angular ledger for a line source."""

    ledger = harmonic_gauge_einstein_ledger(newton_constant, signal_speed)
    cosine = _direction_cosine(direction_cosine)
    sine = sp.sqrt(1 - cosine**2)
    line_tensor = sp.diag(1, 0, 0)
    projected = sp.ImmutableMatrix(
        tt_project_symmetric(line_tensor, [cosine, sine, 0])
    ).applyfunc(sp.simplify)
    norm_squared = sp.simplify(frobenius_norm_squared(projected))
    expected_norm = (1 - cosine**2) ** 2 / 2
    if sp.simplify(norm_squared - expected_norm) != 0:
        raise AssertionError(
            "longitudinal TT angular norm did not close: "
            f"expected {expected_norm}, got {norm_squared}"
        )
    differential = (
        ledger.newton_constant
        * expected_norm
        / (2 * sp.pi * ledger.signal_speed**5)
    )
    azimuth = (
        ledger.newton_constant
        * (1 - cosine**2) ** 2
        / (2 * ledger.signal_speed**5)
    )
    if sp.simplify(2 * sp.pi * differential - azimuth) != 0:
        raise AssertionError(
            "azimuth-integrated power coefficient does not match solid angle"
        )
    return RetardedLineSourceRadiationLedger(
        direction_cosine=cosine,
        projected_longitudinal_tensor=projected,
        projected_norm_squared=expected_norm,
        retarded_far_zone_coefficient=ledger.retarded_far_zone_coefficient,
        differential_solid_angle_power_coefficient=differential,
        azimuth_power_coefficient=azimuth,
    )

"""End-to-end sine-Gordon-to-gravity composition chain.

This module assembles the campaign construction: accepted upstream
sine-Gordon parameters flow through the accepted scale map (C-MED-003), the
spectral induction ledger (:mod:`sg_spectral_induction`), and the accepted
Einstein-scalar and quadrupole machinery (C-STG-001/002, C-GW-001..010), and
produce the static weak-field and radiative predictions with the *derived*
gravitational coupling.

Declared premises (named here, surfaced in every ledger):

1. quantum premise: ``J_scale = hbar`` (one-quantum action identification,
   following the conditional precedent of C-SCL-001), so the normalized
   action increment ``h`` equals ``hbar / J_scale`` and the physical phonon
   mass is ``hbar * omega_0 / c**2``;
2. emergence premise: the spacetime geometry felt by excitations is the
   medium's Gordon geometry (C-GOR-001), so curvature couples to the sG
   vacuum polarization that induces the Einstein-Hilbert term;
3. semiclassical premise: each accepted tower species (phonon,
   kink-antikink, breather levels) contributes to the one-loop induction
   with its normalized mass and declared multiplicity;
4. convention: ``kappa = 8*pi*G/c**4`` is the named C-STG-001 coupling
   convention factor, not a cross-dimensional identification.

Standard imports (role explicit): the Schwarzschild exterior and its
weak-field predictions (Newtonian potential, light-deflection angle
``4*G*M/(c**2*b)``) enter only downstream of the derived coupling; the
quadrupole power law is the accepted C-GW-001/C-GW-004 form
``P = 2*G*<d**2>/15`` with ``d`` the normalized longitudinal moment's third
derivative, converted to physical units below.

This module does not identify the model's free upstream scales with any
observed system and performs no empirical confrontation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .dimensional_sine_gordon import (
    DimensionalSineGordonCoefficients,
    dimensional_sine_gordon_scales,
)
from .sg_spectral_induction import newton_q_normalized, tower_inverse_g_normalized

# Accepted C-SG-010 quadrature bracket for <mu'''(t)^2> at omega = 1/sqrt(2).
MU_THIRD_POWER_LOWER = sp.Float("379.4646380687", 15)
MU_THIRD_POWER_UPPER = sp.Float("379.4646380688", 15)
MU_THIRD_POWER_MID = (MU_THIRD_POWER_LOWER + MU_THIRD_POWER_UPPER) / 2


def _positive_exact(value: Any, name: str) -> sp.Expr:
    expression = sp.sympify(value)
    if expression.is_number:
        if expression.is_real is not True or not float(expression) > 0.0:
            raise ValueError(f"{name} must be real and positive")
    return expression


@dataclass(frozen=True)
class SGGravityChain:
    """The composed upstream-to-gravity ledger.

    Attributes are exact symbolic expressions in the upstream coefficients
    ``(lambda, T, mu)``, the declared action increment ``h``, and the
    induction premises ``xi`` and ``kink_multiplicity``.
    """

    coefficients: DimensionalSineGordonCoefficients
    action_quantum: sp.Expr
    xi: sp.Expr
    kink_multiplicity: int
    signal_speed: sp.Expr
    gap_frequency: sp.Expr
    length: sp.Expr
    energy_scale: sp.Expr
    action_scale: sp.Expr
    inverse_g_normalized: sp.Expr
    q_dimensionless: sp.Expr

    def newton_constant(self, hbar: Any) -> sp.Expr:
        """Derived Newton constant ``G = q * ell**2 * c**3 / hbar`` (C-GRV-001 form)."""

        hb = _positive_exact(hbar, "hbar")
        return self.q_dimensionless * self.length**2 * self.signal_speed**3 / hb

    def einstein_kappa(self, hbar: Any) -> sp.Expr:
        """Derived C-STG-001 coupling ``kappa = 8*pi*G/c**4`` (named convention)."""

        return 8 * sp.pi * self.newton_constant(hbar) / self.signal_speed**4


def assemble_chain(
    inertia: Any,
    gradient: Any,
    onsite: Any,
    action_quantum: Any = 1,
    xi: Any = 0,
    kink_multiplicity: int = 2,
) -> SGGravityChain:
    """Compose upstream sG coefficients into the full gravity ledger."""

    coefficients = DimensionalSineGordonCoefficients(
        inertia=sp.sympify(inertia),
        gradient=sp.sympify(gradient),
        onsite=sp.sympify(onsite),
    )
    scales = dimensional_sine_gordon_scales(coefficients)
    inverse_g = tower_inverse_g_normalized(action_quantum, xi, kink_multiplicity)
    return SGGravityChain(
        coefficients=coefficients,
        action_quantum=sp.sympify(action_quantum),
        xi=sp.sympify(xi),
        kink_multiplicity=kink_multiplicity,
        signal_speed=scales.signal_speed,
        gap_frequency=scales.gap_frequency,
        length=scales.length,
        energy_scale=scales.energy,
        action_scale=scales.action,
        inverse_g_normalized=inverse_g,
        q_dimensionless=newton_q_normalized(action_quantum, xi, kink_multiplicity),
    )


# ---------------------------------------------------------------------------
# Static arm: vacuum exterior of the Einstein-scalar sector with derived G.
# ---------------------------------------------------------------------------


def schwarzschild_areal_radius_factor(mass_parameter: Any, radius: Any) -> sp.Expr:
    """Exact vacuum metric factor ``f = 1 - 2*m/x`` (C-STG-002 conventions).

    In the vacuum region the C-STG-002 reduced equations with vanishing
    scalar amplitude give the Schwarzschild family; ``m`` is the
    dimensionless geometric mass parameter ``mu_gap * M_geo`` and
    ``x = mu_gap * r`` the dimensionless areal radius.
    """

    m = sp.sympify(mass_parameter)
    x = _positive_exact(radius, "radius")
    return 1 - 2 * m / x


def newtonian_potential(mass: Any, radius: Any, newton_g: Any, speed: Any) -> sp.Expr:
    """Weak-field potential ``Phi = -G*M/r`` from the derived coupling."""

    G = _positive_exact(newton_g, "newton_g")
    M = _positive_exact(mass, "mass")
    r = _positive_exact(radius, "radius")
    sp.sympify(speed)  # speed enters the metric form, not the potential
    return -G * M / r


def light_deflection_angle(mass: Any, impact: Any, newton_g: Any, speed: Any) -> sp.Expr:
    """Standard weak-field deflection ``4*G*M/(c**2*b)`` with derived G.

    Imported standard result; its role is the static weak-field
    demonstration and it is evaluated only with the campaign's derived
    coupling.
    """

    G = _positive_exact(newton_g, "newton_g")
    M = _positive_exact(mass, "mass")
    b = _positive_exact(impact, "impact")
    c = _positive_exact(speed, "speed")
    return 4 * G * M / (c**2 * b)


# ---------------------------------------------------------------------------
# Radiative arm: breather quadrupole power with the derived coupling.
# ---------------------------------------------------------------------------


def breather_moment_third_physical_scale(chain: SGGravityChain, hbar: Any) -> sp.Expr:
    """Physical scale of the mass-moment third derivative ``d = mu'''``.

    The accepted normalized moment ``mu(t)`` (C-SG-009) is an energy second
    spatial moment in units ``E_scale * ell**2`` evolving in normalized time
    ``omega_0 * t``.  The C-GW-001 quadrupole convention uses the mass
    moment, so the physical third derivative is

    ``mu'''_phys = (E_scale/c**2) * ell**2 * omega_0**3 * mu'''_normalized``.
    """

    c = chain.signal_speed
    return (chain.energy_scale / c**2) * chain.length**2 * chain.gap_frequency**3


def breather_quadrupole_power(chain: SGGravityChain, hbar: Any) -> sp.Expr:
    """Cycle-averaged breather quadrupole power with the derived coupling.

    Accepted conditional form (C-GW-004): ``<P> = 2*G*<d**2>/15`` with
    ``d = mu'''``; the mass-moment conversion and the restored ``c**5``
    propagation factor of the C-GW-001 physical-units waveform give

    ``<P> = (2*G/(15*c**5)) * mu'''_phys_scale**2 * <mu'''_normalized**2>``

    evaluated at the C-SG-010 bracket midpoint ``379.46463806875``.
    The previously declared input ``A = 2*G`` is now discharged by the
    campaign's derived coupling.
    """

    G = _positive_exact(chain.newton_constant(hbar), "newton_g")
    c = chain.signal_speed
    scale = breather_moment_third_physical_scale(chain, hbar)
    return sp.Rational(2, 15) * G / c**5 * scale**2 * MU_THIRD_POWER_MID


def dimensionless_summary(chain: SGGravityChain) -> dict[str, sp.Expr]:
    """Pure-number summary: q, S, and per-arm dimensionless anchors."""

    return {
        "q": chain.q_dimensionless,
        "inverse_g_normalized": chain.inverse_g_normalized,
        "mu_third_power_midpoint": MU_THIRD_POWER_MID,
    }

"""Einstein dynamics from a declared local causal-horizon entropy density.

This module packages the coefficient interface of the local-horizon
Clausius derivation.  It assumes local Lorentz invariance, an area entropy
law, the Unruh temperature, equilibrium Clausius heat flow, stress
conservation, and the Raychaudhuri equation.  Those hypotheses imply the
nonlinear Einstein equation; they are not a perturbative graviton action.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import sympy as sp

from .exact_symbolic import exact_real as _exact_real
from .exact_symbolic import positive_exact as _positive_exact


@dataclass(frozen=True)
class LocalHorizonEinsteinLedger:
    """Coefficient ledger for the local-equilibrium Einstein equation.

    ``einstein_stress_coupling`` is the coefficient in
    ``G_ab + Lambda*g_ab = coupling*T_ab``.  Entropy is dimensionless, so
    ``entropy_area_density`` has inverse-area units and ``action_scale`` is
    the same action appearing in the Unruh temperature.
    """

    entropy_area_density: sp.Expr
    action_scale: sp.Expr
    signal_speed: sp.Expr
    cosmological_constant: sp.Expr
    unruh_energy_per_inverse_length: sp.Expr
    clausius_stress_coefficient: sp.Expr
    einstein_stress_coupling: sp.Expr
    newton_constant: sp.Expr
    radiative_tensor_polarizations: int


def local_horizon_einstein_ledger(
    entropy_area_density: Any,
    action_scale: Any,
    signal_speed: Any,
    *,
    cosmological_constant: Any = 0,
) -> LocalHorizonEinsteinLedger:
    r"""Derive the Einstein coefficient from local horizon equilibrium.

    For an affinely parameterized local horizon generator ``k``,
    Raychaudhuri gives ``delta A=-integral lambda*R_kk`` while the boost heat
    is ``delta Q=-kappa*integral lambda*T_kk``.  The declared entropy law
    ``delta S=eta*delta A`` and Unruh energy
    ``T_U=J*c*kappa/(2*pi)`` therefore give

    ``R_kk = 2*pi*T_kk/(J*c*eta)``.

    Requiring this for every null generator, then imposing stress
    conservation and the contracted Bianchi identity, yields
    ``G_ab+Lambda*g_ab=2*pi*T_ab/(J*c*eta)``.  Matching the standard
    coefficient ``8*pi*G/c**4`` fixes ``G=c**3/(4*J*eta)``.  ``Lambda`` is
    the integration constant left by the theorem, not a contribution to
    inverse Newton coupling.
    """

    entropy_density = _positive_exact(
        entropy_area_density,
        "entropy_area_density",
    )
    action = _positive_exact(action_scale, "action_scale")
    speed = _positive_exact(signal_speed, "signal_speed")
    cosmological = _exact_real(cosmological_constant, "cosmological_constant")
    unruh_coefficient = sp.simplify(action * speed / (2 * sp.pi))
    clausius_coefficient = sp.simplify(2 * sp.pi / (action * speed))
    einstein_coupling = sp.simplify(clausius_coefficient / entropy_density)
    newton = sp.simplify(speed**3 / (4 * action * entropy_density))
    if sp.simplify(einstein_coupling - 8 * sp.pi * newton / speed**4) != 0:
        raise AssertionError("local-horizon and Einstein coefficients disagree")
    return LocalHorizonEinsteinLedger(
        entropy_area_density=entropy_density,
        action_scale=action,
        signal_speed=speed,
        cosmological_constant=cosmological,
        unruh_energy_per_inverse_length=unruh_coefficient,
        clausius_stress_coefficient=clausius_coefficient,
        einstein_stress_coupling=einstein_coupling,
        newton_constant=newton,
        radiative_tensor_polarizations=2,
    )

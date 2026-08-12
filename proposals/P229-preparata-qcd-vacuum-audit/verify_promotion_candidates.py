#!/usr/bin/env python3
"""Exact verifier for P229's narrow C-CMV-001/C-CMV-003 candidates."""

from __future__ import annotations

import sympy as sp

from substrate_framework.chromomagnetic_background import (
    charged_vector_mass2,
    nielsen_olesen_imaginary_part,
    one_loop_log_coefficient,
    su2_chromomagnetic_beta_bridge,
)
from substrate_framework.chromomagnetic_sectors import (
    sector_heat_kernel_coefficients,
    su3_charged_roots,
)
from substrate_framework.chromomagnetic_two_loop import rg_improved_potential
from substrate_framework.gauge_beta import GaugeFactor, product_gauge_coefficients
from substrate_framework.su3 import fundamental_generators
from substrate_framework.verification import CheckLedger


def run() -> int:
    checks = CheckLedger("P229/narrow-promotion-candidates")
    bridge = su2_chromomagnetic_beta_bridge()
    expected_log = sp.Rational(11, 48) / sp.pi**2
    checks.check(
        "accepted Pauli-half and beta conventions close on the background log",
        bridge.adjoint_casimir == 2
        and bridge.beta_one_loop_gauge == sp.Rational(-22, 3)
        and bridge.background_log_coefficient == expected_log
        and bridge.beta_implied_log_coefficient == expected_log
        and bridge.rg_scale_residual == 0,
    )

    wrong_casimir_beta = product_gauge_coefficients(
        (GaugeFactor("wrong-SU2", sp.Integer(3)),)
    ).one_loop_gauge[0]
    checks.check(
        "wrong adjoint-Casimir mutation breaks the bridge",
        sp.simplify(-wrong_casimir_beta / (32 * sp.pi**2) - expected_log)
        != 0,
    )
    checks.check(
        "wrong logarithm scale weight breaks RG closure",
        sp.simplify(
            bridge.tree_scale_derivative_coefficient
            - bridge.background_log_coefficient
        )
        != 0,
    )

    coeffs = sector_heat_kernel_coefficients()
    checks.check(
        "longitudinal and ghost coefficients cancel exactly",
        coeffs["longitudinal_pair"] == -sp.Rational(1, 3)
        and coeffs["ghost"] == sp.Rational(1, 3)
        and coeffs["longitudinal_pair"] + coeffs["ghost"] == 0,
    )
    checks.check(
        "transverse sectors carry the complete 11/3 logarithm",
        coeffs["transverse_spin_up"] == sp.Rational(11, 6)
        and coeffs["transverse_spin_down"] == sp.Rational(11, 6)
        and sum(coeffs.values()) == sp.Rational(11, 3),
    )

    g_b = sp.Symbol("gB", positive=True)
    checks.check(
        "lowest spin-up mode is tachyonic",
        charged_vector_mass2(g_b, 0, 1) == -g_b,
    )
    checks.check(
        "one-loop imaginary-part magnitude is exact",
        sp.simplify(nielsen_olesen_imaginary_part(g_b) - g_b**2 / (8 * sp.pi))
        == 0,
    )
    checks.check(
        "background coefficient is the exact heat-kernel result",
        sp.simplify(one_loop_log_coefficient() - expected_log) == 0,
    )

    generators = fundamental_generators()
    pairs = ((0, 1), (0, 2), (1, 2))
    t3_roots = tuple(
        sp.simplify(abs(generators[2][i, i] - generators[2][j, j]))
        for i, j in pairs
    )
    t8_roots = tuple(
        sp.simplify(abs(generators[7][i, i] - generators[7][j, j]))
        for i, j in pairs
    )
    checks.check(
        "canonical T3 generator yields the three charged-root magnitudes",
        su3_charged_roots()
        == t3_roots
        == (sp.Integer(1), sp.Rational(1, 2), sp.Rational(1, 2)),
    )
    checks.check(
        "wrong T8 Cartan mutation changes the root pattern",
        t8_roots != t3_roots,
    )

    x = sp.Symbol("x", positive=True)
    derived_two_loop = rg_improved_potential("derived")
    checks.check(
        "two-loop running substitution is singular at the one-loop minimum",
        sp.limit(derived_two_loop, x, 1, dir="-") == sp.oo
        and sp.limit(derived_two_loop, x, 1, dir="+") == -sp.oo,
    )
    return checks.finish()


if __name__ == "__main__":
    run()

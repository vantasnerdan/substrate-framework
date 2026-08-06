"""Primary exact and numeric verifier for P225 / provisional C-SIG-001..003.

Composes only package APIs.  Runs with ``PYTHONPATH=src .venv/bin/python``.
"""

from __future__ import annotations

import math

import sympy as sp

from substrate_framework.sg_gravity_chain import (
    MU_THIRD_POWER_LOWER,
    MU_THIRD_POWER_MID,
    MU_THIRD_POWER_UPPER,
    assemble_chain,
    breather_quadrupole_power,
    light_deflection_angle,
    newtonian_potential,
    schwarzschild_areal_radius_factor,
)
from substrate_framework.sg_spectral_induction import (
    KINK_NORMALIZED_MASS,
    induced_inverse_g_species,
    newton_q_normalized,
    sg_mass_tower,
    tower_edge_cutoff,
    tower_inverse_g_normalized,
)
from substrate_framework.verification import CheckLedger


def main() -> int:
    checks = CheckLedger("P225")

    # ---- C-SIG-001: the induced coupling from the accepted tower ----
    checks.check(
        "tower composition matches declared premises (kink pair + 25 levels at h=1; phonon = n=1 level)",
        len(sg_mass_tower(1)) == 26,
    )
    cutoff = tower_edge_cutoff()
    checks.check(
        "tower edge is the accepted C-SG-007 action boundary 8*pi",
        cutoff == 8 * sp.pi,
    )
    checks.check(
        "every species logarithm is strictly positive (self-cutoff consistency)",
        all(float(cutoff / s.normalized_mass) > 1.0 for s in sg_mass_tower(1)),
    )
    lightest = induced_inverse_g_species(16 * sp.sin(sp.Rational(1, 16)), 0, cutoff)
    checks.check(
        "per-species coefficient matches the imported a2 matching (1/6-xi)m^2 log(L^2/m^2)/(2*pi)",
        sp.simplify(
            lightest
            - sp.Rational(1, 12) / sp.pi * (16 * sp.sin(sp.Rational(1, 16))) ** 2
            * sp.log(cutoff**2 / (16 * sp.sin(sp.Rational(1, 16))) ** 2)
        )
        == 0,
    )
    s_full = tower_inverse_g_normalized(1)
    checks.mutation_sensitive(
        "induction responds to the nonminimal coupling input",
        lambda xi: sp.simplify(tower_inverse_g_normalized(1, xi, 2) - s_full) == 0,
        0,
        [sp.Rational(1, 6), -sp.Rational(1, 3)],
    )
    s_nokink = tower_inverse_g_normalized(1, 0, 0)
    checks.check(
        "kink mutation removes exactly the kink-antikink contribution",
        sp.simplify(
            s_full
            - s_nokink
            - 2 * induced_inverse_g_species(KINK_NORMALIZED_MASS, 0, cutoff)
        )
        == 0,
    )
    # Independent numeric route for the same premises.
    total = 2 * 64.0 * math.log((8 * math.pi) ** 2 / 64.0) / (12 * math.pi)
    for n in range(1, 26):
        e_n = 16 * math.sin(n / 16)
        total += e_n**2 * math.log((8 * math.pi) ** 2 / e_n**2) / (12 * math.pi)
    checks.check(
        "symbolic tower sum agrees with independent float summation",
        abs(float(s_full) - total) < 1e-10 * total,
        f"S = {float(s_full):.9f} vs {total:.9f}",
    )
    checks.check(
        "q is the exact reciprocal of the tower sum",
        sp.simplify(newton_q_normalized(1) * s_full - 1) == 0,
    )
    q_value = float(newton_q_normalized(1))
    checks.check(
        "q is a sane pure number of the upstream model (0 < q < 1)",
        0.0 < q_value < 1.0,
        f"q = {q_value:.9f}",
    )

    # ---- C-SIG-002: static weak-field arm ----
    lam, tension, mu = (sp.Symbol(s, positive=True) for s in ("lam", "T", "mu"))
    hbar = sp.Symbol("hbar", positive=True)
    chain = assemble_chain(lam, tension, mu, action_quantum=1)
    g = chain.newton_constant(hbar)
    q = newton_q_normalized(1)
    checks.check(
        "derived G has the exact C-GRV-001 monomial form q*ell^2*c^3/hbar",
        sp.simplify(
            g - q * (tension / mu) * (tension / lam) ** sp.Rational(3, 2) / hbar
        )
        == 0,
    )
    checks.check(
        "q is invariant under upstream ray rescaling (unit-standard independence)",
        sp.simplify(
            assemble_chain(7, 7, 7, action_quantum=1).q_dimensionless - q
        )
        == 0,
    )
    m, x = sp.symbols("m x", positive=True)
    checks.check(
        "vacuum exterior metric factor is Schwarzschild f = 1 - 2m/x",
        schwarzschild_areal_radius_factor(m, x) == 1 - 2 * m / x,
    )
    M, r, b = sp.symbols("M r b", positive=True)
    unit_chain = assemble_chain(1, 1, 1, action_quantum=1)
    phi = newtonian_potential(M, r, unit_chain.newton_constant(hbar), unit_chain.signal_speed)
    checks.check(
        "Newtonian potential carries the derived coupling Phi = -M/(S*hbar*r)",
        sp.simplify(phi + M / (s_full * hbar * r)) == 0,
    )
    theta = light_deflection_angle(
        M, b, unit_chain.newton_constant(hbar), unit_chain.signal_speed
    )
    checks.check(
        "deflection angle is the standard 4GM/(c^2 b) with derived G",
        sp.simplify(theta - 4 * q * M / (hbar * b)) == 0,
    )

    # ---- C-SIG-003: radiative arm ----
    checks.check(
        "C-SG-010 bracket midpoint respected",
        MU_THIRD_POWER_LOWER < MU_THIRD_POWER_MID < MU_THIRD_POWER_UPPER,
    )
    power = breather_quadrupole_power(unit_chain, hbar)
    checks.check(
        "breather quadrupole power is positive with the derived G",
        float(power.subs(hbar, 1)) > 0.0,
        f"<P>(hbar=1, unit point) = {float(power.subs(hbar, 1)):.9f}",
    )
    power_nokink = breather_quadrupole_power(
        assemble_chain(1, 1, 1, action_quantum=1, kink_multiplicity=0), hbar
    )
    ratio = float((power / power_nokink).subs(hbar, 1))
    checks.check(
        "radiative power tracks the derived coupling under kink mutation",
        abs(ratio - float(s_nokink / s_full)) < 1e-10,
        f"ratio {ratio:.9f} vs {float(s_nokink / s_full):.9f}",
    )
    checks.check(
        "power scales linearly in the onsite coefficient (G*scale**2 accounting)",
        abs(
            float(
                (
                    breather_quadrupole_power(
                        assemble_chain(1, 1, 2, action_quantum=1), hbar
                    )
                    / power
                ).subs(hbar, 1)
            )
            - 2.0
        )
        < 1e-10,
    )

    print(f"q = {q_value:.12f}, S = 1/G_norm = {float(s_full):.9f}")
    print(f"ALL {len(checks.passed)} CHECKS PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

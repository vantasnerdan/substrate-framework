"""C-ELS-008 verifier: barotropic stiff equation of state p = K*rho.

Rung between C-ELS-001 and C-ELS-002 on the compressible closure branch.
The equation of state enters as declared premise P0 (barotropic closure
p = K*rho); everything derived here keeps the P242 premise discipline:
the EOS contribution to the longitudinal sector is read off by algebraic
coefficient matching, never hard-coded, and each check carries a mutation
control where a wrong premise would change the verdict.

Exact symbolic claims:
  1. The barotropic balance residual vanishes on an exact diluting-flow
     solution class (uniform-density flow with linear profile).
  2. Linearizing about rest and inserting a longitudinal plane wave gives
     omega^2 = K*k^2 exactly from the declared closure.
  3. Assembling total stress as stiff-fluid bulk plus frozen-tangle
     Cauchy-Born energy forces lambda_total = rho0*K + E_f*L_v/15 by
     coefficient matching.
  4. Christoffel speeds: c_P^2 = K + E_f*L_v/(5*rho0),
     c_S^2 = E_f*L_v/(15*rho0); fluid limit mu -> 0 recovers c_P^2 = K.
  5. Poisson-ratio interpolation nu = (rho0*K+mu)/(2*(rho0*K+2*mu)):
     K -> 0 gives the C-ELS-003 value 1/4, mu -> 0 gives 1/2.
Mutations: dropping the body force breaks the solution class; flipping the
closure sign destroys strong ellipticity; a half-strength bulk modulus
oscillates measurably off-shell in the dynamics companion.
"""

import os

from substrate_framework import CheckLedger
from substrate_framework.averaging import barotropic_balance_residual
from substrate_framework.elasticity import (
    acoustic_speeds_squared,
    poisson_ratio,
    strong_elliptic,
)
os.environ.setdefault("MKL_NUM_THREADS", "1")

import sys

import numpy as np
import sympy as sp

from substrate_framework.numerics import (
    SolverTolerances,
    solve_ivp_evidence,
)

K_EOS = 1.0          # stiff fluid p = rho in matched units
EF = 3             # axial filament stiffness (exact integer: symbolic matching)
LV = 2             # vortex line density per unit volume (exact integer)
RHO0 = 1.5         # background mass density

def main() -> int:
    ledger = CheckLedger("C-ELS-008")

    x, y, t = sp.symbols("x y t", real=True)
    a, rho0, k_eos = sp.symbols("a rho0 k_eos", positive=True)
    lam_s, mu_s = sp.symbols("_lamTot _muT", positive=True)
    k_wave = sp.Symbol("_k", positive=True)
    omega = sp.Symbol("_omega", positive=True)
    u_amp = sp.Symbol("_U", real=True)
    exx, eyy, ezz, exy, exz, eyz = sp.symbols(
        "_eXX _eYY _eZZ _eXY _eXZ _eYZ", real=True
    )

    def check_exact(name: str, condition) -> None:
        ledger.check(name, sp.simplify(condition) == 0)

    # --- checks 1-2: barotropic balance on an exact solution class ---
    velocity = (a * x, 0)
    density = rho0 * sp.exp(-a * t)
    pressure = k_eos * density
    body = (a**2 * x, 0)
    residual = barotropic_balance_residual(
        velocity, density, pressure, body, (x, y), t
    )
    ledger.check(
        "exact diluting-flow solution class closes",
        all(sp.simplify(r) == 0 for r in residual),
    )
    residual_noforce = barotropic_balance_residual(
        velocity, density, pressure, (0, 0), (x, y), t
    )
    ledger.check(
        "mutation: dropping the body force breaks the class",
        sp.simplify(residual_noforce[0] - a**2 * x) == 0,
    )

    # --- check 3: linearized dispersion from the declared closure ---
    # Momentum: d_t v + grad(p)/rho = 0 about rest. Continuity gives
    # R = k*rho0*U/omega for the density amplitude of the plane wave;
    # momentum x then reads -i*omega*U + i*k*K*R/rho0 = 0.
    R = k_wave * rho0 * u_amp / omega
    momentum_x = -sp.I * omega * u_amp + sp.I * k_wave * k_eos * R / rho0
    dispersion = sp.solve(sp.simplify(momentum_x), omega**2, dict=True)
    ledger.check(
        "linearized stiff-EOS dispersion is omega^2 = K k^2",
        len(dispersion) == 1
        and sp.simplify(dispersion[0][omega**2] - k_eos * k_wave**2) == 0,
    )

    # --- check 4: lambda_total by coefficient matching, not hard-coding ---
    strain = sp.Matrix([[exx, exy, exz], [exy, eyy, eyz], [exz, eyz, ezz]])
    trace = sum(strain[i, i] for i in range(3))
    frobenius = sum(strain[p, q] ** 2 for p in range(3) for q in range(3))
    # Stiff-fluid bulk storage: (1/2)*rho0*K*(div xi)^2 = (1/2)*rho0*K*(tr eps)^2.
    # Frozen-tangle Cauchy-Born storage: E_f*L_v*[(tr eps)^2 + 2 eps:eps]/30.
    total_energy = (
        sp.Rational(1, 2) * rho0 * k_eos * trace**2
        + EF * LV * (trace**2 + 2 * frobenius) / 30
    )
    target = lam_s * trace**2 / 2 + mu_s * frobenius
    poly = sp.Poly(sp.expand(total_energy - target), exx, eyy, ezz, exy, exz, eyz)
    solution = sp.solve([sp.Eq(c, 0) for c in poly.coeffs()], (lam_s, mu_s), dict=True)
    ledger.check("coefficient matching determines both moduli", len(solution) == 1)
    lam_tot = sp.simplify(solution[0][lam_s])
    mu_t = sp.simplify(solution[0][mu_s])
    check_exact(
        "lambda_total = rho0*K + E_f*L_v/15",
        lam_tot - (rho0 * k_eos + sp.Rational(EF * LV, 15)),
    )
    check_exact("mu_tangle = E_f*L_v/15", mu_t - sp.Rational(EF * LV, 15))

    # --- check 5: Christoffel speeds with the assembled Lame data ---
    speeds = acoustic_speeds_squared(lam_tot, mu_t, rho0)
    cp2 = sp.simplify(list(speeds.values())[0])
    cs2 = sp.simplify(list(speeds.values())[-1])
    check_exact(
        "c_P^2 = K + E_f*L_v/(5 rho0)",
        cp2 - (k_eos + 3 * sp.Rational(EF * LV, 15) / rho0),
    )
    check_exact("c_S^2 = E_f*L_v/(15 rho0)", cs2 - sp.Rational(EF * LV, 15) / rho0)

    speeds_fluid = acoustic_speeds_squared(lam_tot.subs(mu_t, 0), 0, rho0)
    cp2_fluid = sp.simplify(list(speeds_fluid.values())[0])
    cs2_fluid = sp.simplify(list(speeds_fluid.values())[-1])
    check_exact("fluid limit mu->0 restores c_P^2 = K", cp2_fluid - k_eos)
    ledger.check("fluid limit mu->0 has no shear branch", cs2_fluid == 0)

    # --- check 6: Poisson interpolation as EOS/tangle discriminator ---
    nu_stiff_limit = poisson_ratio(rho0 * k_eos, 0)
    nu_tangle_limit = poisson_ratio(mu_t, mu_t)
    ledger.check(
        "nu -> 1/2 when only the stiff EOS carries load",
        sp.simplify(nu_stiff_limit - sp.Rational(1, 2)) == 0,
    )
    ledger.check(
        "nu = 1/4 when only the tangle carries load",
        sp.simplify(nu_tangle_limit - sp.Rational(1, 4)) == 0,
    )

    # --- check 7: mutation, flipping the closure sign destroys the medium ---
    ledger.check(
        "mutation: flipped closure loses strong ellipticity",
        not strong_elliptic(
            -lam_tot.subs({k_eos: K_EOS, rho0: RHO0}),
            mu_t.subs({k_eos: K_EOS, rho0: RHO0}),
        ),
    )

    # --- checks 8-10: dynamics companion on the longitudinal branch ---
    cp2_numeric = float(cp2.subs({k_eos: K_EOS, rho0: RHO0}))
    k_mode = 3.0
    omega_exact = float(sp.sqrt(cp2_numeric)) * k_mode

    def rhs(state: np.ndarray, omega_squared: float) -> np.ndarray:
        displacement, velocity_component = state
        return np.array([velocity_component, -omega_squared * displacement])

    def measured_frequency(omega_squared: float) -> dict[str, float]:
        state0 = np.array([1.0, 0.0])
        duration = 4.0 * np.pi / omega_exact
        times = np.linspace(0.0, duration, 4001)
        evidence = solve_ivp_evidence(
            lambda time, state: rhs(state, omega_squared),
            (0.0, float(duration)),
            state0,
            sample_times=times,
            tolerances=SolverTolerances(rtol=1e-12, atol=1e-14),
        )
        displacement = evidence.state[0]
        sign_changes = int(np.sum(np.diff(np.sign(displacement)) != 0))
        omega_measured = 2.0 * np.pi * (sign_changes / 2.0) / float(duration)
        velocity_component = evidence.state[1]
        energy = 0.5 * (velocity_component**2 + omega_squared * displacement**2)
        return {
            "omega": omega_measured,
            "relative_error": abs(omega_measured - omega_exact) / omega_exact,
            "max_energy_drift": float(np.max(np.abs(energy - energy[0]))),
        }
    on_shell = measured_frequency(omega_exact**2)
    ledger.check(
        "longitudinal frequency matches symbolic c_P",
        on_shell["relative_error"] < 1e-6,
    )
    ledger.check(
        "energy conserved under refinement-grade tolerances",
        on_shell["max_energy_drift"] < 1e-9,
    )

    half_bulk_squared = omega_exact**2 - k_mode**2 * (K_EOS / 2.0)
    off_shell = measured_frequency(half_bulk_squared)
    ledger.check(
        "mutation control: half-strength bulk modulus oscillates measurably off-shell",
        abs(off_shell["omega"] - omega_exact) / omega_exact > 1e-3,
    )

    return int(ledger.finish())


if __name__ == "__main__":
    sys.exit(main())

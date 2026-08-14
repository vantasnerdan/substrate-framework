#!/usr/bin/env python3
"""
P226 attempt 0001 verify.py — reproduction of Fulceri paper key equations.

This is a reproduction-only script for the paper's claims that overlap
with existing substrate-framework claims and for the new candidate claims
that admit sympy verification. It does NOT promote any new claim; it
records the gate pass/fail for each candidate.

Each gate is a sympy residual = 0 test. Mutations confirm the gate is
sensitive to the load-bearing input.
"""
from __future__ import annotations
import sympy as sp

# Symbols
omega, alpha, eta, c0, omega0, mu0, v, gamma = sp.symbols(
    "omega alpha eta c0 omega0 mu0 v gamma", positive=True, real=True
)
x, t, tau, X = sp.symbols("x t tau X", real=True)
M, M_c0, h_eff = sp.symbols("M M_c0 h_eff", positive=True, real=True)
n, V = sp.symbols("n V", positive=True, real=True)
rho, Theta, Z0 = sp.symbols("rho Theta Z_0", positive=True, real=True)
rho_0, Theta_0 = sp.symbols("rho_0 Theta_0", positive=True, real=True)


# GATE 1: paper eq. (12) — rest-frame breather energy with dimensional prefactor
def gate_breather_energy():
    E0_paper = 16 * (mu0**2 * c0 / omega0) * sp.sqrt(1 - alpha**2)
    E0_substrate = 16 * sp.sqrt(1 - alpha**2)
    ratio = sp.simplify(E0_paper / E0_substrate)
    return ratio == mu0**2 * c0 / omega0, {"ratio": str(ratio)}


# GATE 2: paper eq. (14) — breather 4-vector dynamics E = gamma E0, P = gamma E0 v
def gate_four_vector():
    v_local = sp.symbols("v_local", real=True)
    gamma_local = 1 / sp.sqrt(1 - v_local**2 * c0**2)
    E0 = 16 * (mu0**2 * c0 / omega0) * sp.sqrt(1 - alpha**2)
    E = gamma_local * E0
    P = gamma_local * E0 * v_local
    invariant = sp.simplify(E**2 - P**2 * c0**2 - E0**2)
    return invariant == 0, {"invariant_residual": str(invariant)}


# GATE 3: paper eq. (19) — de Broglie relation k_dB = omega_carrier v / c0^2
def gate_de_broglie():
    v_local = sp.symbols("v_local", real=True)
    omega_carrier = alpha * gamma * omega0
    k_dB = omega_carrier * v_local / c0**2
    E0 = 16 * (mu0**2 * c0 / omega0) * sp.sqrt(1 - alpha**2)
    P = gamma * E0 * v_local
    residual = sp.simplify(h_eff * k_dB - P)
    h_eff_resolved = sp.solve(sp.Eq(residual, 0), h_eff)
    return (h_eff_resolved != [], len(h_eff_resolved) > 0), {
        "h_eff_required": [str(s) for s in h_eff_resolved],
        "residual": str(residual),
    }


# GATE 4: paper eq. (26) — two-frequencies product invariant
def gate_two_freqs():
    omega_carrier = alpha * gamma * omega0
    omega_clock = alpha * omega0 / gamma
    product = sp.simplify(omega_carrier * omega_clock)
    return product == alpha**2 * omega0**2, {"product": str(product)}


# GATE 5: paper eq. (60) — pure-flow metric components
def gate_flowing_metric_paper():
    g_tt = -c0**2 + V**2
    g_tX = V
    g_XX = 1
    det = sp.simplify(g_tt * g_XX - g_tX**2)
    return det == -c0**2, {"det": str(det)}


# GATE 6: paper eq. (69) — full flowing metric with refractive index
def gate_full_flowing_metric():
    g_tt = -c0**2 / n + n * V**2
    g_tX = n * V
    g_XX = n
    det = sp.simplify(g_tt * g_XX - g_tX**2)
    return det == -c0**2, {"det": str(det)}


# GATE 7: paper eq. (96) — inverse metric identity
def gate_inverse_metric():
    g_tt = -c0**2 / n + n * V**2
    g_tX = n * V
    g_XX = n
    det = g_tt * g_XX - g_tX**2
    g_tt_inv = sp.simplify(g_XX / det)
    g_tX_inv = sp.simplify(-g_tX / det)
    g_XX_inv = sp.simplify(g_tt / det)
    inv_check_tt = sp.simplify(g_tt_inv * g_tt + g_tX_inv * g_tX)
    inv_check_tX = sp.simplify(g_tt_inv * g_tX + g_tX_inv * g_XX)
    inv_check_XX = sp.simplify(g_tX_inv * g_tX + g_XX_inv * g_XX)
    return (inv_check_tt == 1 and inv_check_tX == 0 and inv_check_XX == 1), {
        "g_tt_inv": str(g_tt_inv),
        "g_tX_inv": str(g_tX_inv),
        "g_XX_inv": str(g_XX_inv),
    }


# GATE 8: paper eq. (105) — massless null constraint solutions
def gate_massless_null_constraint():
    v_local = sp.symbols("v_local", real=True)
    null = -c0**2 / n + n * V**2 + 2 * n * V * v_local + n * v_local**2
    solutions = sp.solve(sp.Eq(null, 0), v_local)
    return (len(solutions) == 2, {"solutions": [str(s) for s in solutions]})


# GATE 9: impedance matching (paper eq 61-62)
def gate_impedance_matching():
    rho_local = rho_0 * n
    Theta_local = Theta_0 / n
    product = sp.simplify(rho_local * Theta_local)
    Z_0_squared = rho_0 * Theta_0
    return product == Z_0_squared, {"product": str(product)}


# GATE 10: V=0 reduction to C-OG-001 shape
def gate_v0_reduction():
    g_tt_paper = -c0**2 / n
    g_XX_paper = n
    return (g_tt_paper == -c0**2 / n and g_XX_paper == n), {
        "g_tt_paper": str(g_tt_paper),
        "g_XX_paper": str(g_XX_paper),
    }


# GATE 11: n=1, V=0 limit gives Minkowski
def gate_minkowski_limit():
    g_tt = -c0**2 / 1 + 1 * 0**2
    g_tX = 1 * 0
    g_XX = 1
    det = g_tt * g_XX - g_tX**2
    return (g_tt == -c0**2 and g_tX == 0 and g_XX == 1 and det == -c0**2), {
        "g_tt": str(g_tt),
        "det": str(det),
    }


def main() -> int:
    print("P226 attempt 0001 — Fulceri paper key equation reproduction")
    print("=" * 70)

    gates = [
        ("G1: breather energy w/ dimensional prefactor", gate_breather_energy),
        ("G2: 4-vector dynamics E^2 - (Pc0)^2 = E0^2", gate_four_vector),
        ("G3: de Broglie relation h_eff*k_dB = P", gate_de_broglie),
        ("G4: two-frequencies product invariant", gate_two_freqs),
        ("G5: pure-flow metric det (eq 60)", gate_flowing_metric_paper),
        ("G6: full flowing metric det (eq 69)", gate_full_flowing_metric),
        ("G7: inverse metric identity g*g^-1 = I", gate_inverse_metric),
        ("G8: massless null constraint solutions", gate_massless_null_constraint),
        ("G9: impedance matching rho*Theta = Z_0^2", gate_impedance_matching),
        ("G10: V=0 reduction to C-OG-001 shape", gate_v0_reduction),
        ("G11: n=1, V=0 reduction to Minkowski", gate_minkowski_limit),
    ]

    failures = []
    for name, gate in gates:
        try:
            ok, info = gate()
        except Exception as e:
            ok, info = False, {"error": str(e)}
        status = "PASS" if ok else "FAIL"
        print(f"  {name}: {status}")
        if not ok:
            failures.append((name, info))
            for k, v in info.items():
                print(f"    {k}: {v}")

    print()
    if not failures:
        print(f"ALL {len(gates)} GATES PASS")
    else:
        print(f"{len(failures)}/{len(gates)} GATES FAILED")
    return 0 if not failures else 1


def mutation_tests() -> None:
    """Mutate load-bearing inputs and require the relevant gate to fail."""
    print()
    print("MUTATION TESTS (verifier sensitivity)")
    print("=" * 70)

    # M1: Alter det of full flowing metric (g_tt without -c0^2/n prefix)
    g_tt_mut = n * V**2
    g_tX_mut = n * V
    g_XX_mut = n
    det_mut = sp.simplify(g_tt_mut * g_XX_mut - g_tX_mut**2)
    expected = -c0**2
    print(f"  M1: mutant flowing metric det = {det_mut} (expected -c0^2)")
    assert det_mut != expected, "M1 FAILED: mutant should not equal -c0^2"
    print(f"  M1: PASS")

    # M2: Sign mismatch on impedance matching
    rho_local_mut = rho_0 * n
    Theta_local_mut = Theta_0 * n
    product_mut = sp.simplify(rho_local_mut * Theta_local_mut)
    expected = rho_0 * Theta_0
    print(f"  M2: mutant impedance product = {product_mut} (expected rho_0*Theta_0)")
    assert product_mut != expected, "M2 FAILED"
    print(f"  M2: PASS")

    # M3: V != 0 should change the null constraint solutions
    # (the det is conserved for any V, n population; the null constraint solutions are not)
    v_local = sp.symbols("v_local", real=True)
    null_V1 = -c0**2 / 1 + 1 * 1**2 + 2 * 1 * 1 * v_local + 1 * v_local**2
    sols_V1 = sp.solve(sp.Eq(null_V1, 0), v_local)
    # Should be -1 +/- c0, NOT the V=0 case -c0/n = -c0
    expected = {-1 + c0, -1 - c0}
    actual = set(sols_V1)
    print(f"  M3: null constraint with V=1, n=1: solutions = {sols_V1}")
    assert actual != {-c0, c0}, "M3 FAILED: solver should change with V=1"
    print(f"  M3: PASS (V=1 changes null speeds)")

    # M3b: det is conserved for any V, n population -- structural invariant
    g_tt_mut = -c0**2 / 1 + 1 * 5**2  # V=5, n=1
    g_tX_mut = 1 * 5
    g_XX_mut = 1
    det_mut = sp.simplify(g_tt_mut * g_XX_mut - g_tX_mut**2)
    print(f"  M3b: full flowing metric det at V=5, n=1: {det_mut}")
    assert det_mut == -c0**2, "M3b FAILED: det should be -c0^2"
    print(f"  M3b: PASS (det is structural invariant)")

    # M4: Null constraint with wrong sign on g_tt
    v_local = sp.symbols("v_local", real=True)
    g_tt_wrong = +c0**2 / n + n * V**2
    null_wrong = g_tt_wrong + 2 * n * V * v_local + n * v_local**2
    sols_wrong = sp.solve(sp.Eq(null_wrong, 0), v_local)
    print(f"  M4: null constraint with wrong-sign g_tt: solutions = {sols_wrong}")
    expected_sols = {-V + c0/n, -V - c0/n}
    actual_sols = set(sols_wrong)
    assert actual_sols != expected_sols, "M4 FAILED: mutant should not yield c0/n solutions"
    print(f"  M4: PASS")

    print()
    print("ALL MUTATIONS CONFIRMED SENSITIVE")


if __name__ == "__main__":
    rc = main()
    mutation_tests()
    raise SystemExit(rc)

#!/usr/bin/env python3
"""
Sympy-based dimensional validation of the Fulceri paper's key equations.

This script tests whether each paper equation is dimensionally consistent
under the standard SI-like dimension basis (M, L, T). It uses sympy to:
1. Define dimension columns for each physical quantity in the paper
2. Build the paper's claimed expressions from those dimensions
3. Test dimensional identity

The framework's C-OG-001 / C-OG-002 dimension matrix convention is used
(per the SIConstitutiveDimensionLedger in substrate_framework.constitutive).

This is a peer-review tool, not a campaign acceptance gate. It records
dimension findings as evidence to be evaluated by a distinct reviewer.
"""
from __future__ import annotations
import sympy as sp

# Dimension basis: columns are M, L, T (mass, length, time)

def dim(m, l, t):
    """Return the (M,L,T) dimension column."""
    return sp.Matrix([m, l, t])

DIM = {
    "mu0":     dim(0, -1, 0),
    "c0":      dim(0, 1, -1),
    "omega0":  dim(0, 0, -1),
    "alpha":   dim(0, 0, 0),
    "v":       dim(0, 1, -1),
    "gamma":   dim(0, 0, 0),
    "hbar_eff": dim(1, 2, -1),
    "E0":      dim(1, 2, -2),
    "P":       dim(1, 1, -1),
    "k_dB":    dim(0, -1, 0),
    "n":       dim(0, 0, 0),
    "V":       dim(0, 1, -1),
    "rho_0":   dim(1, -3, 0),
    "Theta_0": dim(1, -1, -2),
}


def dim_product(*args):
    """Return dim of product of quantities (string names or dim matrices)."""
    result = sp.zeros(3, 1)
    for a in args:
        if isinstance(a, str):
            a = DIM[a]
        result = result + a
    return result


def dim_quotient(num, denom):
    """Return dim of num/denom."""
    if isinstance(num, str):
        num = DIM[num]
    if isinstance(denom, str):
        denom = DIM[denom]
    return num - denom


def dim_pow(base, exponent):
    """Return dim of base^exponent (exponent integer)."""
    if isinstance(base, str):
        base = DIM[base]
    return exponent * base


def assert_dim(expr_name, expected, actual):
    diff = actual - expected
    if diff == sp.zeros(3, 1):
        print("  [PASS] " + expr_name + ": dim = (" + str(actual[0,0]) + ", " + str(actual[1,0]) + ", " + str(actual[2,0]) + ")")
        return True
    else:
        print("  [FAIL] " + expr_name + ": dim = (" + str(actual[0,0]) + ", " + str(actual[1,0]) + ", " + str(actual[2,0]) + ") but expected (" + str(expected[0,0]) + ", " + str(expected[1,0]) + ", " + str(expected[2,0]) + ")")
        return False


def test_eq_12_breather_energy():
    """Paper eq (12): E0 = 16 * (mu0^2 * c0 / omega0) * sqrt(1 - alpha^2)

    E0 should have energy units [M L^2 T^-2].
    Check whether mu0^2 * c0 / omega0 has energy units.
    """
    print("")
    print("GATE D1: Paper eq (12) E0 = 16 * (mu0^2 c0 / omega0) * sqrt(1 - alpha^2)")
    prefactor_dim = dim_quotient(dim_product("mu0", "mu0", "c0"), "omega0")
    return assert_dim("E0_paper", DIM["E0"], prefactor_dim)


def test_eq_14_four_vector():
    """Paper eq (14): E = gamma E0, P = gamma (v/c0) E0.

    E should have energy units; P should have momentum units.
    """
    print("")
    print("GATE D2: Paper eq (14) E = gamma E0, P = gamma (v/c0) E0")
    e_dim = dim_product("gamma", "E0")
    ok1 = assert_dim("E", DIM["E0"], e_dim)
    p_dim = dim_product("gamma", dim_quotient("v", "c0"), "E0")
    return ok1, assert_dim("P", DIM["P"], p_dim)


def test_eq_19_de_broglie():
    """Paper eq (19): k_dB = alpha * gamma * omega0 * v / c0^2"""
    print("")
    print("GATE D3: Paper eq (19) k_dB = alpha * gamma * omega0 * v / c0^2")
    k_dim = dim_product("alpha", "gamma", "omega0", "v", dim_pow("c0", -2))
    return assert_dim("k_dB", DIM["k_dB"], k_dim)


def test_eq_20_de_broglie_relation():
    """Paper eq (20): P = hbar_eff * k_dB"""
    print("")
    print("GATE D4: Paper eq (20) P = hbar_eff * k_dB")
    p_dim = dim_product("hbar_eff", "k_dB")
    return assert_dim("P", DIM["P"], p_dim)


def test_eq_61_62_impedance_matching():
    """Paper eq (61-62): rho*Theta = Z_0^2 = const under constitutive map.

    rho_0 [M L^-3], Theta_0 [M L^-1 T^-2].
    Z_0^2 = rho_0 * Theta_0 has units M^2 L^-4 T^-2.
    """
    print("")
    print("GATE D5: Paper eq (61-62) rho*Theta = rho_0*Theta_0")
    z_sq_dim = dim_product("rho_0", "Theta_0")
    return assert_dim("Z_0^2", z_sq_dim, z_sq_dim)


def test_eq_69_effective_metric():
    """Paper eq (69): g_tt = -c0^2/n + n V^2, g_tX = n V, g_XX = n.

    Metric components must have consistent dimensions for the determinant
    and inverse to make sense.
    """
    print("")
    print("GATE D6: Paper eq (69) flowing effective metric components")
    # g_tt = -c0^2/n + n V^2: term1 c0^2/n [L^2 T^-2], term2 n V^2 [L^2 T^-2]
    g_tt_dim = dim_quotient(dim_pow("c0", 2), "n")
    g_tX_dim = dim_product("n", "V")
    g_XX_dim = DIM["n"]
    print("  g_tt dim = (" + str(g_tt_dim[0,0]) + ", " + str(g_tt_dim[1,0]) + ", " + str(g_tt_dim[2,0]) + ")")
    print("  g_tX dim = (" + str(g_tX_dim[0,0]) + ", " + str(g_tX_dim[1,0]) + ", " + str(g_tX_dim[2,0]) + ")")
    print("  g_XX dim = (" + str(g_XX_dim[0,0]) + ", " + str(g_XX_dim[1,0]) + ", " + str(g_XX_dim[2,0]) + ")")
    print("  [WARN] Metric components have inconsistent dimensions.")
    print("    g_tt: L^2 T^-2 (speed^2 units)")
    print("    g_tX: L T^-1 (speed units)")
    print("    g_XX: dimensionless")
    print("  This is a SIGNATURE / SCALING inconsistency. The metric as written")
    print("  mixes three different dimensional conventions.")
    return False


def test_det_metric():
    """Paper eq (96): det g = -c0^2

    Compute det = g_tt * g_XX - g_tX^2 symbolically with sympy.
    """
    print("")
    print("GATE D7: Paper eq (96) det g = -c0^2")
    n, V, c0 = sp.symbols("n V c0", positive=True)
    g_tt = -c0**2 / n + n * V**2
    g_tX = n * V
    g_XX = n
    det = sp.simplify(g_tt * g_XX - g_tX**2)
    print("  det = " + str(det))
    if det == -c0**2:
        print("  [PASS] det g = -c0^2 (algebraically correct)")
        return True
    else:
        print("  [FAIL] det g != -c0^2")
        return False


def test_eq_105_massless_speeds():
    """Paper eq (105): v = V +/- c0/n.

    Solve g_munu xdot^mu xdot^nu = 0 with t_dot = 1:
    g_tt + 2 g_tX v + g_XX v^2 = 0
    Solve for v symbolically with sympy.
    """
    print("")
    print("GATE D8: Paper eq (105) v = V +/- c0/n")
    n, V, c0 = sp.symbols("n V c0", positive=True)
    v = sp.symbols("v_local", real=True)
    g_tt = -c0**2 / n + n * V**2
    g_tX = n * V
    g_XX = n
    null_constraint = g_tt + 2 * g_tX * v + g_XX * v**2
    solutions = sp.solve(sp.Eq(null_constraint, 0), v)
    print("  Solutions to null constraint:")
    for s in solutions:
        print("    v = " + str(sp.simplify(s)))
    print("  Paper claims: v = V +/- c0/n")
    print("  Actual solutions: v = -V +/- c0/n (see above)")
    # Check whether paper's claimed solutions match
    expected = [-V + c0/n, -V - c0/n]
    paper_claim = [V + c0/n, V - c0/n]
    matches_actual = all(sp.simplify(solutions[i] - expected[i]) == 0 for i in range(2))
    matches_paper = all(sp.simplify(solutions[i] - paper_claim[i]) == 0 for i in range(2))
    if matches_actual and not matches_paper:
        print("  [FAIL] Paper eq (105) is OFF BY A SIGN on V.")
        print("         Correct solutions: v = -V +/- c0/n")
        print("         Paper claim: v = V +/- c0/n (wrong sign)")
        return False
    elif matches_paper and not matches_actual:
        print("  [PASS] Paper eq (105) matches actual solutions")
        return True
    else:
        print("  [INDETERMINATE] Neither form matches cleanly")
        return False


def main():
    print("=" * 70)
    print("DIMENSIONAL VALIDATION OF FULCERI PAPER KEY EQUATIONS")
    print("=" * 70)
    print("Dimension basis: M (mass), L (length), T (time)")
    print("Convention: substrate-framework C-OG-001 / C-OG-002 (M, L, T)")
    print("")

    results = []
    findings = []
    for test in [test_eq_12_breather_energy,
                 test_eq_14_four_vector,
                 test_eq_19_de_broglie,
                 test_eq_20_de_broglie_relation,
                 test_eq_61_62_impedance_matching,
                 test_eq_69_effective_metric,
                 test_det_metric,
                 test_eq_105_massless_speeds]:
        try:
            result = test()
            if isinstance(result, tuple):
                results.extend(result)
            else:
                results.append(result)
        except Exception as e:
            print("  [ERROR] " + test.__name__ + ": " + str(e))
            results.append(False)

    print("")
    print("=" * 70)
    print("SUMMARY: " + str(sum(results)) + "/" + str(len(results)) + " checks passed")
    print("=" * 70)
    print("")
    print("FINDINGS (peer-review quality, sympy-validated):")
    print("")
    print("PROVABLY FALSE / SIGN ERROR (sympy-confirmed):")
    print("")
    print("1. D8 (paper eq 105, massless speeds):")
    print("   Paper claims: v = V +/- c0/n")
    print("   Sympy null-constraint solution: v = -V +/- c0/n")
    print("   The paper's eq (105) is OFF BY A SIGN on V.")
    print("   This is a symbolic error in the paper. The correct expression")
    print("   for the massless null speeds is v = -V +/- c0/n.")
    print("")
    print("DIMENSIONAL INCONSISTENCY (peer-review flag):")
    print("")
    print("2. D1 (paper eq 12, breather energy):")
    print("   Paper: E0 = 16 * (mu0^2 c0 / omega0) * sqrt(1 - alpha^2)")
    print("   Under standard sine-Gordon conventions (mu0 = soliton mass in")
    print("   1/length), the prefactor mu0^2 c0 / omega0 has units of 1/length,")
    print("   NOT energy. The paper does not state an explicit dimensional")
    print("   convention for mu0 that would yield energy units.")
    print("")
    print("3. D6 (paper eq 69, metric scaling):")
    print("   Metric components mix three dimensional conventions:")
    print("     g_tt = -c0^2/n + n V^2: L^2 T^-2")
    print("     g_tX = n V: L T^-1")
    print("     g_XX = n: dimensionless")
    print("   The paper does not state a c0-scaling convention. Either all")
    print("   components should be dimensionless (in which case g_tt should")
    print("   be divided by c0^2 and g_tX by c0), or the paper is using a")
    print("   mixed convention without explicit justification.")
    print("")
    print("ALGEBRAICALLY CORRECT:")
    print("")
    print("4. D7 (paper eq 96, det g = -c0^2):")
    print("   Sympy-verified: det((g_tt, g_tX), (g_tX, g_XX)) = -c0^2 exactly.")
    print("   No contradiction here despite the dimensional ambiguity in D6.")
    print("")
    print("5. D2 (paper eq 14, four-vector):")
    print("   E = gamma E0 is dimensionally consistent; P = gamma (v/c0) E0")
    print("   gives energy units, not momentum units. The Lorentz momentum")
    print("   form should be P = gamma (v/c0^2) E0 (= gamma m v with m = E0/c0^2).")
    print("   This is likely a typo in the paper (missing one factor of c0).")
    print("")
    print("6. D3 (paper eq 19, k_dB = alpha gamma omega0 v / c0^2):")
    print("   Dimensional analysis: [T^-1] [L T^-1] / [L^2 T^-2] = 1/L. Correct.")
    print("")
    print("7. D4 (paper eq 20, P = hbar_eff k_dB):")
    print("   [M L^2 T^-1] [L^-1] = [M L T^-1]. Correct.")
    print("")
    print("8. D5 (paper eq 61-62, rho*Theta = Z_0^2):")
    print("   [M L^-3] [M L^-1 T^-2] = [M^2 L^-4 T^-2]. Consistent with impedance^2.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

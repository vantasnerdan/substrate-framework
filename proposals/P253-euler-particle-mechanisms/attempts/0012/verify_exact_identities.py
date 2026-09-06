"""Exact algebra checks for the P253/0012 nearby-pair construction."""

import sympy as sp


def main() -> None:
    S, c = sp.symbols("S c", real=True)
    E10, E20, E1t, E2t = sp.symbols("E10 E20 E1t E2t", real=True)
    X0, Xt, PT = sp.symbols("E12_0 E12_t P_T", real=True)

    delta0 = 2 * S - E10 - E20 + c * PT
    delta_t = 2 * S - E1t - E2t + c * PT
    energy_substitution = {E1t: E10 + E20 + X0 - Xt - E2t}
    assert sp.simplify(delta_t.subs(energy_substitution) - delta0 - Xt + X0) == 0

    P1, P2 = sp.symbols("P1 P2", nonnegative=True)
    product_gap = (P1 + P2) ** 2 / 4 - P1 * P2
    assert sp.simplify(product_gap - (P1 - P2) ** 2 / 4) == 0

    CG, g = sp.symbols("C_G g", positive=True)
    half_gap_bound = CG * PT**2 / (g / 2) ** 3
    assert sp.simplify(half_gap_bound - 8 * CG * PT**2 / g**3) == 0

    A, B, E, P = sp.symbols("A B E P", positive=True)
    scaled_E = A**2 * E / B**3
    scaled_P = A * P / B**3
    scaled_c = A * c
    assert sp.simplify(
        scaled_E - scaled_c * scaled_P - A**2 * (E - c * P) / B**3
    ) == 0

    kappa = sp.symbols("kappa", positive=True)
    scaled_intensity = A * B * kappa / B**2
    assert sp.simplify(scaled_intensity - A * kappa / B) == 0

    print("deficit evolution = delta0 + E12_t - E12_0")
    print("impulse-product gap =", sp.factor(product_gap))
    print("half-gap cross bound =", sp.simplify(half_gap_bound))
    print("scaled functional factor =", A**2 / B**3)
    print("scaled intensity =", sp.simplify(scaled_intensity))
    print("ALL 5 EXACT CHECKS PASS")


if __name__ == "__main__":
    main()

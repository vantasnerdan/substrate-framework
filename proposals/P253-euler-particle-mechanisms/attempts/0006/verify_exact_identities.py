"""Exact algebra checks for the P253/0006 deficit construction."""

import sympy as sp


def main() -> None:
    s, e1, e2, p1, p2, cross, c = sp.symbols(
        "S E1 E2 P1 P2 E12 c", real=True
    )
    d1 = s - (e1 - c * p1)
    d2 = s - (e2 - c * p2)
    total_energy = e1 + e2 + cross
    total_impulse = p1 + p2
    deficit_sum = sp.expand(d1 + d2)
    conserved_form = sp.expand(2 * s - total_energy + cross + c * total_impulse)
    assert sp.simplify(deficit_sum - conserved_form) == 0

    pe, cg, ell = sp.symbols("P_epsilon C_G ell", positive=True)
    product_gap = sp.expand(pe**2 - p1 * p2).subs(p2, 2 * pe - p1)
    assert sp.simplify(product_gap - (pe - p1) ** 2) == 0
    cross_bound = 4 * cg * pe**2 / ell**3
    assert sp.simplify(cross_bound.subs(ell, sp.Symbol("g", positive=True) / 2)) == (
        32 * cg * pe**2 / sp.Symbol("g", positive=True) ** 3
    )

    pi = sp.pi
    cbs = 3 / (2 ** sp.Rational(4, 3) * pi ** sp.Rational(1, 3))
    radius = sp.symbols("b", positive=True)
    omega_inf, omega_l1 = sp.symbols("omega_inf omega_l1", positive=True)
    split_bound = omega_inf * radius + omega_l1 / (4 * pi * radius**2)
    optimizer = (omega_l1 / (2 * pi * omega_inf)) ** sp.Rational(1, 3)
    optimized = sp.simplify(split_bound.subs(radius, optimizer))
    expected = cbs * omega_l1 ** sp.Rational(1, 3) * omega_inf ** sp.Rational(2, 3)
    assert sp.simplify(optimized - expected) == 0

    print("deficit sum =", deficit_sum)
    print("conserved form =", conserved_form)
    print("impulse-product gap =", sp.factor(product_gap))
    print("gap g cross bound =", 32 * cg * pe**2 / sp.Symbol("g", positive=True) ** 3)
    print("Biot-Savart constant =", cbs)
    print("ALL 4 EXACT CHECKS PASS")


if __name__ == "__main__":
    main()

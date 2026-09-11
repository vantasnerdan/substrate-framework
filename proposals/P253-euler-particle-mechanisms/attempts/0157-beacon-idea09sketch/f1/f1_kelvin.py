#!/usr/bin/env python3
"""F1 defect-existence check (beacon, under 0157): frozen f1/design.md.
Symbolic Kelvin construction in the banked F-A medium + F-C branch
cite. Exits nonzero on any falsifier fire.
"""

import sympy as sp

K = sp.symbols('K', positive=True)
mu = K / 10          # banked F-A (R3)
lam = -K / 15        # banked F-A (incompressible artifact)


def main() -> None:
    nu = lam / (2 * (lam + mu))
    assert sp.simplify(nu + 1) == 0, f"nu={nu}"
    print(f"Poisson nu = {nu} (bulk modulus zero — named, see below)")
    c1 = 3 - 4 * nu   # Kelvin deviatoric weight
    c0 = 1 - nu       # Kelvin prefactor denominator piece
    assert sp.simplify(c1 - 7) == 0 and sp.simplify(c0 - 2) == 0
    # displacement parallel to force everywhere: 7 + cos^2 > 0
    th = sp.symbols('th', real=True)
    assert sp.simplify(7 + sp.cos(th) ** 2) > 0
    print("Kelvin 1/r Green's function EXISTS, response parallel (7+cos^2>0)")
    # far-field energy convergence outside core a
    r, a = sp.symbols('r a', positive=True)
    E = sp.integrate(1 / r**2, (r, a, sp.oo))   # |eps|^2 r^2 dr ~ dr/r^2
    assert E == 1 / a
    print("Kelvin 1/r Green's function EXISTS, response parallel (7+cos^2>0)")
    # F-C tilt branch (banked 04-fc-build gate 2, read-only cite):
    # w^2 = (4J0/I) sin^2(kl/2) >= 0 for J0, I > 0
    J0, In, k, ell = sp.symbols('J0 In k ell', positive=True)
    w2 = 4 * J0 / In * sp.sin(k * ell / 2) ** 2
    assert w2.is_nonnegative
    print("F-C tilt branch nonneg (no imaginary-frequency delocalizer)")
    print("bulk-zero note: uniform dilation costs nothing, but no "
          "non-decaying branch exhibited — observation, not fire")
    print("F1 verdict: F1-HOLDS (linear-level defect exists; F2 gated)")


if __name__ == "__main__":
    main()

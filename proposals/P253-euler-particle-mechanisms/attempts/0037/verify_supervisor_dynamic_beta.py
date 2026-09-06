"""Direct symbolic Euler quadratic projection; no modulation theorem claimed."""

import sympy as sp

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P253-0037-supervisor-dynamic-beta")
    r, c, n = sp.symbols("r c N", positive=True)
    f = sp.Function("f")(r)
    b = sp.Function("b")(r)
    d = sp.Function("d")(r)
    e, g = b * d * f / c**2, -d * f / c
    n_eta = 2 * g**2 / r**4 - (f.diff(r) * e - f * e.diff(r)) / r
    n_chi = -(f.diff(r) * g - f * g.diff(r)) / r
    projected = (f * n_eta - b * f * n_chi / c) * r / (2 * n)
    bracket = b.diff(r) * d + 2 * b * d.diff(r) + 2 * d**2 / r**3
    reduced = f**3 * bracket / (2 * n * c**2)
    ledger.check("B5 direct Euler quadratic integrand",
                 sp.simplify(projected - reduced) == 0)
    wrong_dual = (f * n_eta + b * f * n_chi / c) * r / (2 * n)
    ledger.check("B5 opposite dual sign changes the projection",
                 sp.simplify(wrong_dual - reduced) != 0)
    ledger.check("B2 right-vector dual normalization density",
                 sp.simplify(f * e - b * f * g / c - 2 * f * e) == 0)
    ledger.check("B2 opposite-vector annihilation density",
                 sp.simplify(f * e + b * f * g / c) == 0)

    angular = sp.Function("L")(r)
    physical = bracket.subs({b: 2 * angular / r**4,
                            d: angular.diff(r) / r}).doit()
    j = -4 / c**3 * ((angular.diff(r)**2 + angular * angular.diff(r, 2)) / r**4
                    - 3 * angular * angular.diff(r) / r**5)
    ledger.check("B6-B7 exact steady/dynamic source identity",
                 sp.simplify(physical + c**3 * j / r) == 0)
    omega = sp.symbols("Omega", real=True)
    ledger.check("B6 uniform-vorticity core contribution vanishes",
                 sp.simplify(physical.subs(angular, omega * r**2).doit()) == 0)
    ledger.check("B8 projected integrand equals normalized steady source",
                 sp.simplify(f**3 * physical / (2 * n * c**2)
                             + c * j * f**3 / (2 * n * r)) == 0)
    beta = sp.symbols("beta", positive=True)
    sigma = c / (2 * n)
    ledger.check("B8 eigenvalue-slope normalization",
                 sp.simplify(-c * beta / n + 2 * sigma * beta) == 0)
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())

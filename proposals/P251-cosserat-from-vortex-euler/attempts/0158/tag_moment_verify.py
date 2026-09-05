"""Exact leading registered-tag displacement moment versus existing controls."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0158-actual-tag-displacement-moment")
    x = s.Symbol("x", real=True)
    n, jet_order = 8, 7
    laguerre = s.assoc_laguerre(n, 1, x)
    pressure = laguerre-2*s.diff(laguerre, x)
    base = [s.expand(x**j*pressure) for j in range(jet_order+1)]
    # Remove the common exponential ONLY after differentiating it.
    derivatives = [s.expand(-value/2+s.Rational(3, 2)*x*(s.diff(value, x)-value/2))
                   for value in base]
    rows = base+derivatives
    degree = max(s.degree(value, x) for value in rows+[laguerre])

    def coefficients(value):
        return [s.expand(value).coeff(x, power) for power in range(degree+1)]

    matrix = s.Matrix([coefficients(value) for value in rows])
    augmented = matrix.col_join(s.Matrix([coefficients(laguerre)]))
    print("Existing exponential pressure-control rank:", matrix.rank())
    print("Rank after the literal displacement-helicity row L_n:", augmented.rank())
    print("Polynomial degree window:", degree)
    checks.check("the source's sixteen pressure/time/carrier rows remain independent",
                 matrix.rank() == 2*(jet_order+1))
    # This check asks the actual independence proposition, not a rank inferred
    # from the number of candidate controls.
    checks.check("the actual tagged displacement moment has an additional control direction",
                 augmented.rank() == matrix.rank()+1)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

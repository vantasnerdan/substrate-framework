"""Exact derivative bookkeeping for the fixed-cell Euler/Lin diagonal bound."""

import sympy as sp

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0196-diagonal")
    time = sp.Symbol("T", nonnegative=True)
    b1, b2, b3 = sp.symbols("B1 B2 B3", nonnegative=True)
    d0, d1, d2, d3 = sp.symbols("d0 d1 d2 d3", nonnegative=True)
    polynomials = [
        d0,
        d1 + b1 * d0 * time,
        d2 + (2 * b1 * d1 + b2 * d0) * time + b1**2 * d0 * time**2,
        d3 + (3 * b1 * d2 + 3 * b2 * d1 + b3 * d0) * time
        + (3 * b1**2 * d1 + 3 * b1 * b2 * d0) * time**2
        + b1**3 * d0 * time**3,
    ]
    derivative_bounds = [None, b1, b2, b3]
    for order in range(1, 4):
        forcing = sum(
            sp.binomial(order, j) * derivative_bounds[j] * polynomials[order - j]
            for j in range(1, order + 1)
        )
        checks.check(
            f"Duhamel majorant solves the exact order-{order} derivative recurrence",
            sp.expand(sp.diff(polynomials[order], time) - forcing) == 0,
        )
    checks.check(
        "every majorant retains its actual initial derivative bound",
        [p.subs(time, 0) for p in polynomials] == [d0, d1, d2, d3],
    )
    form = sp.Function("F")(time)
    state = sp.Function("y")(time)
    leibniz = sum(
        sp.factorial(3) / (sp.factorial(j) * sp.factorial(left) * sp.factorial(right))
        * sp.diff(form, time, j) * sp.diff(state, time, left) * sp.diff(state, time, right)
        for j in range(4) for left in range(4 - j)
        for right in [3 - j - left]
    )
    checks.check(
        "the complete quadratic action remainder retains all three derivative placements",
        sp.expand(sp.diff(form * state**2, time, 3) - leibniz) == 0,
    )
    # The analytic proof supplies the fixed-cell operator bounds and the
    # genuinely finite control norms; this checks their diagonal bookkeeping.
    accuracy, norm, constant = sp.symbols("epsilon N C", positive=True)
    upper = accuracy / (1 + norm + constant)
    checks.check(
        "one actual scale choice suppresses the Sobolev and normalized cubic remainders",
        sp.simplify(accuracy - upper * norm) == accuracy * (1 + constant) / (1 + norm + constant)
        and sp.simplify(accuracy - upper * constant) == accuracy * (1 + norm) / (1 + norm + constant),
    )
    return int(checks.finish())


if __name__ == "__main__":
    raise SystemExit(main())

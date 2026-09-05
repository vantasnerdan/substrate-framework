"""Exact accumulated boundary charge; no inferred zero initial current."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0235-initial-current")
    t, tau = s.symbols("t tau", real=True)
    i0, ix, sigma, c0, d0, c2, c4, d2, d4 = s.symbols(
        "I0 IX Sigma c0 d0 c2 c4 d2 d4", real=True
    )
    correlation = c0 + c2*t**2/2 + c4*t**4/24
    weighted = d0 + d2*t**2/2 + d4*t**4/24
    trace = 2*(i0+sigma*t**2)+2*(correlation-c0) + 2*(
        t*s.diff(weighted, t)-weighted+d0
    )
    q = (trace-ix)/6
    checks.check("time derivative includes both actual correlation kernels",
                 s.simplify(s.diff(q, t) - (
                     4*sigma*t+2*s.diff(correlation, t)
                     +2*t*s.diff(weighted, t, 2))/6) == 0)
    checks.check("initial current rows retain the full finite-cut virial",
                 q.subs(t, 0) == (2*i0-ix)/6
                 and s.diff(q, t).subs(t, 0) == 0
                 and s.simplify(s.diff(q, t, 2).subs({t: 0, d2: -sigma-c2})
                                -sigma/3) == 0)
    displacement = s.Matrix([1+t+t**3, 2-t**2, t**4-t])
    initial = s.Matrix(s.symbols("G0:3", real=True))
    accumulated = initial + s.integrate(
        (q*displacement.diff(t)).subs(t, tau), (tau, 0, t)
    )
    primitive = (initial+q*displacement-q.subs(t, 0)*displacement.subs(t, 0)
                 -s.integrate((s.diff(q, t)*displacement).subs(t, tau), (tau, 0, t)))
    checks.check("accumulation preserves the independent initial charge and memory",
                 s.simplify(accumulated-primitive) == s.zeros(3, 1)
                 and accumulated.subs(t, 0) == initial)
    checks.check("dropping the derivative-of-q memory changes the measured rate",
                 s.simplify((q*displacement).diff(t)-q*displacement.diff(t))
                 != s.zeros(3, 1))
    fraction = s.symbols("fraction", positive=True)
    scaled = q.subs({a: fraction*a for a in
                    (i0, ix, sigma, c0, d0, c2, c4, d2, d4)}, simultaneous=True)
    checks.check("the current tensor scales linearly with a physical positive tag fraction",
                 s.simplify(scaled-fraction*q) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

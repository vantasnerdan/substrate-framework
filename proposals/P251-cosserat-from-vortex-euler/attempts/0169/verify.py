"""Exposing exact localization checks, not a new global existence oracle."""

import sympy as s

from substrate_framework.euler_acoustic import triangular_euler_array
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0169-pressure-localization")
    x, y = s.symbols("x y", real=True)
    c = s.Symbol("C", positive=True)
    field = triangular_euler_array(1, 1, 1, (x, y))
    p, v = field.pressure, field.velocity
    advected = s.trigsimp(v.dot(s.Matrix([s.diff(p, x), s.diff(p, y)])))
    point = s.simplify(advected.subs({x: s.pi/2, y: s.pi/s.sqrt(3)},
                                    simultaneous=True))
    checks.check("actual triangular pressure is not transported", point != 0)
    print("Exact advected-pressure value at registered point:", point)
    lifted = v.col_join(s.Matrix([s.sqrt(c+field.streamfunction**2)]))
    full_advected = lifted.dot(s.Matrix([s.diff(p, x), s.diff(p, y), 0]))
    checks.check("the actual Bernoulli lift cannot change the pressure residual",
                 s.simplify(full_advected-advected) == 0 and s.diff(full_advected, c) == 0)
    chi = 1+p
    modulated = chi*v
    divergence = s.diff(modulated[0], x)+s.diff(modulated[1], y)
    checks.check("nonconstant pressure modulation exposes the full divergence defect",
                 s.simplify(divergence-advected) == 0)
    # Compare by chain rule before trigonometric expansion: all terms retained.
    gradient = s.Matrix([s.diff(p, x), s.diff(p, y)])
    residual = modulated.jacobian((x, y))*modulated+chi**2*gradient
    base_euler = v.jacobian((x, y))*v+gradient
    # Expand the product rule first, then simplify the lower-degree actual
    # Euler residual. Direct trigsimp of the full high-degree product left
    # an unresolved zero on first execution (diagnosis.stdout).
    chain_difference = residual-chi*v.dot(gradient)*v-chi**2*base_euler
    checks.check("Euler acceleration has the same nonadvected-pressure source",
                 chain_difference.applyfunc(s.expand) == s.zeros(2, 1)
                 and base_euler.applyfunc(s.trigsimp) == s.zeros(2, 1))
    circular = s.Matrix([-y, x])
    pc = (x*x+y*y)/2
    localized = s.exp(-pc)*circular
    localized_pressure = -s.exp(-2*pc)/2
    checks.check("a genuinely advected pressure licenses modulation",
                 s.diff(localized[0], x)+s.diff(localized[1], y) == 0)
    rc = localized.jacobian((x, y))*localized+s.Matrix([
        s.diff(localized_pressure, x), s.diff(localized_pressure, y)])
    checks.check("licensed modulation retains exact full stationary Euler",
                 s.simplify(rc) == s.zeros(2, 1))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

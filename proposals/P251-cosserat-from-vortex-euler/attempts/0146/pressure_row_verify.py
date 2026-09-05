"""Exact adjoint pressure-row reduction and logarithmic separatrix integrability."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0146-pressure-row")
    omega, nu = s.symbols("omega nu", positive=True)
    j = s.Matrix([[0, -1], [1, 0]])
    core = omega*j
    transport = nu*j
    coupling = s.Matrix([[1, 2], [3, 5]])
    row = s.Matrix([[2, -1]])
    h = row*core.inv()
    d = h*coupling*transport.inv()
    f = s.Matrix(s.symbols("f1 f2"))
    b = s.Matrix(s.symbols("b1 b2"))
    derivative = h*(core*f+coupling*b)-d*transport*b
    checks.check("actual core pressure row is the derivative of a core-plus-ambient current",
                 s.simplify(derivative[0]-(row*f)[0]) == 0)
    force_core = s.Matrix(s.symbols("Fc1 Fc2"))
    force_ambient = s.Matrix(s.symbols("Fb1 Fb2"))
    forced = h*(core*f+coupling*b+force_core)+d*(-transport*b+force_ambient)
    checks.check("finite-k response forcing survives the exact row reduction",
                 s.simplify((row*f-forced+h*force_core+d*force_ambient)[0]) == 0)
    zero_row = s.Matrix([[1, -2]])
    d0 = zero_row*coupling*transport.inv()
    checks.check("translation row keeps its ambient conserved contribution",
                 s.simplify((zero_row*coupling*b-d0*transport*b)[0]) == 0)
    theta = s.Symbol("theta", real=True)
    odd_source = s.cos(theta)+2*s.sin(3*theta)
    checks.check("central inversion cancels the complete orbit average",
                 s.simplify(odd_source+odd_source.subs(theta, theta+s.pi)) == 0
                 and s.integrate(odd_source, (theta, 0, 2*s.pi)) == 0)
    level = s.Symbol("a", positive=True)
    logarithmic_mass = s.integrate((-s.log(level))**3, (level, 0, 1))
    checks.check("the unbounded separatrix primitive has finite coarea L2 mass",
                 logarithmic_mass == 6)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

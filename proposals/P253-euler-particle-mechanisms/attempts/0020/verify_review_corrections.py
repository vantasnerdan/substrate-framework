"""Exposing exact tests for the action sign and retained weighted coupling."""

import sympy as sp

t = sp.Symbol("t", real=True)
q, p = sp.Function("q")(t), sp.Function("p")(t)
H = (q*q+p*p)/2
omega_path = q*sp.diff(p, t)-p*sp.diff(q, t)


def euler_lagrange(L, variable):
    return sp.diff(L, variable)-sp.diff(sp.diff(L, sp.diff(variable, t)), t)


wrong = omega_path/2-H
correct = -omega_path/2-H
wrong_rows = [euler_lagrange(wrong, v) for v in (q, p)]
correct_rows = [euler_lagrange(correct, v) for v in (q, p)]
physical = {sp.diff(q, t): p, sp.diff(p, t): -q}
assert [sp.simplify(row.subs(physical)) for row in correct_rows] == [0, 0]
assert [sp.simplify(row.subs(physical)) for row in wrong_rows] == [-2*q, -2*p]

r, xi, kabs = sp.symbols("r xi kabs", positive=True)
b, kz, e = sp.symbols("b kz e", real=True)
a = b/(e*r*r)-2*xi*kz/(e*r**4)
J = sp.I*sp.Matrix([[a, b], [b, 0]])
L = sp.Matrix([[r*r/kabs**2, -e], [-e, 1/r**2]])
weight = sp.diag(1/kabs, 1)
weighted = sp.simplify(weight*J*L*weight.inv()+sp.I*e*b*sp.eye(2))
product = sp.simplify(weighted[0, 1]*weighted[1, 0])
assert product == -2*b*kz*xi/(kabs*kabs*r*r)
N, b0, z0 = sp.symbols("N b0 z0", real=True, nonzero=True)
assert sp.simplify(product.subs({b: N*b0, kz: N*z0, kabs: N})) == -2*b0*z0*xi/r**2
print("Correct quadratic action Euler-Lagrange rows:", correct_rows)
print("Opposite action sign on physical Hamiltonian flow:", [row.subs(physical) for row in wrong_rows])
print("Energy-weighted pressure/shear product:", product)
print("PASS: action sign exposed; pressure product remains order zero at fixed wave direction.")

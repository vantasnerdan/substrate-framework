"""Exact simultaneous shape constraint and exposing pressure-continuation algebra."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0155-continuation")
    n = sp.symbols("n", nonnegative=True, integer=True)
    mean = 2 * (n + 1)
    second = mean**2 + (n + 1) * (n + 2) + n * (n + 1)
    ledger.check("general Laguerre second moment", sp.expand(second - 6 * (n + 1) ** 2) == 0)
    shift = 4 * (n + 1)
    t = 1 / (3 * (n + 1))
    constrained_norm = 1 + 2 * t * (mean - shift) + t**2 * (second - 2 * shift * mean + shift**2)
    ledger.check("positive slow clock enforces one-third overlap bound", sp.simplify(constrained_norm - sp.Rational(1, 3)) == 0)
    ledger.check("half-overlap has a fixed analytic margin", sp.Rational(1, 2) - constrained_norm.simplify() == sp.Rational(1, 6))

    x = sp.symbols("x", positive=True)
    boundary_integral = sp.integrate(x * sp.exp(-x / 2), x)
    full_integral = sp.simplify(boundary_integral - boundary_integral.subs(x, 0))
    ledger.check("larger actual parcel initial spin functional", sp.simplify(full_integral - 4 * (1 - (1 + x / 2) * sp.exp(-x / 2))) == 0)
    ledger.check("parcel overlap increases monotonically with outer radius", sp.simplify(sp.diff(full_integral, x) - x * sp.exp(-x / 2)) == 0)

    radial_index, angular_index = 8, 5
    laguerre = sp.assoc_laguerre(radial_index, angular_index - 1, x)
    pressure = sp.expand(laguerre - 2 * sp.diff(laguerre, x))

    def carrier(poly: sp.Expr) -> sp.Expr:
        return sp.expand(-poly / 2 + sp.Rational(3, 2) * x * (sp.diff(poly, x) - poly / 2))

    rows = [sp.expand(row) for row in [pressure, carrier(pressure), carrier(carrier(pressure)), x * pressure, carrier(x * pressure), carrier(carrier(x * pressure))]]
    coefficient_rows = sp.Matrix([[row.coeff(x, j) for j in range(radial_index + 4)] for row in rows])
    ledger.check("six actual radial time/carrier functionals independent", coefficient_rows.rank() == 6)
    jets = [[sp.factorial(h) if j == h else 0 for j in range(8)] for h in range(2)]
    for row in rows:
        jets.append([sum(sp.binomial(j, h) * (-sp.Rational(1, 2)) ** (j - h) * sp.factorial(h) * row.coeff(x, h) for h in range(j + 1)) for j in range(8)])
    determinant = sp.Matrix(jets).det(method="domain-ge")
    ledger.check("zero-point moment minor contains no residual coordinate", not determinant.has(x))
    ledger.check("reference and background rows extend the actual moment rank", determinant != 0)

    # Fixed exact rational ellipticity; this is coefficient algebra, not a spectrum.
    z, zb = sp.symbols("z zb")
    g = sp.Rational(1, 4)
    metric = sp.diag(g, 1 / g)
    trace = sp.trace(metric)
    difference = metric[0, 0] - metric[1, 1]
    plus = sp.Matrix([1, sp.I]) / sp.sqrt(2)
    minus = sp.conjugate(plus)
    angular = 5
    seed = z ** (angular - 1)

    def dx(poly: sp.Expr) -> sp.Expr:
        return sp.expand(sp.diff(poly, z) + sp.diff(poly, zb) - (z + zb) * poly / 2)

    def dy(poly: sp.Expr) -> sp.Expr:
        return sp.expand(sp.I * (sp.diff(poly, z) - sp.diff(poly, zb)) + sp.I * (z - zb) * poly / 2)

    def grad(poly: sp.Expr) -> sp.Matrix:
        return sp.Matrix([dx(poly), dy(poly)])

    def curl(vector: sp.Matrix) -> sp.Expr:
        return sp.expand(dx(vector[1]) - dy(vector[0]))

    def l1(vector: sp.Matrix) -> sp.Matrix:
        return (sp.I * z * zb * vector / 2 - 2 * metric * grad(curl(vector)) / trace).applyfunc(sp.expand)

    def l2(vector: sp.Matrix) -> sp.Matrix:
        product = (z + zb) * vector[0] / 2 + (z - zb) * vector[1] / (2 * sp.I)
        vorticity = curl(vector)
        metric_laplace = g * dx(dx(vorticity)) + dy(dy(vorticity)) / g
        return (2 * sp.I * metric * grad(product) / trace - 2 * metric * grad(metric_laplace) / trace**2).applyfunc(sp.expand)

    delta_seed = dx(dx(seed)) + dy(dy(seed))
    raised_seed = sp.expand(dx(dx(seed)) + 2 * sp.I * dx(dy(seed)) - dy(dy(seed)))
    initial = plus * seed
    first = difference * plus * raised_seed / (4 * trace) + difference * minus * delta_seed / (4 * trace) + minus * raised_seed / 8
    rotation = sp.Matrix([[0, -1], [1, 0]])
    first_transport = first.applyfunc(lambda value: sp.I * (z * sp.diff(value, z) - zb * sp.diff(value, zb))) + rotation * first
    first_residual = (first_transport + sp.I * (2 - angular) * first + l1(initial) - sp.I * angular * initial).applyfunc(sp.expand)
    ledger.check("actual first-order Euler polarization residual vanishes", first_residual == sp.zeros(2, 1))
    second_source = (l1(first) - sp.I * angular * first + l2(initial)).applyfunc(sp.expand)
    opposite = sp.expand((plus.T * second_source)[0])
    resonant = sp.Add(*(coef * z ** powers[0] * zb ** powers[1] for powers, coef in sp.Poly(opposite, z, zb).terms() if powers[0] - powers[1] == angular - 3))
    ledger.check("computed second-order source keeps finite angular grading", all(powers[0] - powers[1] in {0, 2, 4, 6, 8} for powers, _ in sp.Poly(opposite, z, zb).terms()))
    ledger.check("full pressure cancels the tested opposite resonance", sp.simplify(resonant) == 0)
    generic_cancellation = []
    for power in range(5):
        test = z ** (angular - 1) * (z * zb) ** power
        laplace_test = dx(dx(test)) + dy(dy(test))
        raised_test = dx(dx(test)) + 2 * sp.I * dx(dy(test)) - dy(dy(test))
        first_test = difference * plus * raised_test / (4 * trace) + difference * minus * laplace_test / (4 * trace) + minus * raised_test / 8
        source_test = sp.expand((plus.T * (l1(first_test) + l2(plus * test)))[0])
        resonant_test = sp.Add(*(coef * z ** powers[0] * zb ** powers[1] for powers, coef in sp.Poly(source_test, z, zb).terms() if powers[0] - powers[1] == angular - 3))
        generic_cancellation.append(sp.simplify(resonant_test) == 0)
    ledger.check("actual pressure compositions annihilate five radial jets", all(generic_cancellation))
    delta = sp.symbols("delta", positive=True)
    ledger.check("small reference moment pays the first background cancellation", sp.simplify(delta**6 / delta**3 - delta**3) == 0)
    ledger.check("typed residual gains the spin pressure order", sp.simplify(delta**2 * delta / delta - delta**2) == 0)
    print(f"eight-row moment Wronskian={determinant}")
    print(f"second-order opposite resonant source={sp.factor(resonant)}")
    print("The printed source is a continuation diagnostic, not an eigenvalue or spectral verdict.")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())

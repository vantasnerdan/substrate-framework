"""Exact finite-packet marker and common-circle geometry calculations."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0147")
    x, w = sp.symbols("x w", positive=True)
    delta, omega = sp.symbols("delta Omega", positive=True)
    n, time_order, axial_power, reference_order = 8, 7, 3, 3
    laguerre = sp.assoc_laguerre(n, 1, x)
    poly = sp.Poly(laguerre, x)
    laplace = sum(coef * sp.factorial(degree[0] + 2) * (2 / w) ** (degree[0] + 3) for degree, coef in poly.terms())
    angle = sp.factor(laplace / laplace.subs(w, 1))
    predicted = (-1) ** n * (w - 2) ** (n - 1) * (w - n - 2) / ((n + 1) * w ** (n + 3))
    ledger.check("general excited physical Laplace angle", sp.simplify(angle - predicted) == 0)
    phase_slope = sp.diff(angle, w).subs(w, 1)
    gamma = -2 * omega + ((2 * n + 2) + phase_slope) * sp.sqrt(2) * omega * delta
    ledger.check("n8 physical phase retains negative optical clock", sp.simplify(gamma + 2 * omega + sp.sqrt(2) * omega * delta / 9) == 0)
    leading_square = sp.series(gamma**2, delta, 0, 2).removeO()
    curvature = sp.expand(delta**2 * sp.diff(leading_square, delta, 2) / 4 + 3 * delta * sp.diff(leading_square, delta) / 4)
    ledger.check("n8 natural-scale physical curvature positive", curvature == sp.sqrt(2) * omega**2 * delta / 3)
    norm_poly = sp.Poly(x * laguerre**2, x)
    norm = sum(coef * sp.factorial(degree[0]) for degree, coef in norm_poly.terms())
    ledger.check("n8 KKS Laguerre normalization", norm == n + 1)

    pressure_poly = sp.expand(laguerre - 2 * sp.diff(laguerre, x))
    ledger.check("pressure polynomial has no zero-axis root", pressure_poly.subs(x, 0) == (n + 1) ** 2)
    ledger.check("pressure/log-derivative coprimality", sp.gcd(pressure_poly, x * sp.diff(pressure_poly, x)) == 1)
    rows = []
    for j in range(time_order + 1):
        first = sp.expand(x**j * pressure_poly)
        derivative = sp.expand(-first / 2 + sp.Rational(3, 2) * x * (sp.diff(first, x) - first / 2))
        rows += [first, derivative]
    coefficients = sp.Matrix([[row.coeff(x, j) for j in range(n + time_order + 2)] for row in rows])
    ledger.check("all sixteen slow-time carrier rows independent", coefficients.rank() == 2 * (time_order + 1))
    # At x=0, derivatives of exp(-x/2)P are finite exact binomial sums.
    order = len(rows) + reference_order + 1
    jet_rows = [
        [sp.factorial(h) if j == h else sp.Integer(0) for j in range(order)]
        for h in range(reference_order + 1)
    ]
    for row in rows:
        jet_rows.append([
            sum(sp.binomial(j, h) * (-sp.Rational(1, 2)) ** (j - h) * sp.factorial(h) * row.coeff(x, h) for h in range(j + 1))
            for j in range(order)
        ])
    wronskian = sp.Matrix(jet_rows).det(method="domain-ge")
    ledger.check("twenty actual radial/background rows have nonzero minor", wronskian != 0)
    ledger.check("background moment rows are genuinely additional", sp.Matrix(jet_rows[4:]).rank() == 16)

    c, length, p, p_star, q = sp.symbols("c L p p_star q", positive=True)
    ell_z = c * length
    alpha = length**2 / (length**2 + ell_z**2)
    q_bar = alpha * p + (1 - alpha) * p_star
    gaussian_cost = length**2 * ell_z**2 / (length**2 + ell_z**2)
    square = sp.expand(length**2 * (q - p) ** 2 + ell_z**2 * (q - p_star) ** 2)
    completed = sp.expand((length**2 + ell_z**2) * (q - q_bar) ** 2 + gaussian_cost * (p - p_star) ** 2)
    ledger.check("finite-packet carrier/filter square completion", sp.simplify(square - completed) == 0)
    ledger.check("actual carrier derivative factor", sp.simplify(sp.diff(q_bar, p) - 1 / (1 + c**2)) == 0)
    beta_mass = sp.integrate(length**2 * sp.exp(-length**2 * q**2), (q, -sp.oo, sp.oo))
    ledger.check("finite packet KKS uses full Plancherel mass", beta_mass == sp.sqrt(sp.pi) * length)
    v = sp.symbols("v", positive=True)
    axial_matrix = sp.Matrix([[1, v], [v, 3 * v**2]])
    ledger.check("two physical axial controls independent", sp.factor(axial_matrix.det()) == 2 * v**2)
    z = sp.symbols("z", real=True)
    dual = (sp.Matrix([[1, z**2]]) * axial_matrix.inv()).tolist()[0]
    ledger.check("dual axial controls preserve reference on every sheet", sp.simplify(dual[0] + v * dual[1]) == 1)
    envelope = sp.exp(-gaussian_cost * (p - p_star) ** 2 / 2)
    wrong_difference = sp.simplify((sp.diff(envelope, p, 2) - sp.diff(1 / envelope, p, 2)).subs(p, p_star))
    ledger.check("radial-only packet mismatch exposed", sp.simplify(wrong_difference + 2 * gaussian_cost) == 0)
    loss = sp.simplify(delta ** (time_order + 1) * (delta ** -axial_power) ** 2)
    ledger.check("seven time moments pay full carrier loss", loss == delta**2)
    old_loss = sp.simplify(delta**2 * (delta ** -axial_power) ** 2)
    ledger.check("first time moment alone does not pay loss", old_loss == delta**-4)
    phase_loss = sp.simplify(delta ** (3 * (reference_order + 1) - 1) * (delta ** -axial_power) ** 2)
    ledger.check("background moment repair controls small-Q phase", phase_loss == delta**5)

    # Exact finite representative of the arbitrary-order positive angular rule.
    nodes = [(2 - sp.sqrt(3)) / 4, sp.Rational(1, 2), (2 + sp.sqrt(3)) / 4]
    moments = [sp.simplify(sum(t**h for t in nodes) / 3) for h in range(6)]
    ledger.check("positive angular quadrature moments through degree5", all(moments[h] == sp.binomial(2 * h, h) / 4**h for h in range(6)))
    moment_matrix = sp.Matrix([[t**h for t in nodes] for h in range(3)])
    weights = sp.simplify(moment_matrix.inv() * sp.Matrix([1, sp.Rational(1, 2), sp.Rational(3, 8)]))
    ledger.check("common-circle moment solve preserves positive weights", weights == sp.Matrix([sp.Rational(1, 3)] * 3))
    xx, zz = sp.symbols("xx zz", real=True)
    quartic = sum(t**2 * xx**4 / 24 + t * (1 - t) * xx**2 * zz**2 / 4 + (1 - t)**2 * zz**4 / 24 for t in nodes) / 3
    ledger.check("multi-CK quartic equals actual Lundquist jet", sp.simplify(quartic - (xx**2 + zz**2) ** 2 / 64) == 0)
    tg_quartic = (xx**4 + zz**4) / 96 + xx**2 * zz**2 / 16
    ledger.check("single TG angular anisotropy exposed", sp.expand(tg_quartic - (xx**2 + zz**2) ** 2 / 64 + (xx**4 - 6 * xx**2 * zz**2 + zz**4) / 192) == 0)

    f, fp = sp.symbols("f fp", real=True)
    energy_derivative = sp.expand(2 * f * fp - 2 * fp * (1 + 1 / (4 * x**2)) * f)
    ledger.check("Bessel-zero normalization energy identity", energy_derivative == -f * fp / (2 * x**2))
    inverse_radius, derivative_symbol = sp.symbols("epsilon D", real=True)
    pressure_difference = inverse_radius * derivative_symbol / (1 + inverse_radius * xx) + q**2 * (1 - (1 + inverse_radius * xx) ** -2)
    first_pressure = sp.series(pressure_difference, inverse_radius, 0, 2).removeO()
    ledger.check("full cylindrical pressure leading correction", sp.expand(first_pressure - inverse_radius * (derivative_symbol + 2 * q**2 * xx)) == 0)

    print(f"radial/background Wronskian at zero={wronskian}")
    print(f"finite packet KKS factor={beta_mass}; carrier factor={alpha}")
    print(f"curvature={curvature}; time/carrier remainder={loss}; reference remainder={phase_loss}")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())

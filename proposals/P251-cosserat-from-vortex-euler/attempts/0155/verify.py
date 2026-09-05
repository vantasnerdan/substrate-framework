"""Exact elliptic strain, physical observation and pressure-order identities."""

from __future__ import annotations

import sympy as sp

from substrate_framework.verification import CheckLedger


def main() -> int:
    ledger = CheckLedger("P251-0155")
    s, rho = sp.symbols("s rho", positive=True)
    omega, mass, freq = sp.symbols("Omega M omega", positive=True)
    d = s**4
    strain = sp.Matrix([[0, 1], [-d, 0]])
    scale = sp.diag(s, 1 / s)
    metric = scale**2
    rotation = sp.Matrix([[0, -1], [1, 0]])
    ledger.check("actual strain is elliptic but not rigid", strain**2 == -d * sp.eye(2) and strain + strain.T != sp.zeros(2))
    ledger.check("biorthogonal transport rotation", sp.simplify(scale * strain * scale.inv() + s**2 * rotation) == sp.zeros(2))
    ledger.check("physical pressure metric retained", sp.simplify(scale * rotation * scale.inv() - metric * rotation) == sp.zeros(2))
    ledger.check("Kelvin projected inverse determinant identity", metric * rotation * metric == rotation)
    x, y = sp.symbols("X Y", real=True)
    quartic = sp.expand((d * (x / s) ** 4 + (s * y) ** 4) / 24)
    ledger.check("full axial quartic is anisotropic", quartic == (x**4 + d * y**4) / 24)

    plus = sp.Matrix([1, sp.I]) / sp.sqrt(2)
    minus = sp.conjugate(plus)
    test = x**5 + x * y**3 + sp.I * y**4
    velocity = plus * test
    divergence = sp.diff(velocity[0], x) + sp.diff(velocity[1], y)
    force = metric * sp.Matrix([sp.diff(divergence, x), sp.diff(divergence, y)])
    laplace = sp.diff(test, x, 2) + sp.diff(test, y, 2)
    raise_twice = sp.diff(test, x, 2) + 2 * sp.I * sp.diff(test, x, y) - sp.diff(test, y, 2)
    trace = sp.trace(metric)
    difference = metric[0, 0] - metric[1, 1]
    plus_force = (minus.T * force)[0]
    minus_force = (plus.T * force)[0]
    ledger.check("resonant and raised plus pressure components", sp.simplify(plus_force - trace * laplace / 4 - difference * raise_twice / 4) == 0)
    ledger.check("opposite pressure polarization retained", sp.simplify(minus_force - difference * laplace / 4 - trace * raise_twice / 4) == 0)
    pressure = x**3 * y + x * y**4
    torque = metric[1, 1] * x * sp.diff(pressure, y) - metric[0, 0] * y * sp.diff(pressure, x)
    torque_split = trace * (x * sp.diff(pressure, y) - y * sp.diff(pressure, x)) / 2 - difference * (x * sp.diff(pressure, y) + y * sp.diff(pressure, x)) / 2
    ledger.check("literal physical torque angular splitting", sp.simplify(torque - torque_split) == 0)

    r, z = sp.symbols("r z", positive=True)
    # z lies on the unit circle; complex conjugation is represented by z -> 1/z.
    denominator = 1 + r**2 + r * (z + 1 / z)
    numerator = 1 + r**2 + 2 * r / z
    shear = r * (z - 1 / z) / (sp.I * denominator)
    factor = numerator / denominator
    ledger.check("physical Euclidean observation is exact shear", sp.simplify(factor - (1 - sp.I * shear)) == 0)
    ledger.check("spin-angle observation determinant is constant", sp.simplify((factor + factor.subs(z, 1 / z)) / 2 - 1) == 0)
    a, adot = sp.symbols("a adot", real=True)
    hessian = sp.Matrix([[mass * freq**2, freq * a], [freq * a, (1 + a**2 - adot / freq) / mass]])
    ledger.check("time-dependent physical Hamiltonian determinant", sp.simplify(hessian.det() - freq**2 + freq * adot) == 0)
    theta, momentum, theta_dot = sp.symbols("theta P theta_dot", real=True)
    q = theta + a * momentum / (mass * freq)
    hamiltonian = sp.expand(momentum**2 / (2 * mass) + mass * freq**2 * q**2 / 2 - adot * momentum**2 / (2 * mass * freq))
    ledger.check("observation time connection is included", sp.simplify(hamiltonian - (sp.Matrix([[theta, momentum]]) * hessian * sp.Matrix([theta, momentum]))[0] / 2) == 0)
    b = 1 + a**2 - adot / freq
    stationary_momentum = mass * (theta_dot - freq * a * theta) / b
    reduced = sp.factor((momentum * theta_dot - hamiltonian).subs(momentum, stationary_momentum))
    ledger.check("physical scalar action retains its connection", sp.simplify(reduced - mass * ((theta_dot - freq * a * theta) ** 2 - freq**2 * b * theta**2) / (2 * b)) == 0)
    zz, zzbar = sp.symbols("Z Zbar")
    physical_rotation = metric * rotation * sp.Matrix([(zz + zzbar) / 2, (zz - zzbar) / (2 * sp.I)])
    rotated_z = physical_rotation[0] + sp.I * physical_rotation[1]
    ledger.check("calibrated tensor has physical rigid-rotation response", sp.simplify(rotated_z - sp.I * (trace * zz - difference * zzbar) / 2) == 0)

    # A fixed exact representative, not a numerical spectrum or parameter scan.
    radial, w, delta = sp.symbols("x w delta", positive=True)
    n, m = 8, 5
    polynomial = sp.assoc_laguerre(n, m - 1, radial)
    integral = sum(coef * sp.factorial(degree[0] + m) * (2 / w) ** (degree[0] + m + 1) for degree, coef in sp.Poly(polynomial, radial).terms())
    log_slope = sp.cancel(sp.diff(integral, w) / integral).subs(w, 1)
    expected = -sp.Rational((2 * n + m) ** 2 + m, 2 * n + m)
    ledger.check("actual tagged Laguerre Laplace phase", sp.simplify(log_slope - expected) == 0)
    gamma = 2 * omega + sp.Rational(m, 2 * n + m) * omega * delta
    leading_square = sp.series(gamma**2, delta, 0, 2).removeO()
    curvature = sp.expand(delta**2 * sp.diff(leading_square, delta, 2) / 4 + 3 * delta * sp.diff(leading_square, delta) / 4)
    ledger.check("calibrated physical carrier curvature natural scale", curvature == 3 * omega**2 * sp.Rational(m, 2 * n + m) * delta)
    pressure_polynomial = sp.expand(polynomial - 2 * sp.diff(polynomial, radial))
    ledger.check("actual pressure-spin polynomial has simple roots", sp.gcd(pressure_polynomial, sp.diff(pressure_polynomial, radial)) == 1)
    ledger.check("actual pressure-spin polynomial excludes axis root", pressure_polynomial.subs(radial, 0) != 0)
    norm = sum(coef * sp.factorial(degree[0]) for degree, coef in sp.Poly(radial ** (m - 1) * polynomial**2, radial).terms())
    ledger.check("finite radial KKS normalization", norm == sp.factorial(n + m - 1) / sp.factorial(n))
    print(f"n={n}, m={m}; radial norm={norm}; phase slope={log_slope}")
    print(f"natural curvature={curvature}; physical shear determinant=1")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())

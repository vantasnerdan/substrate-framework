"""Stationary radial ground-mode initial band edge and nonunit-spin map."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0176-stationary-shape")
    p, c, t = s.symbols("P c t", positive=True)
    points = (s.Integer(5), s.Rational(50, 3))
    weights = (s.Rational(32, 35), s.Rational(3, 35))
    mean = sum(w * x for w, x in zip(weights, points, strict=True))
    variance = sum(w * (x - mean)**2 for w, x in zip(weights, points, strict=True))
    third = sum(w * (x - mean)**3 for w, x in zip(weights, points, strict=True))
    checks.check("actual positive two-annulus mean", mean == 6)
    checks.check("actual positive response variance", variance == s.Rational(32, 3))
    checks.check("actual positive-response third moment", third == s.Rational(928, 9))
    # The physical weights are fixed; only the actual radial eigenfunction varies.
    amplitudes = [w * s.exp(-x * (p**s.Rational(3, 2) - 1) / 2) for w, x in zip(weights, points, strict=True)]
    mu = sum(x * a for x, a in zip(points, amplitudes, strict=True)) / sum(amplitudes)
    first = s.simplify(s.diff(mu, p).subs(p, 1))
    second = s.simplify(s.diff(mu, p, 2).subs(p, 1))
    checks.check("fixed physical marker first carrier derivative", first == -3 * variance / 4)
    checks.check("fixed physical marker second carrier derivative", second == -3 * variance / 8 + 9 * third / 16)
    gamma = -2 * c / s.sqrt(p) + c * p * mu / 2
    checks.check("stationary ground-mode physical clock positive", gamma.subs(p, 1) == c)
    checks.check("actual stationary-marker initial carrier slope vanishes", s.simplify(s.diff(gamma, p).subs(p, 1)) == 0)
    checks.check("actual stationary-marker clock curvature positive", s.simplify(s.diff(gamma, p, 2).subs(p, 1)) == 35 * c / 2)
    checks.check("actual squared-clock curvature includes full first derivative", s.simplify(s.diff(gamma**2, p, 2).subs(p, 1)) == 35 * c**2)
    weight = s.symbols("w", real=True)
    gap = points[1] - points[0]
    root_row = 1 + (points[0] + gap * weight) / 2 - 3 * gap**2 * weight * (1 - weight) / 8
    checks.check("smoothing continuation has a nonzero weight Jacobian", s.diff(root_row, weight).subs(weight, weights[1]) != 0)
    # Extract exact time coefficients directly from the two exponentials.
    numerator = sum(w * s.exp(s.I * c * x * t / 2) for w, x in zip(weights, points, strict=True))
    log_jet = s.series(s.log(s.series(numerator, t, 0, 5).removeO()), t, 0, 5).removeO().expand()
    phase_t2 = s.im(s.diff(log_jet, t)).coeff(t, 2)
    amplitude_t1 = s.re(s.diff(log_jet, t)).coeff(t, 1)
    checks.check("actual phase drift retains its third central moment", s.simplify(phase_t2 + c**3 * third / 16) == 0)
    checks.check("actual amplitude drift retains its positive variance", s.simplify(amplitude_t1 + c**2 * variance / 4) == 0)
    current_derivative = 2 * phase_t2 / c + amplitude_t1
    checks.check("natural slow-window current connection is not zero", s.simplify(current_derivative + 140 * c**2 / 9) == 0)

    rho, inertia, lock, eta, h, acoustic = s.symbols("rho j kappa eta h acoustic", positive=True)
    measured = s.Matrix([[1, -eta * inertia * h / (2 * rho)], [h / 2, 1]])
    inverse = measured.inv()
    mass = inverse.T * s.diag(rho, inertia) * inverse
    potential = inverse.T * s.diag(rho * acoustic * h**2, lock) * inverse
    mixed_mass = s.diff(mass[0, 1], h).subs(h, 0)
    mixed_stiffness = s.diff(potential[0, 1], h).subs(h, 0)
    checks.check("nonunit overlap retains the actual mixed kinetic row", s.simplify(mixed_mass - inertia * (eta - 1) / 2) == 0)
    checks.check("same action retains the actual mixed potential row", mixed_stiffness == -lock / 2)
    forcing = s.simplify(-lock * mixed_mass / inertia + mixed_stiffness)
    checks.check("physical optical translational transfer survives positive nonunit overlap", forcing == -eta * lock / 2)
    print(f"mean={mean}; variance={variance}; third={third}")
    print(f"gamma=c, gamma_P=0, gamma_PP=35c/2; current derivative={current_derivative}")
    print(f"physical optical forcing={forcing}")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

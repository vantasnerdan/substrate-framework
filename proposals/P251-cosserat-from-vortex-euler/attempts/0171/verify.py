"""Exact new-field optical pressure, scales, action and observation rows."""

from __future__ import annotations

import sympy as s

from substrate_framework.euler_acoustic import triangular_euler_array
from substrate_framework.verification import CheckLedger


def main() -> int:
    checks = CheckLedger("P251-0171")
    x, y, z = s.symbols("x y z", real=True)
    psi0, amp, lam, rho, cap = s.symbols("psi0 Psi lambda rho C", positive=True)
    om, a, b, p, delta = s.symbols("Omega a b p delta", positive=True)
    field = triangular_euler_array(amp, lam, rho, (x, y))
    psi, v = field.streamfunction, s.Matrix(field.velocity)
    axial = s.sqrt(cap + lam**2 * psi**2)
    full = v.col_join(s.Matrix([axial]))
    curl = s.Matrix([s.diff(axial, y), -s.diff(axial, x), field.vorticity])
    factor = -lam**2 * psi / axial
    checks.check("same imported array has exact variable-factor curl", all(s.simplify(q) == 0 for q in curl - factor * full))
    checks.check("axial velocity preserves the actual streamfunction", s.simplify(v.dot(s.Matrix([s.diff(axial, x), s.diff(axial, y)]))) == 0)
    checks.check("full pressure has its declared Bernoulli constant", s.simplify(field.pressure / rho + full.dot(full) / 2 - cap / 2) == 0)
    checks.check("constant-factor substitution is actually false", s.diff(-lam**2 * psi0 / s.sqrt(cap + lam**2 * psi0**2), psi0) != 0)
    profile = s.sqrt(cap + lam**2 * psi0**2)
    first = s.diff(profile, psi0)
    second = s.diff(profile, psi0, 2)
    checks.check("core axial slope derived from fixed C", s.simplify(first - lam**2 * psi0 / profile) == 0)
    checks.check("core axial second derivative is not erased", s.simplify(second - lam**2 * cap / profile**3) == 0)
    checks.check("same-field core normal Hessian", s.hessian(psi, (x, y)).subs({x: 0, y: 0}) == -3 * amp * lam**2 * s.eye(2) / 2)

    r2 = x**2 + y**2
    j = s.Matrix([[0, -1], [1, 0]])
    pos = s.Matrix([x, y])
    swirl = -om * j * pos
    wquad = -a * om * r2 / 2
    radius = s.symbols("r", real=True)
    deviation = -om * radius**2 / 2 + om * lam**2 * radius**4 / 32
    wjet = s.series(a * deviation + b * deviation**2 / 2, radius, 0, 6).removeO()
    checks.check("nonconstant axial quartic is derived by composition", s.expand(wjet.coeff(radius, 4) - a * om * lam**2 / 32 - b * om**2 / 8) == 0)
    checks.check("positive carrier is trapping", p * wquad == -p * a * om * r2 / 2)
    checks.check("negative carrier mutation reverses the well", -p * wquad != p * wquad)
    ell4 = 2 / (a * p**3)
    checks.check("full pressure and Doppler balance determine width", s.simplify(om / (p**2 * s.sqrt(ell4)) - p * a * om * s.sqrt(ell4) / 2) == 0)
    ell2_scaled = s.simplify(s.sqrt(ell4).subs(p, 2 * a / delta**2))
    checks.check("fixed-profile radial scale", ell2_scaled == delta**3 / (2 * a**2))
    checks.check("extra axial nonlinear term first has order delta four", s.simplify((p * b * om * ell4 / 8).subs(p, 2 * a / delta**2)) == b * om * delta**4 / (16 * a**3))
    checks.check("normal C6 first has order delta six", s.simplify((lam**4 * ell4).subs(p, 2 * a / delta**2)) == lam**4 * delta**6 / (4 * a**4))

    vx, vy, pressure = (s.Function(n)(x, y) for n in ("Vx", "Vy", "P"))
    velocity = s.Matrix([vx, vy])
    div = s.diff(vx, x) + s.diff(vy, y)
    bz = s.I * div / p
    gradp = s.Matrix([s.diff(pressure, x), s.diff(pressure, y)])
    gradw = s.Matrix([s.diff(wquad, x), s.diff(wquad, y)])
    rhs = -velocity.jacobian((x, y)) * swirl - swirl.jacobian((x, y)) * velocity - s.I * p * wquad * velocity - gradp
    normal = s.I * (s.diff(rhs[0], x) + s.diff(rhs[1], y)) / p
    axial_rhs = -swirl[0] * s.diff(bz, x) - swirl[1] * s.diff(bz, y) - s.I * p * wquad * bz - velocity.dot(gradw) - s.I * p * pressure
    expected = p**2 * pressure - s.diff(pressure, x, 2) - s.diff(pressure, y, 2) - 2 * om * (s.diff(vy, x) - s.diff(vx, y)) - 2 * s.I * p * gradw.dot(velocity)
    checks.check("generalized-core pressure follows from all Euler equations", s.simplify((normal - axial_rhs) * p / s.I - expected) == 0)
    checks.check("omitting axial-shear pressure is exposed", s.simplify((normal - axial_rhs) * p / s.I - expected - 2 * s.I * p * gradw.dot(velocity)) != 0)
    checks.check("axial reflection preserves both pressure sources", s.simplify(2 * s.I * (-p) * (-gradw).dot(velocity) - 2 * s.I * p * gradw.dot(velocity)) == 0)

    q, w = s.symbols("q w", positive=True)
    lag = s.assoc_laguerre(2, 1, q)
    radial_mode = radius * s.exp(-radius**2 / 2) * lag.subs(q, radius**2)
    oscillator = -s.diff(radial_mode, radius, 2) - s.diff(radial_mode, radius) / radius + radial_mode / radius**2 + radius**2 * radial_mode
    checks.check("oscillator eigenvalue derived from radial differential equation", s.simplify(oscillator - 12 * radial_mode) == 0)
    real, imag = s.symbols("vr vi", real=True)
    plus = s.Matrix([1, s.I]) / s.sqrt(2)
    vc = plus * (real + s.I * imag)
    xi = s.I * vc / (2 * om)
    checks.check("principal Kelvin inversion is local vorticity, not curl eigenvalue", s.simplify(2 * om * j * xi - vc) == s.zeros(2, 1))
    xi3, vc3 = (obj.col_join(s.Matrix([0])) for obj in (xi, vc))
    u3 = swirl.col_join(s.Matrix([0]))
    physical_spin = (xi3.cross(u3) + s.Matrix([x, y, 0]).cross(vc3 + u3.jacobian((x, y)) * xi))[2]
    checks.check("literal principal material spin cancels pointwise", s.simplify(physical_spin) == 0)
    kks = rho * s.Matrix([0, 0, -2 * om]).dot(s.re(xi3).cross(s.im(xi3)))
    checks.check("complete action leading KKS density sign and factor", s.simplify(kks + rho * (real**2 + imag**2) / (4 * om)) == 0)
    checks.check("axial-shear KKS correction starts at delta squared", s.simplify((a / p).subs(p, 2 * a / delta**2)) == delta**2 / 2)
    norm = sum(coef * s.factorial(degree[0]) for degree, coef in s.Poly(q * lag**2, q).terms())
    checks.check("new mode has full radial norm three pi", norm == 3)
    laplace = sum(coef * s.factorial(degree[0] + 2) * (2 / w) ** (degree[0] + 3) for degree, coef in s.Poly(lag, q).terms())
    slope = s.simplify(-s.diff(laplace, w).subs(w, 1) / laplace.subs(w, 1))
    checks.check("actual material Laplace slope, not eigenfrequency", slope == s.Rational(19, 3))
    gamma = 2 * om + om * delta * (slope - 6)
    sq = s.series(gamma**2, delta, 0, 2).removeO()
    curvature = s.expand(delta**2 * s.diff(sq, delta, 2) / 4 + 3 * delta * s.diff(sq, delta) / 4)
    checks.check("positive fixed-C physical clock curvature", curvature == om**2 * delta)
    pressure_poly = s.expand(lag - 2 * s.diff(lag, q))
    checks.check("physical pressure-spin polynomial independently derived", pressure_poly == q**2 / 2 - 5 * q + 9)
    checks.check("displacement-current row is not spin row", s.simplify(pressure_poly - lag) != 0)

    def derivative(poly):
        return s.expand(-poly / 2 + s.Rational(3, 2) * q * (s.diff(poly, q) - poly / 2))

    rows = []
    for seed in (pressure_poly, q * pressure_poly):
        rows.extend((seed, derivative(seed), derivative(derivative(seed))))
    minor = s.Matrix([[s.expand(poly).coeff(q, h) for h in range(6)] for poly in rows]).det()
    checks.check("six new-field spin time/carrier rows independent", minor != 0)
    jets = [[s.factorial(h) if j0 == h else 0 for j0 in range(8)] for h in range(2)]
    jets.extend([[sum(s.binomial(j0, h) * (-s.Rational(1, 2)) ** (j0 - h) * s.factorial(h) * s.expand(poly).coeff(q, h) for h in range(j0 + 1)) for j0 in range(8)] for poly in rows])
    full_minor = s.Matrix(jets).det(method="domain-ge")
    checks.check("reference and spin rows have independent exact minor", full_minor != 0)
    cutoff = s.expand((q / 2 - slope) * lag)
    cutrows = [cutoff]
    for _ in range(2):
        cutrows.append(s.expand(s.Rational(3, 2) * q * (s.diff(cutrows[-1], q) - cutrows[-1] / 2)))
    checks.check("fixed-tag curvature needs and has three independent rows", s.Matrix([[row.coeff(q, h) for h in range(6)] for row in cutrows]).rank() == 3)
    homogeneous = [sum(coef * s.factorial(deg[0] + 2) * 2 ** (deg[0] + 3) for deg, coef in s.Poly(row, q).terms()) for row in cutrows]
    checks.check("all three ideal physical tag jets vanish", homogeneous == [0, 0, 0])
    integer = s.symbols("n", integer=True, nonzero=True)
    checks.check("full-cell distinct-harmonic cross vanishes exactly", s.simplify(s.integrate(s.exp(s.I * integer * z), (z, 0, 2 * s.pi))) == 0)
    print(f"a={first}; b={s.simplify(second)}; ell4={ell4}")
    print(f"physical gamma={gamma}; p²(gamma²)''={curvature}")
    print(f"six-row minor={minor}; eight-row minor={full_minor}")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

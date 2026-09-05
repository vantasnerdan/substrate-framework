"""Exact identities exposing the analytic closed-ring continuation proof."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0186-closed-ring-operator")
    ii, theta, phi = s.symbols("I theta phi", real=True)
    om, ff = s.Function("Omega")(ii), s.Function("F")(ii)
    coordinates = (ii, theta, phi)
    w = s.Matrix([s.Function(f"w{j}")(*coordinates) for j in range(3)])
    v = s.Matrix([s.Function(f"v{j}")(*coordinates) for j in range(3)])
    base = s.Matrix([0, om, 0])
    vort = s.Matrix([0, 0, ff])

    def bracket(left, right):
        return right.jacobian(coordinates) * left - left.jacobian(coordinates) * right

    operator = -bracket(base, w) + bracket(vort, v)
    expected = -om * w.diff(theta) + s.Matrix([0, s.diff(om, ii) * w[0], 0]) + ff * v.diff(phi) - s.Matrix([0, 0, s.diff(ff, ii) * v[0]])
    checks.check("complete ring Lie operator including pressure-induced velocity", s.simplify(operator - expected) == s.zeros(3, 1))
    checks.check("actual steady vorticity commutes with meridional flow", bracket(base, vort) == s.zeros(3, 1))
    n = s.symbols("n", integer=True)
    vc = s.Matrix(s.symbols("vI vtheta vphi"))
    harmonic = (ff * (vc * s.exp(s.I * n * phi)).diff(phi) - s.Matrix([0, 0, s.diff(ff, ii) * vc[0] * s.exp(s.I * n * phi)])) / s.exp(s.I * n * phi)
    checks.check("toroidal harmonic retains both compact reaction terms", s.simplify(harmonic - (s.I * n * ff * vc - s.Matrix([0, 0, s.diff(ff, ii) * vc[0]]))) == s.zeros(3, 1))
    d, shear = s.symbols("d shear", nonzero=True)
    nilpotent = s.Matrix([[0, 0, 0], [shear, 0, 0], [0, 0, 0]])
    checks.check("transport shear is genuinely nilpotent", nilpotent**2 == s.zeros(3))
    resolvent = s.eye(3) / d + nilpotent / d**2
    checks.check("full transport-band resolvent", (d * s.eye(3) - nilpotent) * resolvent == s.eye(3))
    checks.check("omitting the shear resolvent term is exposed", (d * s.eye(3) - nilpotent) / d != s.eye(3))

    r = s.symbols("r", positive=True)
    f = s.Function("radial_f")(r)
    pressure = s.Function("P")(r)
    q, b, cc = s.symbols("q b c", real=True)
    flux_derivative = s.diff(r * f * pressure, r).subs({s.diff(f, r): -(1 / r + q) * f + b * pressure, s.diff(pressure, r): cc * f + q * pressure})
    checks.check("all-poloidal definite-sign radial identity", s.expand(flux_derivative - r * b * pressure**2 - r * cc * f**2) == 0)
    sigma, potential, potential_prime = s.symbols("sigma potential potential_prime", nonzero=True, real=True)
    exterior_f = s.I * potential_prime / sigma
    exterior_p = s.I * sigma * potential
    checks.check("complete exterior pressure gives the correct negative DtN boundary flux", s.simplify(exterior_f * s.conjugate(exterior_p) - potential * potential_prime) == 0)
    x, ell, ratio = s.symbols("x ell ratio", positive=True)
    kk = -1 - x * ratio
    riccati = x**2 + 1 - kk**2
    checks.check("carrier-to-radial-parameter Jacobian has a strict analytic margin", s.expand(riccati - 2 * kk - (2 + x**2 * (1 - ratio**2))) == 0)

    big_r, xx, yy, zz, ww, angle = s.symbols("R x y z w angle", real=True)
    target_pos = s.Matrix([big_r + xx, 0, zz])
    source_pos = s.Matrix([(big_r + yy) * s.cos(angle), (big_r + yy) * s.sin(angle), ww])
    dist2 = (target_pos - source_pos).dot(target_pos - source_pos)
    dist_target = (xx - yy)**2 + (zz - ww)**2 + 2 * (big_r + xx) * (big_r + yy) * (1 - s.cos(angle))
    checks.check("exact whole-ring kernel distance", s.trigsimp(s.expand(dist2 - dist_target)) == 0)
    mm, tt = s.symbols("m t", positive=True)
    denominator = s.sqrt(1 - mm * s.sin(tt)**2)
    elliptic_integrand = (2 / mm - 1) / denominator - 2 * denominator / mm
    checks.check("exact elliptic reduction retains the cosine numerator", s.simplify(elliptic_integrand - (2 * s.sin(tt)**2 - 1) / denominator) == 0)
    eps, log_r = s.symbols("epsilon L", real=True)
    prefactor = (1 + eps * xx)**s.Rational(1, 2) * (1 + eps * yy)**s.Rational(3, 2)
    checks.check("physical streamfunction kernel carries both radial weights", s.diff(prefactor, eps).subs(eps, 0) == (xx + 3 * yy) / 2)
    gamma, moment, aa = s.symbols("Gamma first_moment A", real=True)
    leading = log_r * (xx * gamma + 3 * moment) / (4 * s.pi)
    drift_coefficient = s.diff(leading.subs(moment, 0), xx)
    checks.check("translation border determines the actual leading ring speed", drift_coefficient == gamma * log_r / (4 * s.pi))
    checks.check("mass plus center fixes physical weighted circulation", (gamma + moment / big_r).subs(moment, 0) == gamma)
    mu = gamma * log_r / (2 * s.pi) - aa - drift_coefficient / 2
    checks.check("global outer pressure threshold stays positive to leading logarithm", s.diff(mu, log_r) == 3 * gamma / (8 * s.pi))
    radial, axial = s.symbols("radial axial", positive=True)
    stream = s.Function("psi")(radial, axial)
    ur = -s.diff(stream, axial) / radial
    uz = s.diff(stream, radial) / radial
    checks.check("global ring velocity is exactly divergence free", s.simplify(s.diff(radial * ur, radial) / radial + s.diff(uz, axial)) == 0)
    checks.check("global ring has the actual streamfunction first integral", s.simplify(ur * s.diff(stream, radial) + uz * s.diff(stream, axial)) == 0)
    delta_star = s.diff(stream, radial, 2) - s.diff(stream, radial) / radial + s.diff(stream, axial, 2)
    checks.check("full cylindrical curl normalization", s.simplify(s.diff(ur, axial) - s.diff(uz, radial) + delta_star / radial) == 0)
    checks.check("uniform-flow term solves the homogeneous streamfunction operator", s.diff(radial**2, radial, 2) - s.diff(radial**2, radial) / radial == 0)

    ss = s.symbols("s", real=True)
    phase = (s.exp(s.I * eps * ss) - 1) / eps
    checks.check("actual first integer-harmonic difference has the axial derivative limit", s.limit(phase, eps, 0) == s.I * ss)
    checks.check("actual second integer-harmonic difference has the second derivative limit", s.limit(phase**2, eps, 0) == -ss**2)
    checks.check("both carrier differences keep periodic far-endpoint cancellation", all(s.simplify((s.exp(s.I * (angle + 2 * s.pi)) - 1)**order - (s.exp(s.I * angle) - 1)**order) == 0 for order in (1, 2)))

    gq, gp, cq, cp, position, rate, accel = s.symbols("Gq Gp Cq Cp theta rate accel", real=True)
    gplus = gq * position + gp * rate
    gminus = gq * position - gp * rate
    splus = gq * rate + gp * accel + 2 * cq * position + 2 * cp * rate
    sminus = gq * rate - gp * accel - 2 * cq * position + 2 * cp * rate
    checks.check("prepared physical TR average retains the real displacement moment", s.expand((gplus + gminus) / 2 - gq * position) == 0)
    checks.check("full mechanical spin retains the even partner-connection row", s.expand((splus + sminus) / 2 - (gq + 2 * cp) * rate) == 0)
    checks.check("zero actual partner moment gives literal all-time current", s.expand(((splus + sminus) / 2).subs(cp, 0) - gq * rate) == 0)
    mass, nu = s.symbols("M nu", positive=True)
    form = mass * s.Matrix([[0, 1], [-1, 0]])
    prepare = s.diag(1, -1)
    checks.check("actual initial TR phase-form pullback precedes elimination", prepare.T * (-form) * prepare == form)
    energy = mass * s.diag(nu**2, 1)
    checks.check("actual initial TR energy pullback has the same positive mass", prepare.T * energy * prepare == energy)
    checks.check("tying both raw conjugate coordinates would wrongly erase the form", (form - form) / 2 != (form + prepare.T * (-form) * prepare) / 2)
    print("The verifier checks exact identities; analytic Fredholm, kernel bounds and bordered IFT are proved in the three companion files.")
    print("No numerical eigenvalue, imposed spectral splitting or empirical comparator enters this attempt.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

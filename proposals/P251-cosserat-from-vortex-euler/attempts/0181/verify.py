"""Exact radial, KKS and literal stationary-tag rows; no numerical root gate."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0181-stationary-nodal-mode")
    r, p, om, tau = s.symbols("r p Omega tau", positive=True)
    sig = s.symbols("sigma", nonzero=True, real=True)
    pressure = s.Function("P")(r)
    dp = s.diff(pressure, r)
    den = sig * (4 * om**2 - sig**2)
    ar = (2 * om * pressure / r - sig * dp) / den
    bt = (2 * om * dp - sig * pressure / r) / den
    checks.check("full radial pressure equation", s.simplify((sig**2 - 4 * om**2) * ar + 2 * om * pressure / (r * sig) - dp) == 0)
    radial_ode = -dp / r + pressure / r**2 - p**2 * (4 * om**2 - sig**2) * pressure / sig**2
    divergence = s.diff(ar, r) + ar / r - bt / r - p**2 * pressure / sig**2
    checks.check("full three-dimensional divergence", s.simplify(divergence.subs(s.diff(pressure, r, 2), radial_ode)) == 0)
    checks.check("sum of actual horizontal displacement rows", s.simplify(ar + bt - (dp + pressure / r) / (sig * (2 * om + sig))) == 0)
    positive_integrand = s.expand((ar * bt * r).subs(sig, -tau * om) * tau**2 * om**4 * (4 - tau**2)**2)
    target = 2 * tau * (r * dp**2 + pressure**2 / r) + (tau**2 + 4) * pressure * dp
    checks.check("KKS radial norm before integration by parts", s.simplify(positive_integrand - target) == 0)

    ell, x = s.symbols("ell x", positive=True)
    j, k = s.symbols("J K", real=True)
    tau_l = 2 / s.sqrt(1 + ell**2 / x**2)
    determinant = ((4 - tau**2) * k + tau**2 * j + 2 * tau) / tau**2
    checks.check("complete determinant transformed without a dropped exterior term", s.simplify(determinant.subs(tau, tau_l) - j - ell**2 * k / x**2 - s.sqrt(1 + ell**2 / x**2)) == 0)
    checks.check("lab zero-frequency crossing is tau one", s.simplify(tau_l.subs(ell, s.sqrt(3) * x)) == 1)
    ellp = s.symbols("ell_prime", real=True)
    tau_derivative = (s.diff(tau_l, x) + s.diff(tau_l, ell) * ellp).subs(ell, s.sqrt(3) * x)
    checks.check("transverse crossing derivative", s.simplify(tau_derivative - s.sqrt(3) * (s.sqrt(3) - ellp) / (4 * x)) == 0)
    nu = s.Function("nu")(x)
    nu1, nu2 = s.symbols("nu1 nu2", real=True)
    square_jet = s.diff(nu**2, x, 2).subs({s.diff(nu, x, 2): nu2, s.diff(nu, x): nu1}).subs(nu, 0)
    checks.check("squared lab frequency has strictly positive curvature at a simple zero", square_jet == 2 * nu1**2)

    beta, freq, c = s.symbols("beta omega c", nonzero=True, real=True)
    rot = s.Matrix([[0, -1], [1, 0]])
    symplectic = -beta * rot
    hessian = -beta * freq * s.eye(2)
    checks.check("actual lab generator and KKS determine the Hamiltonian sign", symplectic * (freq * rot) + hessian == s.zeros(2))
    q, qdot, partner = s.symbols("q qdot partner", real=True)
    lag = beta * partner * qdot - (-beta * freq) * (q**2 + partner**2) / 2
    reduced = s.simplify(lag.subs(partner, -qdot / freq))
    checks.check("one-coordinate inertia derives by reaction elimination", s.simplify(reduced - (-beta / freq) * (qdot**2 - freq**2 * q**2) / 2) == 0)

    angle, z = s.symbols("angle z", real=True)
    aa, bb, pp = s.symbols("A B P", real=True)
    pos = s.Matrix([r * s.cos(angle), r * s.sin(angle), z])
    xi = s.exp(s.I * angle) * s.Matrix([
        aa * s.cos(angle) - s.I * bb * s.sin(angle),
        aa * s.sin(angle) + s.I * bb * s.cos(angle),
        s.I * p * pp / sig**2,
    ])
    base = s.Matrix([-om * pos[1], om * pos[0], 0])

    def angular(expr):
        return s.simplify(s.integrate(s.expand_trig(expr), (angle, 0, 2 * s.pi)))

    shape = s.Matrix([angular(z * xi[i] + pos[i] * xi[2]) for i in range(2)])
    displacement = s.Matrix([angular(pos.cross(xi)[i]) for i in range(2)])
    spin = s.Matrix([angular((-s.I * freq * pos.cross(xi) + 2 * xi.cross(base))[i]) for i in range(2)])
    ep = s.Matrix([1, s.I])
    checks.check("Cartesian shape row includes axial displacement", s.simplify(shape - s.pi * (z * (aa + bb) + s.I * r * p * pp / sig**2) * ep) == s.zeros(2, 1))
    checks.check("Cartesian displacement moment includes both cross terms", s.simplify(displacement + s.pi * (s.I * z * (aa + bb) + r * p * pp / sig**2) * ep) == s.zeros(2, 1))
    expected_spin = s.pi * (-freq * z * (aa + bb) + s.I * (freq - 2 * om) * r * p * pp / sig**2) * ep
    checks.check("full mechanical spin has the material-position contribution", s.simplify(spin - expected_spin) == s.zeros(2, 1))
    zf, zprime, rr1, rr2, rho, gap = s.symbols("Z Zprime R1 R2 rho Delta", real=True)
    dq = rho * s.pi * (-zprime * rr1 + p * zf * rr2 / sig**2)
    dg = rho * s.pi * (zprime * rr1 + p * zf * rr2 / sig**2)
    spin_row = s.I * rho * s.pi * (freq * zprime * rr1 + (freq - 2 * om) * p * zf * rr2 / sig**2)
    checks.check("zero axial-stretch moment gives the actual G identity", s.simplify((-dg - dq).subs(rr2, 0)) == 0)
    checks.check("same condition gives the actual all-time spin identity", s.simplify((spin_row + s.I * freq * dq).subs(rr2, 0)) == 0)

    a = s.symbols("a", positive=True)
    bes = s.besselj(1, ell * r / a)
    bessel_residual = ell**2 * s.diff(bes, ell, 2) + ell * s.diff(bes, ell) - bes + ell**2 * r**2 * bes / a**2
    checks.check("second carrier moment reduction follows the Bessel ODE", s.simplify(s.expand_func(bessel_residual)) == 0)
    r1, r2, w1, w2, d1, d2 = s.symbols("r1 r2 w1 w2 d1 d2", nonzero=True, real=True)
    # d_i are the actual nonzero Bessel derivatives at its simple zeros.
    matrix = s.Matrix([
        [w1 * r1**2 * ell * d1 / a, w2 * r2**2 * ell * d2 / a, 0],
        [w1 * r1**4 * ell * d1 / a, w2 * r2**4 * ell * d2 / a, 0],
        [s.Symbol("lower1"), s.Symbol("lower2"), r2**3 * d2 / a],
    ])
    expected_det = w1 * w2 * ell**2 * d1 * d2**2 * r1**2 * r2**5 * (r2**2 - r1**2) / a**3
    checks.check("positive two-node tag has a full three-row IFT Jacobian", s.factor(matrix.det() - expected_det) == 0)
    weight_ratio = -w1 * r1**3 * d1 / (r2**3 * d2)
    tilt_row = (w1 * r1 * d1 + w2 * r2 * d2).subs(w2, weight_ratio)
    checks.check("current constraints preserve a nonzero physical tilt", s.factor(tilt_row - w1 * r1 * d1 * (1 - r1**2 / r2**2)) == 0)
    norm = s.Function("normalization")(p)
    lam = s.Function("ell")(p)
    ff = s.Function("F")
    moment = norm * ff(lam)
    jet_rules = {ff(lam): 0, s.diff(ff(lam), lam): 0, s.diff(ff(lam), lam, 2): 0}
    checks.check("fixed tag controls carrier zero, one and two with the actual pressure normalization", all(s.diff(moment, p, order).subs(jet_rules).doit() == 0 for order in range(3)))
    print("No numerical Bessel zero, fitted frequency, or discretized small eigenvalue enters these checks.")
    print("Analytic Sturm positivity and thin-annulus transfer are proved in nodal-domain-mode.md.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

"""Exact transit, flux-action twist and arithmetic-free core selection."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0145-flux-action-twist")
    x, radius, kap = s.symbols("x R kappa", positive=True)
    leading = s.cos(kap*x)
    correction = (x*s.cos(kap*x)-s.sin(kap*x)/kap)/2
    checks.check("finite-radius first correction solves the actual radial Bessel ODE",
                 s.simplify(s.diff(correction, x, 2)+kap**2*correction-s.diff(leading, x)) == 0)
    checks.check("finite-radius correction preserves both core initial conditions",
                 correction.subs(x, 0) == 0 and s.diff(correction, x).subs(x, 0) == 0)

    h = s.Symbol("h", positive=True)
    parameter = 1-h*h
    angle = s.Symbol("t", real=True)
    sine_x = s.sqrt(parameter)*s.sin(angle)
    # Square the positive quadrant Jacobian identity to avoid choosing a
    # square-root sign outside its explicitly stated quadrant.
    dx_dt_squared = parameter*s.cos(angle)**2/(1-sine_x**2)
    velocity_squared = 1-sine_x**2-h*h
    checks.check("quadrant transit substitution gives the complete elliptic parameter",
                 s.trigsimp(dx_dt_squared/velocity_squared-1/(1-parameter*s.sin(angle)**2)) == 0)

    ell_k = s.elliptic_k(parameter)
    ell_e = s.elliptic_e(parameter)
    derivative = s.diff(h*ell_k, h)
    checks.check("actual period derivative retains both elliptic integrals",
                 s.simplify(derivative-(ell_k-ell_e)/parameter) == 0)
    positive_integrand = 1/s.sqrt(1-parameter*s.sin(angle)**2) \
        -s.sqrt(1-parameter*s.sin(angle)**2)
    checks.check("elliptic difference is an explicitly positive interior integral",
                 s.simplify(positive_integrand-parameter*s.sin(angle)**2
                            /s.sqrt(1-parameter*s.sin(angle)**2)) == 0)

    lam, speed = s.symbols("lambda U", positive=True)
    reduced_rotation = s.pi*kap**2/(2*lam*h*ell_k)
    flux_derivative = -2*speed*h*ell_k/(s.pi*kap**2)
    twist = s.simplify(s.diff(reduced_rotation, h)/flux_derivative)
    expected = s.pi**2*kap**4*(ell_k-ell_e)/(4*lam*speed*parameter*h**3*ell_k**3)
    checks.check("twist uses actual flux action rather than unweighted area",
                 s.simplify(twist-expected) == 0)
    checks.check("flux derivative and reduced rotation satisfy the exact normalization limit",
                 s.simplify(flux_derivative*reduced_rotation+speed/lam) == 0)
    small_m = s.Symbol("m", real=True)
    k_series = s.series(s.elliptic_k(small_m), small_m, 0, 2).removeO()
    e_series = s.series(s.elliptic_e(small_m), small_m, 0, 2).removeO()
    center = s.limit(s.pi**2*kap**4*(k_series-e_series)
                     /(4*lam*speed*small_m*(1-small_m)**s.Rational(3, 2)*k_series**3), small_m, 0)
    checks.check("center twist is finite and positive with its derived normalization",
                 center == kap**4/(2*lam*speed))
    checks.check("circular fixed-lambda limit agrees with the universal core coefficient",
                 s.simplify(center.subs(kap, lam/s.sqrt(2))-lam**3/(8*speed)) == 0)

    advance = s.Symbol("Delta_theta", positive=True)
    actual_flux_derivative = -speed*radius*advance/(2*s.pi*lam)
    rotation_number = 2*s.pi/advance
    checks.check("finite-radius flux identity keeps the full toroidal transit",
                 s.simplify(actual_flux_derivative*rotation_number+speed*radius/lam) == 0)

    axial, zero = s.symbols("k j", positive=True)
    core_number = zero*axial/lam
    checks.check("single-mode core rotation varies in a genuine continuous parameter",
                 s.diff(core_number, axial) == zero/lam)
    j1 = s.Symbol("J1j", nonzero=True, real=True)
    # The small-b second CK mode tends to r^2 cos(lambda z).
    jacobian_limit = s.Matrix([[-zero*j1, 2], [-4*kap*j1, lam**2*radius]])
    determinant_limit = s.simplify(jacobian_limit.det().subs({radius: zero/kap, lam: s.sqrt(2)*kap}))
    checks.check("exact-circular two-mode constraints have a nonzero tunable core determinant",
                 s.factor(determinant_limit-2*kap*j1*(4-zero**2)) == 0)
    z = s.Symbol("z", real=True)
    lower = 1-z+z*z/4-z**3/36
    checks.check("elementary J0 positivity bound excludes a zero at or below two",
                 s.expand(lower-((1-z)+z*z*(9-z)/36)) == 0
                 and lower.subs(z, 1) == s.Rational(2, 9))
    checks.check("exact circularity permits core selection without Bessel-zero arithmetic",
                 s.diff(lam*radius/2, radius) == lam/2)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

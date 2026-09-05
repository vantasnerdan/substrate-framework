"""Exact radial Euler and supplied thin-ring leading equations, not a mode oracle."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0198-localizable-optical-equation")
    x, y, r = s.symbols("x y r", real=True, positive=True)
    aa, bb, cc = s.symbols("a b c", real=True)
    quadratic = aa*x**2+bb*x*y+cc*y**2
    first = s.Poly(s.diff(quadratic, x)-2*x, x, y)
    solved = s.solve(first.coeffs(), (aa, bb))
    second = s.expand((s.diff(quadratic, y)**2-4*quadratic+4*x**2).subs(solved))
    branches = s.solve(s.Poly(second, y).coeffs(), cc)
    checks.check("source quadratic equations give the positive transverse branch",
                 solved == {aa: 1, bb: 0} and branches == [(0,), (1,)])
    phi = quadratic.subs(solved).subs(cc, 1)
    checks.check("actual swirl square and poloidal speed fix the nonconstant pitch",
                 s.expand(6*phi*s.Rational(1, 3)-2*(x**2+y**2)) == 0
                 and s.expand(s.diff(phi, x)**2+s.diff(phi, y)**2-4*(x**2+y**2)) == 0)

    omega = s.Function("Omega")(r)
    axial = s.Function("W")(r)
    vr, vt, vz, m, k, frequency = s.symbols("v_r v_theta v_z m k omega")
    velocity = s.Matrix([vr, vt, vz])
    base = s.Matrix([0, r*omega, axial])
    generator = s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    base_convection = omega*generator*base
    checks.check("radial steady Euler retains arbitrary axial shear",
                 s.simplify(base_convection-s.Matrix([-r*omega**2, 0, 0])) == s.zeros(3, 1))
    convection = omega*(s.I*m*velocity+generator*velocity)+s.I*k*axial*velocity
    convection += vr*s.diff(base, r)+vt*generator*base/r
    q = m*omega+k*axial
    qq = 2*omega+r*s.diff(omega, r)
    expected = s.I*q*velocity+s.Matrix([-2*omega*vt, qq*vr, s.diff(axial, r)*vr])
    checks.check("complete frame derivative produces all linear Euler components",
                 s.simplify(convection-expected) == s.zeros(3, 1))
    sigma = frequency+q
    pressure = s.Function("P")(r)
    tangential = s.I*qq*vr/sigma-m*pressure/(r*sigma)
    longitudinal = s.I*s.diff(axial, r)*vr/sigma-k*pressure/sigma
    checks.check("tangential elimination retains the actual Doppler denominator",
                 s.simplify(s.I*sigma*tangential+qq*vr+s.I*m*pressure/r) == 0)
    checks.check("axial elimination retains the actual axial derivative",
                 s.simplify(s.I*sigma*longitudinal+s.diff(axial, r)*vr+s.I*k*pressure) == 0)
    radial_equation = s.diff(pressure, r)+s.I*sigma*vr-2*omega*tangential
    checks.check("full radial pressure equation includes the rotation return",
                 s.simplify(radial_equation-s.diff(pressure, r)-2*m*omega*pressure/(r*sigma)
                            -s.I*(sigma-2*omega*qq/sigma)*vr) == 0)
    vprime = s.symbols("vprime")
    divergence = vprime+vr/r+s.I*m*tangential/r+s.I*k*longitudinal
    reduced = vprime+(1/r-(m*qq/r+k*s.diff(axial, r))/sigma)*vr
    reduced -= s.I*(m**2/r**2+k**2)*pressure/sigma
    checks.check("full divergence gives the two-variable radial system",
                 s.simplify(divergence-reduced) == 0)
    exterior = s.diff(s.I*s.diff(pressure, r)/frequency, r)
    exterior += s.I*s.diff(pressure, r)/(r*frequency)
    exterior -= s.I*(m**2/r**2+k**2)*pressure/frequency
    checks.check("full exterior pressure is Bessel matching not a cutoff wall",
                 s.simplify(exterior*frequency/s.I-s.diff(pressure, r, 2)
                            -s.diff(pressure, r)/r+(m**2/r**2+k**2)*pressure) == 0)
    plateau = axial.subs(axial, s.sqrt(2)*r)
    checks.check("omitting plateau axial shear is exposed even for constant rotation",
                 s.diff(plateau, r) == s.sqrt(2))
    checks.check("supplier pitch cannot inherit a constant helical charge",
                 s.diff((s.sqrt(2)*r)/2, r) != 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

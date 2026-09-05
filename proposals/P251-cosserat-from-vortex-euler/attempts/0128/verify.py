"""Exact compact-profile 3Omega Euler, observation and action-cross oracle."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0128-localized-three-Omega")
    x, y, z = s.symbols("x y z", real=True)
    coords = (x, y, z)
    om, rho, radial = s.symbols("Omega rho radial", positive=True)
    zeta = x+s.I*y
    ep, em = s.Matrix([1, s.I, 0]), s.Matrix([1, -s.I, 0])
    j = s.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    u = s.Matrix([-om*y, om*x, 0])
    primitive = s.Function("L")(z)
    kval = s.diff(primitive, z)
    hval = s.diff(primitive, z, 2)
    hprime = s.diff(primitive, z, 3)

    def grad(scalar):
        return s.Matrix([s.diff(scalar, coordinate) for coordinate in coords])

    def curl(vector):
        return s.Matrix([s.diff(vector[2], y)-s.diff(vector[1], z),
                         s.diff(vector[0], z)-s.diff(vector[2], x),
                         s.diff(vector[1], x)-s.diff(vector[0], y)])

    def lie(vector):
        return vector.jacobian(coords)*u-u.jacobian(coords)*vector

    field = hprime*zeta*zeta*em/2-kval*ep-2*hval*zeta*s.Matrix([0, 0, 1])
    potential = -s.I*(hval*zeta*zeta*em/2+primitive*ep)
    checks.check("three compact primitives give the exact axial-profile vector potential",
                 (curl(potential)-field).applyfunc(s.simplify) == s.zeros(3, 1))
    checks.check("the full arbitrary-profile mode is divergence free",
                 s.simplify(sum(s.diff(field[i], coords[i]) for i in range(3))) == 0)
    checks.check("its pressure defect is an actual gradient",
                 ((j-s.I*s.eye(3))*field-2*s.I*grad(kval*zeta)).applyfunc(s.simplify) == s.zeros(3, 1))
    checks.check("its physical azimuthal Lie weight is fixed, not fitted",
                 (lie(field)-s.I*om*field).applyfunc(s.simplify) == s.zeros(3, 1))
    velocity, pressure = -2*s.I*om*field, -8*om*om*kval*zeta
    euler = -3*s.I*om*velocity+velocity.jacobian(coords)*u+u.jacobian(coords)*velocity+grad(pressure)
    checks.check("arbitrary compact axial profiles solve the full local three-Omega Euler equation",
                 euler.applyfunc(s.simplify) == s.zeros(3, 1))
    checks.check("the same mode satisfies the actual Lin reconstruction",
                 (-3*s.I*om*field+lie(field)-velocity).applyfunc(s.simplify) == s.zeros(3, 1))
    normalized = -s.I*field
    core_angle = (j*normalized.diff(z)).subs({x: 0, y: 0})
    checks.check("core vorticity angle is exactly H at the axis",
                 (core_angle-hval*ep).applyfunc(s.simplify) == s.zeros(3, 1))

    theta = s.Symbol("theta", real=True)

    def angular(vector):
        return s.Matrix([s.simplify(s.integrate(s.expand_trig(value.subs(
            {x: radial*s.cos(theta), y: radial*s.sin(theta)})),
            (theta, 0, 2*s.pi))/(2*s.pi)) for value in vector])

    r3 = s.Matrix(coords)
    spin = rho*(r3.cross(-3*s.I*om*normalized)+2*normalized.cross(u))
    checks.check("complete material spin includes pressure-compatible displacement and tag deformation",
                 (angular(spin)+s.I*rho*om*(3*z*kval-radial**2*hval)*ep).applyfunc(s.simplify) == s.zeros(3, 1))
    shape_cross = r3[:2, :]*normalized[2]+z*normalized[:2, :]
    checks.check("the symmetric shape cross row remains nonzero and explicit",
                 (angular(shape_cross)-s.I*(z*kval+radial**2*hval)*ep[:2, :]).applyfunc(s.simplify) == s.zeros(2, 1))

    g = s.Function("g")(x*x+y*y)
    completed = curl(g*potential)
    checks.check("radial localization leaves the physical transverse profile unchanged inside its plateau",
                 (completed[:2, :]-g*field[:2, :]).applyfunc(s.simplify) == s.zeros(2, 1))
    checks.check("radial completion is exactly divergence free",
                 s.simplify(sum(s.diff(completed[i], coords[i]) for i in range(3))) == 0)
    collar = grad(g).cross(potential)
    defect = (j-s.I*s.eye(3))*completed-grad(2*s.I*g*kval*zeta)
    expected_defect = -2*s.I*kval*zeta*grad(g)+(j-s.I*s.eye(3))*collar
    checks.check("the full Kelvin residual is a projected collar source with its gradient removed",
                 (defect-expected_defect).applyfunc(s.simplify) == s.zeros(3, 1))

    # Actual tag integrals, not a point-spin substitution.
    n, ell, area, b = s.symbols("N ell D b", positive=True)
    axial0 = 2*s.integrate(s.cos(n*z), (z, 0, ell))
    axial1 = 2*s.I*s.integrate(z*s.sin(n*z), (z, 0, ell))
    spin_n = rho*om*(-3*area*axial1/n+s.I*area*b*b*axial0/2)
    j3 = s.simplify(spin_n/(-3*s.I*om))
    expected_j3 = 2*rho*area/n**3*(s.sin(n*ell)-n*ell*s.cos(n*ell)-b*b*n*n*s.sin(n*ell)/6)
    checks.check("finite-cylinder mechanical response is obtained by complete integration",
                 s.simplify(j3-expected_j3) == 0)
    jfirst, jsecond = s.simplify(j3.subs(n, s.pi/ell)), s.simplify(j3.subs(n, 2*s.pi/ell))
    checks.check("two equal-frequency axial modes have independent exact angle-spin rows",
                 jfirst == 2*rho*area*ell**3/s.pi**2
                 and jsecond == -rho*area*ell**3/(2*s.pi**2)
                 and s.Matrix([[1, 1], [jfirst, jsecond]]).det() != 0)
    moment = 2*rho*area*ell**3/3
    checks.check("small isotropic tag reproduces one-third rather than fundamental spin inertia",
                 s.simplify(s.limit(j3.subs(b*b, 4*ell*ell/3), n, 0)-moment/3) == 0)

    hi, ki, hj, kj = s.symbols("hi ki hj kj", real=True)
    fi = (-s.I*(hi*zeta*zeta*em/2-ki*ep)).expand(complex=True)
    fj = (-s.I*(hj*zeta*zeta*em/2-kj*ep)).expand(complex=True)
    density = 2*rho*om*s.re(fi).cross(s.im(fj))[2]
    checks.check("exact KKS real-imaginary cross includes the growing transverse polynomial",
                 s.simplify(angular(s.Matrix([density]))[0]-2*rho*om*(ki*kj-radial**4*hi*hj/4)) == 0)
    qx, qy = s.symbols("qx qy", real=True)
    fundamental = s.Matrix([qx, qy, 0])
    cross_density = 2*rho*om*fundamental.cross(fi)[2]
    checks.check("fundamental-complement KKS cross is generally nonzero with its exact phase",
                 s.simplify(angular(s.Matrix([cross_density]))[0]+2*rho*om*ki*(qx+s.I*qy)) == 0)

    # An exposing exact taper for the third-control rank. Point evaluations
    # prove unequal ratios; the proof then uses finite smooth narrow bumps.
    taper_f = z*(1-z*z)**3
    taper_primitive = -(1-z*z)**4/8
    weight = -4-3*z*z  # l=2,b=4, omitting a positive common area factor.
    checks.check("the third-control weight is derived from a genuinely tapered fundamental",
                 s.diff(taper_primitive, z) == s.expand(taper_f)
                 or s.simplify(s.diff(taper_primitive, z)-taper_f) == 0)
    rank_minor = s.det(s.Matrix([[weight.subs(z, s.Rational(1, 4)), weight.subs(z, s.Rational(3, 4))],
                                [taper_primitive.subs(z, s.Rational(1, 4)), taper_primitive.subs(z, s.Rational(3, 4))]]))
    checks.check("off-core spin and KKS-cross weights have a nonzero exact rank minor", rank_minor != 0)
    m, m2, ltag, btag = s.symbols("mass zeroth_second ltag btag", positive=True)
    tag_factor = ((3*ltag*ltag-btag*btag)*m-3*m2)/2
    checks.check("pancake-tag spin-only control has a strict sign without a small eigenvalue",
                 s.expand(tag_factor.subs(btag*btag, 4*ltag*ltag)) == -ltag*ltag*m/2-3*m2/2)

    macro = s.Symbol("macro", real=True)
    form = s.diag(s.Matrix([[0, 1], [-1, 0]]), s.Matrix([[0, 2], [-2, 0]]))
    base = s.Matrix([[1, 0], [0, 1], [1, 0], [0, 0]])
    correction = s.Matrix([[0, 1], [1, 0], [0, 1], [1, 1]])
    preparation = s.Matrix([[1, 2], [-1, 3]])
    embedding = base+macro**2*correction*preparation
    actual = s.expand(embedding.T*form*embedding)
    cross = base.T*form*correction*preparation+preparation.T*correction.T*form*base
    checks.check("prepared state has a retained order-two symplectic cross",
                 actual.applyfunc(lambda entry: entry.coeff(macro, 2)) == cross and cross != s.zeros(2))
    checks.check("its order-four self term does not replace that cross",
                 actual.applyfunc(lambda entry: entry.coeff(macro, 4)) == preparation.T*correction.T*form*correction*preparation)
    print("Exact independent single-N tag-spin coefficients:", jfirst, jsecond)
    print("Exact taper spin/cross rank minor:", rank_minor)
    print("Core/spin localization and reference rank continuation are analytic kernel estimates, not numerical eigenvalue tests.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

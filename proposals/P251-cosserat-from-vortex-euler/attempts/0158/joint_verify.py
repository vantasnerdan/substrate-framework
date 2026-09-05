"""Actual joint mean moments, physical brackets, mass and tagged current."""

import sympy as s

from substrate_framework.euler_fourier import (
    ZERO, add, cross, curl, derivative, inner, leray, mul, trig,
)
from substrate_framework.euler_phase import moving_phase_pullback
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0158-actual-joint-physical-interface")
    rho, beta, j = s.symbols("rho beta j", positive=True)
    m1, m2, n1, n2, b, d, c1, c2 = s.symbols("m1 m2 n1 n2 b d c1 c2", real=True)
    jj = s.Matrix([[0, 1], [-1, 0]])
    cross_block = rho*s.Matrix([[n1, n2], [-m1, -m2]])
    omega = (rho*jj).row_join(cross_block).col_join((-cross_block.T).row_join(beta*jj))
    xrow = s.Matrix([[1, 0, m1, m2]])
    vrow = s.Matrix([[0, 1, n1, n2]])
    theta = s.Matrix([[b, 0, c1, c2]])

    def bracket(left, right):
        return s.simplify(-(left*omega.inv()*right.T)[0])

    checks.check("actual mean displacement and velocity are a canonical initial subpair",
                 bracket(xrow, vrow) == 1/rho)
    checks.check("actual initial material angle commutes with the actual mean position",
                 bracket(xrow, theta) == 0)
    checks.check("omitting the complete optical means fabricates a position bracket",
                 bracket(s.Matrix([[1, 0, 0, 0]]), theta) != 0)
    expected_determinant = rho**2*(beta+rho*(n1*m2-n2*m1))**2
    checks.check("the full joint phase rank uses the actual Schur correction sign",
                 s.factor(omega.det()-expected_determinant) == 0)

    canonical = s.diag(rho*jj, j*jj)
    q = s.Matrix([[1, 0, 0, 0], [b, 0, 1, 0]])
    rate = s.Matrix([[0, 1, 0, 0], [d, b, 0, 1]])
    transform = q.col_join(rate)
    normal = -q*canonical.inv()*rate.T
    mass = s.simplify(normal.inv())
    checks.check("the true initial position/rate chart stays invertible for measured b,d",
                 transform.det() == -1)
    checks.check("the complete physical initial mass is positive with determinant rho*j",
                 mass == s.Matrix([[rho+j*b*b, -j*b], [-j*b, j]])
                 and s.factor(mass.det()) == rho*j)
    pullback = moving_phase_pullback(canonical, s.zeros(4), transform.inv(), s.zeros(4))
    checks.check("the initial physical symplectic chart retains its measured gyroscopic row",
                 pullback.symplectic[:2, :2] == j*d*jj
                 and pullback.symplectic[:2, 2:] == mass
                 and pullback.symplectic[2:, 2:] == s.zeros(2))

    # Full Fourier pressure oracle for the actual potential moments.
    wave = s.Symbol("k", nonzero=True, real=True)
    macro = (3*wave/5, 0, 4*wave/5)
    kv = s.Matrix(macro)
    pk = s.eye(3)-kv*kv.T/kv.dot(kv)
    background = (add(trig(2), trig(1, kind="sin")),
                  add(trig(0), trig(2, kind="sin")),
                  add(trig(1), trig(0, kind="sin")))
    vort = curl(background)
    potential = (add({ZERO: 1}, trig(0, kind="sin"), trig(2)),
                 add({ZERO: 2}, trig(1), trig(2, kind="sin")),
                 add({ZERO: 3}, trig(0), trig(1, kind="sin")))
    lifted = tuple(mul({macro: 1}, component) for component in potential)
    generator = curl(lifted)
    actual_velocity = leray(cross(generator, vort))
    actual_mean = s.Matrix([component.get(macro, 0) for component in generator])
    actual_velocity_mean = s.Matrix([s.cancel(component.get(macro, 0))
                                     for component in actual_velocity])
    potential_mean = s.Matrix([component.get(ZERO, 0) for component in potential])
    derivative_moment = s.Matrix([inner(potential, tuple(derivative(component, axis)
                                                       for component in vort))
                                  for axis in range(3)])
    product_moment = s.Matrix(3, 3, lambda row, col: mul(potential[row], vort[col]).get(ZERO, 0))
    checks.check("actual Bloch curl derives its complete displacement mean",
                 s.simplify(actual_mean-s.I*kv.cross(potential_mean)) == s.zeros(3, 1))
    checks.check("full Euler Leray gives the derivative and product potential moments",
                 s.simplify(actual_velocity_mean-pk*(derivative_moment+s.I*product_moment*kv))
                 == s.zeros(3, 1))
    checks.check("the product-moment pressure term is genuinely exercised",
                 s.simplify(pk*product_moment*kv) != s.zeros(3, 1))

    # A nondegenerate CK cone actually spans the three helical components.
    z, angle = s.symbols("z angle", real=True)
    radial = s.sqrt(1-z*z)
    normal_cone = s.Matrix([radial*s.cos(angle), radial*s.sin(angle), z])
    helical = s.Matrix([-s.sin(angle)-s.I*z*s.cos(angle),
                       s.cos(angle)-s.I*z*s.sin(angle), s.I*radial])
    checks.check("the actual cone polarization has the correct curl helicity",
                 s.simplify(s.I*normal_cone.cross(helical)-helical) == s.zeros(3, 1))
    cone_basis = s.Matrix.hstack(*(helical.subs(angle, value) for value in (0, s.pi, s.pi/2)))
    cone_det = s.factor(cone_basis.det())
    checks.check("three actual cone polarizations span the control value space",
                 s.simplify(cone_det/(2*s.I*radial*(1-z*z))) in (1, -1))
    k1, k2 = s.symbols("k1 k2", positive=True)
    roots = [k1, -k1, k2, -k2]
    vandermonde = s.Matrix([[value**power for power in range(4)] for value in roots])
    checks.check("two distinct signed CK cones exceed a cubic eigenvalue relation",
                 s.factor(vandermonde.det()/(k1*k2*(k1-k2)**2*(k1+k2)**2)) in (4, -4))

    # Literal tag integrals expose the initial displacement/spin distinction.
    radius, profile, eps, marker, rotation = s.symbols("r f eps btag Omega", real=True)
    xx, yy = radius*s.cos(angle), radius*s.sin(angle)
    xi1 = s.Matrix([xx*profile, -yy*profile, 0])/s.sqrt(2)
    xi2 = s.Matrix([yy*profile, xx*profile, 0])/s.sqrt(2)
    position = s.Matrix([xx, yy, 0])
    base_velocity = rotation*s.Matrix([-yy, xx, 0])
    fraction = 1+eps*marker*s.cos(2*angle)
    g1 = s.integrate(fraction*position.cross(xi1)[2], (angle, 0, 2*s.pi))
    g2 = s.integrate(fraction*position.cross(xi2)[2], (angle, 0, 2*s.pi))
    gyroscopic = s.integrate(fraction*xi1.cross(base_velocity)[2], (angle, 0, 2*s.pi))
    checks.check("the actual initial G row is the independent displacement-helicity moment",
                 g1 == 0 and s.simplify(g2-s.pi*eps*marker*radius**2*profile/s.sqrt(2)) == 0)
    checks.check("the literal moving spin retains its order-one background gyroscopic term",
                 s.simplify(gyroscopic-rotation*g2) == 0 and gyroscopic != 0)
    print("Derived complete joint Schur determinant:", s.factor(omega.det()))
    print("Derived physical initial mass:", mass)
    print("Derived actual CK cone determinant:", cone_det)
    print("Scope: actual initial joint chart, compact reaction and physical current;")
    print("full-time configuration bracket and optical spatial closure remain explicit.")
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

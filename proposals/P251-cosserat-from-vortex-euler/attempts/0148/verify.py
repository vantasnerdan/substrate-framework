"""Exact material-angle variation and complete conditional field-chart algebra."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0148-material-frame")
    aa, bb = s.symbols("A B", positive=True)
    ex, ey, exy, eyx = s.symbols("ex ey exy eyx", real=True)
    inertia = s.diag(aa, bb)
    strain = s.Matrix([[ex, exy], [eyx, ey]])
    variation = strain*inertia+inertia*strain.T
    angle = variation[0, 1]/(aa-bb)
    omega = (eyx-exy)/2
    symmetric = (eyx+exy)/2
    checks.check("direct central quadrupole includes the symmetric strain row",
                 s.factor(angle-omega-(aa+bb)*symmetric/(aa-bb)) == 0)
    checks.check("rigid rotation has unit physical angle response",
                 s.simplify(angle.subs(exy, -eyx)-eyx) == 0)
    checks.check("constant core translation has zero central-angle response",
                 angle.subs({exy: 0, eyx: 0}) == 0)
    checks.check("symmetric shear is a counterexample to an unaveraged polar-only row",
                 angle.subs({exy: 1, eyx: 1}) != 0)
    t = s.symbols("t", real=True)
    rotation = s.Matrix([[s.cos(t), -s.sin(t)], [s.sin(t), s.cos(t)]])
    e0 = s.Matrix([[ex, symmetric], [symmetric, -ex]])
    local = rotation.T*e0*rotation
    checks.check("whole in-plane pose averaging removes the linear STF angle row",
                 s.integrate(s.expand_trig(local[0, 1]), (t, 0, 2*s.pi)) == 0)

    # Coordinate components in (r,theta,z), not orthonormal components:
    # both fields have zero radial component, so their exact Lie bracket
    # vanishes for arbitrary radial angular/axial profiles.
    r, theta, z = s.symbols("r theta z", positive=True)
    u = s.Matrix([0, s.Function("O")(r), s.Function("W")(r)])
    eta = s.Matrix([0, s.Function("g")(r), 0])
    bracket = eta.jacobian((r, theta, z))*u-u.jacobian((r, theta, z))*eta
    checks.check("arbitrary radial passive-label rotation has zero exact Lin velocity",
                 bracket == s.zeros(3, 1))

    rho, j, kap, mu, d = s.symbols("rho j kappa mu d", positive=True)
    b, h, frequency = s.symbols("b h frequency", real=True)
    c = j/(2*rho)
    observation = s.Matrix([[1, -c*h], [b*h, 1]])
    pullback = observation.inv()
    mass0 = s.diag(rho, j)
    potential0 = s.diag(mu*h*h, kap+d*h*h)
    mass = s.simplify(pullback.T*mass0*pullback)
    potential = s.simplify(pullback.T*potential0*pullback)
    den = 1+b*c*h*h
    checks.check("complete observed mass retains both gradient diagonal terms",
                 s.factor(mass[0, 0]-(rho+j*b*b*h*h)/den**2) == 0
                 and s.factor(mass[1, 1]-(j+rho*c*c*h*h)/den**2) == 0)
    checks.check("mixed inertia is derived rather than omitted",
                 s.factor(mass[0, 1]-(rho*c-j*b)*h/den**2) == 0)
    checks.check("potential transformation retains shear and optical gradients",
                 s.factor(potential[0, 1]
                          -(mu*c*h**3-b*h*(kap+d*h*h))/den**2) == 0)
    checks.check("one-half frame removes the leading mixed inertia",
                 mass[0, 1].subs(b, s.Rational(1, 2)) == 0)
    invariant = s.diff(potential[0, 1]-kap*mass[0, 1]/j, h).subs(h, 0)
    checks.check("optical first-gradient dynamical coupling is frame independent",
                 s.simplify(invariant+kap/2) == 0)
    checks.check("discarding mixed inertia at zero frame loses the nonzero coupling",
                 potential[0, 1].subs(b, 0) == mu*c*h**3
                 and s.simplify(invariant) != 0)
    optical = observation*s.Matrix([0, 1])
    checks.check("actual hybrid optical displacement ratio is independent of frame",
                 s.factor(optical[0]/optical[1]+j*h/(2*rho)) == 0)
    pencil = potential-frequency**2*mass
    checks.check("invertible field chart preserves the complete dispersion pencil",
                 s.factor(pencil.det()
                          -(potential0-frequency**2*mass0).det()/den**2) == 0)
    print("Scope: exact material variations and conditional full-action field map;")
    print("actual common acoustic/optical autonomous closure remains a separate construction")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()

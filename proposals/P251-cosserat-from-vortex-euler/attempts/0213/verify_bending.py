"""Exact curved material bending, rigid Ward rows and core-strain order."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0213-exact-curved-bending")
    r, z, p, radius = s.symbols("r z phi R", positive=True)
    q, capq, vertical, caph = (s.Function(n)(p) for n in ("q", "Q", "Z", "H"))
    rules = {s.diff(capq, p): q, s.diff(caph, p): vertical}
    xx = r-radius
    aa = s.diff(q, p)+capq
    potential = s.Matrix([
        0,
        (-radius*q*z+(r**3-radius**3)*vertical/(3*radius))/r,
        xx*capq+xx**2*aa/(2*radius)+z*xx*s.diff(vertical, p)/radius-z*caph,
    ])
    curl = s.Matrix([
        s.diff(potential[2], p)/r-s.diff(potential[1], z),
        s.diff(potential[0], z)-s.diff(potential[2], r),
        (s.diff(r*potential[1], r)-s.diff(potential[0], p))/r,
    ]).subs(rules)
    xi = s.Matrix([
        q-z*vertical/radius+(xx**2*(s.diff(q, p, 2)+q)/2
                            +z*xx*(s.diff(vertical, p, 2)+vertical))/(radius*r),
        -capq-xx*aa/radius-z*s.diff(vertical, p)/radius,
        r*vertical/radius,
    ])
    checks.check("independent cylindrical curl constructs the exact bending lift",
                 s.simplify(curl-xi) == s.zeros(3, 1))
    divergence = s.diff(r*xi[0], r)/r+s.diff(xi[1], p)/r+s.diff(xi[2], z)
    checks.check("full curved displacement preserves volume exactly",
                 s.simplify(divergence.subs(rules)) == 0)
    gradient = s.Matrix.hstack(
        xi.diff(r),
        (xi.diff(p)+s.Matrix([-xi[1], xi[0], 0]))/r,
        xi.diff(z),
    ).subs(rules)
    strain = s.simplify((gradient+gradient.T)/2)
    checks.check("cross-section center has zero symmetric strain at finite radius",
                 s.simplify(strain.subs({r: radius, z: 0})) == s.zeros(3))
    a, b = s.symbols("a b", real=True)
    translation = {q: a*s.cos(p), capq: a*s.sin(p), vertical: 0, caph: 0}
    rotation = {q: 0, capq: 0, vertical: b*s.cos(p), caph: b*s.sin(p)}
    checks.check("first radial harmonic is exact Cartesian translation",
                 s.simplify(xi.subs(translation).doit()
                            -s.Matrix([a*s.cos(p), -a*s.sin(p), 0])) == s.zeros(3, 1))
    checks.check("first vertical harmonic is exact Euclidean rotation germ",
                 s.simplify(xi.subs(rotation).doit()
                            -s.Matrix([-b*z*s.cos(p)/radius,
                                       b*z*s.sin(p)/radius,
                                       b*r*s.cos(p)/radius])) == s.zeros(3, 1))
    for name, replacement in (("translation", translation), ("rotation", rotation)):
        checks.check(f"{name} retains identically zero full strain",
                     s.simplify(strain.subs(replacement).doit()) == s.zeros(3))
    eps, localx = s.symbols("eps x", real=True)
    scaled = strain.subs(r, radius+localx).subs(radius, 1/eps)
    checks.check("all fixed-core strain entries have vanishing constant and first curvature jets",
                 all(s.simplify(s.limit(entry, eps, 0)) == 0
                     and s.simplify(s.limit(entry/eps, eps, 0)) == 0
                     for entry in scaled))
    strain2 = scaled.applyfunc(lambda entry: s.simplify(s.limit(entry/eps**2, eps, 0)))
    checks.check("second harmonic leaves a real nonrigid second-curvature strain",
                 strain2.subs({q: s.cos(2*p), capq: s.sin(2*p)/2,
                               vertical: 0, caph: 0}).doit() != s.zeros(3))
    naive = s.Matrix([q-z*vertical/radius, -capq+z*caph/radius, r*vertical/radius])
    naive_gradient = s.Matrix.hstack(
        naive.diff(r),
        (naive.diff(p)+s.Matrix([-naive[1], naive[0], 0]))/r,
        naive.diff(z),
    ).subs(rules)
    naive_strain = (naive_gradient+naive_gradient.T)/2
    checks.check("omitting the cross-section rotation restores the unwanted first-order shear",
                 s.simplify(naive_strain[0, 1].subs({r: radius, z: 0})
                            -(s.diff(q, p)+capq)/(2*radius)) == 0)
    coords = s.symbols("x1 x2 x3", real=True)
    vel = s.Matrix([coords[1]*coords[2], coords[0]*coords[2], coords[0]*coords[1]])
    disp = s.Matrix([coords[1]**2, coords[2]**2, coords[0]**2])
    vort = s.Matrix([s.diff(vel[2], coords[1])-s.diff(vel[1], coords[2]),
                     s.diff(vel[0], coords[2])-s.diff(vel[2], coords[0]),
                     s.diff(vel[1], coords[0])-s.diff(vel[0], coords[1])])
    bracket = disp.jacobian(coords)*vel-vel.jacobian(coords)*disp
    grad_scalar = s.Matrix([s.diff(disp.dot(vel), c) for c in coords])
    checks.check("Cartesian Kelvin offset is minus twice strain modulo its exact pressure gradient",
                 s.expand(disp.cross(vort)-bracket-grad_scalar
                          +(disp.jacobian(coords)+disp.jacobian(coords).T)*vel)
                 == s.zeros(3, 1))
    print("Derived R^-2 strain coefficient:")
    print(strain2)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

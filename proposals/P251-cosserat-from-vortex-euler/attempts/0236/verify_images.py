"""Exact full-periodic local Green structure and circular image forcing."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0236-periodic-image-forcing")
    x, y, z, radius, theta = s.symbols("x y z R theta", real=True)
    coords = (x, y, z)
    rr = x*x+y*y+z*z
    harmonic = x**4+y**4+z**4-s.Rational(3, 5)*rr**2
    checks.check("the full periodic zero-mode subtraction fixes its quadratic regular term",
                 sum(s.diff(rr/6, q, 2) for q in coords) == 1)
    checks.check("the first cubic lattice anisotropy is the exact harmonic quartic",
                 s.simplify(sum(s.diff(harmonic, q, 2) for q in coords)) == 0)
    checks.check("its retained planar Laplacian is fixed by the full three-dimensional kernel",
                 s.simplify((s.diff(harmonic, x, 2)+s.diff(harmonic, y, 2)).subs(z, 0)
                            -s.Rational(12, 5)*(x*x+y*y)) == 0)
    source = s.Matrix([radius*s.cos(theta), radius*s.sin(theta), 0])
    tangent = source.diff(theta)
    shifted = s.Matrix([s.diff(harmonic, q).subs(
        {x: x-source[0], y: y-source[1], z: z}, simultaneous=True) for q in coords])
    image = shifted.cross(tangent)
    integrated = s.Matrix([s.simplify(s.integrate(s.expand_trig(v.subs(z, 0)),
                                                 (theta, 0, 2*s.pi))) for v in image])
    expected = s.Matrix([0, 0, -12*s.pi*(radius**2*(x*x+y*y)+radius**4/2)/5])
    checks.check("actual circular vorticity integration retains the complete quartic image velocity",
                 s.simplify(integrated-expected) == s.zeros(3, 1))
    on_core = integrated.subs({x: radius*s.cos(theta), y: radius*s.sin(theta)})
    checks.check("the first quartic image is a uniform axial drift on the circular core",
                 s.simplify(on_core-s.Matrix([0, 0, -18*s.pi*radius**4/5])) == s.zeros(3, 1))
    monopole = s.Matrix([s.integrate(v, (theta, 0, 2*s.pi)) for v in tangent])
    moment = s.Matrix([s.integrate(v, (theta, 0, 2*s.pi))
                       for v in source.cross(tangent)])
    checks.check("actual closed-ring vorticity has zero monopole and nonzero impulse",
                 monopole == s.zeros(3, 1) and moment == s.Matrix([0, 0, 2*s.pi*radius**2]))
    quadratic = (s.Matrix(coords)-source).cross(tangent)/3
    q_image = s.Matrix([s.simplify(s.integrate(v, (theta, 0, 2*s.pi))) for v in quadratic])
    checks.check("the leading periodic image drift has the derived negative impulse sign",
                 q_image == -moment/3)
    checks.check("treating the harmonic quartic as a zero image misses a real velocity",
                 integrated != s.zeros(3, 1))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

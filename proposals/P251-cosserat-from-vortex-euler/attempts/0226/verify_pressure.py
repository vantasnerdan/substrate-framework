"""Exact compact pressure-moment and first toroidal correction identities."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0226-compact-pressure-recursion")
    r = s.symbols("r", positive=True)
    m, j = s.symbols("m j", integer=True, nonnegative=True)
    test = r**(m+2*j)
    radial_lap = s.diff(test, r, 2)+s.diff(test, r)/r-m**2*test/r**2
    checks.check("radial Green integration has its exact polynomial-moment coefficient",
                 s.simplify(radial_lap-4*j*(m+j)*r**(m+2*j-2)) == 0)
    f = s.Function("f")(r)
    dm = s.diff(f, r, 2)+s.diff(f, r)/r-m**2*f/r**2
    flux = r**(m+1)*s.diff(f, r)-m*r**m*f
    checks.check("the first exterior harmonic moment is exactly a compact boundary flux",
                 s.simplify(r**(m+1)*dm-s.diff(flux, r)) == 0)
    aa = s.Function("A")(r)
    g0 = -(s.diff(aa, r)+aa/r)/2
    g2 = -(s.diff(aa, r)-aa/r)/2
    checks.check("first toroidal correction has identically zero exterior monopole",
                 s.simplify(r*g0+s.diff(r*aa, r)/2) == 0)
    checks.check("its quadrupole tail is the existing compact m1 return moment",
                 s.simplify(r**3*g2-2*r**2*aa+s.diff(r**3*aa, r)/2) == 0)
    x, z, h, nn = s.symbols("x z h n", real=True)
    pp = s.Function("p")(x, z)
    ring_lap = s.diff(pp, x, 2)+s.diff(pp, z, 2)
    ring_lap += h*s.diff(pp, x)/(1+h*x)-nn**2*h**2*pp/(1+h*x)**2
    checks.check("actual toroidal Laplacian supplies the retained first curvature operator",
                 s.diff(ring_lap, h).subs(h, 0) == s.diff(pp, x))
    checks.check("second curvature operator retains both metric and global harmonic",
                 s.diff(ring_lap, h, 2).subs(h, 0)/2
                 == -x*s.diff(pp, x)-nn**2*pp)
    om = s.symbols("Omega", real=True)
    psi = s.Function("psi")(x, z)
    vel = s.Matrix([-s.diff(psi, z), s.diff(psi, x)])
    jj = s.Matrix([[0, -1], [1, 0]])
    checks.check("uniform-core vector rotation retains its exact compact pressure",
                 2*om*jj*vel+s.Matrix([s.diff(2*om*psi, x),
                                      s.diff(2*om*psi, z)]) == s.zeros(2, 1))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

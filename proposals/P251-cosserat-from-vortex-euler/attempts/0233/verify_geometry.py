"""Actual centered ring, off-resonant pressure and physical centroid anchors."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0233-geometric-physical-current")
    r, lam, seed, om, w, t = s.symbols("r lambda seed Omega W t", positive=True)
    phi0 = seed*s.besselj(0, lam*r)
    phi1 = -seed*r*s.besselj(2, lam*r)/2
    lap1 = s.diff(phi1, r, 2)+s.diff(phi1, r)/r-phi1/r**2
    residual = s.expand_func(lap1+lam**2*phi1-s.diff(phi0, r))
    checks.check("centered inner coefficient solves the actual affine Helmholtz source",
                 s.simplify(residual) == 0)
    checks.check("the actual centered coefficient does not translate the core again",
                 s.limit(phi1/r, r, 0) == 0)
    x, y = s.symbols("x y", real=True)
    local_phi1 = -om*x*(x*x+y*y)/8
    jj = s.Matrix([[0, -1], [1, 0]])
    u0 = om*jj*s.Matrix([x, y])
    u1 = -jj*s.Matrix([s.diff(local_phi1, x), s.diff(local_phi1, y)])-x*u0
    checks.check("first geometric velocity includes the actual volume-divergence term",
                 s.simplify(s.diff(u1[0], x)+s.diff(u1[1], y)+u0[0]) == 0)
    checks.check("large axial Euler connection cancels its actual field-gradient partner",
                 s.simplify(-s.diff(-w*x, x)-w) == 0)

    a = s.Function("A")(r)
    plus = s.diff(a, r)+a/r
    minus = s.diff(a, r)-a/r
    z0 = 2*s.I*om*plus*s.exp(s.I*om*t)
    z2 = -2*s.I*om*minus*(2*s.exp(2*s.I*om*t)-s.exp(s.I*om*t))
    checks.check("actual generated monopole solves the complete forced Euler row",
                 s.simplify(s.diff(z0, t)+2*om**2*plus*s.exp(s.I*om*t)) == 0)
    checks.check("the full pressure retains a nonzero off-resonant double-frequency row",
                 s.simplify(s.diff(z2, t)-2*s.I*om*z2
                            -2*om**2*minus*s.exp(s.I*om*t)) == 0
                 and s.simplify(-4*s.I*om*minus) != 0)

    hh, ff = s.Function("h")(r), s.Function("F")(r)
    df = s.diff(ff, r)+ff/r
    cc = (r**2*s.diff(hh, r)+3*r*hh)/16-w*plus/om+w*lam*df/om
    axial = cc*(s.exp(s.I*om*t)-1)-s.I*w*lam*t*df*s.exp(s.I*om*t)
    forcing = (s.I*om*(r**2*s.diff(hh, r)+3*r*hh)/16
               -s.I*w*plus+w*lam*om*t*df)*s.exp(s.I*om*t)
    checks.check("actual axial Lin row retains its constant and secular polar pieces",
                 s.simplify(s.diff(axial, t)-forcing) == 0
                 and axial.subs(t, 0) == 0)
    flux = r**3*hh/16-w*r*a/om+w*lam*r*ff/om
    checks.check("centroid repair has the exact zero moment needed for a compact m0 inverse",
                 s.simplify(r*cc-s.diff(flux, r)) == 0)
    repaired = axial+cc
    checks.check("the actual initial toroidal return removes its constant centroid offset",
                 s.simplify(repaired-(cc-s.I*w*lam*t*df)*s.exp(s.I*om*t)) == 0)

    radius, vx, vy, vz = s.symbols("R xi_x xi_y xi_z", real=True)
    qrow = -y*vx-(radius+x)*vy-s.I*y*vz
    grow = y*vz+s.I*((radius+x)*vy-y*vx)
    checks.check("full Euclidean displacement moment retains the toroidal tag term",
                 s.simplify(grow+s.I*qrow-2*y*(vz-s.I*vx)) == 0)
    ang = s.symbols("theta", real=True)
    checks.check("leading geometric sidebands are absent from angle only by actual selection",
                 all(s.integrate(s.exp(s.I*j*ang)*s.sin(ang), (ang, 0, 2*s.pi)) == 0
                     for j in (0, -2)))
    pot = s.Function("potential")(x, y)
    source = s.diff(pot, x, 2)+s.diff(pot, y, 2)
    k, kap = s.symbols("k kappa", real=True)
    control = s.Matrix([-s.I*k*kap*s.diff(pot, x),
                        -s.I*k*kap*s.diff(pot, y), kap*source])
    checks.check("centroid return includes its physically necessary second-order divergence completion",
                 s.simplify(s.diff(control[0], x)+s.diff(control[1], y)
                            +s.I*k*control[2]) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

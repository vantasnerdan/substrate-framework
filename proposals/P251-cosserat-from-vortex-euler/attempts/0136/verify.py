"""Exact shaped force-free column and toroidal Euler construction."""

import sympy as s

from substrate_framework.euler_forcefree import planar_bernoulli_lift
from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0136-shaped-forcefree-core")
    r, a, omega, axial, rho = s.symbols("r a Omega U0 rho", positive=True)
    v, w = s.Function("V")(r), s.Function("W")(r)
    zeta = s.diff(v, r)+v/r
    factor = zeta/w
    axial_rule = {s.diff(w, r): -zeta*v/w}
    checks.check("derived axial velocity supplies the azimuthal force-free equation",
                 s.simplify((-s.diff(w, r)-factor*v).subs(axial_rule)) == 0)
    checks.check("the axial curl supplies the same nonconstant factor",
                 s.simplify(s.diff(r*v, r)/r-factor*w) == 0)
    p = -rho*(v**2+w**2)/2
    checks.check("full force-free pressure equals the centrifugal Euler pressure",
                 s.simplify((s.diff(p, r)-rho*v**2/r).subs(axial_rule)) == 0)
    checks.check("integrated axial deficit follows from the swirl rather than a fitted constant",
                 s.simplify(s.diff(v**2, r)+2*v**2/r-2*zeta*v) == 0)

    profile = omega*r/(1+r**2/a**2)
    axial_squared = axial**2-omega**2*a**2+omega**2*a**2/(1+r**2/a**2)**2
    core_zeta = s.factor(s.diff(r*profile, r)/r)
    checks.check("explicit smooth core vorticity is derived from circulation",
                 s.simplify(core_zeta-2*omega/(1+r**2/a**2)**2) == 0)
    checks.check("explicit axial square solves the full force-free ODE",
                 s.simplify(s.diff(axial_squared, r)+2*core_zeta*profile) == 0)
    checks.check("axis velocity and vorticity have regular nonzero limits",
                 s.limit(profile/r, r, 0) == omega
                 and s.limit(core_zeta, r, 0) == 2*omega
                 and s.limit(axial_squared, r, 0) == axial**2)
    checks.check("the explicit core has a real positive axial floor under its stated hypothesis",
                 s.simplify(axial_squared-(axial**2-omega**2*a**2))
                 == omega**2*a**6/(a**2+r**2)**2)
    checks.check("actual Galilean-frame axial correction is inverse-speed small",
                 s.simplify(s.limit((s.sqrt(axial_squared)-axial)*axial,
                                     axial, s.oo)
                            +(omega**2*a**2-omega**2*a**2/(1+r**2/a**2)**2)/2) == 0)
    checks.check("an imposed constant axial velocity would not solve the shaped force-free equation",
                 s.simplify(core_zeta*profile/axial) != 0)

    # Toroidal cylindrical coordinates: x is r_toroidal - R, theta suppressed.
    x, z, major = s.symbols("x z R", real=True)
    radius = major+x
    phi = s.Function("phi")(x, z)
    c = s.Function("C")(phi)
    cp = s.diff(c, phi)
    ur = -major*s.diff(phi, z)/radius
    ut = major*c/radius
    uz = major*s.diff(phi, x)/radius
    field = (ur, ut, uz)
    div = s.diff(radius*ur, x)/radius+s.diff(uz, z)
    curl = (-s.diff(ut, z), s.diff(ur, z)-s.diff(uz, x),
            s.diff(radius*ut, x)/radius)
    pde = {s.diff(phi, z, 2): -s.diff(phi, x, 2)+s.diff(phi, x)/radius-c*cp}
    checks.check("actual toroidal velocity is fully divergence free", s.simplify(div) == 0)
    for i in range(3):
        checks.check(f"toroidal curl component {i} has the derived common factor",
                     s.simplify((curl[i]-cp*field[i]).subs(pde)) == 0)
    pressure = -rho*sum(component**2 for component in field)/2
    euler = (
        ur*s.diff(ur, x)+uz*s.diff(ur, z)-ut**2/radius+s.diff(pressure, x)/rho,
        ur*s.diff(ut, x)+uz*s.diff(ut, z)+ur*ut/radius,
        ur*s.diff(uz, x)+uz*s.diff(uz, z)+s.diff(pressure, z)/rho,
    )
    checks.check("the complete toroidal Euler pressure and all convective terms balance",
                 all(s.simplify(value.subs(pde)) == 0 for value in euler))
    checks.check("the constructed torus surface is an exact material surface",
                 s.simplify(ur*s.diff(phi, x)+uz*s.diff(phi, z)) == 0)

    h, q = s.Function("h")(r), s.Function("Q")(r)
    lhs = r*s.diff(h, r)**2+(1/r-r*q)*h**2
    square = r*v**2*s.diff(h/v, r)**2
    boundary = s.diff(r*s.diff(v, r)*h**2/v, r)
    radial_rule = {s.diff(v, r, 2): -s.diff(v, r)/r-(q-1/r**2)*v}
    checks.check("Dirichlet angular-sector ground-state identity retains its boundary derivative",
                 s.simplify((lhs-square-boundary).subs(radial_rule)) == 0)
    aa, bb = s.symbols("A B", real=True)
    outer = aa+bb*s.log(r)
    checks.check("the only radial zero-mode candidate is harmonic outside the core",
                 s.simplify(s.diff(outer, r, 2)+s.diff(outer, r)/r) == 0)
    circulation, far_axial = s.symbols("Gamma W_infinity", positive=True)
    rotation = circulation/(r**2*far_axial)
    checks.check("actual outer-section rotation has nonzero flux-action twist",
                 s.simplify(s.diff(rotation, r)/(r*far_axial)
                            +2*circulation/(r**4*far_axial**2)) == 0)
    checks.check("the shaped squared-velocity deficit is independent of the free throughflow speed",
                 s.diff(axial_squared-axial**2, axial) == 0)
    # Derive the global periodic example through the reusable API, then
    # independently differentiate its complete Cartesian field.
    y = s.Symbol("y", real=True)
    planar = [s.sin(x)*s.cos(y), -s.cos(x)*s.sin(y)]
    planar_p = rho*(s.cos(2*x)+s.cos(2*y))/4
    lift = planar_bernoulli_lift(planar, planar_p, rho, (x, y), 1)
    uu = lift.velocity
    actual_curl = s.Matrix([s.diff(uu[2], y), -s.diff(uu[2], x),
                           s.diff(uu[1], x)-s.diff(uu[0], y)])
    checks.check("global periodic Bernoulli lift satisfies the full force-free equation",
                 s.simplify(actual_curl-lift.curl_factor*uu) == s.zeros(3, 1))
    checks.check("the full pressure has exactly constant Bernoulli after the lift",
                 s.simplify(planar_p/rho+uu.dot(uu)/2-1) == 0)
    checks.check("the planar lift really needs its nonconstant axial flow",
                 s.simplify(s.diff(uu[2], x)) != 0)
    xx, xy, yx, yy, wx, wy, zz = s.symbols("xi_x xi_y eta_x eta_y w_x w_y zeta")
    full_omega = s.Matrix([wy, -wx, zz])
    planar_cross = s.Matrix([xx, xy, 0]).cross(s.Matrix([yx, yy, 0]))
    checks.check("full KKS contraction on the invariant planar orbit keeps precisely planar vorticity",
                 s.expand(full_omega.dot(planar_cross)-zz*(xx*yy-xy*yx)) == 0)
    print("Scope: exact shaped core/local torus/global planar lift; finite-time transfer is analytic;")
    print("global EPS embedding and a physical optical eigenmode are not asserted by these checks")
    raise SystemExit(checks.finish())


if __name__ == "__main__":
    main()

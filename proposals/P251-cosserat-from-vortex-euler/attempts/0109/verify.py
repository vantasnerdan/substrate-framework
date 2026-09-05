"""Exact pressure, Lin, and full material moment construction for 0109."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0109-actual-optical-parcel")
    x, y, z, t = s.symbols("x y z t", real=True)
    k, om, a, radius, mass, zc, delta, lam = s.symbols(
        "k Omega a R M zc delta lambda", positive=True)
    coords = (x, y, z)
    r = s.Matrix(coords)
    ez = s.Matrix([0, 0, 1])
    b = s.Matrix([s.cos(om*t), -s.sin(om*t), 0])
    jb = ez.cross(b)
    f = s.Function("F")((x*x+y*y)/a**2)
    primitive = s.Function("H")(k*z)
    h = s.diff(primitive, z)/k

    def curl(v):
        return s.Matrix([s.diff(v[2], y)-s.diff(v[1], z),
                         s.diff(v[0], z)-s.diff(v[2], x),
                         s.diff(v[1], x)-s.diff(v[0], y)])

    def grad(q):
        return s.Matrix([s.diff(q, c) for c in coords])

    def zero(v):
        return all(s.simplify(entry) == 0 for entry in v)

    def ball_average(poly):
        terms = s.Poly(s.expand(poly), x, y, z).terms()
        result = s.S.Zero
        for powers, coefficient in terms:
            if powers == (0, 0, 0):
                result += coefficient
            elif powers in ((2, 0, 0), (0, 2, 0), (0, 0, 2)):
                result += coefficient*radius**2/5
            elif any(power % 2 for power in powers):
                continue
            else:
                raise ValueError(f"Unprovided exact ball moment {powers}")
        return s.simplify(result)

    packet = curl(jb*f*h/k)
    u = om*ez.cross(r)
    residual = packet.diff(t)+packet.jacobian(coords)*u+u.jacobian(coords)*packet
    pressure = 2*om*jb.dot(grad(f))*primitive/k**2
    improved = 2*om*grad(jb.dot(grad(f)))*primitive/k**2
    ledger.check("the full pressure gradient cancels the leading axial return residual",
                 zero(residual+grad(pressure)-improved))
    ledger.check("improved residual remains nonzero and purely transverse",
                 not zero(improved) and s.simplify(improved[2]) == 0)

    w = -b*f*h/k+ez*b.dot(grad(f))*primitive/k**2
    xi = curl(w)/(2*om)
    lin_residual = xi.diff(t)+xi.jacobian(coords)*u-u.jacobian(coords)*xi-packet
    lin_target = ez.cross(s.hessian(f, coords)*jb)*primitive/k**2
    ledger.check("the independent full Lin residual has the exact second-order return",
                 zero(lin_residual-lin_target))
    ledger.check("the new envelope retains exact compact Kelvin preparation",
                 zero(curl(xi.cross(2*om*ez))-curl(packet))
                 and s.simplify(s.trace(xi.jacobian(coords))) == 0)

    zz = s.symbols("zeta", real=True)
    inner_h = s.diff(zz-zz**3/6, zz)
    ledger.check("the new polynomial inner envelope gives an exact affine physical tilt",
                 s.diff(inner_h, zz) == -zz and s.diff(inner_h, zz, 2) == -1)

    q = k*jb/(2*om)
    phi = ez.cross(q)
    inner_xi = q*(z+zc)
    inner_v = b*k*(z+zc)
    ledger.check("actual inner Lin reconstruction is not an arbitrary angle-rate ansatz",
                 zero(inner_xi.diff(t)+inner_xi.jacobian(coords)*u
                      -u.jacobian(coords)*inner_xi-inner_v))
    ledger.check("physical rotation vector is determined by actual core vorticity",
                 zero(curl(inner_v)/(2*om)-q) and zero(phi.cross(ez)-q))
    ledger.check("physical tilt and rotation vector follow the same genuine optical frequency",
                 zero(q.diff(t)+om*ez.cross(q))
                 and zero(phi.diff(t, 2)+om**2*phi))

    centroid = inner_xi.applyfunc(ball_average)
    material_velocity = inner_xi.diff(t)+inner_xi.jacobian(coords)*u
    momentum = mass*material_velocity.applyfunc(ball_average)
    spin_density = (inner_xi-centroid).cross(u)+r.cross(
        material_velocity-centroid.diff(t))
    spin = mass*spin_density.applyfunc(ball_average)
    inertia = mass*radius**2/5
    ledger.check("actual material centroid is nonzero and geometrically tied to the tilt",
                 zero(centroid-zc*q) and s.simplify(centroid.dot(centroid)) > 0)
    ledger.check("actual parcel momentum equals mass times the transported centroid rate",
                 zero(momentum-mass*centroid.diff(t))
                 and zero(momentum+mass*zc*ez.cross(phi.diff(t))))
    ledger.check("full moving-boundary material spin gives a positive derived inertia",
                 zero(spin-inertia*phi.diff(t))
                 and inertia.is_positive
                 and s.simplify(spin.dot(spin)) > 0)

    bulk = mass*r.cross(inner_v).applyfunc(ball_average)
    ledger.check("fixed-domain velocity spin overcounts the actual result by exactly two",
                 zero(bulk-2*spin) and zero((spin-bulk)+inertia*phi.diff(t)))
    ledger.check("selected fluid inertia differs from the ball's transverse rigid inertia",
                 s.simplify(2*mass*radius**2/5-inertia) == inertia)

    hessian_p = s.diag(om**2, om**2, 0)
    gradient_p = hessian_p*r
    force = -mass*(hessian_p*inner_xi).applyfunc(ball_average)
    torque = -mass*((inner_xi-centroid).cross(gradient_p)
                   +r.cross(hessian_p*inner_xi)).applyfunc(ball_average)
    ledger.check("ambient background-pressure force supplies the exact optical translation",
                 zero(force+mass*om**2*centroid) and zero(force-momentum.diff(t)))
    ledger.check("ambient moving-boundary pressure torque supplies the exact spin evolution",
                 zero(torque+inertia*om**2*phi) and zero(torque-spin.diff(t)))

    ell = s.Rational(1, 8)
    c = s.Rational(1, 4)
    ledger.check("the explicit noncentral material ball stays strictly inside the inner slab",
                 c-ell > 0 and c+ell < s.Rational(1, 2))

    i_f, i_h1, i_f2, i_h0 = s.symbols("IF Ih1 IF2 IH", positive=True)
    n_horizontal_sq = a**2*i_f*i_h1/k
    return_sq = 4*om**2*i_f2*i_h0/(a**2*k**5)
    ledger.check("the pressure-corrected global residual gains the required second aspect power",
                 s.simplify(return_sq/n_horizontal_sq
                            -4*om**2*i_f2*i_h0/(i_f*i_h1*(k*a)**4)) == 0)
    chosen = {a: delta/lam, k: lam*delta**s.Rational(-4, 3)}
    error_bound = 1/(k*a)**2+delta*k*a+delta
    ledger.check("finite carrier selection makes the complete velocity error order delta two-thirds",
                 s.simplify(error_bound.subs(chosen)-2*delta**s.Rational(2, 3)-delta) == 0)
    ledger.check("the actual parcel-moment sensitivity still leaves a vanishing error",
                 s.simplify((k*a*error_bound).subs(chosen)
                            -2*delta**s.Rational(1, 3)-delta**s.Rational(2, 3)) == 0)
    ledger.check("actual near-axis parcel transport occurs on the smaller derived scale",
                 s.simplify((lam/k).subs(chosen)-delta**s.Rational(4, 3)) == 0)

    lx, ly = s.Function("Lx")(x, y), s.Function("Ly")(x, y)
    div_l = s.diff(lx, x)+s.diff(ly, y)
    ledger.check("the compact net-spin proof uses an exact first-moment divergence identity",
                 s.expand(s.diff(x*lx, x)+s.diff(x*ly, y)-lx-x*div_l) == 0
                 and s.expand(s.diff(y*lx, x)+s.diff(y*ly, y)-ly-y*div_l) == 0)

    print("EXACT parcel inertia:", inertia)
    print("EXACT centroid:", centroid.T)
    print("EXACT material spin:", spin.T)
    print("EXACT moving-boundary contribution: minus one half of the bulk velocity spin")
    print("EXACT residual scales: velocity delta^(2/3), physical parcel moments delta^(1/3)")
    print("SCOPE: actual transported core parcel; distant EPS knot and independent continuum remain distinct")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())

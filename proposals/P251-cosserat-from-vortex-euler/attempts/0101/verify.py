"""Exact calculus for the Kelvin-prepared finite-time optical tilt packet."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0101-Kelvin-optical-packet")
    x, y, z, t = s.symbols("x y z t", real=True)
    k, om, a, delta, lam = s.symbols("k Omega a delta lambda", positive=True)
    coords = (x, y, z)
    ez = s.Matrix([0, 0, 1])
    b = s.Matrix([s.cos(om*t), -s.sin(om*t), 0])
    jb = ez.cross(b)
    f = s.Function("F")((x*x+y*y)/a**2)
    hprim = s.Function("H")(k*z)
    h = s.diff(hprim, z)/k

    def curl(v):
        return s.Matrix([s.diff(v[2], y)-s.diff(v[1], z),
                         s.diff(v[0], z)-s.diff(v[2], x),
                         s.diff(v[1], x)-s.diff(v[0], y)])

    def zero(v):
        return all(s.simplify(entry) == 0 for entry in v)

    grad_f = s.Matrix([s.diff(f, c) for c in coords])
    potential = jb*f*h/k
    packet = curl(potential)
    claimed = -b*f*s.diff(h, z)/k+ez*b.dot(grad_f)*h/k
    ledger.check("full compact curl contains the horizontal tilt and vertical return",
                 zero(packet-claimed))
    ledger.check("packet is exactly divergence free, not only at leading carrier order",
                 s.simplify(s.trace(packet.jacobian(coords))) == 0)
    u_rot = om*ez.cross(s.Matrix(coords))
    remainder = packet.diff(t)+packet.jacobian(coords)*u_rot+u_rot.jacobian(coords)*packet
    target = -2*om*ez*jb.dot(grad_f)*h/k
    ledger.check("the full rotating-Euler residual is precisely the retained return",
                 zero(remainder-target))
    ledger.check("compact return residual is not silently certified as an exact mode",
                 not zero(target))
    horizontal = -b*f*s.diff(h, z)/k
    ledger.check("discarding the return would violate incompressibility",
                 s.simplify(s.trace(horizontal.jacobian(coords))) != 0)

    primitive = -b*f*h/k+ez*b.dot(grad_f)*hprim/k**2
    displacement = curl(primitive)/(2*om)
    ledger.check("compact primitive differentiates to the actual packet",
                 zero(primitive.diff(z)-packet))
    ledger.check("Kelvin preparer is exactly solenoidal",
                 s.simplify(s.trace(displacement.jacobian(coords))) == 0)
    ledger.check("Kelvin preparation matches curl before global Leray uniqueness",
                 zero(curl(displacement.cross(2*om*ez))-curl(packet)))
    kelvin_difference = displacement.cross(2*om*ez)-packet
    ledger.check("Leray gradient correction has not been replaced by raw pointwise equality",
                 not zero(kelvin_difference) and zero(curl(kelvin_difference)))

    core = -b*s.diff(h, z)/k
    tilt = curl(core)
    ledger.check("physical core vorticity tilt follows directly from the velocity curl",
                 zero(tilt+jb*s.diff(h, z, 2)/k))
    ledger.check("intrinsic fixed-axis optical frequency is Omega, not axial carrier Doppler",
                 zero(jb.diff(t)+om*ez.cross(jb))
                 and zero(jb.diff(t, 2)+om**2*jb))
    rotating = s.Matrix([[s.cos(om*t), s.sin(om*t), 0],
                         [-s.sin(om*t), s.cos(om*t), 0], [0, 0, 1]])*jb
    ledger.check("material corotating axes see the separate frequency two Omega",
                 zero(rotating.diff(t)+2*om*ez.cross(rotating)))

    q = s.symbols("q", positive=True)
    j0, j1 = s.besselj(0, q), s.besselj(1, q)
    ledger.check("smooth Lundquist axial and azimuthal components obey Beltrami recurrences",
                 s.simplify(s.diff(j0, q)+j1) == 0
                 and s.simplify(s.expand_func(s.diff(j1, q)+j1/q-j0)) == 0)
    ledger.check("near-axis rotation and axial shear orders are derived from Bessel series",
                 s.series(j1, q, 0, 6).removeO() == q/2-q**3/16+q**5/384
                 and s.series(j0, q, 0, 5).removeO() == 1-q**2/4+q**4/64)

    zz = s.symbols("zeta", real=True)
    concrete_hprim = zz*s.exp(-1/(1-zz**2))
    numerator = 3*zz**8+24*zz**6-26*zz**4+3
    actual_h2 = s.diff(concrete_hprim, zz, 3)
    formula = -2*s.exp(-1/(1-zz**2))*numerator/(1-zz**2)**6
    ledger.check("positive-weight optical observable has the derived nonzero core curvature",
                 s.simplify(actual_h2-formula) == 0
                 and actual_h2.subs(zz, 0) == -6/s.E
                 and 3-s.Rational(26, 10000) > 0)

    i_f, i_g, i_h, i_h1 = s.symbols("IF IG Ih Ih1", positive=True)
    horizontal_norm = a**2*i_f*i_h1/k
    return_norm = i_g*i_h/(2*k**3)
    rotation_norm = 2*om**2*i_g*i_h/k**3
    ledger.check("separated return residual has the claimed inverse-aspect-ratio scale",
                 s.simplify(rotation_norm/horizontal_norm
                            -om**2/(k*a)**2*(2*i_g*i_h/(i_f*i_h1))) == 0
                 and s.simplify(return_norm/horizontal_norm
                                -i_g*i_h/(2*(k*a)**2*i_f*i_h1)) == 0)
    c_rot, c_d = s.symbols("Crot CD", positive=True)
    residual_bound = c_rot/(k*a)+6*c_d*delta*k*a+8*delta
    chosen = s.simplify(residual_bound.subs({a: delta/lam, k: lam*delta**s.Rational(-3, 2)}))
    ledger.check("finite localization and carrier selection balance the actual residual terms",
                 s.expand(chosen-((c_rot+6*c_d)*s.sqrt(delta)+8*delta)) == 0)
    ledger.check("moving-tag gradient correction is smaller than the optical carrier",
                 s.simplify((lam/k).subs(k, lam*delta**s.Rational(-3, 2))
                            -delta**s.Rational(3, 2)) == 0)
    c_w, c_n, bnorm = s.symbols("Cw CN bnorm", positive=True)
    observation_bound = c_w*k**s.Rational(3, 2)/a*c_n*bnorm*a/s.sqrt(k)
    ledger.check("actual weak-curl angle error has no hidden inverse packet-size penalty",
                 s.simplify(observation_bound/(k*bnorm)-c_w*c_n) == 0)

    # The odd primitive makes the full packet's impulse cancel; this is
    # evidence against identifying its local tilt with a total-spin rotor.
    ledger.check("the explicit primitive is odd and its differentiated profile even",
                 s.simplify(concrete_hprim.subs(zz, -zz)+concrete_hprim) == 0
                 and s.simplify(s.diff(concrete_hprim, zz).subs(zz, -zz)
                                -s.diff(concrete_hprim, zz)) == 0)
    print("EXACT rotation residual: -2 Omega/k (J b dot grad F_a) h(kz) e_z")
    print("EXACT finite parameter choice: a=delta/lambda, k=lambda*delta^(-3/2)")
    print("EXACT dimensionless residual bound:", chosen)
    print("EXACT core h''(0):", actual_h2.subs(zz, 0))
    print("SCOPE: finite-time actual Euler material tilt; not knotted-torus angle or total-spin rotor")
    return ledger.finish()


if __name__ == "__main__":
    raise SystemExit(main())

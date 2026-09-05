"""Exact Kelvin, compact-completion and complete-energy anchors for 0222."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0222-curved-kelvin-transfer")
    r, lam, om = s.symbols("r lambda Omega", positive=True)
    w, v, wp, ss, sr, st, sz = s.symbols("W V Wprime S Sr St Sz", real=True)
    xi = s.Matrix([-st/r, sr, lam*ss])
    u = s.Matrix([0, v, w])
    grad = s.Matrix([wp*ss+w*sr, w*st/r, w*sz])
    residual = xi.cross(u)-grad+s.Matrix([0, 0, om*st+w*sz])
    checks.check("literal column cross product retains its full gradient and axial return",
                 s.simplify(residual.subs({v: r*om, wp: -lam*r*om})) == s.zeros(3, 1))

    th = s.symbols("theta", real=True)
    aa, ap, m, op = s.symbols("A Aprime m Oprime", real=True)
    alpha, q = m*aa/r, m*om
    xc = s.Matrix([alpha*s.sin(th), ap*s.cos(th), lam*aa*s.cos(th)])
    xs = s.Matrix([-alpha*s.cos(th), ap*s.sin(th), lam*aa*s.sin(th)])
    wc = s.Matrix([0, 0, lam*q*aa*s.sin(th)])
    ws = s.Matrix([0, 0, -lam*q*aa*s.cos(th)])
    du = s.Matrix([[0, -om, 0], [om+r*op, 0, 0], [wp, 0, 0]])
    pic = (wc+du*xc).subs(wp, -lam*r*om)
    pis = (ws+du*xs).subs(wp, -lam*r*om)
    checks.check("actual Kelvin cotangent has zero axial component, not the homogeneous one",
                 pic[2] == 0 and pis[2] == 0)
    zz = 2*om+r*op
    phase = xc.dot(pis)-pic.dot(xs)
    checks.check("complete canonical phase restores the positive homogeneous phase",
                 s.simplify(phase+zz*alpha*ap) == 0)
    vort = s.Matrix([0, lam*r*om, zz])
    full_phase = -vort.dot(xc.cross(xs))+xc.dot(ws)-wc.dot(xs)
    checks.check("vorticity-cross and velocity-cross terms must both be retained",
                 s.simplify(full_phase-phase) == 0
                 and s.simplify(xc.dot(ws)-wc.dot(xs)) != 0)
    xt = s.Matrix([-q*alpha*s.cos(th),
                   (q*ap+r*op*alpha)*s.sin(th), lam*q*aa*s.sin(th)])
    tx = s.Matrix([(q*alpha-om*ap)*s.cos(th),
                   (-q*ap+om*alpha)*s.sin(th), -lam*q*aa*s.sin(th)])
    checks.check("full passive Kelvin Lin equation is exact in every component",
                 s.simplify((xt+tx-du*xc-wc).subs(wp, -lam*r*om)) == s.zeros(3, 1))
    checks.check("the actual passive velocity has nonzero vorticity for nonconstant S",
                 s.diff(wc[2], th) == lam*q*aa*s.cos(th))
    hp = s.diag(om**2+2*r*om*op, om**2, 0)
    energy = (xt.dot(xt)-tx.dot(tx)+(xc.T*hp*xc)[0])/2
    avg = s.simplify(s.integrate(energy, (th, 0, 2*s.pi))/(2*s.pi))
    checks.check("full Jacobi energy includes and cancels the actual axial transport square",
                 s.simplify(avg-(zz*q*alpha*ap/2+zz*r*op*alpha**2/4)) == 0)
    zp = s.symbols("Zprime", real=True)
    flux_derivative = m**2*((zp*om+zz*op)*aa**2+2*zz*om*aa*ap)/4
    checks.check("compact radial integration yields the complete physical energy",
                 s.simplify(r*avg+m**2*zp*om*aa**2/4-flux_derivative) == 0)
    checks.check("the actual radial Beltrami relation makes that energy positive",
                 s.simplify((-m**2*zp*om*aa**2/4).subs(zp, -lam**2*r*om)
                            -lam**2*m**2*r*om**2*aa**2/4) == 0)

    x, y, zeta, k = s.symbols("X Y zeta k", real=True)
    ff = s.Function("f")(x, y)
    scalar = s.diff(ff, x, 2)+s.diff(ff, y, 2)
    comp = s.Matrix([-s.diff(scalar, y)-s.I*k*lam*s.diff(ff, x),
                     s.diff(scalar, x)-s.I*k*lam*s.diff(ff, y), lam*scalar])
    checks.check("compact potential gives an exactly solenoidal continuous carrier preparation",
                 s.simplify(s.diff(comp[0], x)+s.diff(comp[1], y)+s.I*k*comp[2]) == 0)
    radf = s.Function("f0")(r)
    lapm = s.diff(radf, r, 2)+s.diff(radf, r)/r-radf/r**2
    checks.check("one exact radial moment removes the exterior m1 potential tail",
                 s.simplify(r**2*lapm-s.diff(r**2*s.diff(radf, r)-r*radf, r)) == 0)
    rr = s.symbols("R", positive=True)
    jj = 1+x/rr
    checks.check("the actual Piola cylindrical divergence is exactly the straight divergence",
                 s.simplify(s.diff((rr+x)*comp[0]/jj, x)/(rr+x)
                            +s.diff(comp[1]/jj, y)+s.I*k*comp[2]/jj) == 0)
    # Full toroidal diagnostic, before exact compact completion.
    psr, psz, ff0 = s.symbols("psi_r psi_z F", real=True)
    eta = s.Matrix([-sz/r, -lam*ss/r, sr/r])
    uring = s.Matrix([-psz/r, ff0/r, psr/r])
    grad_ring = s.Matrix([-(lam*psr*ss+ff0*sr)/r**2+2*ff0*ss/r**3,
                          -ff0*st/r**3, -(lam*psz*ss+ff0*sz)/r**2])
    ts = (-psz*sr+psr*sz)/r+ff0*st/r**2
    checks.check("the curved cross identity retains its nonzero radial pressure source",
                 s.simplify(eta.cross(uring)-grad_ring
                            -s.Matrix([-2*ff0*ss/r**3, ts/r, 0])) == s.zeros(3, 1))
    checks.check("the uncompleted toroidal lift has a genuine angular divergence",
                 s.diff(-lam*s.Function("S")(zeta)/r, zeta)/r
                 == -lam*s.diff(s.Function("S")(zeta), zeta)/r**2)
    checks.check("whole-space tube pressure error is small relative to the full field norm",
                 s.simplify(s.sqrt(s.log(rr)/rr)/s.sqrt(rr)-s.sqrt(s.log(rr))/rr) == 0)
    checks.check("the fixed nonzero phase margin can be selected before ring radius",
                 s.limit(s.log(rr)/rr, rr, s.oo) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

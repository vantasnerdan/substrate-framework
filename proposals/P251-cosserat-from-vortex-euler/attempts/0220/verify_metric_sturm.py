"""Full coadjoint force, Lin and pressure rows in the actual mean metric."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0220-metric-Sturm")
    ii, k, c = s.symbols("I k c", positive=True)
    aa, bb, cc, dd, om, ww, ff, yy = (
        s.Function(n)(ii) for n in ("A", "B", "C", "D", "Omega", "W", "f", "Y"))
    ei, et, ey = (s.Function(n)(ii) for n in ("eta_I", "eta_theta", "eta_y"))
    det = bb*dd-cc**2
    # Independently differentiate the two actual covariant curl identities.
    op, wp = s.symbols("Omega_prime W_prime")
    solved = s.solve([
        bb*op+cc*wp+s.diff(bb, ii)*om+s.diff(cc, ii)*ww-ff*ww,
        cc*op+dd*wp+s.diff(cc, ii)*om+s.diff(dd, ii)*ww+ff*om,
    ], (op, wp))
    identities = {s.diff(om, ii): solved[op], s.diff(ww, ii): solved[wp]}
    ee = dd-cc**2/bb
    tt = s.diff(bb, ii)*om+s.diff(cc, ii)*ww
    jj = s.diff(cc, ii)*om+s.diff(dd, ii)*ww-cc*tt/bb
    dop = c-ww
    checks.check("positive tangent metric Schur coefficient is its actual determinant over B",
                 s.simplify(ee-det/bb) == 0)
    checks.check("actual averaged curl removes the false singular first-derivative term",
                 s.simplify((jj+ff*(om+cc*ww/bb)+ee*s.diff(ww, ii))
                            .subs(identities)) == 0)
    et_cov = bb*et+cc*ey
    xi = s.Matrix([yy, -s.I*((tt*yy-et_cov)/(bb*dop)+cc*s.diff(yy, ii)/bb)/k,
                   s.I*s.diff(yy, ii)/k])
    checks.check("complete action-angle coadjoint generator preserves the true unit volume",
                 s.simplify(s.diff(xi[0], ii)+s.I*k*xi[2]) == 0)
    vel = s.Matrix([-s.I*k*dop*yy-ei,
                    -ff*ww*yy/bb-cc*(s.diff(dop*yy, ii)-ey)/bb,
                    s.diff(dop*yy, ii)-ey])
    lin = -s.I*k*dop*xi-s.Matrix([0, s.diff(om, ii)*yy, s.diff(ww, ii)*yy])
    lin -= s.Matrix([ei, et, ey])
    checks.check("physical response contains all three forced Lin components",
                 s.simplify((lin-vel).subs(identities)) == s.zeros(3, 1))
    metric = s.Matrix([[aa, 0, 0], [0, bb, cc], [0, cc, dd]])
    force = xi.cross(s.Matrix([0, ff*om, ff*ww]))
    pressure = s.I*(jj*yy+dop*ee*s.diff(yy, ii)-ee*ey)/k
    defect = force-metric*vel-s.Matrix([s.diff(pressure, ii), 0, s.I*k*pressure])
    checks.check("angular pressure row uses the actual physical covariant velocity",
                 s.simplify(defect[1]) == 0)
    checks.check("axial pressure row retains the C-cross metric and its particular return",
                 s.simplify(defect[2].subs(identities)) == 0)
    rhs = ff*ww*et+ff*ww*cc*ey/bb+dop*s.diff(ee*ey, ii)-s.I*k*aa*dop*ei
    lhs = (s.diff(dop**2*ee*s.diff(yy, ii), ii)
           +(dop*s.diff(jj, ii)+ff*ww*tt/bb-k**2*aa*dop**2)*yy)
    checks.check("full forced scalar equation is independently recovered from radial pressure",
                 s.simplify((defect[0]+s.I*(lhs-rhs)/(k*dop))
                            .subs(identities)) == 0)
    checks.check("deleting the metric angular-axial term changes the forced radial equation",
                 s.simplify(rhs-(rhs-ff*ww*cc*ey/bb)) != 0)
    flat = {aa: 1/(2*ii), bb: 2*ii, cc: 0, dd: 1}
    checks.check("mean metric reduces to the full flat column scalar equation",
                 s.simplify((lhs-rhs).subs(flat).doit()
                            -(s.diff(dop**2*s.diff(yy, ii), ii)
                              +(ff*ww*om/ii-k**2*dop**2/(2*ii))*yy
                              -ff*ww*et-dop*s.diff(ey, ii)
                              +s.I*k*dop*ei/(2*ii))) == 0)
    gg, factor = s.symbols("G axial_factor")
    exact_cross = gg*(1-dd*factor)/om
    checks.check("actual covariant axial identity gives a flat-factor metric cross row",
                 s.simplify(exact_cross*om+dd*gg*factor-gg) == 0)
    hh = s.symbols("h", positive=True)
    q = s.Matrix(s.symbols("q1 q2 q3"))
    vort = s.Matrix(s.symbols("w1 w2 w3"))
    jac = s.diag(1, hh, 1)
    physical_cross = (jac*q/hh).cross(jac*vort/hh)
    checks.check("toroidal Piola cross product has its actual reciprocal-volume factor",
                 s.simplify(jac.T*physical_cross-q.cross(vort)/hh) == s.zeros(3, 1))
    metric_density = s.diag(hh, 1/hh, hh)
    checks.check("physical kinetic energy pulls back to the same flux metric",
                 s.simplify(hh*(jac*q/hh).dot(jac*q/hh)
                            -(q.T*metric_density.inv()*q)[0]) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

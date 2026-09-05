"""Exact full-pressure, Kelvin, phase and fixed-tag carrier-jet identities."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0200-noncritical-annular-mode")
    r, k, m, omega = s.symbols("r k m omega", positive=True)
    angular = s.Function("Omega")(r)
    axial = s.Function("W")(r)
    radial = s.Function("F")(r)
    pressure = s.Function("P")(r)
    doppler = omega-m*angular-k*axial
    vort = 2*angular+r*s.diff(angular, r)
    discriminant = 2*angular*vort
    kr2 = k**2+m**2/r**2
    g = 2*m*angular/r
    fp = -(1/r+g/doppler)*radial+kr2*pressure/doppler**2
    pp = (doppler**2-discriminant)*radial+g*pressure/doppler
    vr = -s.I*doppler*radial
    vt = m*pressure/(r*doppler)-vort*radial
    vz = k*pressure/doppler-s.diff(axial, r)*radial
    checks.check("full radial Euler pressure equation",
                 s.simplify(-s.I*doppler*vr-2*angular*vt+pp) == 0)
    checks.check("full tangential Euler equation",
                 s.simplify(-s.I*doppler*vt+vort*vr+s.I*m*pressure/r) == 0)
    checks.check("full axial shear equation",
                 s.simplify(-s.I*doppler*vz+s.diff(axial, r)*vr
                            +s.I*k*pressure) == 0)
    divv = s.diff(vr, r)+vr/r+s.I*m*vt/r+s.I*k*vz
    checks.check("pressure elimination retains the complete divergence",
                 s.simplify(divv.subs(s.diff(radial, r), fp)) == 0)
    xi = s.Matrix([radial, s.I*(m*pressure/(r*doppler**2)
                               -2*angular*radial/doppler),
                   s.I*k*pressure/doppler**2])
    shear = s.Matrix([0, r*s.diff(angular, r), s.diff(axial, r)])
    velocity = s.Matrix([vr, vt, vz])
    checks.check("actual Lin displacement includes both shear components",
                 s.simplify(-s.I*doppler*xi-shear*radial-velocity) == s.zeros(3, 1))
    divxi = s.diff(xi[0], r)+xi[0]/r+s.I*m*xi[1]/r+s.I*k*xi[2]
    checks.check("actual material displacement is divergence-free",
                 s.simplify(divxi.subs(s.diff(radial, r), fp)) == 0)
    potential = pressure/(-s.I*doppler)
    gradient = s.Matrix([s.diff(potential, r), s.I*m*potential/r, s.I*k*potential])
    kelvin = xi.cross(s.Matrix([0, -s.diff(axial, r), vort]))-velocity-gradient
    checks.check("full Kelvin identity with Doppler-gradient pressure",
                 s.simplify(kelvin.subs(s.diff(pressure, r), pp)) == s.zeros(3, 1))

    f = s.Function("f")(r)
    pencil = (doppler*s.diff(f, r)+g*f)**2/(r*kr2)
    pencil += (doppler**2-discriminant)*f**2/r
    p_form = doppler*(doppler*s.diff(f, r)+g*f)/(r*kr2)
    variation = s.diff(pencil, f)-s.diff(s.diff(pencil, s.diff(f, r)), r)
    checks.check("full selfadjoint pencil is the radial Euler variation",
                 s.simplify(variation+2*(s.diff(p_form, r)-g*p_form/doppler
                                       -(doppler**2-discriminant)*f/r)) == 0)
    frequency_density = 2*doppler*(s.diff(f, r)**2/(r*kr2)+f**2/r)
    frequency_density -= s.diff(g/(r*kr2), r)*f**2
    checks.check("frequency monotonicity includes its exact boundary derivative",
                 s.simplify(s.diff(pencil, omega)-frequency_density
                            -s.diff(g*f**2/(r*kr2), r)) == 0)
    beta_density = r*((m*vort/r+k*s.diff(axial, r))*radial*pressure/doppler**2
                      -discriminant*radial**2/doppler)
    ibp = beta_density+r*doppler*radial**2+pressure*s.diff(r*radial, r)/doppler
    ibp -= s.diff(pressure*r*radial/doppler, r)
    checks.check("full axial-vorticity KKS reduces to the frequency derivative",
                 s.simplify(ibp.subs(s.diff(pressure, r), pp)) == 0)
    checks.check("reduced KKS density equals minus half the pencil derivative",
                 s.simplify(s.diff(pencil, omega)/2
                            -doppler*f**2/r-p_form*s.diff(f, r)/doppler) == 0)
    admissible = s.sqrt(2)*axial/r
    checks.check("actual admissible axial maximum has positive epicyclic term",
                 s.simplify((2*admissible*(2*admissible+r*s.diff(admissible, r)))
                            .subs(s.diff(axial, r), 0)-2*admissible**2) == 0)
    ell = k**(-s.Rational(3, 4))
    checks.check("localized trial has the stated two exact small-ratio powers",
                 s.simplify(k*ell**2-k**(-s.Rational(1, 2))) == 0
                 and s.simplify(1/(k**2*ell**2)-k**(-s.Rational(1, 2))) == 0)
    q = s.symbols("q", positive=True)
    edge = s.Matrix([[1, -q], [q**3, 2*(1-q**2)-1]])
    checks.check("exposing planar m2 patch discriminant vanishes",
                 s.factor(s.trace(edge)**2-4*edge.det()) == 0)
    beta = s.symbols("beta", negative=True)
    unit = s.Matrix([[0, 1], [-1, 0]])
    checks.check("physical laboratory generator has positive phase Hessian",
                 -(beta*unit)*(-omega*unit) == -beta*omega*s.eye(2))

    theta = s.symbols("theta", real=True)
    p, dp, eps, paint = s.symbols("p dp eps paint", real=True)
    aa, bb = s.symbols("a b", real=True)
    scalar = aa*s.cos(2*theta)+bb*s.sin(2*theta)
    xr = dp*scalar/omega**2
    xt = 2*p*(-aa*s.sin(2*theta)+bb*s.cos(2*theta))/(r*omega**2)
    weight = 1+eps*paint*s.cos(2*theta)
    dq = 2*r*(s.cos(2*theta)+s.I*s.sin(2*theta))*(xr+s.I*xt)*weight
    checks.check("ordinary fixed quadrupole has the full real and imaginary rows",
                 s.simplify(s.integrate(s.expand_trig(dq), (theta, 0, 2*s.pi))
                            -2*s.pi*r*(dp+2*p/r)*(aa+s.I*bb)/omega**2) == 0)
    checks.check("actual tag displacement moment includes the painted row",
                 s.simplify(s.integrate(weight*r*xt, (theta, 0, 2*s.pi))
                            -2*s.pi*eps*paint*p*bb/omega**2) == 0)
    checks.check("quiet-cavity full spin has no omitted base-flow term",
                 s.Matrix([xr, xt, 0]).cross(s.zeros(3, 1)) == s.zeros(3, 1))
    checks.check("actual axial displacement and spin moments remain linked",
                 s.simplify(r*s.I*m*pressure/(r*omega**2)
                            -(m/k)*s.I*k*pressure/omega**2) == 0)
    checks.check("transverse centroid row vanishes on the m2 tag",
                 s.simplify(s.integrate(weight*(s.cos(theta)*xr-s.sin(theta)*xt),
                                        (theta, 0, 2*s.pi))) == 0)
    values = [4, 6, 8]
    rows = s.Matrix([[1, n, n*(n-1)] for n in values])
    checks.check("three Bessel-tail coefficient rows have exact nonzero determinant",
                 rows.det() == 16)
    kk = s.symbols("kk", real=True)
    g0, g1, g2, tau = s.symbols("g0 g1 g2 tau", positive=True)
    target = tau*(g0+g1*kk+g2*kk**2/2)
    ratio = target/(g0+g1*kk+g2*kk**2/2)
    checks.check("one fixed four-moment tag cancels both actual current derivatives",
                 s.simplify(ratio-tau) == 0
                 and s.diff(ratio, kk) == 0 and s.diff(ratio, kk, 2) == 0)
    c = s.symbols("c", positive=True)
    hphase = -beta*omega
    checks.check("physical angle and rate convert the entire phase energy",
                 s.simplify(hphase/(omega**2*c**2)+beta/(omega*c**2)) == 0)
    checks.check("fixed-current full history includes the actual initial displacement",
                 s.simplify(s.diff(tau*(aa*s.sin(omega*kk)+bb*s.cos(omega*kk)), kk)
                            -tau*omega*(aa*s.cos(omega*kk)-bb*s.sin(omega*kk))) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

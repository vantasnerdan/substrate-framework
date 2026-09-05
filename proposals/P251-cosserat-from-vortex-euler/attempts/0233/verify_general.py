"""Full general-axial first-carrier pressure/Lin and measured polynomial row."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0233-general-axial-pressure-row")
    r, t, lam, om = s.symbols("r t lambda Omega", positive=True)
    w = 2*om/lam-lam*om*r**2/2
    p, q = s.Function("p")(r), s.Function("q")(r)

    def lap(value):
        return s.diff(value, r, 2)+s.diff(value, r)/r-value/r**2

    def grad(value):
        return s.Matrix([s.diff(value, r), -s.I*value/r])

    def jgrad(value):
        return s.Matrix([s.I*value/r, s.diff(value, r)])

    f, g = lap(p), lap(q)
    aa, hh = lap(f), lap(g)
    big_h, ff = lam*p-q, lam*f-g
    source = lam*aa-hh
    u_inverse = w*ff+2*lam*om*(r*s.diff(big_h, r)-big_h)
    dilation_inverse = r*s.diff(big_h, r)-2*big_h
    zinitial_potential = s.I*lam*(w*g+lam*om*r*s.diff(q, r))
    zinitial = s.I*lam*w*hh-s.I*lam**2*r*om*s.diff(g, r)
    checks.check("actual initial Kelvin curl has its complete compact pressure potential",
                 s.simplify(lap(zinitial_potential)-zinitial) == 0)
    checks.check("retained axial shear changes the inverse of the carrier forcing",
                 s.simplify(lap(u_inverse)-w*source) == 0)
    checks.check("both compact potentials represent the actual prescribed axial lift",
                 s.simplify(lap(big_h)-ff) == 0
                 and s.simplify(lap(ff)-source) == 0)

    transverse0 = s.I*lam*om*jgrad(ff)
    axial0 = s.I*lam*om*aa-lam**2*om**2*t*ff
    vort0 = s.I*lam*om*source
    vort1 = zinitial+t*(lam*om*w*source-2*lam*om**2*aa
                       -lam**2*om**2*r*s.diff(ff, r))
    vort1 -= s.I*lam**2*om**3*t**2*ff
    true_forcing = -s.I*w*vort0+2*s.I*om*axial0-s.I*s.diff(w, r)*transverse0[1]
    checks.check("general first carrier uses the full axial-shear vorticity forcing",
                 s.simplify(s.diff(vort1, t)-true_forcing) == 0)
    stream1 = zinitial_potential+t*(lam*om*u_inverse-2*lam*om**2*f
                                   -lam**2*om**2*dilation_inverse)
    stream1 -= s.I*lam**2*om**3*t**2*big_h
    grad1 = lam*om*f+s.I*lam**2*om**2*t*big_h
    checks.check("full pressure reconstruction has the computed velocity curl and divergence",
                 s.simplify(lap(stream1)-vort1) == 0
                 and s.simplify(lap(grad1)+s.I*axial0) == 0)

    velocity1 = jgrad(stream1)+grad(grad1)
    xi0 = jgrad(aa)+t*transverse0
    xi1 = -s.I*grad(g)+t*(jgrad(zinitial_potential)+lam*om*grad(f))
    xi1 += t**2*(jgrad(lam*om*u_inverse-2*lam*om**2*f
                      -lam**2*om**2*dilation_inverse)
                  +s.I*lam**2*om**2*grad(big_h))/2
    xi1 -= s.I*lam**2*om**3*t**3*jgrad(big_h)/3
    xi1 -= s.I*w*t*jgrad(aa)
    xi1 += lam*om*w*t**2*jgrad(ff)/2
    checks.check("general corrected physical lift solves the complete first Lin row",
                 s.simplify(s.diff(xi1, t)-velocity1+s.I*w*xi0) == s.zeros(2, 1))
    div_xi1 = s.diff(r*xi1[0], r)/r-s.I*xi1[1]/r
    checks.check("general corrected lift remains solenoidal with its actual axial component",
                 s.simplify(div_xi1+s.I*hh) == 0)

    chi = s.Function("chi")(r)
    bp = r*s.diff(chi, r)
    dp = r**2*s.diff(chi, r)
    row = -s.I*dp*s.diff(g, r)
    row += t*(bp*w*aa-lam*bp*w*g-lam**2*om*dp*s.diff(q, r)
              +lam*om*dp*s.diff(f, r))
    row += s.I*lam*om*t**2*(bp*w*ff+lam*om*dp*s.diff(big_h, r)-om*bp*f)
    row += lam**2*om**3*t**3*bp*big_h/3
    checks.check("measured angle polynomial comes from the actual radial displacement",
                 s.simplify(dp*xi1[0]-row) == 0)

    raw_current = r*chi*(w*(s.diff(aa, r)+aa/r
                              +s.I*lam*om*t*(s.diff(ff, r)+ff/r))-om*r*hh)
    predicted = -(bp*w*aa+om*r**2*chi*(hh-lam*aa)
                  +s.I*lam*om*t*(bp*w*ff-lam*om*r**2*chi*ff))
    flux = r*chi*w*(aa+s.I*lam*om*t*ff)
    checks.check("complete moving-spin defect equals the two retained physical current rows",
                 s.simplify(raw_current-predicted-s.diff(flux, r)) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

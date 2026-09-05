"""Exact first-carrier pressure/Lin, physical row-rank and Kelvin-form checks."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0233-first-curved-physical-matrix")
    r, t, lam, om = s.symbols("r t lambda Omega", positive=True)
    w0 = 2*om/lam
    w = w0-lam*om*r**2/2
    f = s.Function("f")(r)
    aa = s.diff(f, r, 2)+s.diff(f, r)/r-f/r**2
    gf = s.Matrix([s.diff(f, r), -s.I*f/r])
    jf = s.Matrix([s.I*f/r, s.diff(f, r)])
    ja = s.Matrix([s.I*aa/r, s.diff(aa, r)])
    xi1 = lam*(-s.I+om*t)*gf+(2*s.I*lam*om*t-lam*om**2*t**2)*jf
    xi1 -= s.I*w*t*ja
    vel1 = lam*om*gf+(2*s.I*lam*om-2*lam*om**2*t)*jf
    checks.check("first carrier displacement solves the complete leading Lin equation",
                 s.simplify(s.diff(xi1, t)-vel1+s.I*w*ja) == s.zeros(2, 1))
    div_xi1 = s.diff(r*xi1[0], r)/r-s.I*xi1[1]/r
    checks.check("retained axial shear cancels the apparent displacement divergence",
                 s.simplify(div_xi1+s.I*lam*aa) == 0)
    div_v1 = s.diff(r*vel1[0], r)/r-s.I*vel1[1]/r
    curl_v1 = s.diff(r*vel1[1], r)/r+s.I*vel1[0]/r
    checks.check("full pressure gives the actual nonzero horizontal velocity divergence",
                 s.simplify(div_v1-lam*om*aa) == 0)
    checks.check("the first velocity curl retains its resonant time forcing",
                 s.simplify(curl_v1-(2*s.I*lam*om-2*lam*om**2*t)*aa) == 0)
    wrong_xi = xi1+s.I*(w-w0)*t*ja
    wrong_div = s.diff(r*wrong_xi[0], r)/r-s.I*wrong_xi[1]/r+s.I*lam*aa
    checks.check("dropping the small relative axial shear produces a real Lin error",
                 s.simplify(wrong_div-lam*om*t*aa) == 0
                 and s.simplify(wrong_div) != 0)

    chi = s.Function("chi")(r)
    bweight = r*s.diff(chi, r)
    dweight = -s.diff(r**2*s.diff(chi, r), r)
    aweight = s.diff(bweight, r, 2)-s.diff(bweight, r)/r
    nweight = s.diff(r**2*chi, r, 2)-s.diff(r**2*chi, r)/r
    checks.check("direct toroidal current is B(f)-D(f), not both rows separately",
                 s.simplify(nweight-bweight+dweight) == 0)
    exp_chi = {chi: s.exp(-r)}
    weights3 = [s.simplify(v.subs(exp_chi).doit()/s.exp(-r))
                for v in (bweight, dweight, aweight)]
    matrix3 = s.Matrix([[v.subs(r, point) for point in (1, 2, 3)]
                       for v in weights3])
    det3 = s.factor(matrix3.det())
    print("three physical radial rows:", matrix3, "determinant:", det3)
    checks.check("actual three-bump tag matrix has a nonzero determinant",
                 det3 == s.Rational(17, 3))
    gaussian = {chi: s.exp(-r**2)}
    bad_weights = [s.simplify(v.subs(gaussian).doit()/s.exp(-r**2))
                   for v in (bweight, dweight, aweight)]
    bad_matrix = s.Matrix([[v.subs(r, point) for point in (1, 2, 3)]
                          for v in bad_weights])
    checks.check("a Gaussian tag exposes why generic rank cannot be presumed",
                 bad_matrix.det() == 0)

    eweight = dweight-(s.diff(r**3*s.diff(chi, r), r, 2)
                       -s.diff(r**3*s.diff(chi, r), r)/r)/2
    ddweight = s.diff(dweight, r, 2)-s.diff(dweight, r)/r
    weights4 = [s.simplify(v.subs(exp_chi).doit()/s.exp(-r))
                for v in (bweight, aweight, eweight, ddweight)]
    matrix4 = s.Matrix([[v.subs(r, point) for point in (1, 2, 3, 4)]
                       for v in weights4])
    det4 = s.factor(matrix4.det())
    print("four axial-return rows:", matrix4, "determinant:", det4)
    checks.check("actual axial-return matrix has four independent physical rows",
                 det4 != 0)
    hh = s.Function("H")(r)
    lap_h = s.diff(hh, r, 2)+s.diff(hh, r)/r-hh/r**2
    dil_h = r*s.diff(hh, r)-2*hh
    lap_dil = s.diff(dil_h, r, 2)+s.diff(dil_h, r)/r-dil_h/r**2
    checks.check("compact pressure dilation identity retains its inverse source",
                 s.simplify(lap_dil-r*s.diff(lap_h, r)) == 0)

    bh, dh, bdelh, brrdelh, bf, bw_a, df = s.symbols(
        "BH DH BDeltaH Br2DeltaH Bf BWA Df", real=True)
    bw_f = w0*bdelh-lam*om*brrdelh/2
    conditions = {bh: 0, bdelh: 0, brrdelh: 2*dh, bf: 0}
    c1 = bw_f-lam*om*(bh-dh)
    quadratic = bw_f+lam*om*dh-om*bf
    checks.check("computed axial-return rows cancel the full linear-time spin defect",
                 s.simplify(c1.subs(conditions)) == 0)
    checks.check("the same actual matrix cancels both unwanted first-carrier time rows",
                 s.simplify(quadratic.subs(conditions)) == 0
                 and (lam**2*om**3*bh/3).subs(conditions) == 0)
    dg = lam*df+bw_a/om
    checks.check("constant spin current fixes the initial physical amplitude row",
                 s.simplify(bw_a+om*(lam*df-dg)) == 0)

    radial_om = s.Function("Omega_r")(r)
    radial_w = s.Function("W_r")(r)
    ar = s.Function("A_r")(r)
    fr = s.Function("f_r")(r)
    m = s.symbols("m", integer=True, nonzero=True)
    raw = r*radial_om*radial_w*ar**2+lam*r**2*radial_om**2*ar*s.diff(fr, r)
    raw -= m**2*radial_om*radial_w*ar*fr/r
    raw -= r*radial_w*(s.diff(radial_om, r)*ar+radial_om*s.diff(ar, r))*s.diff(fr, r)
    flux = r*radial_om*radial_w*ar*s.diff(fr, r)
    energy_remainder = raw-2*r*radial_om*radial_w*ar**2+s.diff(flux, r)
    energy_remainder = energy_remainder.subs(s.diff(radial_w, r), -lam*r*radial_om)
    energy_remainder = energy_remainder.subs(
        s.diff(fr, r, 2), ar-s.diff(fr, r)/r+m**2*fr/r**2)
    checks.check("actual first-carrier energy has its full radial boundary cancellation",
                 s.simplify(energy_remainder) == 0)
    cc, energy, remainder = s.symbols("c energy remainder", real=True)
    full_energy = energy+(1-cc)**2*remainder-2*(1-cc)*energy
    checks.check("mixed Kelvin signature keeps the positive horizontal kinetic remainder",
                 s.expand(full_energy-((2*cc-1)*energy+(1-cc)**2*remainder)) == 0)
    nup, nun, amp = s.symbols("nu_positive nu_negative amplitude", positive=True)
    checks.check("zero-energy return supplies its actual nonzero phase gap",
                 s.simplify(amp/nup-amp/nun-amp*(nun-nup)/(nun*nup)) == 0)
    checks.check("zero-phase return supplies a separately adjustable energy",
                 s.simplify(amp/nup-(amp*nun/nup)/nun) == 0
                 and s.simplify(amp-amp*nun/nup-amp*(1-nun/nup)) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

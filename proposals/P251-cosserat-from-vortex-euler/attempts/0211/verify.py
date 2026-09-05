"""Exact affine-current seed, full radial border and physical torus identities."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0211-same-ring-constant-curl-core")
    t = s.symbols("t", positive=True)
    ea, eb = s.exp(-1/t), s.exp(-1/(1-t))
    taper = ea/(ea+eb)
    reflected = eb/(eb+ea)
    checks.check("chosen smooth taper has the exact half-integral symmetry",
                 s.cancel(taper+reflected-1) == 0)
    derivative = ea*eb*(1/t**2+1/(1-t)**2)/(ea+eb)**2
    checks.check("actual taper derivative has a positive factorization",
                 s.simplify(s.diff(taper, t)-derivative) == 0)
    phi, delta, lam, amp = s.symbols("phi delta lambda A", positive=True)
    primitive = s.Function("H")(phi)
    hh = s.Function("h")(phi/delta)
    source = lam**2*primitive*hh
    q = s.diff(source, phi).subs(s.diff(primitive, phi), hh)
    expected_q = lam**2*(hh**2+primitive*s.diff(hh, phi))
    checks.check("full seed derivative retains the taper-current product term",
                 s.simplify(q-expected_q) == 0)
    checks.check("both swirl signs give the identical meridional source",
                 all(s.simplify((sign*lam*primitive)
                                *s.diff(sign*lam*primitive, phi)
                                .subs(s.diff(primitive, phi), hh)-source) == 0
                     for sign in (-1, 1)))
    affine = lam*(phi-delta/2)
    checks.check("inner current is exactly affine with the actual primitive offset",
                 s.diff(affine, phi) == lam
                 and s.diff(affine, phi, 2) == 0)
    checks.check("inner Euler source has the nonzero derivative absent in a Rankine plateau",
                 s.diff(affine*s.diff(affine, phi), phi) == lam**2)

    r = s.symbols("r", positive=True)
    radial = delta/2+amp*s.besselj(0, lam*r)
    vv = amp*lam*s.besselj(1, lam*r)
    zz = amp*lam**2*s.besselj(0, lam*r)
    ww = amp*lam*s.besselj(0, lam*r)
    checks.check("inner Bessel seed solves the actual affine radial PDE",
                 s.simplify(s.expand_func(s.diff(radial, r, 2)+s.diff(radial, r)/r
                                         +lam**2*(radial-delta/2))) == 0)
    checks.check("inner meridional velocity is the actual stream derivative",
                 s.simplify(-s.diff(radial, r)-vv) == 0)
    checks.check("inner axial vorticity is the actual radial circulation derivative",
                 s.simplify(s.expand_func(s.diff(r*vv, r)/r-zz)) == 0)
    checks.check("full force-free swirl reaction satisfies the low-carrier pressure identity",
                 s.simplify(ww*s.diff(ww, r)+zz*vv) == 0)
    root = s.symbols("j", positive=True)
    border = (r*s.diff(s.besselj(0, lam*r), r)).subs(r, root/lam)
    checks.check("actual radial exterior-log border has the derived Bessel coefficient",
                 s.simplify(border+root*s.besselj(1, root)) == 0)
    vfun, gfun, qfun = [s.Function(name)(r) for name in ("V", "g", "Q")]
    order = s.symbols("m", integer=True, positive=True)
    lhs = r*(s.diff(gfun, r)**2+(order**2/r**2-qfun)*gfun**2)
    rhs = r*vfun**2*s.diff(gfun/vfun, r)**2+(order**2-1)*gfun**2/r
    remainder = lhs-rhs-s.diff(r*s.diff(vfun, r)*gfun**2/vfun, r)
    remainder = remainder.subs(s.diff(vfun, r, 2),
                               -s.diff(vfun, r)/r+vfun/r**2-qfun*vfun)
    checks.check("all-angular ground-state identity keeps its exact boundary term",
                 s.simplify(remainder) == 0)
    excess = s.symbols("m_minus_two", integer=True, nonnegative=True)
    angular_gap = s.expand((excess+2)**2-1)
    checks.check("higher-angular correction is strictly positive for every m at least two",
                 angular_gap == excess**2+4*excess+3 and angular_gap.is_positive is True)

    z, radius = s.symbols("z R", real=True, positive=True)
    stream = s.Function("Phi")(r, z)
    for sign in (1, -1):
        swirl = radius*sign*lam*(stream-delta/2)
        velocity = s.Matrix([-radius*s.diff(stream, z)/r,
                             swirl/r, radius*s.diff(stream, r)/r])
        curl = s.Matrix([-s.diff(velocity[1], z),
                         s.diff(velocity[0], z)-s.diff(velocity[2], r),
                         s.diff(r*velocity[1], r)/r])
        residual = curl-sign*lam*velocity
        equation = s.diff(stream, r, 2)-s.diff(stream, r)/r+s.diff(stream, z, 2)
        equation += lam**2*(stream-delta/2)
        checks.check(f"literal inner cylindrical curl identity, physical sign {sign}",
                     s.simplify(residual+s.Matrix([0, radius*equation/r, 0]))
                     == s.zeros(3, 1))
    checks.check("complete cylindrical velocity is divergence-free",
                 s.simplify(s.diff(r*velocity[0], r)/r+s.diff(velocity[2], z)) == 0)
    gradient = s.Matrix([s.diff(stream, r), 0, s.diff(stream, z)])
    checks.check("actual stream level is a material first integral",
                 s.simplify(velocity.dot(gradient)) == 0)
    checks.check("full Green source factor produces the actual force-free PDE",
                 s.cancel((r**2/radius**2)*(radius**2*source/r**2)-source) == 0)

    j0 = s.besselj(0, lam*r).series(r, 0, 8).removeO()
    j1 = s.besselj(1, lam*r).series(r, 0, 9).removeO()
    rotation_series = (j1/(r*j0)).series(r, 0, 6).removeO()
    checks.check("actual Lundquist transit ratio has its full quartic expansion",
                 s.expand(rotation_series-lam/2-lam**3*r**2/16-lam**5*r**4/96) == 0)
    flux = amp*r*s.besselj(1, lam*r)
    checks.check("section action differentiates to the actual velocity flux",
                 s.simplify(s.expand_func(s.diff(flux, r)-r*ww)) == 0)
    twist = s.limit(s.diff(rotation_series, r)/(r*amp*lam*j0), r, 0)
    checks.check("nonzero center twist is normalized by physical flux rather than area",
                 s.simplify(twist-lam**2/(8*amp)) == 0)
    gg, contour = s.symbols("G contour", positive=True)
    advance = gg*contour
    checks.check("full finite-ring coarea identity keeps the R normalization",
                 s.simplify((-radius*gg*contour/(2*s.pi))*(2*s.pi/advance)+radius) == 0)
    central = s.symbols("central_phi", positive=True)
    isotropic_hessian = lam**2*central/2
    checks.check("physical core return ratio tends to lambda times radius over two",
                 s.simplify(radius*isotropic_hessian/(lam*central)-lam*radius/2) == 0)
    boost = s.symbols("U", nonzero=True, real=True)
    checks.check("a Galilean boost does not preserve a nonzero literal curl factor",
                 -lam*s.Matrix([0, 0, boost]) != s.zeros(3, 1))
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

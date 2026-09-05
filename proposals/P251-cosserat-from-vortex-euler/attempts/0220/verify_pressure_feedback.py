"""Exact pressure, gauge-moment and full phase-Schur algebra."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0220-pressure-feedback")
    eps = s.symbols("epsilon", real=True)
    d11, d12, d13, d22, d23, d33 = s.symbols("d11 d12 d13 d22 d23 d33")
    ee = s.diag(*s.symbols("e1 e2 e3"))
    dd = s.Matrix([[d11, d12, d13], [d12, d22, d23], [d13, d23, d33]])
    nn = s.eye(3)+eps*dd+eps**2*ee
    grad = s.Matrix([0, 0, 1])
    pp = s.eye(3)-grad*grad.T
    exact = nn-nn*grad*grad.T*nn/(grad.T*nn*grad)[0]
    checks.check("exact pressure solve annihilates its unchanged gradient subspace",
                 s.simplify(exact*grad) == s.zeros(3, 1))
    checks.check("first pressure coefficient is P D P with both projections",
                 s.simplify(exact.diff(eps).subs(eps, 0)-pp*dd*pp) == s.zeros(3))
    checks.check("second pressure coefficient retains the eliminated gradient reaction",
                 s.simplify(exact.diff(eps, 2).subs(eps, 0)/2
                            -pp*ee*pp+pp*dd*(s.eye(3)-pp)*dd*pp) == s.zeros(3))
    checks.check("deleting the inner projection creates a false gradient response",
                 s.simplify((pp*dd)*grad) != s.zeros(3, 1))
    a, b = s.symbols("a b", positive=True)
    # Same divergence constraint, two actual positive kinetic metrics.
    basis = s.Matrix([[1, 0], [0, 1], [1, 1]])
    metric_a = s.diag(a, 1, 1)
    metric_b = s.diag(1, b, 1)
    sa = basis*(basis.T*metric_a*basis).inv()*basis.T
    sb = basis*(basis.T*metric_b*basis).inv()*basis.T
    checks.check("finite constrained metric inverse has the exact physical-velocity feedback identity",
                 s.simplify(sa-sb+sa*(metric_a-metric_b)*sb) == s.zeros(3))
    n, xx, c0, log_coefficient = s.symbols("n x C0 L")
    fy, fx, xfy = s.symbols("force_y_mass force_x_mass force_y_xmoment")
    kk = n*eps
    mass = s.I*kk*fy
    first = -fx+s.I*kk*xfy
    pressure = c0*mass+eps*log_coefficient*(xx*mass+first)
    row = s.Matrix([s.diff(pressure, xx), s.I*kk*pressure])
    checks.check("actual divergence moments suppress both constant and affine-log return through first order",
                 row.subs(eps, 0) == s.zeros(2, 1)
                 and row.diff(eps).subs(eps, 0) == s.zeros(2, 1))
    fake_mass = s.symbols("independent_monopole")
    fake_pressure = c0*fake_mass+eps*log_coefficient*(xx*fake_mass+first)
    checks.check("an invented independent pressure monopole would invalidate the return order",
                 s.diff(fake_pressure, xx, eps).subs(eps, 0) != 0)
    ii, cc, wave = s.symbols("I c k", positive=True)
    ww, yy = (s.Function(nm)(ii) for nm in ("w", "Y"))
    raw = ((cc+ww)**2*s.diff(yy, ii)**2
           +wave**2*(cc+ww)**2*yy**2/(2*ii)+ww*s.diff(ww, ii)*yy**2/ii)
    squares = (ww**2*((s.diff(yy, ii)-yy/(2*ii))**2+yy**2/(4*ii**2)
                      +wave**2*yy**2/(2*ii))
               +(cc**2+2*cc*ww)*(s.diff(yy, ii)**2+wave**2*yy**2/(2*ii)))
    checks.check("actual-I Hardy decomposition includes its complete boundary derivative",
                 s.simplify(raw-squares-s.diff(ww**2*yy**2/(2*ii), ii)) == 0)
    jj = s.Function("J")(ii)
    checks.check("metric J-prime is paired before estimating the flat edge",
                 s.simplify(-(cc+ww)*s.diff(jj, ii)*yy**2
                            -jj*(s.diff(ww, ii)*yy**2+2*(cc+ww)*yy*s.diff(yy, ii))
                            +s.diff((cc+ww)*jj*yy**2, ii)) == 0)
    j, r1, r2, h1, h2, e1, e2 = s.symbols("j r1 r2 h1 h2 eta1 eta2")
    response = s.Matrix([[a, 1], [0, b]])
    pairing = s.Matrix([[r1, r2]])
    column = s.Matrix([h1, h2])
    forcing = s.Matrix([e1, e2])
    correction = (pairing*response*column)[0]
    phase = -(pairing*response*forcing)[0]/(j+correction)
    fast = response*(forcing+column*phase)
    checks.check("full finite phase correction enforces exact geometric KKS orthogonality",
                 s.simplify(j*phase+(pairing*fast)[0]) == 0)
    checks.check("the same finite correction satisfies the unabridged fast equation",
                 s.simplify(response.inv()*fast-forcing-column*phase) == s.zeros(2, 1))
    freq, p1, p2, k1, k2, b1, b2, h0 = s.symbols("omega p1 p2 k1 k2 b1 b2 h0")
    fast_form = s.diag(p1-freq*k1, p2-freq*k2)
    coupling = s.Matrix([b1, b2])
    mode_fast = -fast_form.inv()*coupling
    schur = h0-freq*j-(coupling.T*fast_form.inv()*coupling)[0]
    inherited = j+(mode_fast.T*s.diag(k1, k2)*mode_fast)[0]
    checks.check("Schur frequency derivative is the inherited whole-mode Krein form",
                 s.simplify(-s.diff(schur, freq)-inherited) == 0)
    checks.check("omitting the reaction derivative loses a real fast-phase action row",
                 s.simplify(inherited-j) != 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

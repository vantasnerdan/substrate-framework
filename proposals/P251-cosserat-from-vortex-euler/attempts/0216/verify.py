"""Exact actual-cell orbit and stationary-tag harmonic coefficients."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    checks = CheckLedger("P251-0216-fixed-tag-harmonics")
    e, w = s.symbols("epsilon Omega", positive=True)
    t = s.symbols("theta", real=True)
    a3 = (3*w**2-1)/192
    b1 = w*(1-w**2)/16
    b3 = w*(3-w**2)/192
    w2 = -w*(1+w**2)/16
    a = e*s.cos(t)+e**3*a3*s.cos(3*t)
    b = -w*e*s.sin(t)+e**3*(b1*s.sin(t)+b3*s.sin(3*t))
    frequency = w+e**2*w2
    adot = frequency*s.diff(a, t)
    bdot = frequency*s.diff(b, t)
    ar = s.series(adot-b+b**3/6, e, 0, 5).removeO()
    br = s.series(bdot+w**2*(a-a**3/6), e, 0, 5).removeO()
    checks.check("actual a orbit equation through cubic order", s.trigsimp(s.expand_trig(ar)) == 0)
    checks.check("actual b orbit equation through cubic order", s.trigsimp(s.expand_trig(br)) == 0)
    energy = s.series((w**2*a**2+b**2)/2-(w**2*a**4+b**4)/24, e, 0, 6).removeO()
    checks.check("orbit energy is independent of the physical time angle", s.trigsimp(s.expand_trig(s.diff(energy, t))) == 0)
    covariance = s.series(a*b, e, 0, 6).removeO().expand()
    dilation = s.series(-w**2*a**2-b**2+(w**2*a**4+b**4)/6, e, 0, 6).removeO().expand()

    def harmonic(expr, order, trig):
        return s.simplify(s.integrate(s.expand_trig(expr*trig(order*t)), (t, -s.pi, s.pi))/s.pi)

    cov4 = harmonic(covariance, 4, s.sin)
    dil4 = harmonic(dilation, 4, s.cos)
    dil2 = harmonic(dilation, 2, s.cos)
    checks.check("actual covariance fourth harmonic includes orbital distortion", s.simplify(cov4-e**4*w*(1-w**2)/96) == 0)
    checks.check("actual dilation fourth harmonic includes orbital distortion", s.simplify(dil4-e**4*w**2*(1+w**2)/96) == 0)
    checks.check("second harmonic dilation has its distinct nonlinear order", s.simplify(dil2-e**4*w**2*(1-w**2)/24) == 0)
    rho, q = s.symbols("rho Q", positive=True)
    spin_ratio4 = s.cancel(rho*q*dil4/(4*w*cov4))
    checks.check("fourth harmonic physical spin normalization", s.cancel(spin_ratio4-rho*q*(1+w**2)/(4*(1-w**2))) == 0)
    checks.check("actual elliptic cell has positive fourth harmonic inertia", spin_ratio4.subs(w, s.Rational(1, 10)).is_positive is True)
    spin_ratio2 = s.cancel(rho*q*dil2/(2*w*(-w*e**2/2)))
    checks.check("second harmonic instead has negative leading inertia", spin_ratio2.subs(w, s.Rational(1, 10)).is_negative is True)
    wrong_dil4 = harmonic((w**2*(e*s.cos(t))**4+(-w*e*s.sin(t))**4)/6, 4, s.cos)
    checks.check("omitting exact orbit distortion changes the spin coefficient", s.simplify(wrong_dil4-dil4) != 0)
    checks.check("affine orbit alone has no fourth covariance observation", harmonic(-w*e**2*s.cos(t)*s.sin(t), 4, s.sin) == 0)
    phase = s.symbols("phase", real=True)
    checks.check("transported quadrature gives the observed angle sign", s.simplify(s.integrate(s.sin(4*t)*s.cos(4*t-phase), (t, -s.pi, s.pi))/s.pi-s.sin(phase)) == 0)
    return checks.finish()


if __name__ == "__main__":
    raise SystemExit(main())

"""Exact cell and full point-street symbol identities; no spectral numerics."""

import sympy as s

from substrate_framework.verification import CheckLedger


def main():
    ledger = CheckLedger("P251-0134-street")
    x, y, aa, bb = s.symbols("x y aa bb", real=True)
    u = s.Matrix([-bb*s.sin(y), aa*s.sin(x), 0])
    w = s.Matrix([0, -bb*s.sin(y), 0])
    tangent = -u.diff(y)
    pressure = -aa*bb*s.cos(x)*s.cos(y)
    pi = aa*bb*s.sin(x)*s.sin(y)
    jac_u = u.jacobian([x, y, s.Symbol("z")])
    jac_w = w.jacobian([x, y, s.Symbol("z")])
    residual = jac_w*u+jac_u*w+s.Matrix([s.diff(pi,x),s.diff(pi,y),0])+u[0]*tangent+s.Matrix([-s.diff(pressure,y),0,0])
    ledger.check("planar corrector retains full pressure equation", all(s.simplify(v)==0 for v in residual))
    ledger.check("planar corrector exact divergence", s.simplify(s.diff(w[0],x)+s.diff(w[1],y)+tangent[0])==0)
    stress = s.integrate(u[0]*w[1],(y,-s.pi,s.pi))/(2*s.pi)
    ledger.check("planar physical response negative", s.simplify(-stress+bb**2/2)==0)

    a, b, gamma, g = s.symbols("a b gamma g", positive=True)
    d = s.symbols("d", real=True)
    energy = gamma*s.log(a*s.sqrt(s.cosh(2*s.pi*b/a)-s.cos(2*s.pi*d/a)))
    hbb = s.simplify(s.diff(energy,b,2).subs(d,a/2))
    hdd = s.simplify(s.diff(energy,d,2).subs(d,a/2))
    ledger.check("relative-width Hessian positive formula", s.simplify(hbb-gamma*s.pi**2/(a**2*s.cosh(s.pi*b/a)**2))==0)
    ledger.check("relative-offset Hessian opposite sign", s.simplify(hbb+hdd)==0)
    h = gamma*s.log(a*s.cosh(s.pi*b/a))
    det = s.diff(h,a,2)*s.diff(h,b,2)-s.diff(h,a,b)**2
    ledger.check("fixed-frame strain-momentum determinant negative", s.simplify(det+gamma**2*s.pi**2/a**4)==0)
    ledger.check("positive fluid impulse mass", s.simplify(g**2/hbb-g**2*a**2*s.cosh(s.pi*b/a)**2/(gamma*s.pi**2))==0)

    # Canonical physical coordinates (X,b,Y,d), original (x+,y+,x-,y-).
    transform = s.Matrix([[1,0,0,-s.Rational(1,2)], [0,s.Rational(1,2),1,0], [1,0,0,s.Rational(1,2)], [0,-s.Rational(1,2),1,0]])
    j = s.Matrix([[0,1],[-1,0]])
    omega = s.diag(j,-j)
    ledger.check("actual vortex form gives impulse pairs", transform.T*omega*transform == s.diag(j,j))

    k,t = s.symbols("k t", real=True)
    ff = s.pi/s.cosh(t)*((s.pi-k)*s.cosh((s.pi-k)*t/s.pi)-s.pi*s.tanh(t)*s.sinh((s.pi-k)*t/s.pi))
    gg = s.pi/s.cosh(t)*((s.pi-k)*s.sinh((s.pi-k)*t/s.pi)-s.pi*s.tanh(t)*s.cosh((s.pi-k)*t/s.pi))
    dd = s.pi*k-k**2/2-s.pi**2/s.cosh(t)**2
    cap_a = -1+2*t*s.tanh(t)+t**2/s.cosh(t)**2
    ledger.check("cross-row constant", s.simplify(ff.subs(k,0)-s.pi**2/s.cosh(t)**2)==0)
    ledger.check("neutrality cancels constant stiffness", s.simplify((ff+dd).subs(k,0))==0)
    ledger.check("neutrality cancels first-order stiffness", s.simplify(s.diff(ff+dd,k).subs(k,0))==0)
    ledger.check("full lattice derives longitudinal stiffness", s.simplify(s.diff(ff+dd,k,2).subs(k,0)-cap_a)==0)
    ledger.check("full lattice derives convective cross", s.simplify(s.diff(gg,k).subs(k,0)+s.pi*(s.tanh(t)+t/s.cosh(t)**2))==0)
    ledger.check("cross-row hyperbolic difference identity", s.simplify(ff**2-gg**2-s.pi**2/s.cosh(t)**2*((s.pi-k)**2-s.pi**2*s.tanh(t)**2))==0)
    f0 = s.pi**2/s.cosh(t)**2
    ledger.check("exact acoustic discriminant leading coefficient", s.simplify(s.diff(ff**2-dd**2,k,2).subs(k,0)/2-f0*cap_a)==0)
    ledger.check("fixed-frame drift exceeds relative acoustic speed", s.simplify((s.tanh(t)+t/s.cosh(t)**2)**2-cap_a/s.cosh(t)**2-1)==0)
    f,gf,df=s.symbols("F G D",real=True)
    z=s.diag(1,-1)
    cross=s.Matrix([[f,s.I*gf],[s.I*gf,-f]])
    full=(df*z).row_join(cross).col_join(cross.conjugate().T.row_join(df*z))
    transformed=s.simplify(transform.T*full*transform)
    expected=s.Matrix([[2*(df+f),-s.I*gf,0,0],[s.I*gf,(f-df)/2,0,0],[0,0,-2*(df+f),s.I*gf],[0,0,-s.I*gf,(df-f)/2]])
    ledger.check("full row Hessian separates longitudinal/transverse sectors", transformed==expected)
    freq=s.symbols("freq",real=True)
    block=expected[:2,:2]
    generator=s.Matrix([[0,1],[-1,0]])*block
    ledger.check("exact physical-clock dispersion", s.expand((generator+s.I*freq*s.eye(2)).det()-(f**2-df**2-(freq+gf)**2))==0)
    raise SystemExit(ledger.finish())


if __name__ == "__main__":
    main()

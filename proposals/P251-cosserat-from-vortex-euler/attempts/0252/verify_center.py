"""Exact center jets of radial Euler and its first integrals; no numerics."""
import sympy as s
from substrate_framework.verification import CheckLedger


def main():
    checks=CheckLedger('P251-0252-localizable-center')
    x,z=s.symbols('x z',real=True)
    r0=s.Symbol('r0',positive=True)
    a,b,c,d,p0,pr,pz,w0,wr,wz=s.symbols('a b c d p0 pr pz w0 wr wz',real=True)
    r=r0+x
    v=s.Matrix([a*x+b*z,c*x+d*z])
    A=v.jacobian([x,z])
    p=p0+pr*x+pz*z
    w=w0+wr*x+wz*z
    center={x:0,z:0}
    grad=lambda f:s.Matrix([s.diff(f,x),s.diff(f,z)])
    derivative=lambda f:grad(v.dot(grad(f))).subs(center)
    checks.check('actual localizability derivative is transpose poloidal Jacobian times pressure gradient',
                 s.simplify(derivative(p)-A.T*s.Matrix([pr,pz]))==s.zeros(2,1))
    checks.check('nondegenerate derivative makes the pressure gradient vanish',
                 s.simplify(A.T.inv()*derivative(p))==s.Matrix([pr,pz]))
    radial=(v.dot(grad(v[0]))-w*w/r+s.diff(p,x)).subs(center)
    checks.check('physical radial Euler retains the centripetal core pressure',
                 s.simplify(radial-(pr-w0*w0/r0))==0)
    checks.check('vanishing pressure gradient requires zero real core speed',
                 s.solve(radial.subs(pr,0),w0)==[0])
    angular=r*w
    speed=v.dot(v)+w*w
    angular_gradient=s.simplify(A.T.inv()*derivative(angular))
    speed_gradient=s.simplify(A.T.inv()*derivative(speed))
    conditions=s.solve(list(angular_gradient),(wr,wz))
    checks.check('actual angular momentum first integral fixes radial swirl derivative',
                 conditions=={wr:-w0/r0,wz:0})
    checks.check('independent speed first integral exposes the same nonzero-core obstruction',
                 s.simplify(speed_gradient.subs(conditions))==s.Matrix([-2*w0*w0/r0,0]))
    checks.check('removing cylindrical centripetal term would falsely allow arbitrary core speed',
                 (v.dot(grad(v[0]))+s.diff(p,x)).subs(center).subs(pr,0)==0
                 and radial.subs(pr,0)!=0)
    # A genuine localizable pure-swirl solution on r>0 has degenerate v.
    swirl=s.Symbol('W',nonzero=True,real=True)
    pure_pressure=swirl**2*s.log(r)
    checks.check('degenerate-center counterexample is exact steady pure swirl',
                 s.simplify(s.diff(pure_pressure,x)-swirl**2/r)==0
                 and s.zeros(2).det()==0)
    return checks.finish()


if __name__=='__main__':
    raise SystemExit(main())

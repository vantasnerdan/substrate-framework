"""Exact exposing examples for nonlinear inference, not a new Euler closure theorem."""
import sympy as s

# Same measured linear restoring coefficient, distinct periodic finite-angle laws.
theta,k,beta=s.symbols('theta k beta', real=True)
V1=k*(1-s.cos(theta))
V2=V1+beta*(1-s.cos(theta))**2
assert s.diff(V1,theta,2).subs(theta,0)==s.diff(V2,theta,2).subs(theta,0)==k
extra_torque=-s.expand_trig(s.diff(V2-V1,theta))  # mechanical torque = -dV/dtheta
assert s.simplify(extra_torque+2*beta*(1-s.cos(theta))*s.sin(theta))==0
print('Same linear stiffness:',s.diff(V1,theta,2).subs(theta,0))
print('Finite-angle torque difference:',s.factor(extra_torque))
print('First nonlinear separation:',s.series(extra_torque,theta,0,5))

# Exact incompressible steady Euler strain changes material covariance eigenvalues.
x,y,z,a,t=s.symbols('x y z a t', real=True)
coords=s.Matrix([x,y,z]);u=s.Matrix([a*x,-a*y,0]);p=-a*a*(x*x+y*y)/2
assert s.trace(u.jacobian(coords))==0
residual=s.simplify(u.jacobian(coords)*u+s.Matrix([s.diff(p,q) for q in coords]))
assert residual==s.zeros(3,1)
flow=s.diag(s.exp(a*t),s.exp(-a*t),1)
C0=s.diag(1,2,3);Ct=flow*C0*flow.T
assert flow.det()==1
assert s.simplify(Ct.det()-C0.det())==0
assert s.diff(s.trace(Ct),t).subs(t,0)==-2*a
print('Exact strain Euler residual:',residual.T)
print('Advected material covariance:',Ct)
print('Covariance eigenvalues change despite incompressibility; trace derivative at0:',s.diff(s.trace(Ct),t).subs(t,0))

# Finite rotations require composition, not additive infinitesimal angles.
Rx=s.Matrix([[1,0,0],[0,0,-1],[0,1,0]])
Ry=s.Matrix([[0,0,1],[0,1,0],[-1,0,0]])
assert Rx.T*Rx==Ry.T*Ry==s.eye(3)
assert Rx*Ry!=Ry*Rx
print('Noncommuting right-angle rotations commutator:',Rx*Ry-Ry*Rx)
print('All exact separation checks pass. Affine strain is unbounded and nonperiodic; no compactness or nonlinear closure claim.')

"""Exact generic-column and free-exterior checks; no spectral numerics."""
import sympy as s
r=s.Symbol("r", positive=True)
c=s.Symbol("c", nonzero=True, real=True)
L=s.Function("L")(r)
W=s.Function("W")(r)
q=W-c
D=lambda a:s.diff(a,r)/(r*q)
H=L*s.diff(L,r)/(r*q)
Bp=s.diff(W,r)/r+H/r**2
Q=s.simplify(D(H)-r**2*D(Bp))
assert s.simplify(Q-(2*L*s.diff(L,r)/(r**3*q**2)-(s.diff(W,r,2)-s.diff(W,r)/r)/q))==0
print("PASS full axial-shear column linearization")
D0=lambda a:-s.diff(a,r)/(r*c)
H0=-L*s.diff(L,r)/(r*c)
Q0=s.simplify(D0(H0)-r*r*D0(H0/r**2))
J=s.simplify(D0(D0(H0))-r*r*D0(D0(H0/r**2)))
assert s.simplify(J+2*s.diff(Q0,r)/(r*c))==0
print("PASS quadratic nonlinearity from actual background labels")
f=s.Function("f")(r)
Qf=(-s.diff(f,r,2)+s.diff(f,r)/r)/f
h=r*s.diff(f,r)/f
boundary=Qf*f**3/r**2+s.Rational(3,2)*f*s.diff(f,r)**2/r**2-2*f*f*s.diff(f,r)/r**3
integrand=s.diff(Qf,r)*f**3/r**2
negative=-f**3*h*(3*h*h-8*h+8)/(2*r**5)
assert s.simplify(integrand-s.diff(boundary,r)-negative)==0
x=s.Symbol("h", real=True)
assert s.expand(3*x*x-8*x+8-(3*(x-s.Rational(4,3))**2+s.Rational(8,3)))==0
print("PASS projected nonlinear sign with exact endpoint term")
k,R=s.symbols("k R", positive=True)
u=r*s.besselk(1,k*r)
assert s.simplify(s.expand_func(s.diff(u,r)+k*r*s.besselk(0,k*r)))==0
assert s.simplify(s.expand_func(s.diff(u,r,2)-s.diff(u,r)/r-k*k*u))==0
trace=u/(R*s.besselk(1,k*R))
DtN=-k*s.besselk(0,k*R)/s.besselk(1,k*R)
assert s.simplify(s.expand_func(s.diff(trace,r).subs(r,R)-DtN))==0
# Exterior quadratic energy boundary formula, before integration at infinity.
flux=s.diff(u*s.diff(u,r)/r,r)
assert s.simplify(s.expand_func(flux-(s.diff(u,r)**2+k*k*u*u)/r))==0
print("PASS exact exterior equation, DtN and kinetic boundary flux")
inside=r*s.besselj(1,k*r)
assert s.simplify(s.expand_func(s.diff(inside,r)-k*r*s.besselj(0,k*r)))==0
assert s.simplify(s.expand_func(s.diff(inside,r,2)-s.diff(inside,r)/r+k*k*inside))==0
print("PASS Rankine interior and distinct J0 versus J1 boundary rows")
# Five-dimensional desingularization derives axis regularity explicitly.
phi=s.Function("phi")(r)
assert s.expand((s.diff(r*r*phi,r,2)-s.diff(r*r*phi,r)/r)/r**2-s.diff(phi,r,2)-3*s.diff(phi,r)/r)==0
print("PASS axis-regular five-dimensional radial conjugacy")

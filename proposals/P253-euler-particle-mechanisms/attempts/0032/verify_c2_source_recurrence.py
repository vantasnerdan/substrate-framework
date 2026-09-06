#!/usr/bin/env python3
"""Exact source-equation recurrence and returned BAS jet for P253/0032."""
import sympy as s

e,q,X,Z=s.symbols('e q X Z', real=True)

def es(v,n=3):
    return s.series(v,e,0,n).removeO().expand()

def ep(v,n=3):
    v=s.expand(v)
    return s.expand(sum(v.coeff(e,j)*e**j for j in range(n)))

def mp(v,n=3):
    p=s.Poly(s.expand(v),X,Z)
    return s.Add(*[coef*X**ij[0]*Z**ij[1] for ij,coef in p.terms() if sum(ij)<n])

# Baldi (4.29)--(4.33), with P_n obtained from Gavrilov alpha_2 (4.14).
P3=4*s.sin(q)**3-s.sin(q)
P4=s.Rational(11,8)+s.sin(q)**2
W0=1/s.sqrt(2)
W1=-P3/8
W2=s.factor(-(2*W1**2+3*P3*W0**2*W1+P4*W0**4)/(4*W0))
assert s.trigsimp(4*W0*W2+2*W1**2+3*P3*W0**2*W1+P4*W0**4)==0

# F_c(g_c)=sigma.  f=theta+e*f1+e^2*f2 and its inverse are solved recursively.
Q=s.Rational(4,3)*s.cos(q)**3-3*s.cos(q)+s.Rational(5,3)
f1=-s.Rational(3,4)*Q
f2=-2*s.sin(q)**5*s.cos(q)-s.Rational(15,8)*s.sin(q)*s.cos(q)+s.sin(q)*s.cos(q)*(s.cos(q)**2-s.sin(q)**2)/2
assert s.trigsimp(s.diff(f2,q)-(3*P3**2-4*P4)/4)==0
d1=-f1
d2=f1*s.diff(f1,q)-f2
theta=q+e*d1+e**2*d2

# sqrt(2 gamma)=sqrt(c)w, c=h(I)=2e^2+O(e^6), I=e^2/2.
r2=-P3/4
r3=2*s.sqrt(2)*W2+d1*s.diff(r2,q)
r=e+e**2*r2+e**3*r3
sn=s.sin(q)+e*d1*s.cos(q)+e**2*(d2*s.cos(q)-d1**2*s.sin(q)/2)
cs=s.cos(q)-e*d1*s.sin(q)+e**2*(-d2*s.sin(q)-d1**2*s.cos(q)/2)
xx=ep(r*sn,4)
rc=ep(r*cs,4)
zz=ep(rc*(1-xx+xx**2),4)

# Gavrilov (4.4), (4.10), and Baldi (4.14): exact local field jet.
a=2*X**2+2*Z**2+3*X**3-X*Z**2+s.Rational(19,8)*X**4+s.Rational(15,4)*X**2*Z**2+s.Rational(11,8)*Z**4
invr=1-X+X**2-X**3
ur=s.expand(s.diff(a,Z)*invr/4)
uz=s.expand(-s.diff(a,X)*invr/4)
c=2*e**2
sqrtH=2*s.sqrt(2)*e*(1-s.Rational(21,8)*e**2)
Hp=4-21*c+s.Rational(117,32)*c**2
up=sqrtH*invr/4
ratio=s.sqrt(2)/e*(1-s.Rational(63,8)*e**2) # H'(c)/sqrt(H(c))
upr=ratio*s.diff(a,X)*invr/8-sqrtH*(1-2*X+3*X**2)/4
upz=ratio*s.diff(a,Z)*invr/8
sub={X:xx,Z:zz}
Lraw=s.Matrix([[s.diff(ur,X),-up*invr,s.diff(ur,Z)],
               [upr,ur*invr,upz],
               [s.diff(uz,X),0,s.diff(uz,Z)]])
# The H'(c)/sqrt(H(c)) factor in upr is O(e^-1).  Consequently a spatial
# cubic in Lraw contributes at order e^2 after X,Z=O(e) are substituted.
# mp(...,4), rather than mp(...,3), is therefore required for C2.
L=Lraw.applyfunc(lambda z:ep(mp(z,4).subs(sub),3))
L_spatial_quadratic=Lraw.applyfunc(lambda z:ep(mp(z,3).subs(sub),3))

# k is the returned sigma covector: grad sigma || (z_I,-rho_I).
ki=s.Matrix([s.diff(zz,e),0,-s.diff(xx,e)])
n2=ep((ki.T*ki)[0],3)
n1=n2.coeff(e,1); n22=n2.coeff(e,2)
invkn=1-e*n1/2+e**2*(-n22/2+3*n1**2/8)
kh=ki.applyfunc(lambda z:ep(z*invkn,3))
p=s.Matrix([-kh[2],0,kh[0]])
E=s.Matrix.hstack(p,s.Matrix([0,1,0]))
phidot=ep(up.subs(sub)*(1-xx+xx**2),3)
Oc=s.Matrix([[0,-phidot,0],[phidot,0,0],[0,0,0]])
Edot=E.diff(q) # Omega_1=1+O(e^4) on the flat cutoff plateau

def pulled_coefficient(Ljet):
    M=-Ljet+2*kh*(kh.T*Ljet)
    return (E.T*(M-Oc)*E-E.T*Edot).applyfunc(lambda z:ep(z,3))

C=pulled_coefficient(L)
C_spatial_quadratic=pulled_coefficient(L_spatial_quadratic)
C0=C.applyfunc(lambda z:s.simplify(z.coeff(e,0)))
C1=C.applyfunc(lambda z:s.trigsimp(z.coeff(e,1),method='fu'))
C2=C.applyfunc(lambda z:s.trigsimp(z.coeff(e,2),method='fu'))
C2_spatial_quadratic=C_spatial_quadratic.applyfunc(
    lambda z:s.trigsimp(z.coeff(e,2),method='fu'))
C210_increment=s.sqrt(2)*(s.Rational(3,2)*s.cos(q)**4
                          -s.Rational(13,8)*s.cos(q)**2
                          -s.Rational(9,16))
assert s.trigsimp(C2[1,0]-C2_spatial_quadratic[1,0]
                  -C210_increment,method='fu')==0
print('W2 =',W2)
print('g1 =',d1)
print('g2 =',d2)
print('C0 =',C0)
print('C2 =',C2)

U=lambda x:s.eye(2)+x*C0
T=2*s.pi
C1s=s.Matrix([[s.cos(q)-s.Rational(3,2)*s.cos(3*q),s.sqrt(2)*s.sin(q)],
              [s.sqrt(2)*s.sin(q)/8+3*s.sqrt(2)*s.sin(3*q)/8,-s.cos(q)]])
assert all(s.simplify(s.expand_trig(x))==0 for x in C1-C1s)
print('C1 =',C1s)
print('C2[1,0]_cubic_spatial_increment =',C210_increment)
single=s.integrate(C2[0,0]+C2[1,1]-T/s.sqrt(2)*C2[0,1],(q,0,T))
u=s.symbols('u', real=True)
C1u=C1s.subs(q,u)
double=s.integrate(s.integrate(s.expand_trig(s.expand(s.trace(U(T-q)*C1s*U(q-u)*C1u*U(u)))),(u,0,q)),(q,0,T))
total=s.trigsimp(single+double)
assert s.simplify(single-9*s.pi**2)==0
assert s.simplify(double-2*s.pi**2)==0
assert s.simplify(total-11*s.pi**2)==0
print('single_C2_trace =',s.trigsimp(single))
print('double_C1_trace =',s.trigsimp(double))
print('epsilon2_trace =',total)
print('I_trace =',s.simplify(2*total))
print('I_discriminant =',s.simplify(4*2*total))
print('ALL 8 SOURCE-RECURRENCE CHECKS PASSED')

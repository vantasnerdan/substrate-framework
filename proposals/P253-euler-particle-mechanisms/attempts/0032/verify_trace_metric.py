#!/usr/bin/env python3
import sympy as sp

def p(msg):
    print("PASS " + msg)

l11,l12,l13,l21,l22,l23,l31,l32=sp.symbols('l11 l12 l13 l21 l22 l23 l31 l32', real=True)
L=sp.Matrix([[l11,l12,l13],[l21,l22,l23],[l31,l32,-l11-l22]])
k1,k2,k3,A1,A2=sp.symbols('k1 k2 k3 A1 A2', nonzero=True, real=True)
k=sp.Matrix([k1,k2,k3])
A=sp.Matrix([A1,A2,-(k1*A1+k2*A2)/k3])
kd=-L.T*k
M=-L+2*k*(k.T*L)/(k.dot(k))
Ad=M*A
omega=sp.Matrix([l32-l23,l13-l31,l21-l12])
c=k.cross(A)
cd=kd.cross(A)+k.cross(Ad)
assert sp.simplify(cd-L*c-omega.dot(k)*A)==sp.zeros(3,1)
assert sp.simplify(kd.dot(A)+k.dot(Ad))==0
p('full-pressure c identity and transversality')

v1=sp.Matrix([k2,-k1,0])
v2=k.cross(v1)
area=sp.expand((M*v1).cross(v2).dot(k)+v1.cross(M*v2).dot(k)+(v1.cross(v2)).dot(kd))
assert sp.simplify(area)==0
p('polarization area conservation implies det(M_return)=1')

rhs_u=-(-k)+2*k*(k.dot(-k))/k.dot(k)
assert sp.simplify(rhs_u+k)==sp.zeros(3,1)
p('pressure-normal neutral amplitude')

s=sp.symbols('s', nonzero=True, real=True)
Mh=sp.diag(sp.exp(s),sp.exp(-s))
Mp=sp.Matrix([[1,1],[0,1]])
Mm=-sp.eye(2)
assert Mh.det()==1 and Mp.det()==1 and Mm.det()==1
assert sp.simplify(Mh.trace()**2-4)!=0
assert Mp.trace()**2-4==0 and (Mp-sp.eye(2)).rank()==1
assert Mm.trace()**2-4==0 and (Mm+sp.eye(2)).rank()==0
p('hyperbolic, nontrivial Jordan, and -I cases are separated')

h11,h12,h22=sp.symbols('h11 h12 h22', real=True)
H=sp.Matrix([[h11,h12],[h12,h22]])
assert sp.simplify((Mh.T*H*Mh-H)[0,0])==h11*(sp.exp(2*s)-1)
assert sp.simplify((Mp.T*H*Mp-H)[0,1])==h11
theta=sp.symbols('theta', real=True)
R=sp.Matrix([[sp.cos(theta),-sp.sin(theta)],[sp.sin(theta),sp.cos(theta)]])
assert sp.simplify(R.T*R-sp.eye(2))==sp.zeros(2)
p('positive invariant metric exists for rotation, not hyperbolic/Jordan')

print('ALL 5 EXACT CHECKS PASSED')

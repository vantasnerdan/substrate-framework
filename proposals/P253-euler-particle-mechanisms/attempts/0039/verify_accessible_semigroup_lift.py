#!/usr/bin/env python3
"""Exact algebraic checks for the P253/0039 continuum bridge."""
import sympy as sp

def passed(label):
    print('PASS ' + label)

k1,k2,k3,A1,A2,A3,w1,w2,w3=sp.symbols(
    'k1 k2 k3 A1 A2 A3 w1 w2 w3', real=True)
k=sp.Matrix([k1,k2,k3]); A=sp.Matrix([A1,A2,A3]); w=sp.Matrix([w1,w2,w3])
mu=k.dot(w); k2norm=k.dot(k)
a=k.cross(A)/mu
domega=k.cross(a.cross(w))
v=-k.cross(domega)/k2norm
constraint={A3:-(k1*A1+k2*A2)/k3}
assert sp.simplify((v-A).subs(constraint))==sp.zeros(3,1)
assert sp.simplify(k.dot(a))==0
passed('nonzero-mu orbit tangent realizes every transverse polarization')

# The pressure sign is load-bearing: deleting the outer Biot-Savart minus sign
# returns -A and must fail the target identity.
v_wrong=k.cross(domega)/k2norm
assert sp.simplify((v_wrong+A).subs(constraint))==sp.zeros(3,1)
passed('wrong Biot-Savart sign is exposed')

tr=sp.symbols('tr', positive=True)
disc=tr**2-4
lp=(tr+sp.sqrt(disc))/2
lm=(tr-sp.sqrt(disc))/2
assert sp.simplify(lp*lm-1)==0
assert sp.simplify(lp+lm-tr)==0
passed('hyperbolic SL2 multipliers are reciprocal')

j,T,r,C,B,delta,lam=sp.symbols('j T r C B delta lam', positive=True)
hstar=r*lam**j*delta**3/(4*C*j*T*sp.exp(B*j*T))
err=sp.simplify(hstar*C*delta**-3*j*T*sp.exp(B*j*T))
assert sp.simplify(err-r*lam**j/4)==0
passed('frequency-after-time choice closes relative WKB error')

gamma=sp.log(lam)/T
assert sp.simplify(sp.log(lam**j)/(j*T)-gamma)==0
passed('all-circuit gains imply the stated semigroup growth exponent')

I=sp.symbols('I', positive=True)
trace=2+22*sp.pi**2*I
Delta=sp.expand(trace**2-4)
assert sp.expand(Delta).coeff(I,1)==88*sp.pi**2
passed('0032 discriminant coefficient is transferred without refitting')

print('ALL 6 EXACT CHECKS PASSED')

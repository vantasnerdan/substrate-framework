#!/usr/bin/env python3
"""Exact pullback checks for the activated resonant coefficient."""
import sympy as sp

rho, rs, rI, es, eI, zs, zI = sp.symbols('rho rs rI es eI zs zI', nonzero=True, real=True)
J = sp.Matrix([[rs, 0, rI], [rho*es, rho, rho*eI], [zs, 0, zI]])
g = sp.simplify(J.T*J)
assert sp.simplify(g[0,0]-(rs**2+rho**2*es**2+zs**2)) == 0
assert sp.simplify(g[0,1]-rho**2*es) == 0
assert sp.simplify(g[0,2]-(rs*rI+rho**2*es*eI+zs*zI)) == 0
assert sp.simplify(g[1,1]-rho**2) == 0
assert sp.simplify(g[1,2]-rho**2*eI) == 0
assert sp.simplify(g[2,2]-(rI**2+rho**2*eI**2+zI**2)) == 0
assert sp.simplify(J.det()-rho*(rs*zI-zs*rI)) == 0
print('PASS Baldi metric entries and det(J)')

n, m = sp.symbols('n m', nonzero=True, real=True)
kappa = sp.Matrix([n,-m,0])
e = sp.Matrix([m,n,0])
f = sp.Matrix([0,0,1])
g11,g12,g13,g22,g23,g33 = sp.symbols('g11 g12 g13 g22 g23 g33', real=True)
G = sp.Matrix([[g11,g12,g13],[g12,g22,g23],[g13,g23,g33]])
h = G.inv()*kappa
den = sp.simplify(kappa.dot(h))
def R(v):
    return G.inv()*(v.cross(h))/den
Re = R(e)
Rf = R(f)
assert sp.simplify(kappa.dot(Re)) == 0
assert sp.simplify(kappa.dot(Rf)) == 0
q = n**2+m**2
Eee = sp.simplify(e.dot(Re)/q)
Eef = sp.simplify(e.dot(Rf)/q)
Efe = sp.simplify(f.dot(Re))
Eff = sp.simplify(f.dot(Rf))
assert sp.simplify(R(e)-e*Eee-f*Efe) == sp.zeros(3,1)
assert sp.simplify(R(f)-e*Eef-f*Eff) == sp.zeros(3,1)
print('PASS metric cross operator closes on kappa-perp scalar basis')

a, b = sp.symbols('a b', nonzero=True, real=True)
B = sp.Matrix([[0,0,a],[0,0,b],[0,0,0]])
assert sp.simplify(B*f-(a/m)*e).subs(b,n*a/m) == sp.zeros(3,1)
print('PASS resonance reduces D(Omega) to the scalar upper entry a/m')

print('ALL 3 EXACT CHECKS PASSED')

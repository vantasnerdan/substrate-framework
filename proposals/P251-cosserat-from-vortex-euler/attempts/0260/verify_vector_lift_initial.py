"""Exact coordinate exposure of0260; the physical chart is a hypothesis."""

import sympy as s

I = s.symbols("I", real=True)
N, m, n, lam = s.symbols("N m n lam", real=True, nonzero=True)
J, A, O1, O2, b1, b2 = [s.Function(name)(I) for name in
                         ("J", "A", "O1", "O2", "b1", "b2")]
nu = m*O1+n*O2
f = A*s.exp(s.I*N*I)

# Single harmonic in an orthogonal-I coordinate example; no Euclidean
# Beltrami existence is inferred from arbitrary metric coefficient inputs.
alpha1, alpha2 = J*b2*f, -J*b1*f
wI = (s.I*m*alpha2-s.I*n*alpha1)/J
w1 = -s.diff(alpha2, I)/J
w2 = s.diff(alpha1, I)/J
xiI = wI/(s.I*lam*nu)
xi1 = w1/(s.I*lam*nu)+xiI*s.diff(O1, I)/(s.I*nu)
xi2 = w2/(s.I*lam*nu)+xiI*s.diff(O2, I)/(s.I*nu)
assert s.simplify(lam*s.I*nu*xiI-wI) == 0
assert s.simplify(lam*(s.I*nu*xi1-xiI*s.diff(O1, I))-w1) == 0
assert s.simplify(lam*(s.I*nu*xi2-xiI*s.diff(O2, I))-w2) == 0
divergence = (s.diff(J*xiI, I)+J*s.I*m*xi1+J*s.I*n*xi2)/J
assert s.simplify(divergence) == 0
wrong_xi1 = w1/(s.I*lam*nu)
wrong_shear_row = s.simplify(
    lam*(s.I*nu*wrong_xi1-xiI*s.diff(O1, I))-w1
)
assert wrong_shear_row != 0
normal = s.simplify(xiI/f)
assert s.simplify(normal+(m*b1+n*b2)/(lam*nu)) == 0
assert s.simplify(normal.subs({b1: O1, b2: O2})+1/lam) == 0
assert s.simplify(normal.subs({b1: -n, b2: m})) == 0

# General coordinate triple product retains the physical Jacobian J.
p, q, o1, o2, j = s.symbols("p q o1 o2 j", real=True)
freq = m*o1+n*o2
omega = s.Matrix([0, lam*o1, lam*o2])
source_u = s.Matrix([-p/lam, N*o1*p/(lam*freq), N*o2*p/(lam*freq)])
source_perp = s.Matrix([0, -N*n*q/(lam*freq), N*m*q/(lam*freq)])
triple = s.factor(j*omega.dot(source_u.cross(source_perp)))
assert s.simplify(triple-j*N*p*q/lam) == 0

print("exact normal coefficient:", normal)
print("divergence:", s.simplify(divergence))
print("omitted shear residual:", wrong_shear_row)
print("distinct-polarization KKS integrand:", triple)
print("PASS: vector bracket, shear, divergence, two normal rows and KKS sign")

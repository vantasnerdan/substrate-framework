#!/usr/bin/env python3
"""Exact exposing checks for P253/0025's moving-fibre reduction."""

import sympy as sp


def ok(label: str) -> None:
    print(f"PASS {label}")


I, t = sp.symbols("I t", positive=True, real=True)
a, b = sp.symbols("a b", nonzero=True, real=True)
n, m = sp.symbols("n m", nonzero=True, real=True)

# Baldi source-series consequence: a nonconstant derivative-frequency ratio.
cK = sp.Rational(1065, 1024)
K = I + cK * I**3
Omega1 = sp.diff(K, I)
Omega2 = sp.sqrt(I) * (1 + sp.Rational(7, 4) * I) * Omega1
aprime = sp.diff(Omega1, I)
bprime = sp.diff(Omega2, I)
D = sp.cancel(bprime / aprime)
assert sp.limit(D * I ** sp.Rational(3, 2), I, 0, dir="+") == sp.Rational(256, 3195)
assert sp.limit(sp.diff(D, I) * I ** sp.Rational(5, 2), I, 0, dir="+") == -sp.Rational(128, 1065)
ok("Baldi-series derivative ratio and decreasing leading slope")

# Exact action-angle deformation and covector resonance.
F = sp.Matrix([[1, 0, t * a], [0, 1, t * b], [0, 0, 1]])
assert F.det() == 1
k0 = sp.Matrix([n, -m, 0])
kt = sp.simplify(F.inv().T * k0)
expected_kt = sp.Matrix([n, -m, -t * (a * n - b * m)])
assert sp.simplify(kt - expected_kt) == sp.zeros(3, 1)
assert sp.simplify(kt.subs(n, m * b / a) - k0.subs(n, m * b / a)) == sp.zeros(3, 1)
ok("cotangent shear and exact derivative resonance")

# Six-dimensional cotangent lift is genuinely canonical symplectic.
Z = sp.zeros(3)
J6 = Z.row_join(sp.eye(3)).col_join((-sp.eye(3)).row_join(Z))
C = F.row_join(Z).col_join(Z.row_join(F.inv().T))
assert sp.simplify(C.T * J6 * C - J6) == sp.zeros(6)
ok("moving base/covector map preserves the canonical ray symplectic form")

# Full Kelvin--Leray algebra. Impose k.A=0 by eliminating A3.
l11, l12, l13, l21, l22, l23, l31, l32 = sp.symbols(
    "l11 l12 l13 l21 l22 l23 l31 l32", real=True
)
L = sp.Matrix(
    [[l11, l12, l13], [l21, l22, l23], [l31, l32, -l11 - l22]]
)
k1, k2, k3, A1, A2 = sp.symbols("k1 k2 k3 A1 A2", nonzero=True, real=True)
kv = sp.Matrix([k1, k2, k3])
Av = sp.Matrix([A1, A2, -(k1 * A1 + k2 * A2) / k3])
kdot = -L.T * kv
mu_pressure = kv.dot(L * Av) / kv.dot(kv)
Adot = -L * Av + 2 * kv * mu_pressure
assert sp.simplify(kdot.dot(Av) + kv.dot(Adot)) == 0
bad_constraint_dot = sp.simplify(kdot.dot(Av) + kv.dot(-L * Av))
assert sp.simplify(bad_constraint_dot + 2 * kv.dot(L * Av)) == 0
assert bad_constraint_dot != 0
ok("full Leray coefficient preserves transversality and omission is detected")

# The failed shortcut is exposed: c=kxA has an unavoidable vorticity term.
omega = sp.Matrix([l32 - l23, l13 - l31, l21 - l12])
cvec = kv.cross(Av)
cdot = sp.simplify(kdot.cross(Av) + kv.cross(Adot))
correct_cdot = sp.simplify(L * cvec + omega.dot(kv) * Av)
assert sp.simplify(cdot - correct_cdot) == sp.zeros(3, 1)
shortcut_residual = sp.simplify(cdot - L * cvec)
assert sp.simplify(shortcut_residual - omega.dot(kv) * Av) == sp.zeros(3, 1)
assert shortcut_residual != sp.zeros(3, 1)
ok("Cauchy-vector shortcut fails by exactly (vorticity dot k)A")

# Polarization-area invariant, the source of the SL(2,R) returned cocycle.
v1 = sp.Matrix([k2, -k1, 0])
v2 = kv.cross(v1)


def amp_dot(v: sp.Matrix) -> sp.Matrix:
    return -L * v + 2 * kv * (kv.dot(L * v)) / kv.dot(kv)


area_dot = sp.expand(
    (amp_dot(v1).cross(v2) + v1.cross(amp_dot(v2))).dot(kv)
    + (v1.cross(v2)).dot(kdot)
)
assert sp.simplify(area_dot) == 0
ok("Kelvin polarization triple area is conserved")

# For k=grad P, steady Euler gives L*u=-k and k.u=0. Thus A=u is exact.
q2 = kv.dot(kv)
rhs_from_Lu = -(-kv) + 2 * kv * kv.dot(-kv) / q2
assert sp.simplify(rhs_from_Lu - (-kv)) == sp.zeros(3, 1)
ok("pressure-gradient covector admits the exact neutral amplitude A=u")

# Check the normalized complement in a local orthonormal frame.  Take
# u=e1, k=e3, Lu=-k and grad|u|^2=2L^T u parallel to k.
alpha, ell22, ell23, ell32 = sp.symbols("alpha ell22 ell23 ell32", real=True)
Lnormal = sp.Matrix([[0, 0, alpha], [0, ell22, ell23], [-1, ell32, -ell22]])
unormal = sp.Matrix([1, 0, 0])
knormal = sp.Matrix([0, 0, 1])
Bnormal = knormal.cross(unormal)
kdot_normal = -Lnormal.T * knormal
udot_normal = Lnormal * unormal
h_dot = kdot_normal.cross(unormal) + knormal.cross(udot_normal)
q_dot = 2 * knormal.dot(kdot_normal)
s_dot = 2 * unormal.dot(udot_normal)
Bdot_normal = sp.simplify(h_dot - Bnormal * (q_dot + s_dot))
MB_normal = sp.simplify(
    -Lnormal * Bnormal
    + 2 * knormal * knormal.dot(Lnormal * Bnormal)
)
assert sp.simplify(Bdot_normal - MB_normal) == sp.zeros(3, 1)
ok("localizability makes the normalized pressure-phase complement exact")

# Dynamically accessible leading-vorticity triple product retains both terms.
d1, d2, d3, w1, w2, w3 = sp.symbols("d1 d2 d3 w1 w2 w3", real=True)
dv = sp.Matrix([d1, d2, d3])
wv = sp.Matrix([w1, w2, w3])
triple = sp.simplify(kv.cross(dv.cross(wv)))
expected = sp.simplify(dv * kv.dot(wv) - wv * kv.dot(dv))
assert sp.simplify(triple - expected) == sp.zeros(3, 1)
ok("coadjoint displacement symbol retains both accessibility terms")

print("ALL 9 EXACT CHECKS PASSED")

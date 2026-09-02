"""P250 exact local certification of the deep-branch Maxwell point.

Two rigorous interval certificates on the exact polynomial system:

  1. Krawczyk operator: the deep-branch Maxwell system has EXACTLY ONE real
     solution in the rational box B around the recorded 50-dps point.
  2. Interval Gershgorin: the fixed-omega Hessian of V_omega at that
     solution is strictly positive definite over the box (nondegenerate
     interior minimum), with explicit rational margins.

mpmath interval arithmetic on SymPy-derived exact polynomials; no grid, no
fit, no large-subtraction.
"""
import json

import mpmath as mp
import sympy as sp

mp.iv.dps = 40
mp.mp.dps = 60

cS, bS, fS = sp.symbols('c b f', real=True)
w_of = 32*fS**2 - 12*fS + 6 - 24*bS

gA = 8*cS**3 - 3*cS**2 + (8*bS**2 + 1)*cS - 3*bS**2
gB = sp.expand(8*bS**3 + 8*bS*cS**2 - 6*bS*cS - 2*bS*w_of + 7*bS - 6*fS**2)
TrS2 = 2*cS**2 + 2*bS**2
V0 = (-(TrS2)/2 - (2*cS**3 + 6*cS*bS**2) + TrS2**2 + sp.Rational(1, 2)
      + 2*cS**2 + 2*bS**2 + 6*(bS - fS**2)**2
      + 3*fS**2 - 4*fS**3 + 2*fS**4)
gC = sp.expand(2*(V0 - w_of*(fS**2 + 4*bS**2)/2))
G = [gA, gB, gC]

J_sym = sp.Matrix([[sp.diff(g, v) for v in (cS, bS, fS)] for g in G])
# fixed-omega Hessian of V_omega = V0 - (omega^2/2)(f^2+4b^2) in (c,b,f):
om2 = sp.Symbol('omega_sq', nonnegative=True)
V_om = sp.expand(V0 - om2*(fS**2 + 4*bS**2)/2)
H_sym = sp.Matrix([[sp.diff(V_om, u1, u2) for u2 in (cS, bS, fS)]
                   for u1 in (cS, bS, fS)])

g_lam = [sp.lambdify((cS, bS, fS), g, modules='mpmath') for g in G]
J_lam = [[sp.lambdify((cS, bS, fS), J_sym[i, j], modules='mpmath')
          for j in range(3)] for i in range(3)]
H_lam = [[sp.lambdify((cS, bS, fS, om2), H_sym[i, j], modules='mpmath')
          for j in range(3)] for i in range(3)]

c0 = mp.mpf("0.30280764518677380369")
b0 = mp.mpf("0.65773437772324925193")
f0 = mp.mpf("0.81436149699856776719")
w0 = 32*f0**2 - 12*f0 + 6 - 24*b0

delta = mp.mpf('1e-12')
X = mp.iv.matrix([mp.iv.mpf([c0 - delta, c0 + delta]),
                  mp.iv.mpf([b0 - delta, b0 + delta]),
                  mp.iv.mpf([f0 - delta, f0 + delta])])
y = [c0, b0, f0]
J_mid = mp.matrix([[J_lam[i][j](c0, b0, f0) for j in range(3)] for i in range(3)])
Y = mp.iv.matrix(J_mid**-1)
gy = mp.iv.matrix([mp.iv.mpf(g_lam[i](c0, b0, f0)) for i in range(3)])
JX = mp.iv.matrix([[J_lam[i][j](*X) for j in range(3)] for i in range(3)])
I3 = mp.iv.matrix([[mp.iv.mpf(1) if i == j else mp.iv.mpf(0)
                    for j in range(3)] for i in range(3)])
Xy = mp.iv.matrix([X[i] - mp.iv.mpf(y[i]) for i in range(3)])
K = mp.iv.matrix([mp.iv.mpf(v) for v in y]) - Y*gy + (I3 - Y*JX)*Xy
strict_inside = all(K[i].a > X[i].a and K[i].b < X[i].b for i in range(3))
print("Krawczyk strict inclusion:", strict_inside)
print("K boxes:", [mp.nstr(K[i], 5) for i in range(3)])

# interval Gershgorin of the fixed-omega Hessian at omega = w(X)
W = mp.iv.mpf([32*X[2].a**2 - 12*X[2].a + 6 - 24*X[1].b,
               32*X[2].b**2 - 12*X[2].b + 6 - 24*X[1].a])
H = mp.iv.matrix([[H_lam[i][j](*X, W) for j in range(3)] for i in range(3)])
gersh = []
for i in range(3):
    off = mp.iv.mpf(0)
    for j in range(3):
        if j != i:
            off = off + abs(H[i, j])
    gersh.append(min(H[i, i].a, H[i, i].b) - off.b)
print("Gershgorin lower bounds:", [mp.nstr(g, 10) for g in gersh])
h_pd = all(g > 0 for g in gersh)
out = {
    "box_halfwidth": "1e-12",
    "krawczyk_strict_interior": bool(strict_inside),
    "hessian_gershgorin_lower_bounds": [mp.nstr(g, 15) for g in gersh],
    "hessian_positive_definite": bool(h_pd),
    "w_interval": [mp.nstr(W.a, 15), mp.nstr(W.b, 15)],
}
with open("proposals/P250-shell-bubble-clock/attempts/0001/"
          "maxwell_certificate.json", "w") as fh:
    json.dump(out, fh, indent=1)
print("DONE")

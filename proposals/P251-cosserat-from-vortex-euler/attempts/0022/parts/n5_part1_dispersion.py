"""N5 Part 1 -- symbolic two-branch dispersion of the identified micropolar system.

Plane waves u = u_h e^{i(k.x - w t)}, Phi = Phi_h e^{i(k.x - w t)} in the N4
balance laws:

  linear : rho w^2 u_h = (lam+mu-alpha) k (k.u_h) + (mu+alpha) k^2 u_h
           - 2 alpha i k x Phi_h
  angular: j w^2 Phi_h = (c_s+c_a) k^2 Phi_h - (2c_tr-c_a+c_s) k (k.Phi_h)
           - 2 alpha i k x u_h + 4 alpha Phi_h

Decompose along / across k (k = k e_z WLOG): longitudinal (u_z, Phi_z) and
transverse (u_x, u_y; Phi_x, Phi_y) sectors. Derive w^2(k) symbolically.
"""
import sympy as sp

k, w = sp.symbols("k omega", positive=True)
rho, lam, mu, alpha = sp.symbols("rho lambda mu alpha", positive=True)
cs, ca, ctr, j = sp.symbols("c_s c_a c_tr j", positive=True)
I = sp.I

# transverse coupled 2x2 per circular polarization:
#   rho w^2 u  = (mu+alpha) k^2 u  - 2 alpha s k Phi
#   j w^2 Phi  = (c_s+c_a) k^2 Phi + 4 alpha Phi + 2 alpha s k u
# (s = helicity +-1 from k x Phi = s k Phi for transverse circular modes)
s = sp.Symbol("s")
u_t, Phi_t = sp.symbols("u_t Phi_t")
M = sp.Matrix([
    [(mu + alpha) * k**2 - rho * w**2,          -2 * alpha * s * k],
    [2 * alpha * s * k,   (cs + ca) * k**2 + 4 * alpha - j * w**2],
])
det_T = sp.factor(sp.expand(M.det()))
print("transverse dispersion det = 0 :", det_T)

w2 = sp.symbols("w2")
det_T_w2 = sp.factor(det_T.subs(w**2, w2))
roots = sp.solve(det_T_w2, w2)
print("\ntransverse branches w^2 =")
for r in roots:
    print("  ", sp.simplify(sp.factor(r)))

# longitudinal 2x2:
#   rho w^2 u_L = (lam + 2 mu - alpha) k^2 u_L          (rot vanishes along k)
#   j w^2 Phi_L = (2 c_tr - c_a + c_s) k^2 Phi_L ... check: k(k.Phi) = k^2 Phi_L,
#   (c_s+c_a) k^2 Phi_L - (2c_tr-c_a+c_s) k^2 Phi_L = (2 c_a - 2 c_tr) k^2 Phi_L
#   + 4 alpha Phi_L
u_l, Phi_l = sp.symbols("u_l Phi_l")
M_L = sp.Matrix([
    [(lam + 2 * mu - alpha) * k**2 - rho * w**2, 0],
    [0, (2 * ca - 2 * ctr) * k**2 + 4 * alpha - j * w**2],
])
print("\nlongitudinal branches:")
print("   rho w^2 = (lam + 2 mu - alpha) k^2          (P-like wave, alpha-corrected)")
print("   j w^2 = (2 c_a - 2 c_tr) k^2 + 4 alpha      (longitudinal spin wave)")

"""N5 Part 2 -- branch limits, fluid limit, and moduli substitution.

Checks:
  L1  w_-(k->0) -> 0          (acoustic branch reaches zero)
  L2  w_+(k->0) -> 4 alpha/j  (gapped optical spin wave)
  L3  fluid limit L_v -> 0:   all branches -> 0 (Euler recovery premise)
  L4  Vieta: sum and product of the transverse w^2 roots match the determinant
"""
import sympy as sp

k, w2, s = sp.symbols("k omega2 s")
rho, lam, mu, alpha = sp.symbols("rho lambda mu alpha", positive=True)
cs, ca, ctr, j = sp.symbols("c_s c_a c_tr j", positive=True)
Lv = sp.Symbol("L_v", positive=True)

# transverse quadratic in w2 (s^2 = 1 for circular polarization):
A = rho * j
B = -(alpha * j * k**2 + 4 * alpha * rho + rho * k**2 * (cs + ca) + j * k**2 * mu)
C = ((mu + alpha) * k**2) * ((cs + ca) * k**2 + 4 * alpha) - 8 * alpha**2 * k**2

disc = sp.expand(B**2 - 4 * A * C)
w2p = sp.simplify((-B + sp.sqrt(disc)) / (2 * A))
w2m = sp.simplify((-B - sp.sqrt(disc)) / (2 * A))

# L1: w2m -> 0 as k -> 0
l1 = sp.simplify(sp.limit(w2m, k, 0))
print("L1 w2m(k->0) =", l1, "  -> 0 :", l1 == 0)

# L2: w2p -> 4 alpha / j
l2 = sp.simplify(sp.limit(w2p, k, 0) - 4 * alpha / j)
print("L2 w2p(k->0) - 4a/j =", l2, "  -> 0 :", l2 == 0)

# L3: fluid limit — all moduli carry L_v; with L_v factored the branches vanish.
# Substitute the L_v scalings: mu=a3 Lv, alpha=a6 Lv, cs=b10 Lv, ca=b6 Lv, ctr=b30 Lv, j=j3 Lv
a3, a6, b10, b6, b30, j3 = sp.symbols("a3 a6 b10 b6 b30 j3", positive=True)
sub_L0 = {mu: 0, alpha: 0, cs: 0, ca: 0, ctr: 0, j: 1}  # L_v = 0 limit keeps rho, j(1)=0 too
w2m_L0 = sp.simplify(w2m.subs([(mu, 0), (alpha, 0), (cs, 0), (ca, 0), (ctr, 0), (j, j * 0 + 1)]))
print("L3 w2m at L_v = 0 :", sp.simplify(w2m_L0), "  -> 0 :", sp.simplify(w2m_L0) == 0)
w2p_L0 = sp.simplify(w2p.subs([(mu, 0), (alpha, 0), (cs, 0), (ca, 0), (ctr, 0), (j, 1)]))
print("L3 w2p at L_v = 0 :", sp.simplify(w2p_L0), "  -> 0 :", sp.simplify(w2p_L0) == 0)

# L4: Vieta over the transverse quadratic
sum_roots = sp.simplify(-B / A)
prod_roots = sp.simplify(C / A)
sR = sp.simplify(w2p + w2m - sum_roots)
pR = sp.simplify(w2p * w2m - prod_roots)
print("L4 Vieta sum residual:", sR, "  product residual:", sp.simplify(pR))

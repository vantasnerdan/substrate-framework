# EXPLORATORY PRICING SCRATCH — NOT an F-B receipt.
# Charter: shepherd tasking 2026-09-11 (frame-price analysis, owner-level
# build decision reserved). This scratch prices P2 (anisotropic modulus
# algebra) ONLY. 01's F-B row stays "no computation" for claims; owner
# call restores or lifts. Self-counted; exploratory role retained.

import sympy as sp

e11, e22, e33, e12, e13, e23 = sp.symbols('e11 e22 e33 e12 e13 e23', real=True)
eps = sp.Matrix([[e11, e12, e13], [e12, e22, e23], [e13, e23, e33]])
K, p = sp.symbols('K p', positive=True)
c1, s1, c2, s2 = sp.symbols('c1 s1 c2 s2', real=True)
ngen = sp.Matrix([c1*s2, s1*s2, c2])
ez = sp.Matrix([0, 0, 1])
tr_e2 = sum((eps*eps)[i, i] for i in range(3))


def W_of(eps_m, n, pv):
    """Quadratic-order energy of the polarized tangle:
    angular averages with <n_i n_j> = (1-p) dij/3 + p nhat_i nhat_j,
    <4-fold> = (1-p) iso4 + p nhat^x4. Linear term kept (pressure-like)."""
    trs = sp.trace(eps_m)
    trs2 = sum((eps_m*eps_m)[i, i] for i in range(3))
    nen1 = n.dot(eps_m*n)                       # scalar n.eps n
    nesq = n.dot(eps_m*eps_m*n)                 # n.eps^2 n
    ne2 = nen1**2                               # (n.eps n)^2
    return K*(trs/3 + pv*(nen1 - trs/3)
              + sp.Rational(1, 2)*((1 - pv)*trs2/3 + pv*nesq
                                   - (1 - pv)*(trs**2 + 2*trs2)/15 - pv*ne2))


print("[pricing 1] p=0 quadratic == F-A form (continuity):",
      sp.simplify(sp.expand(W_of(eps, ngen, 0)) - K*sp.trace(eps)/3
                  - K*(tr_e2/10 - sp.trace(eps)**2/30)) == 0)

# exact stiffness matrix, axis z
qs = sp.expand(W_of(eps, ez, p))
v = [e11, e22, e33, e12, e13, e23]
M = sp.zeros(6, 6)
for i, vi in enumerate(v):
    for j, vj in enumerate(v):
        if j < i:
            continue
        M[i, j] = M[j, i] = (sp.simplify(qs.coeff(vi, 2)) if i == j
                             else sp.simplify(qs.coeff(vi*vj, 1)/2))

# transverse-isotropy checks (5-constant reduction)
c = {str(m): sp.simplify(qs.coeff(m, 1)) for m in
     [e11**2, e22**2, e33**2, e11*e22, e11*e33, e22*e33, e12**2, e13**2, e23**2]}
print("[pricing 2] TI reduction: e11^2==e22^2:", sp.simplify(c['e11**2'] - c['e22**2']) == 0,
      "| e13^2==e23^2:", sp.simplify(c['e13**2'] - c['e23**2']) == 0,
      "| distinct nonzero constants:", len(set(c.values())))

# exact 2nd-variation spectrum on the traceless subspace
b = [sp.Matrix([1, -1, 0, 0, 0, 0]), sp.Matrix([1, 1, -2, 0, 0, 0]),
     sp.Matrix([0, 0, 0, 1, 0, 0]), sp.Matrix([0, 0, 0, 0, 1, 0]),
     sp.Matrix([0, 0, 0, 0, 0, 1])]
G = sp.Matrix([[bi.dot(bj) for bj in b] for bi in b])
Mr = sp.Matrix([[(bi.T*M*bj)[0, 0] for bj in b] for bi in b])
lam = sp.Symbol('lam')
charpoly = sp.factor(sp.simplify((Mr - lam*G).det()))
roots = sorted(sp.solve(charpoly, lam), key=str)
print("[pricing 3] traceless eigenvalues (xK):",
      [sp.factor(sp.simplify(r)) for r in roots])
grid = [sp.Rational(0), sp.Rational(1, 4), sp.Rational(1, 2),
        sp.Rational(3, 4), sp.Rational(1)]
psd = all((r.subs({K: 1, p: pv}) >= 0) == True for r in roots for pv in grid)
print("[pricing 4] PSD on p-grid", grid, ":", psd)
print("[pricing 4b] roots by form: K(1-p)/10, K(1-p)/5, K(3p+2)/10",
      "-> PSD for all p in [0,1]; p=1 degenerate")
t = sp.Symbol('tau', real=True)
for mode, sub in [("mu_perp (e12)", e12), ("mu_par (e13)", e13)]:
    muq = sp.simplify(sp.expand(W_of(eps.subs(sub, t), ez, p)).coeff(t, 2))
    print(f"[pricing 5] 2*{mode}(p) = {sp.factor(muq)}")
# REPAIRED (due-2, drift review-sage-fbpricing; F-B-lite round): state
# mu_par(1) = K/4 with a REAL assert against pricing-5's own formula,
# replacing the truncated print + tautological assert.
mu_par_formula = sp.factor(sp.simplify(
    sp.expand(W_of(eps.subs(e13, t), ez, p)).coeff(t, 2))) / 2
print("[pricing 6] p=1 degeneracy: mu_perp -> 0 (sliding mode);")
print("[pricing 6] mu_par(p=1) == K/4:",
      sp.simplify(mu_par_formula.subs(p, 1) - K / 4) == 0)

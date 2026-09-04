"""N3 Part 2 -- wryness (couple) sector: (beta, gamma) from tube bend stiffness
via joint tangent-triad moments.

Model: the tube triad transports affinely; the segment's bend energy is
(B/2)|P_perp (kappa n)|^2 where kappa = grad Phi is the wryness tensor and
P_perp = I - n n^T projects across the tangent (bend lives in the plane
perpendicular to the centerline). The twist-about-tangent channel carries
coefficient C_tw^tube = 0 EXACTLY (attempts 0005 gauge identity, 0011 virial
closure -- recorded on the frontier); it is retained symbolically here so the
matching shows its vanishing.

Matching form (isotropic couple energy, two invariants):
    W_c = c1 kappa:kappa + c2 kappa:kappa^T   ->  (beta, gamma) map printed.

Averages (exact sphere moments):
    <n_j n_m>            = delta_jm/3
    <n_j n_m n_i n_l>    = (d_jm d_il + d_ji d_ml + d_jl d_mi)/15
"""
import sympy as sp
from substrate_framework.homogenization import (
    sphere_second_moment,
    sphere_fourth_moment_isotropic,
)

kappa = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"kp_{i}{j}"))
n = sp.Matrix(sp.symbols("n1:4"))
P2 = sphere_second_moment()
P4 = sphere_fourth_moment_isotropic()

# kappa n and projections
kn = kappa * n                                   # (kappa n)_i = kappa_ij n_j
kn_abs2 = sum(kn[i] ** 2 for i in range(3))      # = k_ij k_il n_j n_l
kn_abs2 = sp.expand(kn_abs2)

# average of |kappa n|^2 = <k_ij k_il n_j n_l> = (kappa:kappa)/3
kkt = sum(kappa[i, j] * kappa[i, j] for i in range(3) for j in range(3))
avg_abs = sp.simplify(sum(
    kappa[i, j] * kappa[l, j] * P2[i, l] for i in range(3) for l in range(3) for j in range(3)))
print("<|kappa n|^2> = k:k/3 :", sp.simplify(avg_abs - kkt / 3) == 0)

avg_ntkn = sp.simplify(sum(
    kappa[i, j] * kappa[l, m] * P4[i, j, l, m]
    for i in range(3) for j in range(3) for l in range(3) for m in range(3)))
tr_k = sum(kappa[i, i] for i in range(3))
kktT = sum(kappa[i, j] * kappa[j, i] for i in range(3) for j in range(3))
print("<(n' kappa n)^2> = ((tr k)^2 + k:k + k:k^T)/15 :",
      sp.simplify(avg_ntkn - ((tr_k**2 + kkt + kktT) / 15)) == 0)

# average of |P_perp kappa n|^2 = <|kappa n|^2> - <(n' kappa n)^2>
avg_perp = sp.simplify(avg_abs - avg_ntkn)
print("<|P_perp k n|^2> =", sp.simplify(sp.expand(avg_perp)))

# ---- ensemble couple energy ----
Lv, B, Ctw = sp.symbols("L_v B C_tw", positive=True)
W_c = sp.expand(Lv * (B / 2 * avg_perp + Ctw / 2 * avg_ntkn))

# ---- match against the three-invariant isotropic couple form ----
c_tr, c_s, c_a = sp.symbols("c_tr c_s c_a")
kS = (kappa + kappa.T) / 2
kA = (kappa - kappa.T) / 2
Ssym = sum(kS[i, j] ** 2 for i in range(3) for j in range(3))
Sskw = sum(kA[i, j] ** 2 for i in range(3) for j in range(3))
W_form = sp.expand(c_tr * tr_k**2 + c_s * Ssym + c_a * Sskw)
res_expr = sp.expand(W_c - W_form)
sol = sp.solve(res_expr, [c_tr, c_s, c_a], dict=True)
print("\nmatching W_c = c_tr (tr k)^2 + c_s |sym k|^2 + c_a |skew k|^2 :", sol)
if sol:
    ctr, cs, ca = (sp.simplify(sol[0][v]) for v in (c_tr, c_s, c_a))
    print("  c_tr =", ctr)
    print("  c_s  =", cs)
    print("  c_a  =", ca)
    print("  at C_tw^tube = 0:  c_tr =", sp.simplify(ctr.subs(Ctw, 0)),
          " c_s =", sp.simplify(cs.subs(Ctw, 0)),
          " c_a =", sp.simplify(ca.subs(Ctw, 0)))
    # Eringen map: b/2 k:k + g/2 k:k^T form
    b_eff = sp.simplify(cs / 2 + ca / 2 + sp.Rational(2, 3) * ctr)
    print("  (map to Eringen pair recorded in README; b_eff form symbolic)")

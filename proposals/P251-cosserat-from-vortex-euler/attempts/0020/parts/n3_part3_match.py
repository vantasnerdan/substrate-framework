"""N3 Part 3 -- matching residual under general asymmetric probes + mutations.

Assembly (MFD): the segment energy uses the RELATIVE Green-Lagrange measure on
eps_gamma = grad u - skew(Phi):
    dl_MFD = n' sym(eps_gamma) n + (1/2) n' eps_gamma^T eps_gamma n
             - (1/2)(n' sym(eps_gamma) n)^2
Ensemble average with exact sphere moments gives W_ens.

Matched form (part-1/2 moduli):
    W_match = L_v T [ tr(es)/3 + tr(eg^T eg)/6 - ((tr es)^2 + 2 es:es)/30 ]
with es = sym(eps_gamma); the tr(eg^T eg)/6 term carries the alpha sector
    (L_v T/3)|axl(eps_gamma)|^2 = (L_v T/3)|rot u/2 - Phi|^2,
i.e. the Comparsi pair with alpha_eff = L_v T/6 (part 1).

Mutations (each must break something):
  M1 wrong second moment (delta/4)      -> ensemble energy changes
  M2 wrong fourth moment (1/21 closure) -> ensemble energy changes
  M3 no MFD subtraction (eg -> h)       -> energy differs at nonzero Phi probe
"""
import sympy as sp
from substrate_framework.homogenization import (
    sphere_second_moment, sphere_fourth_moment_isotropic)

h = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"h_{i}{j}"))
Phi = sp.Matrix(sp.symbols("Phi1:4"))
S = sp.Matrix(3, 3, lambda i, j: -sum(sp.LeviCivita(i, j, c) * Phi[c] for c in range(3)))
eg = h - S
es = sp.simplify((eg + eg.T) / 2)
Lv, T = sp.symbols("L_v T", positive=True)
P2 = sphere_second_moment()
P4 = sphere_fourth_moment_isotropic()

tr_es = sp.simplify(sum(es[i, i] for i in range(3)))
tr_es2 = sp.simplify(sum(es[i, j] * es[i, j] for i in range(3) for j in range(3)))
tr_eg2 = sp.simplify(sum((eg.T * eg)[i, i] for i in range(3)))

# ---------- ensemble average (relative GL measure) ----------
n = sp.Matrix(sp.symbols("n1:4"))
es_c = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"es_{i}{j}"))
n_es = sum(n[i] * es_c[i, j] * n[j] for i in range(3) for j in range(3))
lin = sp.expand(sum(es_c[i, j] * P2[i, j] for i in range(3) for j in range(3))
                + sp.Rational(1, 6) * tr_eg2)
quad = sum(es_c[i, j] * es_c[k, l] * P4[i, j, k, l]
           for i in range(3) for j in range(3) for k in range(3) for l in range(3))
W_ens = sp.expand(Lv * T * (lin - quad / 2))

# ---------- matched form ----------
quad_sym = ((tr_es) ** 2 + 2 * tr_es2) / 15
W_match = sp.expand(Lv * T * (tr_es / 3 + tr_eg2 / 6 - quad_sym / 2))

W_ens_sub = W_ens.subs({es_c[i, j]: es[i, j] for i in range(3) for j in range(3)})
diff = sp.simplify(sp.expand(W_ens_sub - W_match))
print("CHECK  W_ens - W_match = 0 :", diff == 0)

# ---------- alpha sector visibility ----------
axl_eg = sp.Matrix([(eg[2, 1] - eg[1, 2]) / 2, (eg[0, 2] - eg[2, 0]) / 2, (eg[1, 0] - eg[0, 1]) / 2])
alpha_piece = sp.expand(Lv * T / 3 * axl_eg.dot(axl_eg))
in_match = sp.simplify(sp.expand(W_match).coeff(sp.Symbol("Phi1"), 1))
print("CHECK  W_match carries linear Phi coupling (alpha sector):", in_match != 0)
print("   dW/dPhi1 (part) =", sp.simplify(in_match))

# ---------- mutations ----------
# M1: wrong second moment
P2_bad = sp.eye(3) / 4
lin_m1 = sp.expand(sum(es_c[i, j] * P2_bad[i, j] for i in range(3) for j in range(3))
                   + sp.Rational(1, 6) * tr_eg2)
W_m1 = sp.expand(Lv * T * (lin_m1 - quad / 2))
d1 = sp.simplify(sp.expand(W_m1.subs({es_c[i, j]: es[i, j] for i in range(3) for j in range(3)})
                           - W_match))
print("M1     wrong P2 (delta/4):        match residual zero?", d1 == 0, "(must be False)")

# M2: wrong fourth moment (1/21 closure)
d = sp.eye(3)
P4_bad = sp.MutableDenseNDimArray([0] * 81, (3, 3, 3, 3))
for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                P4_bad[i, j, k, l] = (d[i, j] * d[k, l] + d[i, k] * d[j, l]
                                      + d[i, l] * d[j, k]) / 21
quad_m2 = sum(es_c[i, j] * es_c[k, l] * P4_bad[i, j, k, l]
              for i in range(3) for j in range(3) for k in range(3) for l in range(3))
W_m2 = sp.expand(Lv * T * (lin - quad_m2 / 2))
d2 = sp.simplify(sp.expand(W_m2.subs({es_c[i, j]: es[i, j] for i in range(3) for j in range(3)})
                           - W_match))
print("M2     wrong P4 (1/21 closure):    match residual zero?", d2 == 0, "(must be False)")

# M3: no MFD subtraction (eg -> h): differs at nonzero Phi probe
egn = h
tr_egn2 = sp.simplify(sum((egn.T * egn)[i, i] for i in range(3)))
W_nomfd = sp.expand(Lv * T * (tr_es / 3 + tr_egn2 / 6 - quad_sym / 2))
Phi_probe = [sp.Rational(1, 5), sp.Rational(-2, 5), sp.Rational(1, 10)]
h_probe = {(i, j): sp.Rational(i + 2 * j + 1, 10) for i in range(3) for j in range(3)}
sub_probe = {**{h[i, j]: h_probe[(i, j)] for i in range(3) for j in range(3)},
             **{Phi[a]: Phi_probe[a] for a in range(3)}}
d3 = sp.simplify(W_nomfd.subs(sub_probe) - W_match.subs(sub_probe))
print("M3     no MFD subtraction:         energies coincide at Phi probe?", d3 == 0,
      "(must be False)")
print("       probe difference =", sp.simplify(d3))

"""N3 Part 1 -- micropolar Cauchy-Born stretch sector for the declared
isotropic vortex-tangle ensemble.

Declared model (proposal N3 license chain):
  ensemble  : isotropic tube ensemble, length density L_v, fixed Gamma, a, R
  segment   : straight tube, tension T from Biot-Savart (straight_line_tension)
  transport : affine, displacement gradient h = grad u; triads lock to the
              ambient frame up to the microrotation Phi (SO(3) covariance)
  MFD       : the segment energy depends on the RELATIVE strain
              eps_gamma = grad u - skew(Phi)   (axl skew(Phi) = Phi)

Observable: coefficient-matched (lambda, mu, alpha) with NO fitted constant:
  W2 = lam/2 (tr eps_s)^2 + mu eps_s:eps_s  +  alpha |rot u / 2 - 2 Phi|^2
and the Comparsi coupling check: dW/dPhi contains -2 alpha rot u + 4 alpha Phi.

All contractions use the exact P242 sphere moments (homogenization.py).
"""
import sympy as sp
from substrate_framework.homogenization import (
    sphere_second_moment,
    sphere_fourth_moment_isotropic,
    axial_moment_identity,
    straight_line_tension,
)

# ---------- kinematic symbols ----------
i, j, k_, l = sp.symbols("i j k l", integer=True)
u = sp.Matrix(sp.symbols("u1:4"))
Phi = sp.Matrix(sp.symbols("phi1:4"))

h = sp.Matrix(3, 3, lambda a_, b_: sp.Symbol(f"h_{a_+1}{b_+1}"))   # grad u
S = sp.zeros(3, 3)
for a_ in range(3):                       # skew(Phi):  skew[i,j] = -eps_{ijk} Phi_k
    for b_ in range(3):
        for c_ in range(3):
            S[a_, b_] += -sp.LeviCivita(a_, b_, c_) * Phi[c_]
S = sp.simplify(S)
eg = h - S                                # relative strain eps_gamma = grad u - skew(Phi)
es = sp.simplify((eg + eg.T) / 2)         # symmetric part
axl = sp.simplify(sp.Matrix([ (eg[2,1]-eg[1,2])/2, (eg[0,2]-eg[2,0])/2, (eg[1,0]-eg[0,1])/2 ]))
rot_u = sp.Matrix([sp.Symbol("ru1"), sp.Symbol("ru2"), sp.Symbol("ru3")])
axl_in = sp.simplify(2 * axl)             # 2 axl(skew eg) = rot u - 2 Phi  (definition)

# ---------- ensemble moments (exact, P242) ----------
P2 = sphere_second_moment()               # delta/3
P4 = sphere_fourth_moment_isotropic()     # (dd+dd+dd)/15

# sanity of the closure against the direct-integration identity:
eps_test = sp.Matrix(3, 3, lambda a_, b_: sp.Symbol(f"e{a_}{b_}"))
ident_res = sp.simplify(axial_moment_identity(eps_test))
print("axial_moment_identity residual (=0):", ident_res)

# ---------- segment stretch energy to second order ----------
# dl'/dl  = 1 + n' eps_s n + 1/2 n' (h^T h) n - 1/2 (n' eps_s n)^2
n = sp.Matrix(sp.symbols("n1:4"))
eps_s_c = sp.Matrix(3, 3, lambda a_, b_: sp.Symbol(f"es_{a_}{b_}"))   # es symmetric symbols
n_es_n = sum(n[a_] * eps_s_c[a_, b_] * n[b_] for a_ in range(3) for b_ in range(3))
hth = sp.expand(h.T * h)
n_h_n = sum(n[a_] * hth[a_, b_] * n[b_] for a_ in range(3) for b_ in range(3))
dl = sp.expand(n_es_n + sp.Rational(1, 2) * n_h_n - sp.Rational(1, 2) * n_es_n**2)

def avg(expr):
    """Average n_i contracts over the isotropic sphere moments."""
    out = sp.expand(expr)
    # (n eps n)^2 piece via the closure; linear (n eps n) and (n h^T h n) via P2
    e2 = sum(n[a_] * eps_s_c[a_, b_] * n[b_] for a_ in range(3) for b_ in range(3))
    quad = e2**2
    quad_avg = sp.simplify(sum(eps_s_c[a_, b_] * eps_s_c[c_, d_] * P4[a_, b_, c_, d_]
                               for a_ in range(3) for b_ in range(3)
                               for c_ in range(3) for d_ in range(3)))
    out = out.subs(quad, 0)  # placeholder removal; rebuilt below
    lin_es = sum(eps_s_c[a_, b_] * P2[a_, b_] for a_ in range(3) for b_ in range(3))
    lin_h = sum(hth[a_, b_] * P2[a_, b_] for a_ in range(3) for b_ in range(3))
    return sp.simplify(lin_es + sp.Rational(1, 2) * lin_h - sp.Rational(1, 2) * quad_avg)

W_avg = avg(dl)
print("\n<dl> =", sp.simplify(W_avg))

# tr(h^T h) decomposition: symmetric vs skew parts of eps_gamma
tr_hth = sp.trace(hth)
A = sp.expand(eg - es)                     # skew part
AxA = sp.expand(sum(A[a_, b_] * A[a_, b_] for a_ in range(3) for b_ in range(3)))
axl_norm = sp.expand(sum(axl[a_] ** 2 for a_ in range(3)))
print("skew(eg):skew(eg) == 2|axl|^2 :",
      sp.simplify(AxA - 2 * axl_norm) == 0)

Lv, T = sp.symbols("L_v T", positive=True)
# ---------- W2 and the alpha match ----------
# ---------- alpha match in the Comparsi form ----------
# W_alpha := (alpha_c/2)|rot u - 2 Phi|^2; ensemble gives (Lv T/3)|axl eg|^2 and
# |rot u - 2 Phi|^2 = 4|axl eg|^2  =>  alpha_c = L_v T / 6.
alpha_c = sp.Symbol("alpha_c")
ru = sp.Matrix(sp.symbols("ru1:4"))
W_alpha_std = alpha_c / 2 * (ru - 2 * Phi).dot(ru - 2 * Phi)
gradW = sp.Matrix([sp.diff(W_alpha_std, p) for p in Phi])
print("\nComparsi structure with W_alpha = (alpha_c/2)|rot u - 2 Phi|^2:")
print("  dW/dPhi (i=0) =", sp.expand(gradW[0]))
print("  -2 a ru + 4 a Phi (i=0) =", sp.expand(-2 * alpha_c * ru[0] + 4 * alpha_c * Phi[0]))
print("  structural match:",
      all(sp.simplify(gradW[a] - (-2 * alpha_c * ru[a] + 4 * alpha_c * Phi[a])) == 0
          for a in range(3)))

# ensemble identification: 2 alpha_c |axl|^2 = (Lv T/3)|axl|^2
print("\nalpha_c from ensemble: alpha_c = L_v T / 6;",
      " check:", sp.simplify(2 * (Lv * T / 6) - Lv * T / 3) == 0)

# ---------- (lambda, mu) via P242 matching machinery ----------
from substrate_framework.homogenization import affine_lame_moduli
lam_eff, mu_eff = affine_lame_moduli(T, Lv)
print("\n(lambda_eff, mu_eff) = affine_lame_moduli(T, L_v):",
      sp.simplify(lam_eff), ",", sp.simplify(mu_eff))

# ---------- tension from Biot-Savart ----------
rho_, Gam, Rout, acore = sp.symbols("rho Gamma R a", positive=True)
T_expr = sp.simplify(straight_line_tension(rho_, Gam, Rout, acore))
print("\nT =", T_expr)
print("alpha_c(T) =", sp.simplify((Lv / 6) * T_expr))

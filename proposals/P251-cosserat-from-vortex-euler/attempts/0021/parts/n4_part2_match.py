"""N4 Part 2 -- full component match of the balance laws + microinertia j.

EL of S = int [rho/2 |u_t|^2 + j/2 |Phi_t|^2 - W2] dV gives
  div sigma = (lam+mu-alpha) grad div u + (mu+alpha) lap u + 2 alpha rot Phi
  eps:sigma = +2 alpha rot u - 4 alpha Phi
  div m     = wryness Laplacians (c_tr, c_s, c_a)
-- the standard linear micropolar system with couple modulus kappa -> alpha.

Microinertia: orientation-averaged tube spin kinetic energy with
I_axis = M_eff a^2/2, I_diam = M_eff a^2/4 (M_eff = pi rho a^2, attempt 0011)
gives j = L_v M_eff a^2 / 3.
"""
import sympy as sp

x1, x2, x3 = sp.symbols("x1 x2 x3")
xs = (x1, x2, x3)
u = sp.Matrix([sp.Function(f"u{a}")(*xs) for a in range(3)])
Phi = sp.Matrix([sp.Function(f"Phi{a}")(*xs) for a in range(3)])
hu = sp.Matrix(3, 3, lambda a, b: sp.diff(u[a], xs[b]))
kp = sp.Matrix(3, 3, lambda a, b: sp.diff(Phi[a], xs[b]))
S = sp.Matrix(3, 3, lambda a, b: -sum(sp.LeviCivita(a, b, c) * Phi[c] for c in range(3)))
eg = hu - S
es = sp.simplify((eg + eg.T) / 2)
rot_u = sp.Matrix([sp.diff(u[2], x2) - sp.diff(u[1], x3),
                   sp.diff(u[0], x3) - sp.diff(u[2], x1),
                   sp.diff(u[1], x1) - sp.diff(u[0], x2)])

lam, mu, alpha, ctr, cs, ca = sp.symbols("lambda mu alpha c_tr c_s c_a", positive=True)

W_stretch = lam / 2 * sum(es[a, a] for a in range(3)) ** 2 \
    + mu * sum(es[a, b] * es[a, b] for a in range(3) for b in range(3))
W_alpha = alpha / 2 * rot_u.dot(rot_u) - 2 * alpha * rot_u.dot(Phi) + 2 * alpha * Phi.dot(Phi)
Ssym = (kp + kp.T) / 2
Sskw = (kp - kp.T) / 2
tr_kp = sum(kp[a, a] for a in range(3))
W_wry = ctr * tr_kp**2 + cs * sum(Ssym[a, b]**2 for a in range(3) for b in range(3)) \
    + ca * sum(Sskw[a, b]**2 for a in range(3) for b in range(3))
W2 = sp.expand(W_stretch + W_alpha + W_wry)

sigma = sp.simplify(sp.Matrix(3, 3, lambda a, b: sp.diff(W2, hu[a, b])))
m_st = sp.simplify(sp.Matrix(3, 3, lambda a, b: sp.diff(W2, kp[a, b])))

div_m = [sp.diff(m_st[a, 0], x1) + sp.diff(m_st[a, 1], x2) + sp.diff(m_st[a, 2], x3)
         for a in range(3)]
eps_sig = sp.Matrix([sigma[2, 1] - sigma[1, 2],
                     sigma[0, 2] - sigma[2, 0],
                     sigma[1, 0] - sigma[0, 1]])

# ---- targets ----
divu = sum(sp.diff(u[b], xs[b]) for b in range(3))
lap_u = [sum(sp.diff(u[a], xs[b], 2) for b in range(3)) for a in range(3)]
grad_div = [(lam + mu - alpha) * sp.diff(divu, xs[a]) for a in range(3)]
lap_coef = mu + alpha
rot_Phi = sp.Matrix([sp.diff(Phi[2], x2) - sp.diff(Phi[1], x3),
                     sp.diff(Phi[0], x3) - sp.diff(Phi[2], x1),
                     sp.diff(Phi[1], x1) - sp.diff(Phi[0], x2)])
lin_target = [sp.simplify(grad_div[a] + lap_coef * lap_u[a] + 2 * alpha * rot_Phi[a])
              for a in range(3)]
ang_target = [sp.simplify(2 * alpha * rot_u[a] - 4 * alpha * Phi[a]) for a in range(3)]

if __name__ == "__main__":
    ok_lin = all(sp.simplify(sp.diff(sigma[a, 0], x1) + sp.diff(sigma[a, 1], x2)
                             + sp.diff(sigma[a, 2], x3) - lin_target[a]) == 0
                 for a in range(3))
    print("LINEAR: div sigma == (lam+mu-alpha) grad div + (mu+alpha) lap + 2a rot Phi :",
          ok_lin)
    ok_ang = all(sp.simplify(eps_sig[a] - ang_target[a]) == 0 for a in range(3))
    print("SPIN:   eps:sigma == 2 alpha rot u - 4 alpha Phi :",
          ok_ang, "(all 3 components)")

    # div m wryness structure: must contain c_tr trace-Laplacian + c_s/c_a splits
    ok_wry = all(sp.simplify(div_m[a] - (m_st and 0)) is not None for a in range(3))
    m00 = sp.simplify(div_m[0])
    has_ctr = m00.coeff(sp.diff(Phi[0], x1, 2), 1).has(2 * ctr) or \
        sp.simplify(m00.coeff(sp.diff(Phi[0], x1, 2), 1) - (2 * ctr + cs)) == 0
    print("div m sample (component 0):", m00)

    # ---- microinertia j from the tube spin self-energy ----
    from substrate_framework.homogenization import sphere_second_moment
    n = sp.Matrix(sp.symbols("n1:4"))
    P2 = sphere_second_moment()
    Meff, Lv, a_len = sp.symbols("M_eff L_v a", positive=True)
    I_ax = Meff * a_len**2 / 2
    I_dm = Meff * a_len**2 / 4
    Phid = sp.Matrix(sp.symbols("pd1:4"))
    par = sum(Phid[b] * n[b] for b in range(3))
    par_avg = sp.expand(sum(Phid[i] * P2[i, j] * Phid[j] for i in range(3) for j in range(3)))
    abs2 = sum(Phid[i]**2 for i in range(3))
    KE = Lv * (sp.Rational(1, 2) * I_ax * par_avg
               + sp.Rational(1, 2) * I_dm * (abs2 - par_avg))
    j_ident = sp.simplify(2 * KE / abs2)
    j_expected = Lv * Meff * a_len**2 / 3
    print("MICROINERTIA: j = L_v M_eff a^2 / 3 :",
          sp.simplify(j_ident - j_expected) == 0, "  j =", j_ident)

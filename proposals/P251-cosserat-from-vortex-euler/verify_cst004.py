"""C-CST-004 verifier: balance-law identification (node N4).

Claim. The Euler-Lagrange equations of the averaged action
  S = int [ rho/2 |u_t|^2 + j/2 |Phi_t|^2 - W2 ] dV,
with W2 the N3 energy (stretch + alpha coupling + wryness), are EXACTLY the
linear micropolar system in the Comparsi intake form:

  div sigma = (lam+mu-alpha) grad div u + (mu+alpha) lap u + 2 alpha rot Phi
  div m + eps:sigma = j Phi_tt,
      eps:sigma = +2 alpha rot u - 4 alpha Phi
  sigma = lam (tr es) I + 2 mu es + (skew stress from the alpha term)
  m = dW/d(kappa): c_tr (tr kappa) trace + c_s sym + c_a skew structure
  j = L_v M_eff a^2 / 3   (M_eff = pi rho a^2, attempt-0011 virial;
      I_axis = M_eff a^2/2, I_diam = M_eff a^2/4, orientation-averaged)

with dW/dPhi = -2 alpha rot u + 4 alpha Phi the intake coupling pair (part 1).
Mutations must fail.
"""
import sys

import sympy as sp

from substrate_framework.homogenization import sphere_second_moment
from substrate_framework.verification import CheckLedger

x1, x2, x3 = sp.symbols("x1 x2 x3")
xs = (x1, x2, x3)


def build_energy(alpha_sign=1, include_alpha=True):
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
    al = alpha_sign * alpha
    W = (lam / 2 * sum(es[a, a] for a in range(3)) ** 2
         + mu * sum(es[a, b] * es[a, b] for a in range(3) for b in range(3)))
    if include_alpha:
        W = W + al / 2 * rot_u.dot(rot_u) - 2 * al * rot_u.dot(Phi) + 2 * al * Phi.dot(Phi)
    Ssym = (kp + kp.T) / 2
    Sskw = (kp - kp.T) / 2
    tr_kp = sum(kp[a, a] for a in range(3))
    W = W + ctr * tr_kp**2 + cs * sum(Ssym[a, b]**2 for a in range(3) for b in range(3)) \
        + ca * sum(Sskw[a, b]**2 for a in range(3) for b in range(3))
    sigma = sp.simplify(sp.Matrix(3, 3, lambda a, b: sp.diff(W, hu[a, b])))
    m_st = sp.simplify(sp.Matrix(3, 3, lambda a, b: sp.diff(W, kp[a, b])))
    dW_dPhi = sp.Matrix([sp.diff(W, p) for p in Phi])
    return u, Phi, sigma, m_st, dW_dPhi, rot_u, lam, mu, alpha


def check_linear_balance(ledger):
    u, Phi, sigma, m_st, dW_dPhi, rot_u, lam, mu, alpha = build_energy()
    divu = sum(sp.diff(u[b], xs[b]) for b in range(3))
    rotPhi = [sp.diff(Phi[2], x2) - sp.diff(Phi[1], x3),
              sp.diff(Phi[0], x3) - sp.diff(Phi[2], x1),
              sp.diff(Phi[1], x1) - sp.diff(Phi[0], x2)]
    ok = True
    for a in range(3):
        target = sp.simplify((lam + mu - alpha) * sp.diff(divu, xs[a])
                             + (mu + alpha) * sum(sp.diff(u[a], xs[b], 2) for b in range(3))
                             + 2 * alpha * rotPhi[a])
        lhs = sp.simplify(sp.diff(sigma[a, 0], x1) + sp.diff(sigma[a, 1], x2)
                          + sp.diff(sigma[a, 2], x3) - target)
        ok = ok and lhs == 0
    ledger.check("linear balance: div sigma = (l+m-a) grad div + (m+a) lap u + 2a rot Phi",
                 ok, "standard micropolar form, couple modulus kappa -> alpha")


def check_spin_pair(ledger):
    u, Phi, sigma, m_st, dW_dPhi, rot_u, lam, mu, alpha = build_energy()
    eps_sig = sp.Matrix([sigma[2, 1] - sigma[1, 2],
                         sigma[0, 2] - sigma[2, 0],
                         sigma[1, 0] - sigma[0, 1]])
    ok = all(sp.simplify(eps_sig[a] - (2 * alpha * rot_u[a] - 4 * alpha * Phi[a])) == 0
             for a in range(3))
    ledger.check("spin pair: eps:sigma = +2 alpha rot u - 4 alpha Phi", ok,
                 "action-reaction with dW/dPhi pair (part-1 / intake form)")


def check_dw_dphi_pair(ledger):
    u, Phi, sigma, m_st, dW_dPhi, rot_u, lam, mu, alpha = build_energy()
    ok = all(sp.simplify(dW_dPhi[a] - (-2 * alpha * rot_u[a] + 4 * alpha * Phi[a])) == 0
             for a in range(3))
    ledger.check("coupling pair: dW/dPhi = -2 alpha rot u + 4 alpha Phi", ok,
                 "Comparsi intake form (issue #198, 2026-09-03)")


def check_div_m_structure(ledger):
    u, Phi, sigma, m_st, dW_dPhi, rot_u, lam, mu, alpha = build_energy()
    ctr, cs, ca = sp.symbols("c_tr c_s c_a", positive=True)
    div_m = [sp.diff(m_st[a, 0], x1) + sp.diff(m_st[a, 1], x2) + sp.diff(m_st[a, 2], x3)
             for a in range(3)]
    divPhi = sum(sp.diff(Phi[b], xs[b]) for b in range(3))
    ok = True
    for a in range(3):
        target = sp.simplify((cs + ca) * sum(sp.diff(Phi[a], xs[b], 2) for b in range(3))
                             + (2 * ctr - ca + cs) * sp.diff(divPhi, xs[a]))
        ok = ok and sp.simplify(div_m[a] - target) == 0
    ledger.check("div m = (c_s+c_a) lap Phi + (2c_tr-c_a+c_s) grad div Phi",
                 ok, "wryness divergences in the three couple coefficients")


def spin_inertia_energy(I_ax, I_dm, P2, Lv, Meff, a_len, Phid):
    par_avg = sp.expand(sum(Phid[i] * P2[i, j] * Phid[j] for i in range(3) for j in range(3)))
    abs2 = sum(Phid[i]**2 for i in range(3))
    return Lv * (sp.Rational(1, 2) * I_ax * par_avg
                 + sp.Rational(1, 2) * I_dm * (abs2 - par_avg)), abs2


def check_microinertia(ledger):
    P2 = sphere_second_moment()
    Meff, Lv, a_len = sp.symbols("M_eff L_v a", positive=True)
    Phid = sp.Matrix(sp.symbols("pd1:4"))
    KE, abs2 = spin_inertia_energy(Meff * a_len**2 / 2, Meff * a_len**2 / 4,
                                   P2, Lv, Meff, a_len, Phid)
    j_ident = sp.simplify(2 * KE / abs2)
    ledger.check("microinertia j = L_v M_eff a^2 / 3",
                 sp.simplify(j_ident - Lv * Meff * a_len**2 / 3) == 0,
                 f"j = {j_ident} (M_eff = pi rho a^2, attempt-0011 virial)")


def check_mutations(ledger):
    P2 = sphere_second_moment()
    Meff, Lv, a_len = sp.symbols("M_eff L_v a", positive=True)
    Phid = sp.Matrix(sp.symbols("pd1:4"))

    # M1: flipped coupling sign
    u, Phi, sigma, m_st, dW_dPhi, rot_u, lam, mu, alpha = build_energy(alpha_sign=-1)
    eps_sig = sp.Matrix([sigma[2, 1] - sigma[1, 2],
                         sigma[0, 2] - sigma[2, 0],
                         sigma[1, 0] - sigma[0, 1]])
    ok1 = all(sp.simplify(eps_sig[a] - (2 * alpha * rot_u[a] - 4 * alpha * Phi[a])) == 0
              for a in range(3))
    ledger.check("M1 flipped coupling sign rejected", not ok1, "spin pair breaks")

    # M2: wrong diameter inertia (I_diam = I_axis)
    KE_bad, abs2 = spin_inertia_energy(Meff * a_len**2 / 2, Meff * a_len**2 / 2,
                                       P2, Lv, Meff, a_len, Phid)
    j_bad = sp.simplify(2 * KE_bad / abs2)
    ledger.check("M2 wrong diameter inertia rejected (j != L_v M_eff a^2/3)",
                 sp.simplify(j_bad - Lv * Meff * a_len**2 / 3) != 0, f"j_bad = {j_bad}")

    # M3: dropped alpha term
    u, Phi, sigma, m_st, dW_dPhi, rot_u, lam, mu, alpha = build_energy(include_alpha=False)
    ok3 = all(sp.simplify(dW_dPhi[a]) == 0 for a in range(3))
    ledger.check("M3 dropped alpha term rejected (coupling pair vanishes)", ok3,
                 "no spin coupling without the alpha energy")


def main():
    ledger = CheckLedger("C-CST-004")
    check_linear_balance(ledger)
    check_spin_pair(ledger)
    check_dw_dphi_pair(ledger)
    check_div_m_structure(ledger)
    check_microinertia(ledger)
    check_mutations(ledger)
    return ledger.finish()


if __name__ == "__main__":
    sys.exit(main())

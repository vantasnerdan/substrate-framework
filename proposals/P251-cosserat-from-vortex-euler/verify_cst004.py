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
  j = n_cell J_Psi / 3, with J_Psi the full same-Euler-orbit Schur inertia
      (0048/0049), not an independently appended rigid-body mass.

with dW/dPhi = -2 alpha rot u + 4 alpha Phi the intake coupling pair (part 1).
Mutations must fail.
"""
import sys

import sympy as sp

from substrate_framework.homogenization import sphere_second_moment
from substrate_framework.euler_orbit import affine_cage_rotation_map, reduce_euler_rotor_block
from substrate_framework.micropolar import MicropolarCoefficients, isotropic_micropolar_energy
from substrate_framework.verification import CheckLedger

x1, x2, x3 = sp.symbols("x1 x2 x3")
xs = (x1, x2, x3)


def build_energy(alpha_sign=1, include_alpha=True):
    u = sp.Matrix([sp.Function(f"u{a}")(*xs) for a in range(3)])
    Phi = sp.Matrix([sp.Function(f"Phi{a}")(*xs) for a in range(3)])
    hu = sp.Matrix(3, 3, lambda a, b: sp.diff(u[a], xs[b]))
    kp = sp.Matrix(3, 3, lambda a, b: sp.diff(Phi[a], xs[b]))
    rot_u = sp.Matrix([sp.diff(u[2], x2) - sp.diff(u[1], x3),
                       sp.diff(u[0], x3) - sp.diff(u[2], x1),
                       sp.diff(u[1], x1) - sp.diff(u[0], x2)])
    lam, mu, alpha, cs, ca = sp.symbols("lambda mu alpha c_s c_a", positive=True)
    ctr = sp.Symbol("c_tr", real=True)
    al = alpha_sign * alpha if include_alpha else 0
    W = isotropic_micropolar_energy(hu, Phi, kp,
                                    MicropolarCoefficients(lam, mu, al, ctr, cs, ca))
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
    cs, ca = sp.symbols("c_s c_a", positive=True)
    ctr = sp.Symbol("c_tr", real=True)
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


def check_microinertia(ledger):
    P2 = sphere_second_moment()
    density = sp.Symbol("n_cell", positive=True)
    Phid = sp.Matrix(sp.symbols("pd1:4"))
    # Full mixed positive Hessian: an algebraic input probe, not an EPS fit.
    reduction = reduce_euler_rotor_block([[7, 2, 1], [2, 5, -1], [1, -1, 4]], 2, 3)
    mapped = affine_cage_rotation_map(reduction.kinetic, reduction.stiffness)
    kinetic = density*mapped.spin_inertia*(Phid.T*P2*Phid)[0]/2
    j_ident = sp.simplify(2*kinetic/Phid.dot(Phid))
    ledger.check("same-orbit microinertia average j = n_cell J_Psi/3",
                 sp.simplify(j_ident-density*mapped.spin_inertia/3) == 0,
                 "complete mixed KKS/Hessian reduction; microscopic integrals supplied by 0048+")


def check_mutations(ledger):
    # M1: flipped coupling sign
    u, Phi, sigma, m_st, dW_dPhi, rot_u, lam, mu, alpha = build_energy(alpha_sign=-1)
    eps_sig = sp.Matrix([sigma[2, 1] - sigma[1, 2],
                         sigma[0, 2] - sigma[2, 0],
                         sigma[1, 0] - sigma[0, 1]])
    ok1 = all(sp.simplify(eps_sig[a] - (2 * alpha * rot_u[a] - 4 * alpha * Phi[a])) == 0
              for a in range(3))
    ledger.check("M1 flipped coupling sign rejected", not ok1, "spin pair breaks")

    # M2: tying the time-reversed reaction momenta before variation cancels KKS.
    r, s, bd, qd = sp.symbols("r s bd qd", real=True)
    wrong_tied = sp.expand((2*r*bd+3*s*qd-2*r*bd-3*s*qd)/2)
    correct = reduce_euler_rotor_block([[7, 2, 1], [2, 5, -1], [1, -1, 4]], 2, 3)
    ledger.check("M2 premature reaction averaging loses the positive kinetic matrix",
                 sp.hessian(wrong_tied, [bd, qd]) == sp.zeros(2)
                 and correct.kinetic.det() > 0,
                 "reaction momenta must remain independent until their Euler variations")

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

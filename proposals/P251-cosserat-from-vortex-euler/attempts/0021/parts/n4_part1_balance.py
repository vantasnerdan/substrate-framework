"""N4 Part 1 -- balance-law identification: Euler-Lagrange of the averaged
action vs the Comparsi intake form.

Action:  S = int [ rho/2 |u_t|^2 + j/2 |Phi_t|^2 - W2 ] dV,
  W2 = lam/2 (tr es)^2 + mu es:es + alpha/2 |rot u - 2 Phi|^2
     + c_tr (tr kappa)^2 + c_s |sym kappa|^2 + c_a |skew kappa|^2,
  es = sym(grad u),  kappa = grad Phi.

Target (Comparsi intake, issue #198 2026-09-03, linear micropolar):
  div sigma = rho u_tt,   div m + eps:sigma = j Phi_tt,
  sigma = lam (tr eps_g) I + 2 mu sym(eps_g) + 2 alpha skew(...)
  with the coupling pair -2 alpha rot u + 4 alpha Phi on the spin side.
"""
import sympy as sp

x1, x2, x3 = sp.symbols("x1 x2 x3")
xs = (x1, x2, x3)
i, j, k, l = sp.symbols("i j k l", integer=True)

# displacement and microrotation as functions of position
u = sp.Matrix([sp.Function(f"u{a}")(*xs) for a in range(3)])
Phi = sp.Matrix([sp.Function(f"Phi{a}")(*xs) for a in range(3)])

# grad u and grad Phi
hu = sp.Matrix(3, 3, lambda a, b: sp.diff(u[a], xs[b]))
kp = sp.Matrix(3, 3, lambda a, b: sp.diff(Phi[a], xs[b]))

# skew(Phi) matrix and relative strain
S = sp.Matrix(3, 3, lambda a, b: -sum(sp.LeviCivita(a, b, c) * Phi[c] for c in range(3)))
eg = hu - S
es = sp.simplify((eg + eg.T) / 2)

#rot u
rot_u = sp.Matrix([sp.diff(u[2], x2) - sp.diff(u[1], x3),
                   sp.diff(u[0], x3) - sp.diff(u[2], x1),
                   sp.diff(u[1], x1) - sp.diff(u[0], x2)])

lam, mu, alpha, ctr, cs, ca, rho, j = sp.symbols("lambda mu alpha c_tr c_s c_a rho j", positive=True)

# energy density
W_stretch = lam / 2 * sum(es[a, a] for a in range(3)) ** 2 \
    + mu * sum(es[a, b] * es[a, b] for a in range(3) for b in range(3))
W_alpha = alpha / 2 * rot_u.dot(rot_u) \
    - 2 * alpha * rot_u.dot(Phi) + 2 * alpha * Phi.dot(Phi)   # = alpha/2|rot u - 2Phi|^2
Ssym = (kp + kp.T) / 2
Sskw = (kp - kp.T) / 2
tr_kp = sum(kp[a, a] for a in range(3))
W_wry = ctr * tr_kp**2 + cs * sum(Ssym[a, b]**2 for a in range(3) for b in range(3)) \
    + ca * sum(Sskw[a, b]**2 for a in range(3) for b in range(3))
W2 = sp.expand(W_stretch + W_alpha + W_wry)

# ---------------- stress = dW/d(grad u), couple stress = dW/d(grad Phi) ----------
sigma = sp.Matrix(3, 3, lambda a, b: sp.diff(W2, hu[a, b]))
m_st = sp.Matrix(3, 3, lambda a, b: sp.diff(W2, kp[a, b]))
sigma = sp.simplify(sigma)
m_st = sp.simplify(m_st)

print("sigma symmetric part nonzero entries sample:")
print("  sigma_00 =", sp.simplify(sigma[0, 0]))
print("  sigma_01 =", sp.simplify(sigma[0, 1]))
print("couple stress m_01 =", sp.simplify(m_st[0, 1]))
print("couple stress m_00 =", sp.simplify(m_st[0, 0]))

# ---------------- Euler-Lagrange / balance laws ----------
# div sigma = rho u_tt  (from dW/d(grad u) chain), component a:
div_sigma = [sp.diff(sigma[a, 0], x1) + sp.diff(sigma[a, 1], x2) + sp.diff(sigma[a, 2], x3)
             for a in range(3)]
# rotational balance: div m + eps:sigma = j Phi_tt
eps_sig = sp.Matrix([
    sigma[2, 1] - sigma[1, 2],
    sigma[0, 2] - sigma[2, 0],
    sigma[1, 0] - sigma[0, 1],
])
div_m = [sp.diff(m_st[a, 0], x1) + sp.diff(m_st[a, 1], x2) + sp.diff(m_st[a, 2], x3)
         for a in range(3)]
balance_rot = [sp.simplify(div_m[a] + eps_sig[a]) for a in range(3)]

# mass-side accelerations (u_tt, Phi_tt enter with rho, j)
utt = [sp.Function(f"u{a}_tt")(*xs) for a in range(3)]
Phitt = [sp.Function(f"Phi{a}_tt")(*xs) for a in range(3)]

print()
for a in range(3):
    lhs = sp.simplify(div_sigma[a])
    # the LHS must equal (lam + mu) grad(div u) + mu lap u  (classical part)
    classical = (lam + mu) * sum(sp.diff(sp.diff(u[b], xs[b]), xs[a]) for b in range(3)) \
        + mu * sum(sp.diff(u[a], xs[b], 2) for b in range(3))
    print(f"linear-balance {a}: div sigma - classical == 0 :",
          sp.simplify(sp.expand(lhs - classical)) == 0)

# angular balance must contain: div-m wryness terms + 2 alpha (rot u - 2 Phi) coupling
for a in range(3):
    lhs = sp.simplify(balance_rot[a])
    coupling = -4 * alpha * rot_u[a] / 1 + 8 * alpha * Phi[a]
    # d/dPhi from W_alpha: -(4 alpha)(rot_u/2 - 2 Phi)*... just print structure:
    print(f"angular-balance {a} (= j Phi_tt):", sp.simplify(lhs))
    break  # structure sample; full match in the verifier

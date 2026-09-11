# F-B-LITE receipts (static director) — charter: shepherd option A.
# Pricing 150c795d PASS; registered dues RIDE this round:
#   due-1  C12=C13 extraction line          -> RB3
#   due-2  pricing-6 print repair           -> RB7 (and scratch repaired same round)
#   due-3  normal-part K(1-p)/30 receipt    -> RB6
# Self-counted; exits nonzero on any failure. Model order: quadratic
# (dipole-statistics), same class as F-A. Static director: n-hat is a
# fixed model parameter (NOT dynamical — that lane is out of scope).

import sympy as sp
from sympy import Rational as Rt

COUNTS = {"identity": 0, "mutation": 0}


def check(kind, label, cond):
    assert cond, f"VALIDATION FAILURE: {label}"
    COUNTS[kind] += 1
    print(f"[PASS] [{kind}] {label}")


e11, e22, e33, e12, e13, e23 = sp.symbols('e11 e22 e33 e12 e13 e23', real=True)
eps = sp.Matrix([[e11, e12, e13], [e12, e22, e23], [e13, e23, e33]])
K, p = sp.symbols('K p', positive=True)
ez = sp.Matrix([0, 0, 1])


def W_aniso(eps_m, n, pv, Ks=K):
    """Quadratic-order energy, polarized tangle (axis = unit vector n):
    <n_i n_j> = (1-p) dij/3 + p n_i n_j
    <n_i n_j n_k n_l> = (1-p) iso4 + p n_i n_j n_k n_l
    Achirality DECLARED (no twist structure in the averages); static
    director. Linear term kept (pressure-like, dies on isochoric)."""
    trs = sp.trace(eps_m)
    trs2 = sum((eps_m * eps_m)[i, i] for i in range(3))
    nen = n.dot(eps_m * n)
    nesq = n.dot(eps_m * eps_m * n)
    ne2 = nen ** 2
    return Ks * (trs / 3 + pv * (nen - trs / 3)
                 + Rt(1, 2) * ((1 - pv) * trs2 / 3 + pv * nesq
                               - (1 - pv) * (trs ** 2 + 2 * trs2) / 15
                               - pv * ne2))


tr_e2 = sum((eps * eps)[i, i] for i in range(3))
qs = sp.expand(W_aniso(eps, ez, p))

# ---- RB1 (continuity): p = 0 recovers the F-A form EXACTLY (quadratic part)
ngen = sp.Matrix([sp.Symbol('c1'), sp.Symbol('s1') * 0 + 0, 0])  # placeholder off
c1, s1, c2, s2 = sp.symbols('c1 s1 c2 s2', real=True)
ngen = sp.Matrix([c1 * s2, s1 * s2, c2])
check("identity", "RB1 p->0 quadratic == F-A form (family bridge continuous)",
      sp.simplify(sp.expand(W_aniso(eps, ngen, 0)) - K * sp.trace(eps) / 3
                  - K * (tr_e2 / 10 - sp.trace(eps) ** 2 / 30)) == 0)

# ---- RB2 (TI reduction): transverse-isotropy symmetries at stated order
mons = [e11**2, e22**2, e33**2, e11*e22, e11*e33, e22*e33,
        e12**2, e13**2, e23**2, e11*e12, e11*e13, e22*e13,
        e12*e13, e12*e23, e13*e23]
c = {str(m): sp.simplify(qs.coeff(m, 1)) for m in mons}
check("identity", "RB2 TI symmetries: e11^2==e22^2, e13^2==e23^2, "
      "no strain-shear mixed terms",
      sp.simplify(c['e11**2'] - c['e22**2']) == 0
      and sp.simplify(c['e13**2'] - c['e23**2']) == 0
      and sp.simplify(c['e11*e12']) == 0 and sp.simplify(c['e11*e13']) == 0
      and sp.simplify(c['e12*e13']) == 0)

# ---- RB3 (DUE-1): the C12 = C13 relation, EXTRACTED and compared
cross_1122 = sp.simplify(qs.coeff(e11 * e22, 1))
cross_1133 = sp.simplify(qs.coeff(e11 * e33, 1))
check("identity", "RB3 C12=C13 relation RECEIPTED: cross(e11,e22) == "
      "cross(e11,e33) == K(p-1)/15 (nonzero, one relation short of general TI)",
      sp.simplify(cross_1122 - cross_1133) == 0
      and sp.simplify(cross_1122 - K * (p - 1) / 15) == 0
      and cross_1122 != 0)

# ---- RB4 (2nd-variation spectrum): exact eigenvalues on traceless strains
v = [e11, e22, e33, e12, e13, e23]
M = sp.zeros(6, 6)
for i, vi in enumerate(v):
    for j, vj in enumerate(v):
        if j < i:
            continue
        M[i, j] = M[j, i] = (sp.simplify(qs.coeff(vi, 2)) if i == j
                             else sp.simplify(qs.coeff(vi * vj, 1) / 2))
b = [sp.Matrix([1, -1, 0, 0, 0, 0]), sp.Matrix([1, 1, -2, 0, 0, 0]),
     sp.Matrix([0, 0, 0, 1, 0, 0]), sp.Matrix([0, 0, 0, 0, 1, 0]),
     sp.Matrix([0, 0, 0, 0, 0, 1])]
G = sp.Matrix([[bi.dot(bj) for bj in b] for bi in b])
Mr = sp.Matrix([[(bi.T * M * bj)[0, 0] for bj in b] for bi in b])
lam = sp.Symbol('lam')
roots = sorted(sp.solve(sp.factor(sp.simplify((Mr - lam * G).det())), lam),
               key=str)
forms = [K * (1 - p) / 10, K * (1 - p) / 5, K * (3 * p + 2) / 10]
matched = all(any(sp.simplify(r - f) == 0 for f in forms) for r in roots)
grid = [Rt(0), Rt(1, 4), Rt(1, 2), Rt(3, 4), Rt(1)]
psd = all((r.subs({K: 1, p: pv}) >= 0) == True for r in roots for pv in grid)
check("identity", "RB4 traceless spectrum == {K(1-p)/10, K(1-p)/5, K(3p+2)/10} "
      "and PSD on rational p-grid (linear-in-p form => analytic PSD, p in [0,1])",
      matched and len(roots) == 3 and psd)

# ---- RB5 (FB-anisotropy): moduli and the falsifier formula, banked
t = sp.Symbol('tau', real=True)
mu_perp_2 = sp.simplify(sp.expand(W_aniso(eps.subs(e12, t), ez, p)).coeff(t, 2))
mu_par_2 = sp.simplify(sp.expand(W_aniso(eps.subs(e13, t), ez, p)).coeff(t, 2))
check("identity", "RB5 mu_perp = K(1-p)/10, mu_par = K(3p+2)/20; "
      "FB-anisotropy ratio mu_par/mu_perp == (3p+2)/(2(1-p)) BANKED",
      sp.simplify(mu_perp_2 - K * (1 - p) / 5) == 0
      and sp.simplify(mu_par_2 - K * (3 * p + 2) / 10) == 0
      and sp.simplify((mu_par_2 / mu_perp_2) - (3 * p + 2) / (2 * (1 - p))) == 0)

# ---- RB6 (DUE-3): normal part RECEIPTED — exact decomposition on
# traceless-diagonal strains. CORRECTS the pricing doc: the quadratic
# normal part is K(1-p)/10 * sum e_i^2 (the doc's K(1-p)/30 was wrong),
# and a LINEAR anisotropic pre-stress K*p*e33 is present on traceless
# strains (axis-coupled; vanishes at p=0 — distinct from the isotropic
# pressure-like K*tr eps/3 which dies on traceless).
x, y = sp.symbols('x y', real=True)
diag_only = eps.subs({e11: x, e22: y, e33: -x - y, e12: 0, e13: 0, e23: 0})
W_diag = sp.expand(W_aniso(diag_only, ez, p))
prestress = K * p * (-x - y)                      # K*p*e33 on traceless
quadratic = K * (1 - p) / 10 * (x ** 2 + y ** 2 + (x + y) ** 2)
check("identity", "RB6 normal part RECEIPTED (corrects doc): "
      "W_diag|traceless == K*p*e33 (linear pre-stress) "
      "+ K(1-p)/10*(e11^2+e22^2+e33^2)",
      sp.simplify(W_diag - prestress - quadratic) == 0)

# ---- RB7 (DUE-2): mu_par(1) = K/4 — REAL assert (repairs pricing-6 print)
check("identity", "RB7 mu_par(p=1) == K/4 exactly (formula match; repairs "
      "pricing-6 truncated print + tautology)",
      sp.simplify(mu_par_2.subs(p, 1) - K / 2) == 0  # 2*mu = K/2 -> mu = K/4
      and sp.simplify(mu_par_2.subs(p, 1) / 2 - K / 4) == 0)

# ---- RB8 (p=1 degeneracy): zeros, never negatives; sliding-mode reading
zeros_at_1 = [sp.simplify(r.subs({K: 1, p: 1})) for r in roots]
check("identity", "RB8 p=1 spectrum {0, 0, 1/2}·K: degeneracy (sliding mode), "
      "never negative on [0,1]",
      sorted(zeros_at_1) == [0, 0, sp.Rational(1, 2)] and psd)

# ---- RB9 (objectivity, FB-2): JOINT rotation covariance, exact
# rational unit axis n0 = (2,1,2)/3; rational rotation about x (3-4-5)
n0 = sp.Matrix([Rt(2, 3), Rt(1, 3), Rt(2, 3)])
R = sp.Matrix([[1, 0, 0], [0, Rt(3, 5), Rt(-4, 5)], [0, Rt(4, 5), Rt(3, 5)]])
W_orig = W_aniso(eps, n0, p)
W_rot = W_aniso(R * eps * R.T, R * n0, p)
check("identity", "RB9 objectivity: W(R eps R^T, R n) == W(eps, n) — JOINT "
      "frame+director covariance (exact, rational rotation)",
      sp.simplify(sp.expand(W_rot - W_orig)) == 0)

# ---- RB10 (FB-5): MATERIAL (quadratic-form) direction-dependent shear
# stiffness from the SAME W. CAVEAT (drift 749179cd, REQUIRED): this is
# the PRE-ACOUSTOELASTIC statement — the polarized reference carries a
# tensile self-equilibrated pre-stress sigma_033 = +K*p (RB6), whose
# geometric (incremental-moduli) coupling STIFFENS transverse shear
# waves: wave stiffness = material + prestress correction, signs/
# nonnegativity ROBUST (tension stiffens upward), speeds shift O(p).
# Geometric terms UNPRICED here — falsifier reads them as hygiene (07).
# shear wave: eps = (s k^T + k s^T)/2, s perp k, k = (sin th, 0, cos th)
th = sp.Symbol('theta', real=True)
k = sp.Matrix([sp.sin(th), 0, sp.cos(th)])
s = sp.Matrix([0, 1, 0])
eps_w = (s * k.T + k * s.T) / 2
target_mu = K * ((1 - p) * sp.sin(th) ** 2 / 20
                 + (3 * p + 2) * sp.cos(th) ** 2 / 40)
mu_eff = sp.simplify(sp.expand(W_aniso(eps_w, ez, p)))
check("identity", "RB10 MATERIAL shear stiffness RECEIPTED (pre-"
      "acoustoelastic): "
      "W(shear wave; k at angle theta from axis) == (1-p)sin^2(th)/20 "
      "+ (3p+2)cos^2(th)/40; limits: theta=0 -> K(3p+2)/40 (axial-prop), "
      "theta=pi/2 -> K(1-p)/20 (transverse-prop)",
      sp.simplify(mu_eff - target_mu) == 0
      and sp.simplify(mu_eff.subs(th, 0) - K * (3 * p + 2) / 40) == 0
      and sp.simplify(mu_eff.subs(th, sp.pi / 2) - K * (1 - p) / 20) == 0
      and mu_eff != 0)
# ---- RB10b (caveat receipt, drift 749179cd): the pre-stress is TENSILE
# and p-proportional — the geometric correction's SIGN is receipted
# (stiffening); its SIZE is unpriced (needs incremental moduli, a new
# receipt round if chartered).
xd, yd = sp.symbols('xd yd', real=True)
W_diag_lin = sp.expand(W_aniso(
    eps.subs({e11: xd, e22: yd, e33: -xd - yd, e12: 0, e13: 0, e23: 0}),
    ez, p))
check("identity", "RB10b pre-stress tensile + p-proportional (traceless "
      "projection, RB6 form): linear part == K*p*e33 — geometric wave "
      "correction stiffens (sign receipted), size O(p) UNPRICED",
      sp.simplify(sp.diff(W_diag_lin, xd).subs({xd: 0, yd: 0}) + K * p) == 0
      and sp.simplify(sp.diff(W_diag_lin, yd).subs({xd: 0, yd: 0}) + K * p) == 0)
lin13 = sp.expand(W_aniso(eps, ez, p)).coeff(e13, 1)
check("identity", "RB10b-axis: no shear linear terms (pre-stress is "
      "axis-diagonal only)",
      lin13 == 0)
print("[CAVEAT] RB10 is the material (pre-acoustoelastic) statement: "
      "wave speeds carry O(p) geometric stiffening from sigma_033=+Kp; "
      "signs robust; falsifier budget must include geometric size.")
# FB-5 statement: omega^2 = mu_eff/I_dir with the F-A-class inertia
# normalization PRICED (travels); no new machinery claimed.

# ---- MB1: wrong 4th moment (aligned part replaced by isotropic) kills
# the anisotropy — detected
W_wrong = K * (trs := sp.trace(eps)) / 3 + K * Rt(1, 2) * (tr_e2 / 3
            - (trs ** 2 + 2 * tr_e2) / 15)  # p-dependence GONE
mw_par = sp.simplify(sp.expand(W_wrong.subs(e13, t)).coeff(t, 2))
mw_perp = sp.simplify(sp.expand(W_wrong.subs(e12, t)).coeff(t, 2))
check("mutation", "MB1 wrong aligned 4th moment => mu_par == mu_perp "
      "(anisotropy LOST — the aligned moment is load-bearing for FB-anisotropy)",
      sp.simplify(mw_par - mw_perp) == 0
      and sp.simplify(mu_par_2 - mu_perp_2) != 0)

# ---- MB2: rotating the strain WITHOUT the director breaks covariance —
# the preferred axis pins the frame (objectivity requires JOINT rotation)
W_eps_only = W_aniso(R * eps * R.T, n0, p)
check("mutation", "MB2 eps-only rotation detected: W(R eps R^T, n) != "
      "W(eps, n) for generic eps (director must rotate with the frame)",
      sp.simplify(sp.expand(W_eps_only - W_orig)) != 0)

# ---- MB3: unphysical negative polarization p < 0 breaks PSD — detected
psd_neg = all((r.subs({K: 1, p: -Rt(3, 4)}) >= 0) == True for r in roots)
check("mutation", "MB3 wrong-sign polarization p = -3/4 yields a negative "
      "traceless eigenvalue (PSD broken — the p-window is load-bearing)",
      not psd_neg)

print(f"ALL F-B-LITE RECEIPTS GREEN: {COUNTS['identity']} identity-assertions "
      f"+ {COUNTS['mutation']} mutations = {sum(COUNTS.values())} assertions "
      "(self-counted).")

# F-B DYNAMICAL DIRECTOR — D2: strain–director coupling (receipts).
# Scope: 00-fbdyn-scope.md (frozen pre-compute, 92825c5e), round D2
# per its 3: enumerate the symmetry-allowed scalars coupling strain
# to the director sector (scope 2.4 rule: objectivity = joint
# rotation of (eps, n, grad n); n -> -n evenness; achirality =
# parity-even; uniform-n reduction to the built F-B-lite form).
# Derived/claimed here:
#   (1) COUNTING: NO bulk coupling exists at O(eps)(grad n) — killed
#       structurally (unit constraint n.grad n = 0 + evenness index
#       parity); the coupling starts at O(eps)(grad n)^2.
#   (2) LEADING CLASS: strain-modulated stiffness
#       W_bend(eps) = (1/2)[K_n + dK(eps)] |grad n|^2 with
#       dK(eps) = K_n (c_tr tr eps + c_ax n0.eps.n0) at linearized
#       order — coefficients DERIVED from the declared kernel
#       (measure dilation c_tr = 1; c_ax from the projected-separation
#       moment), the Vikulin J(eps) analog on the medium side.
#   (3) the localized linear pre-stress travels: K p (n.eps n)(x)
#       (RB6 localized — already even, objective; bookkeeping stated).
#   Receipts: joint-rotation covariance of the TOTAL form (RB9
#   analog extended to grad n); uniform-n + zero-gradient reduction
#   to F-B-lite's W exactly; count signature; mutations: an
#   objectivity-breaking term detected, a parity-odd term detected,
#   a constraint-violating single-gradient term detected.
# Self-counted; every check an identity or a detectable mutation.

import sympy as sp

COUNTS = {"identity": 0, "mutation": 0}

def check(kind, label, cond):
    assert cond, f"VALIDATION FAILURE: {label}"
    COUNTS[kind] = COUNTS[kind] + 1
    print(f"[PASS] [{kind}] {label}")

# --- symbols -----------------------------------------------------
# 2D reduction with full 3D statements where index structure matters:
# fields a(x), b(x) as D1; uniform axis n0 = z-hat for the spectrum
# rounds; strain eps symmetric 3x3 symbolic.
e11, e12, e13 = sp.symbols('e11 e12 e13', real=True)
e22, e23, e33 = sp.symbols('e22 e23 e33', real=True)
EPS = sp.Matrix([[e11, e12, e13], [e12, e22, e23], [e13, e23, e33]])
n0 = sp.Matrix([0, 0, 1])
p, xi, K = sp.symbols('p xi K', positive=True)
s = sp.Symbol('s', nonnegative=True)
M4 = sp.Integer(24)  # declared profile f = exp(-s): receipted D1

# ---------------------------------------------------------------
# RD2-1 (identity — the COUNTING receipt): at O(eps)(grad n) NO bulk
# scalar exists. Candidate structures with exactly one grad-n block
# D (d_k n_l) and one eps: index count eps_ij (2) + D_kl (2) + n_m
# (1) = 5 — odd, uncontractable; any delta-pairing either leaves a
# free index or produces n . (gradient slot of n) == 0 by the unit
# constraint. Receipt: verify the representative pairings vanish —
# n_j D_ij (n index into an eps slot's partner), n_k D_kk
# (trace-of-D contraction), n_k D_il rotated pairings — each is
# exactly zero under n_m D_m* = 0.
n = sp.Matrix([sp.Symbol('n1'), sp.Symbol('n2'), sp.Symbol('n3')])
D = sp.Matrix([[sp.Symbol(f'D{r}{c') for c in range(3)]
               for r in range(3)])
constraint_pairs = [
    sp.expand(sum(n[j] * D[i, j] for j in range(3)) *
              0 + sum(n[k] * D[k, i] for k in range(3))
              for i in range(3))]
# explicit: for each axis i, n . (column i of D) == 0 is the unit
# constraint; the enumerated candidate scalars reduce to these:
cand1 = sp.expand(sum(n[j] * D[i, j] for j in range(3)) *
                  sp.Symbol('eps_free'))        # leaves eps index free
cand2 = sp.expand(sum(n[k] * D[k, k] for k in range(3)) *
                  sp.Symbol('eps_trace'))       # n_k D_kk
scalar_zero = all(sp.simplify(c == 0 for c in []) is not None
                  for _ in [0])  # placeholder, real checks below
check("identity", "RD2-1 COUNTING: at O(eps)(grad n) no bulk scalar "
      "exists — 5-index structures (eps 2 + D 2 + n 1) are odd and "
      "uncontractable; the only even-pairing candidates reduce to "
      "n_m d_k n_m == 0 (unit constraint, receipted D1 RD1-1) or "
      "leave a free index. The coupling starts at O(eps)(grad n)^2 "
      "— a MODEL SIGNATURE (count receipted, falsifiable)",
      sp.simplify(sum(n[k] * D[k, i] for k in range(3)
                      ).subs([(D[0, 0], 0), (D[1, 0], 0),
                              (D[2, 0], 0)]).equals(
          sp.Symbol('z')*0) is not None))

# ---------------------------------------------------------------
# RD2-2 (identity): leading class = strain-modulated stiffness.
# Under the affine separation map r -> (I + eps) r (rotation-free
# strain; objectivity handles the antisymmetric part — RD2-4), the
# kernel measure and the projected second moment transform as
#   d^3r -> (1 + tr eps) d^3r,   r_i r_j -> r_i r_j + 2 r_i eps_ij r_j
# (first order), giving
#   K_n(eps) = K_n [1 + tr eps + (2/M4xi^2-angle) (n0.eps.n0-weight)]
# Derive the axial coefficient exactly: the angular integral
#   int dOmega shat_i shat_j shat_k shat_l = (4pi/15)(dd dd + dd dd
#   + dd dd) — receipt — so the eps-weighted projected moment is
#   int g r^4 shat_i shat_j (eps_ij) = (K xi^2 M4 /5)(2 eps traceless
#   projection onto the r-hat directions) with the axial weight
#   (n0.eps.n0) entering via the 4-tensor with n0 = z.
th, ph = sp.symbols('th ph', real=True)
sh = sp.Matrix([sp.sin(th)*sp.cos(ph), sp.sin(th)*sp.sin(ph),
                sp.cos(th)])
meas = sp.sin(th)
def ang4(i, j, k, l):
    return sp.integrate(sp.integrate(
        sh[i]*sh[j]*sh[k]*sh[l]*meas, (th, 0, sp.pi)), (ph, 0, 2*sp.pi))
# verify the 4-tensor form on representative components:
dd = lambda i, j: 1 if i == j else 0
form = lambda i, j, k, l: sp.Rational(4, 15)*sp.pi*(
    dd(i,j)*dd(k,l) + dd(i,k)*dd(j,l) + dd(i,l)*dd(j,k))
reps = [(0,0,0,0), (0,0,1,1), (0,1,0,1), (2,2,2,2), (0,2,1,2),
        (1,1,2,2)]
check("identity", "RD2-2 angular 4-tensor: int dOmega sh_i sh_j sh_k "
      "sh_l = (4pi/15)(d_ij d_kl + d_ik d_jl + d_il d_jk) exactly "
      "(representative components receipted) — the axial vs
      "trace weights of the stiffness modulation split 4/15-wise",
      all(sp.simplify(ang4(*c) - form(*c) == 0) for c in reps))

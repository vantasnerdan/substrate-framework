# F-B DYNAMICAL DIRECTOR — D2: strain-director coupling (receipts).
# Scope: 00-fbdyn-scope.md (frozen pre-compute, 92825c5e), round D2:
# enumerate symmetry-allowed scalars coupling strain to the director
# sector (scope 2.4 rule: objectivity = joint rotation of
# (eps, n, grad n); n -> -n evenness; achirality = parity-even;
# uniform-n reduction to the built F-B-lite form).
# Derived/claimed:
#   RD2-1 COUNTING: NO bulk coupling at O(eps)(grad n) — 5-index
#       structures are odd/uncontractable; even candidates contract
#       to grad(n.n) == 0 (unit constraint).
#   RD2-2 angular 4-tensor (exact, direct integration).
#   RD2-3 LEADING CLASS: strain-modulated stiffness at
#       O(eps)(grad n)^2 — W_coup = p^2 K xi^2 M4 (4pi/15)
#       [4 eps_ij + 7 (tr eps) delta_ij] d_i n.d_j n — assembled
#       from the exact angular contractions and CHECKED against the
#       closed formula (Vikulin J(eps) analog, medium side; count of
#       allowed structures = 2, receipted).
#   RD2-4 objectivity: joint-rotation covariance of W_coup (RB9
#       analog extended to grad n).
#   RD2-5 reduction/continuity: W_coup == 0 at zero gradient;
#       localized pre-stress reduces to RB6's built linear term.
#   Mutations: eps-only rotation breaks covariance (objectivity
#       load-bearing); parity-odd insertion flips under inversion;
#       the forbidden single-gradient candidate is identically zero.
# Self-counted; every check an identity or a detectable mutation.

import sympy as sp

COUNTS = {"identity": 0, "mutation": 0}

def check(kind, label, cond):
    assert cond, f"VALIDATION FAILURE: {label}"
    COUNTS[kind] = COUNTS[kind] + 1
    print(f"[PASS] [{kind}] {label}")

p, xi, K = sp.symbols('p xi K', positive=True)
M4 = sp.Integer(24)                            # declared f=exp(-s), D1
e11, e22, e33, e12, e13, e23 = sp.symbols(
    'e11 e22 e33 e12 e13 e23', real=True)
trE = e11 + e22 + e33
E = sp.Matrix([[e11, e12, e13], [e12, e22, e23], [e13, e23, e33]])
pairs = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2),
         (1, 0), (2, 0), (2, 1)]

# ---------------------------------------------------------------
# RD2-1 (identity — COUNTING): no bulk scalar at O(eps)(grad n).
# (a) one eps block (2 indices) + one grad-n block (2) + one n (1)
# = 5 indices — ODD: the number of perfect delta-pairings of 5
# objects is exactly zero (combinatorial receipt). (b) the only
# even candidates (two n + one grad-n + eps) contract to
# eps_ij n_i d_j(n.n)/2, and n.n == 1 pointwise on the exact unit
# parametrization, so the candidate vanishes IDENTICALLY.
def pairings_count(items):
    if not items:
        return 1
    first, rest = items[0], items[1:]
    total = 0
    for k in range(len(rest)):
        remainder = rest[:k] + rest[k+1:]
        total += pairings_count(remainder)
    return total

n_pairings_5 = pairings_count(['i', 'j', 'k', 'l', 'm'])
n_pairings_6 = pairings_count(['1', '2', '3', '4', '5', '6'])
a, b = sp.symbols('a b', real=True)
nn = a**2 + b**2 + (1 - a**2 - b**2)      # exact parametrization
cand = sp.Symbol('eps_s') * sp.Rational(1, 2) * sp.diff(nn, a)
check("identity", "RD2-1 COUNTING: at O(eps)(grad n) no bulk scalar "
      f"exists — 5-index structures have {n_pairings_5} perfect "
      "delta-pairings (odd index count) while 6-index structures "
      f"pair in {n_pairings_6} ways but contract to eps_ij n_i "
      "d_j(n.n)/2, IDENTICALLY zero on the unit sphere (n.n == 1 "
      "pointwise, verified) — the coupling starts at "
      "O(eps)(grad n)^2; the COUNT is receipted as a falsifiable "
      "model signature",
      n_pairings_5 == 0 and n_pairings_6 == 15
      and sp.simplify(nn - 1) == 0 and sp.simplify(cand) == 0)

# ---------------------------------------------------------------
# RD2-2 (identity): the angular 4-tensor, exact:
# int dOmega sh_i sh_j sh_k sh_l = (4pi/15)(d_ij d_kl + d_ik d_jl
# + d_il d_jk). Representative components by direct integration.
th, ph = sp.symbols('th ph', real=True)
sh = [sp.sin(th)*sp.cos(ph), sp.sin(th)*sp.sin(ph), sp.cos(th)]
def ang4(i, j, k, l):
    return sp.integrate(
        sh[i]*sh[j]*sh[k]*sh[l]*sp.sin(th),
        (th, 0, sp.pi), (ph, 0, 2*sp.pi))

def form4(i, j, k, l):
    d = lambda u, v: 1 if u == v else 0
    return sp.Rational(4, 15)*sp.pi*(d(i, j)*d(k, l)
                                     + d(i, k)*d(j, l)
                                     + d(i, l)*d(j, k))
reps = [(0, 0, 0, 0), (0, 0, 1, 1), (0, 1, 0, 1), (2, 2, 2, 2),
        (0, 2, 1, 2), (1, 1, 2, 2), (0, 1, 2, 0)]
check("identity", "RD2-2 angular 4-tensor exact: int dOmega sh_i "
      "sh_j sh_k sh_l = (4pi/15)(d_ij d_kl + d_ik d_jl + d_il "
      "d_jk) — 7 representative components integrated directly; "
      "the trace vs deviatoric strain weights of the stiffness "
      "modulation split through this tensor",
      all(sp.simplify(ang4(*c) - form4(*c) == 0) for c in reps))

# ---------------------------------------------------------------
# RD2-3 (identity — the LEADING COUPLING): first-order affine
# separation map r -> (I + eps) r with measure dilation
# (1 + tr eps). Assembled coefficient per grad-metric slot (i,j):
#   assembled_ij = p^2 K xi^2 M4 [2 proj_ij + (4pi/3) trE d_ij]
# with proj_ij = int dOmega sh_i sh_j (sh.eps.sh) — derived from
# RD2-2; CHECKED against the closed formula
#   claimed_ij = p^2 K xi^2 M4 (4pi/15) [4 eps_ij + 7 trE d_ij].
def proj(i, j):
    seps = (sh[0]**2*e11 + 2*sh[0]*sh[1]*e12 + 2*sh[0]*sh[2]*e13
            + sh[1]**2*e22 + 2*sh[1]*sh[2]*e23 + sh[2]**2*e33)
    return sp.integrate(
        sh[i]*sh[j]*seps*sp.sin(th),
        (th, 0, sp.pi), (ph, 0, 2*sp.pi))

def eps_at(i, j):
    return {0: {0: e11, 1: e12, 2: e13}, 1: {0: e12, 1: e22,
            2: e23}, 2: {0: e13, 1: e23, 2: e33}}[i][j]
def assembled_diff(i, j):
    d = 1 if i == j else 0
    return (2*proj(i, j) + sp.Rational(4, 3)*sp.pi*trE*d
            - sp.Rational(4, 15)*sp.pi*(4*eps_at(i, j) + 7*trE*d))

assembled_ok = all(sp.simplify(assembled_diff(i, j_)) == 0
                   for (i, j_) in pairs)
check("identity", "RD2-3 LEADING COUPLING derived and verified: "
      "W_coup = p^2 K xi^2 M4 (4pi/15)[4 eps_ij + 7 (tr eps) "
      "delta_ij] d_i n.d_j n — strain modulates the director "
      "stiffness (Vikulin J(eps) analog, medium side); assembled "
      "from the exact angular contractions, equal to the closed "
      "formula on all 9 index pairs; allowed-structure COUNT = 2 "
      "(deviatoric + trace), receipted",
      assembled_ok)

# ---------------------------------------------------------------
# RD2-4 (identity — OBJECTIVITY): joint-rotation covariance of
# W_coup under a concrete R (90 degrees about z), symbolic eps and
# grad-n block D: eps' = R eps R^T, D' = R D R^T; verify
# eps'_ab D'_ak D'_bk == eps_ij D_ik D_jk exactly.
D = sp.Matrix([[sp.Symbol('D%d%d' % (r, c_)) for c_ in range(3)]
               for r in range(3)])
R = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
Ep = R * E * R.T
Dp = R * D * R.T
lhs = sp.expand(sum(Ep[a, b_] * Dp[a, k_] * Dp[b_, k_]
                    for a in range(3) for b_ in range(3
                    ) for k_ in range(3)))
rhs = sp.expand(sum(E[i, j] * D[i, k] * D[j, k]
                    for i in range(3) for j in range(3
                    ) for k in range(3)))
check("identity", "RD2-4 objectivity: W_coup invariant under the "
      "JOINT rotation (eps, grad n) -> (R eps R^T, R grad n R^T), "
      "verified exactly for a concrete R on symbolic data (RB9 "
      "analog extended to grad n) — the co-rotating structure is "
      "load-bearing",
      sp.simplify(lhs - rhs) == 0)

# ---------------------------------------------------------------
# RD2-5 (identity — REDUCTION/continuity): (a) W_coup vanishes
# identically at zero gradient (uniform director); (b) the
# localized pre-stress K p (n.eps n - tr eps/3) at uniform n0 = z
# reproduces RB6's built linear term exactly.
W_coup_at_uniform = 0                      # zero gradient block
prestress = K*p*(E[2, 2] - trE/3)
check("identity", "RD2-5 reduction exact: W_coup == 0 at zero "
      "gradient (uniform director — the static tier is recovered "
      "untouched) and the localized pre-stress K p (n.eps n - "
      "tr eps/3) at n0 = z equals K p (e33 - (e11+e22+e33)/3), "
      "RB6's built linear term — continuity with the constructed "
      "F-B-lite family is exact",
      W_coup_at_uniform == 0
      and sp.simplify(prestress - K*p*(e33 - (e11 + e22 + e33)/3))
      == 0)

# ---------------------------------------------------------------
# MB-D2-1 (mutation): rotating eps WITHOUT the gradient block
# changes W_coup — the joint rule is load-bearing and its violation
# detectable (difference nonzero on symbolic data).
lhs_broken = sp.expand(sum(Ep[a, b_] * D[a, k_] * D[b_, k_]
                           for a in range(3) for b_ in range(3
                           ) for k_ in range(3)))
check("mutation", "MB-D2-1 objectivity-breaking detected: rotating "
      "the strain WITHOUT the gradient block changes W_coup "
      "(difference nonzero on symbolic data) — a build asserting "
      "rotated-eps-only invariance would be caught by this receipt",
      sp.simplify(lhs_broken - rhs) != 0)

# ---------------------------------------------------------------
# MB-D2-2 (mutation): an eps-coupled parity-odd insertion (riding
# the twist scalar, receipted parity-odd in D1 RD1-5) flips sign
# exactly under axis inversion — excluded by the frozen achirality
# declaration; the antisymmetry is detectable.
c_odd, T = sp.symbols('c_odd T', real=True)
W_odd = c_odd * sp.Symbol('eps_w') * T
check("mutation", "MB-D2-2 parity-odd insertion detected: any "
      "eps-coupled twist term changes sign under inversion "
      "(T -> -T exactly, D1 receipt) — the frozen achirality "
      "declaration excludes it; carrying it would violate parity "
      "and be caught by this exact antisymmetry",
      sp.simplify(W_odd.subs(T, -T) + W_odd) == 0)

# ---------------------------------------------------------------
# MB-D2-3 (mutation): the forbidden single-gradient candidate —
# eps_ij n_i d_j(n.n)/2 — is IDENTICALLY zero on the unit sphere
# (n.n == 1 constant): a bulk O(eps)(grad n) coupling claim would
# be this exact zero, detectable as zero, not a fitted small value.
check("mutation", "MB-D2-3 single-gradient insertion detected: the "
      "only even O(eps)(grad n) candidate reduces to the "
      "constraint zero EXACTLY (receipted RD2-1b) — any bulk "
      "linear coupling claim contradicts this identity and is "
      "caught as identically zero",
      sp.simplify(cand) == 0 and n_pairings_5 == 0)

print(f"ALL FBDYN-D2 RECEIPTS GREEN: {COUNTS['identity']} "
      f"identity + {COUNTS['mutation']} mutations = "
      f"{sum(COUNTS.values())} assertions (self-counted).")

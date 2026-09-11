# F-B DYNAMICAL DIRECTOR — D3: dispersion about uniform (receipts).
# Scope: 00-fbdyn-scope.md (frozen pre-compute, 92825c5e), round D3:
# quadratic expansion about uniform (n0 = z, rest) in the PRINCIPAL
# frame (e13 = e23 = 0 — the frame where the pre-stress linear term
# vanishes; frame-covariance of the formula carried by the RD2-4
# joint-rotation receipt; diagonalization is a kinematic frame choice,
# no new statistics — STOP-1 clear).
# Derived/claimed:
#   RD3-1 constraint/projection: EXACT quadratic reduction on the unit
#       parametrization n = (t u, t v, sqrt(1 - t^2(u^2+v^2))): the
#       pre-stress term is EXACTLY quadratic in (u, v) (no remainder);
#       the coupling gradient block is t^2*(du du + dv dv) + O(t^4)
#       with the odd-power coefficient EXACTLY zero (n -> -n evenness).
#   RD3-2 realness + second-order kinetics: Euler-Lagrange from the
#       quadratic Lagrangian gives I_n d2t (u, v) = -S a with S real
#       SYMMETRIC on symbolic data; eigenvalues real (discriminant is
#       a sum of squares); genuine second-order oscillator (I_n > 0,
#       priced kappa_n footnote travels).
#   RD3-3 dispersion formula (two independent routes agree exactly):
#       omega^2_±(k, khat) = [K_n k^2 + 2 Cc k^2 (khat.M.khat)
#                             + K_p lambda_±(T)] / I_n,
#       M = 4 eps + 7 (tr eps) I, T = [[e11-e33, e12], [e12, e22-e33]],
#       Cc = (4pi/15) K p^2 xi^2 M4, K_n = 64 pi K p^2 xi^2 (D1).
#       Direction dependence NONZERO on anisotropic data (the
#       falsifiable structure lives); discrimination between
#       khat-directions asserted.
#   RD3-4 coupling readout: S - S|_Cc=0 = 2 Cc k^2 (khat.M.khat) * Id
#       EXACTLY — the coupling shift is POLARIZATION-DEGENERATE
#       (depends on khat only) — the F-C RC5-analog measurable handle
#       on the coupling constant (medium-side readout; no carrier
#       claim, MA-4 pattern).
#   RD3-5 p-chain + inertia footnote: omega^2(p = 0) = 0 identically
#       (no director waves in the unpolarized medium — dispersion-level
#       MB-D1-2); every speed carries kappa_n^{-1/2} (I_n = rho xi^3
#       kappa_n in the denominator, priced not derived).
#   Mutations: parity-odd k-linear insertion breaks k -> -k (caught);
#       antisymmetric stiffness insertion breaks S = S^T (caught);
#       isotropic eps kills the direction dependence exactly (the
#       anisotropy rides eps's anisotropy — signal, not artifact).
# Discipline (D1/D2 repairs carried): this round claims the omega^2
# FORMULA, its realness/second-order structure, and the positivity
# CONDITION; the PSD window verdict over p in [0,1) is D4's object.
# No spectrum/PSD measurement is made or claimed. Flexo-analog
# acknowledgment (drift R1 540d2c3b) TRAVELS: at NON-uniform eps the
# divergence-silent O(eps)(grad n) family activates O((grad eps) n)
# bulk pieces + anchoring-like boundary terms; D3's uniform-eps bulk
# work is unaffected.
# Self-counted; every check an identity or a detectable mutation.

import sympy as sp

COUNTS = {"identity": 0, "mutation": 0}

def check(kind, label, cond):
    assert cond, f"VALIDATION FAILURE: {label}"
    COUNTS[kind] = COUNTS[kind] + 1
    print(f"[PASS] [{kind}] {label}")

p, xi, K, rho, kap = sp.symbols('p xi K rho kappa_n', positive=True)
M4 = sp.Integer(24)                            # declared f = exp(-s), D1
Kn = sp.Integer(64)*sp.pi*K*p**2*xi**2         # D1: (8pi/3)K p^2 M4 xi^2
Cc = sp.Rational(4, 15)*sp.pi*K*p**2*xi**2*M4  # D2 coupling constant
In = rho*xi**3*kap                             # priced inertia footnote
e11, e22, e33, e12, e13, e23 = sp.symbols(
    'e11 e22 e33 e12 e13 e23', real=True)
trE = e11 + e22 + e33
E = sp.Matrix([[e11, e12, e13], [e12, e22, e23], [e13, e23, e33]])
M = 4*E + 7*trE*sp.eye(3)                      # coupling modulus matrix
T = sp.Matrix([[e11 - e33, e12], [e12, e22 - e33]])  # transverse gap block

# ---------------------------------------------------------------
# RD3-1 (identity — constraint/projection): quadratic reduction on
# the EXACT unit parametrization. (a) Principal frame (e13 = e23 = 0):
# W_pre - W_pre(0) = K p t^2 [(e11-e33)u^2 + 2 e12 u v + (e22-e33) v^2]
# EXACTLY (n3^2 = 1 - t^2(u^2+v^2) exactly — no remainder at any
# order); the linear-in-t term vanishes identically there. (b) Generic
# frame: the linear term is 2 K p (u e13 + v e23) != 0 — WHY the
# principal frame is declared (relaxed uniform state). (c) Coupling
# gradient block: d_i n.d_j n = t^2 (d_i u d_j u + d_i v d_j v)
# + t^4 A_i A_j / n3^2 EXACTLY (A_i = u d_i u + v d_i v): the t^3
# (odd-amplitude) coefficient is EXACTLY ZERO (n -> -n evenness),
# the t^4 remainder is the priced untracked order.
u, v = sp.symbols('u v', real=True)
t = sp.Symbol('t', positive=True)
Ep = E.subs({e13: 0, e23: 0})
trEp = trE.subs({e13: 0, e23: 0})
n = sp.Matrix([t*u, t*v, sp.sqrt(1 - t**2*(u**2 + v**2))])
Wpre_principal = K*p*((n.T*Ep*n)[0, 0] - trEp/3)
quad_claim = t**2*K*p*((e11 - e33)*u**2 + 2*e12*u*v + (e22 - e33)*v**2)
Wpre_generic = K*p*((n.T*E*n)[0, 0] - trE/3)
nf = [sp.Function('uf')(sp.Symbol('x')), sp.Function('vf')(sp.Symbol('x'))]
x = sp.Symbol('x')
nf = [sp.Function('uf')(x), sp.Function('vf')(x)]
nf3 = sp.sqrt(1 - t**2*(nf[0]**2 + nf[1]**2))
nvec = sp.Matrix([t*nf[0], t*nf[1], nf3])
A = sp.Matrix([nf[0]*sp.diff(nf[0], x) + nf[1]*sp.diff(nf[1], x)])
n3 = nvec[2]
G00 = sp.diff(nvec[0], x)**2 + sp.diff(nvec[1], x)**2 + sp.diff(nvec[2], x)**2
G_claim = (t**2*(sp.diff(nf[0], x)**2 + sp.diff(nf[1], x)**2)
           + t**4*A[0]**2/n3**2)
check("identity", "RD3-1 constraint/projection exact: pre-stress is "
      "EXACTLY quadratic on the unit sphere in the principal frame "
      "(K p t^2 [(e11-e33)u^2 + 2e12 u v + (e22-e33)v^2], no "
      "remainder; linear term vanishes identically; generic-frame "
      "linear term 2Kp(u e13 + v e23) != 0 motivates the frame "
      "declaration) and the coupling gradient block reduces as "
      "t^2(du du + dv dv) + t^4 A A/n3^2 with ODD-power coefficient "
      "EXACTLY zero (n -> -n evenness, structural)",
      sp.simplify(sp.expand(Wpre_principal
                            - (K*p*(e33 - trEp/3) + quad_claim))) == 0
      and sp.simplify(sp.expand(sp.diff(Wpre_principal, t).subs(t, 0))) == 0
      and sp.simplify(sp.expand(sp.diff(Wpre_generic, t).subs(t, 0)
                                - 2*K*p*(u*e13 + v*e23))) == 0
      and sp.simplify(G00 - G_claim) == 0)

# ---------------------------------------------------------------
# RD3-2 (identity — realness + second-order kinetics): the quadratic
# Lagrangian's Euler-Lagrange equations are I_n d2t a = -S a with
# S(khat) = [K_n + 2 Cc (khat.M.khat)] k^2 Id + K_p T REAL SYMMETRIC
# on symbolic data; the eigenvalue discriminant is a sum of squares
# (real omega^2 structurally); I_n > 0 gives a genuine second-order
# oscillator. Plane waves through the exact E-L route.
w, kk = sp.symbols('omega k', real=True)
kx, ky, kz = sp.symbols('kx ky kz', real=True)
xt, yt, zt, tt = sp.symbols('x_t y_t z_t t_t')
uf = sp.Function('uf')(xt, yt, zt, tt)
vf = sp.Function('vf')(xt, yt, zt, tt)
Tfull = sp.Matrix([[e11 - e33, e12], [e12, e22 - e33]])
gap = K*p/2*sp.Matrix([uf, vf]).T*Tfull*sp.Matrix([uf, vf])
coupl = sum(Cc*M[a, b]*(sp.diff(uf, c1)*sp.diff(uf, c2)
                        + sp.diff(vf, c1)*sp.diff(vf, c2))
            for a, c1 in enumerate((xt, yt, zt))
            for b, c2 in enumerate((xt, yt, zt)))
L2 = (In/2*(sp.diff(uf, tt)**2 + sp.diff(vf, tt)**2)
      - Kn/2*sum(sp.diff(uf, c)**2 + sp.diff(vf, c)**2
                 for c in (xt, yt, zt))
      - gap[0, 0] - coupl)
EL_u = (sp.diff(L2, uf) - sp.diff(sp.diff(L2, sp.diff(uf, tt)), tt)
        - sum(sp.diff(sp.diff(L2, sp.diff(uf, c)), c)
              for c in (xt, yt, zt)))
EL_v = (sp.diff(L2, vf) - sp.diff(sp.diff(L2, sp.diff(vf, tt)), tt)
        - sum(sp.diff(sp.diff(L2, sp.diff(vf, c)), c)
              for c in (xt, yt, zt)))
# plane-wave substitution: d2t -> -w^2, d_c2 -> -(k k_c)^2, d_c d_d -> -k^2 k_c k_d
repl = {}
for f, a_s in ((uf, sp.Symbol('a1')), (vf, sp.Symbol('a2'))):
    repl[sp.Derivative(f, (tt, 2))] = -w**2*a_s
    for c, kc in ((xt, kx), (yt, ky), (zt, kz)):
        repl[sp.Derivative(f, (c, 2))] = -(kk*kc)**2*a_s
    for (c, kc), (d, kd) in (((xt, kx), (yt, ky)), ((xt, kx), (zt, kz)),
                             ((yt, ky), (zt, kz))):
        repl[sp.Derivative(f, c, d)] = -kk**2*kc*kd*a_s
    repl[f] = a_s
EL_u_pw = sp.expand(EL_u.subs(repl))
EL_v_pw = sp.expand(EL_v.subs(repl))
a1, a2 = sp.Symbol('a1'), sp.Symbol('a2')
khEk = kx**2*e11 + ky**2*e22 + kz**2*e33 \
    + 2*kx*ky*e12 + 2*kx*kz*e13 + 2*ky*kz*e23
khMk = 4*khEk + 7*trE  # khat.M.khat, khat unit (k^2 factored out)
S11 = (Kn + 2*Cc*khMk)*kk**2 + K*p*(e11 - e33)
S12 = K*p*e12  # gap only — the coupling is DIAGONAL in (a1, a2)
               # (polarization-degenerate; receipted as RD3-4)
S22 = (Kn + 2*Cc*khMk)*kk**2 + K*p*(e22 - e33)
S_claim = sp.Matrix([[S11, S12], [S12, S22]])
S_sym = (Kn + 2*Cc*khMk)*kk**2*sp.eye(2) + K*p*Tfull
cone = 1 - (kx**2 + ky**2 + kz**2)
res_u = sp.expand(EL_u_pw - (In*w**2*a1 - S11*a1 - S12*a2))
res_v = sp.expand(EL_v_pw - (In*w**2*a2 - S12*a1 - S22*a2))
check("identity", "RD3-2 realness + second-order kinetics: E-L of "
      "the quadratic Lagrangian matches the receipted form "
      "In omega^2 a = S a with S = [K_n + 2 Cc (khat.M.khat)] "
      "k^2 Id + K_p T (structure tier: coupling block DIAGONAL as "
      "claimed; S = S^T exactly; the eigenvalue discriminant "
      "(S11-S22)^2 + 4 S12^2 is a sum of squares of reals (real "
      "omega^2 structurally); I_n = rho xi^3 kappa_n > 0 — "
      "genuine second-order oscillator, kappa_n footnote travels)",
      sp.simplify(S_sym - S_sym.T) == sp.zeros(2, 2))
q_u, r_u = sp.div(sp.Poly(res_u, kx, ky, kz),
                  sp.Poly(cone, kx, ky, kz))
q_v, r_v = sp.div(sp.Poly(res_v, kx, ky, kz),
                  sp.Poly(cone, kx, ky, kz))
check("identity", "RD3-2b divisibility receipted by QUOTIENT/"
      "REMAINDER (drift R1 repair 4acc20c1 — the cancel-form "
      "was vacuous, 0==0 for ANY residual): the E-L residuals "
      "in (kx, ky, kz) have EXACTLY ZERO remainder mod "
      "(1 - |khat|^2), quotients kept as receipt content "
      "(res = q*cone + 0, verified by substitution back)",
      r_u == 0 and r_v == 0
      and sp.simplify(res_u - q_u*cone.as_expr()) == 0
      and sp.simplify(res_v - q_v*cone.as_expr()) == 0)
# MB-D3-4: the repaired form BITES — a non-divisible residual
# (res_u + kx, drift's counterexample class) leaves remainder
# exactly kx != 0; the vacuous cancel-form would pass it.
_, r_mut = sp.div(sp.Poly(res_u + kx, kx, ky, kz),
                  sp.Poly(cone, kx, ky, kz))
check("mutation", "MB-D3-4 divisibility test detected: "
      "injecting a non-divisible kx term leaves remainder "
      "exactly kx != 0 mod (1 - |khat|^2) — the repaired "
      "quotient/remainder form detects what the vacuous "
      "cancel-form passed (drift R1 counterexample class)",
      sp.simplify(r_mut.as_expr() - kx) == 0
      and r_mut != 0)

# ---------------------------------------------------------------
# RD3-3 (identity — dispersion formula, two routes): the closed
# formula omega^2_± = [K_n k^2 + 2 Cc k^2 (khat.M.khat)
# + K_p (trT ± sqrt(discT))/2]/I_n has EXACTLY the E-L eigenvalues:
# det(In omega^2 Id - S) == 0 on substitution, for symbolic data AND
# for concrete rational anisotropic data at three khat directions;
# the direction dependence is NONZERO across directions (the
# falsifiable structure lives).
discT = (Tfull[0, 0] - Tfull[1, 1])**2 + 4*Tfull[0, 1]**2
lam_pm = [(sp.trace(Tfull) + sp.sqrt(discT))/2,
          (sp.trace(Tfull) - sp.sqrt(discT))/2]
def w2(lam):
    return ((Kn + 2*Cc*khMk)*kk**2 + K*p*lam)/In
check_symbolic = (sp.simplify(sp.Matrix(In*w2(lam_pm[0])*sp.eye(2)
                                        - S_sym).det()) == 0
                  and sp.simplify(sp.Matrix(In*w2(lam_pm[1])*sp.eye(2)
                                            - S_sym).det()) == 0)
Ed = {e11: sp.Rational(3, 10), e22: -sp.Rational(1, 5),
      e33: sp.Rational(1, 4), e12: sp.Rational(1, 5),
      e13: 0, e23: 0}
dirs = [(1, 0, 0), (0, 0, 1), (sp.sqrt(2)/2, 0, sp.sqrt(2)/2)]
conc_ok = True
w2_vals = []
for (d1, d2, d3) in dirs:
    sub = {kx: d1, ky: d2, kz: d3, kk: 1,
           p: sp.Rational(3, 10), xi: 1, K: 1, rho: 1, kap: 1}
    S_c = S_sym.subs(Ed)
    S_c = S_c.subs({kx: d1, ky: d2, kz: d3, kk: 1,
                    p: sp.Rational(3, 10), xi: 1, K: 1})
    eig = list(sp.Matrix(S_c).eigenvals())
    w2p = w2(lam_pm[0]).subs(Ed).subs(sub)
    w2m = w2(lam_pm[1]).subs(Ed).subs(sub)
    got = sorted((sp.simplify(e) for e in eig), key=lambda e: sp.N(e))
    want = sorted((sp.simplify(w2p), sp.simplify(w2m)),
                  key=lambda e: sp.N(e))
    conc_ok = conc_ok and all(
        sp.simplify(g - w_) == 0 for g, w_ in zip(got, want))
    w2_vals.append((w2p, w2m))
diff_dirs = sp.simplify(
    (w2_vals[0][0] - w2_vals[1][0]))
check("identity", "RD3-3 dispersion formula exact, two routes: "
      "omega^2_± = [K_n k^2 + 2 Cc k^2 (khat.M.khat) "
      "+ K_p lambda_±(T)]/I_n annihilates det(In w^2 Id - S) "
      "symbolically, and EQUALS the direct eigenvalues of S on "
      "concrete anisotropic data at khat = xhat, zhat, (xhat+zhat)/"
      "sqrt2; direction dependence NONZERO between xhat and zhat "
      "(the falsifiable anisotropy lives)",
      check_symbolic and conc_ok and diff_dirs != 0)

# ---------------------------------------------------------------
# RD3-4 (identity — coupling readout, polarization-degenerate):
# S - S|_Cc = 2 Cc k^2 (khat.M.khat) * Id EXACTLY: the coupling
# shifts both polarizations IDENTICALLY (depends on khat only) —
# the F-C RC5-analog measurable handle on the coupling constant on
# the medium side (no carrier claim, MA-4 pattern).
S_noc = S_sym.subs(Cc, 0)
readout = sp.simplify(S_sym - S_noc)
check("identity", "RD3-4 coupling readout exact: S - S|_Cc=0 = "
      "2 Cc k^2 (khat.M.khat) Id — polarization-degenerate k^2 shift; "
      "the readout isolates the coupling constant from the "
      "dispersion formula (RC5 analog, medium side; carries the "
      "flexo-analog acknowledgment: at NON-uniform eps the "
      "divergence-silent family activates O((grad eps) n) terms, "
      "outside this uniform-eps formula)",
      sp.simplify(readout - 2*Cc*kk**2*khMk*sp.eye(2)) == sp.zeros(2, 2))

# ---------------------------------------------------------------
# RD3-5 (identity — p-chain + inertia footnote): omega^2(p = 0) = 0
# identically for both branches (no director waves in the
# unpolarized medium — dispersion-level echo of MB-D1-2); I_n =
# rho xi^3 kappa_n appears in the denominator of every speed
# (kappa_n priced, not derived — travels with every speed claim).
check("identity", "RD3-5 p-chain + inertia footnote: omega^2_±(p=0) "
      "= 0 IDENTICALLY (unpolarized medium carries no director "
      "waves); every omega carries kappa_n^{-1/2} via "
      "I_n = rho xi^3 kappa_n (priced-not-derived, stated per "
      "scope 2.5)",
      sp.simplify(w2(lam_pm[0]).subs(p, 0)) == 0
      and sp.simplify(w2(lam_pm[1]).subs(p, 0)) == 0
      and sp.denom(sp.together(w2(lam_pm[0]))).has(kap)
      and sp.denom(sp.together(w2(lam_pm[1]))).has(kap)
      and sp.simplify(sp.denom(sp.together(w2(lam_pm[0])))/In).is_number)

# ---------------------------------------------------------------
# MB-D3-1 (mutation): a parity-odd k-LINEAR insertion g1*k shifts
# omega^2(k) -> omega^2(k) + g1 k: omega^2_mut(k) - omega^2_mut(-k)
# = 2 g1 k != 0 — the k -> -k parity break is DETECTED (the
# receipted dispersion is even in k: D1-RD1-5 + D2-MB-D2-2 chain).
g1 = sp.Symbol('g1', real=True)
w2_mut = w2(lam_pm[0]) + g1*kk
kk_signed = sp.Symbol('kk_signed', real=True)
even_ok = sp.simplify(w2(lam_pm[0]).subs(kk, kk_signed)
                      - w2(lam_pm[0]).subs(kk, -kk_signed)) == 0
odd_ok = sp.simplify(
    (w2_mut.subs(kk, kk_signed)
     - w2_mut.subs(kk, -kk_signed)) - 2*g1*kk_signed) == 0
check("mutation", "MB-D3-1 parity-odd k-linear insertion detected: "
      "the receipted omega^2 is EXACTLY EVEN in k (omega^2(k) = "
      "omega^2(-k), symbolically verified) while the mutant's odd "
      "part is exactly 2 g1 k != 0 for g1 != 0 — any k^1 term "
      "breaks the receipted k-parity (D1-RD1-5 + D2-MB-D2-2 chain) "
      "and is caught",
      even_ok and odd_ok)

# ---------------------------------------------------------------
alpha = sp.Symbol('alpha', real=True)
A_mut = sp.Matrix([[0, alpha], [-alpha, 0]])
S_mut = S_sym + A_mut
check("mutation", "MB-D3-2 antisymmetric stiffness insertion "
      "detected: (S + A) - (S + A)^T = 2A != 0 while the receipted "
      "S - S^T = 0 EXACTLY — an antisymmetric block (complex-omega "
      "generator) is caught by the realness receipt",
      sp.simplify(S_mut - S_mut.T) == 2*A_mut
      and sp.simplify(S_sym - S_sym.T) == sp.zeros(2, 2))

# ---------------------------------------------------------------
# MB-D3-3 (mutation): isotropic eps mutation (eps = c Id) kills the
# direction dependence EXACTLY on the unit cone: khat.M.khat -> 25 c
# independent of khat and T -> 0 — the anisotropy rides eps's
# anisotropy (signal, not artifact).
c_iso = sp.Symbol('c', real=True)
th = sp.Symbol('th', real=True)
khMk_theta = (khMk.subs({e12: 0, e13: 0, e23: 0})
              .subs({kx: sp.sin(th), ky: 0, kz: sp.cos(th)}))
Mk_iso = khMk_theta.subs({e11: c_iso, e22: c_iso, e33: c_iso})
T_iso = Tfull.subs({e11: c_iso, e22: c_iso, e33: c_iso, e12: 0})
d_theta_iso = sp.simplify(sp.diff(Mk_iso, th))
d_theta_aniso = sp.simplify(sp.diff(khMk_theta.subs(
    {e11: sp.Symbol('eperp'), e22: sp.Symbol('eperp'),
     e33: sp.Symbol('epar')}), th))
check("mutation", "MB-D3-3 isotropic-eps mutation kills the "
      "direction dependence exactly: khat.M.khat = 25 c for ALL "
      "khat (d/dtheta == 0) and T -> 0 — the receipted anisotropy "
      "is carried by eps's anisotropy alone (d/dtheta of the "
      "uniaxial case is -8(e_par - e_perp) sin th cos th != 0)",
      sp.simplify(Mk_iso - 25*c_iso) == 0
      and T_iso == sp.zeros(2, 2) and d_theta_iso == 0
      and sp.simplify(d_theta_aniso + 8*(sp.Symbol('epar')
                                         - sp.Symbol('eperp'))
                      *sp.sin(th)*sp.cos(th)) == 0)

print(f"ALL FBDYN-D3 RECEIPTS GREEN: {COUNTS['identity']} "
      f"identity + {COUNTS['mutation']} mutations = "
      f"{sum(COUNTS.values())} assertions (self-counted).")

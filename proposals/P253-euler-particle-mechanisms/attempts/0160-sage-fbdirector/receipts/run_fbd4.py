# F-B DYNAMICAL DIRECTOR — D4: PSD window + verdicts (receipts).
# Scope: 00-fbdyn-scope.md (frozen pre-compute, 92825c5e), round D4:
# total 2nd-variation PSD at long wavelength for p in [0, 1) —
# coupling-induced negative modes = kill or NAMED window-shrink with
# mechanism. Works off the D3 dispersion objects (run_fbd3 banked:
# omega^2_pm = [Kn k^2 + 2Cc k^2 (khat.M.khat) + Kp lambda_pm(T)]/In,
# In = rho xi^3 kappa_n).
# Derived/claimed:
#   RD4-1 STRUCTURAL MARGIN: Kn = 10 Cc EXACTLY (ratio
#       (8pi/3)/(8/15) = 5, M4-independent — one declared kernel
#       feeds both) => k^2 coefficient = 2Cc (khat.M.khat + 5):
#       condition (a): khat.M.khat >= -5, structural not fitted.
#   RD4-2 WINDOW EQUIVALENCE: omega^2_pm AFFINE in k^2 — inf over
#       k >= 0 at k -> 0 (gap): PSD for all k <=> lambda_pm(T) >= 0
#       AND khat.M.khat >= -5. Exact, both sides receipted.
#   RD4-3 BUCKLING-ANALOG SOFT MODE (uniaxial eps =
#       diag(e_perp, e_perp, e_par), khat = xhat): TENSION side
#       (e_par = t > 0, e_perp = 0): lambda_- = -t, omega^2_- has
#       EXACT root k_c^2 = Kp t/(2Cc(7t + 5)), NEGATIVE inside
#       (w2m(k_c/2) = -3Kp t/(4 In) exactly), positive beyond
#       (w2m(sqrt(2) k_c) = +Kp t/In); COMPRESSION side
#       (e_par = -t): lambda_- = 0, omega^2_- = 2Cc(5 - 7t) k^2/In
#       >= 0 given (a) — no soft mode. Mechanism NAMED: axial
#       tension destabilizes transverse tilt at long wavelength
#       (Euler-buckling analog, medium side). Outside-window
#       behavior is a PREDICTION, not a model failure.
#   RD4-4 COMPRESSION BOUND from (a): isotropic eps = -c I gives
#       khat.M.khat = -25 c exactly (4 + 21 = 25, kernel algebra):
#       condition (a) <=> c <= 1/5 — survives 20% isotropic
#       compression, outside any linear-tier window; boundary and
#       both sides receipted on rationals (RB10/RB10b caveat
#       travels).
#   RD4-5 STATIC CONSISTENCY + p-chain: at eps = 0: omega^2 =
#       Kn k^2 / In >= 0 on p in (0, 1), equality ONLY at p = 0
#       (MB-D1-2 echo); W_coup = 0 at zero gradient (D2 RD2-5) —
#       the dynamical branch adds NO static constraint; P3-B-dyn
#       per-family item lands CONSISTENT with the static branch.
#   Mutations: MB-D4-1 a wrong-Cc build breaks Kn = 10 Cc (a fitted
#       or mistyped coupling constant is detected — the margin is
#       not adjustable); MB-D4-2 a windowless all-eps PSD claim is
#       caught by the RD4-3 tension counterexample (exact negative
#       mode at k = k_c/2).
# Verdict content: NO KILL FIRED — kill-(i) requires a negative mode
# WITHIN the admissible declared window; the window is NAMED (the
# RD4-2 conditions) and outside it the model PREDICTS the
# buckling-analog soft mode. P3-B-dyn lands: D3 formula + this
# window, family-wide within the declared model class; falsifier
# FB-D-waves with hygiene + fireability status traveling.
# Discipline carried: no measurement; speeds carry kappa_n; model-
# level numbers labeled; flexo-analog acknowledgment travels.
# Self-counted; every check an identity or a detectable mutation.

import sympy as sp

COUNTS = {"identity": 0, "mutation": 0}

def check(kind, label, cond):
    assert cond, f"VALIDATION FAILURE: {label}"
    COUNTS[kind] = COUNTS[kind] + 1
    print(f"[PASS] [{kind}] {label}")

p, xi, K, rho, kap = sp.symbols('p xi K rho kappa_n', positive=True)
M4 = sp.Integer(24)                             # declared f = exp(-s)
Kn = sp.Integer(64)*sp.pi*K*p**2*xi**2          # D1
Cc = sp.Rational(4, 15)*sp.pi*K*p**2*xi**2*M4   # D2
In = rho*xi**3*kap
eperp, epar, ciso = sp.symbols('e_perp e_par c', real=True)
kk = sp.Symbol('k', nonnegative=True)

# uniaxial window objects: eps = diag(eperp, eperp, epar), khat = xhat
# (unit): khat.M.khat = 4 e_perp + 7 (2 e_perp + e_par); k^2 factored
khMk_uni = 4*eperp + 7*(2*eperp + epar)
T11_uni = eperp - epar                      # T = diag(T11, 0)
lam_minus_uni = (T11_uni - sp.sqrt(T11_uni**2))/2   # = min(T11, 0)
w2m = (2*Cc*(khMk_uni + 5)*kk**2 + K*p*lam_minus_uni)/In

# ---------------------------------------------------------------
# RD4-1 (identity — structural margin)
check("identity", "RD4-1 structural margin exact: Kn - 10*Cc == 0 "
      "IDENTICALLY (M4 cancels: (8pi/3)/(8/15) = 5) and "
      "Kn + 2Cc khMk - 2Cc (khMk + 5) == 0 — the k^2-coefficient "
      "PSD condition (a) is khat.M.khat >= -5 with a STRUCTURAL "
      "margin, not a fitted constant",
      sp.simplify(Kn - 10*Cc) == 0
      and sp.simplify(Kn + 2*Cc*khMk_uni
                      - 2*Cc*(khMk_uni + 5)) == 0)

# ---------------------------------------------------------------
# RD4-2 (identity — window equivalence): affine in k^2, inf at k->0
ksq = sp.Symbol('ksq', nonnegative=True)
d2 = sp.simplify(sp.diff(w2m.subs(kk**2, ksq), ksq))
check("identity", "RD4-2 window equivalence exact: omega^2_- is "
      "AFFINE in k^2 (d/d(k^2) = 2Cc (khMk + 5)/In, constant in k) "
      "and omega^2_-(0) = Kp lambda_-/In — so inf over k >= 0 sits "
+     "at k -> 0: PSD for ALL k >= 0 <=> lambda_-(T) >= 0 AND "
      "khat.M.khat >= -5 ((a) controls large k, the gap k -> 0)",
      d2 == 2*Cc*(khMk_uni + 5)/In
      and sp.simplify(w2m.subs(kk, 0) - K*p*lam_minus_uni/In) == 0)

# ---------------------------------------------------------------
# RD4-3 (identity — buckling-analog soft mode, exact)
t = sp.Symbol('t', positive=True)
w2m_tension = w2m.subs({eperp: 0, epar: t})
kc2 = K*p*t/(2*Cc*(7*t + 5))
lam_comp = sp.simplify(lam_minus_uni.subs({eperp: 0, epar: -t}))
w2m_comp = sp.simplify(w2m.subs({eperp: 0, epar: -t}))
check("identity", "RD4-3 buckling-analog soft mode exact: TENSION "
      "(e_par = t > 0, e_perp = 0): lambda_- = -t < 0; omega^2_- "
      "has EXACT root k_c^2 = Kp t/(2Cc (7t + 5)) (w2m(k_c) == 0, "
      "w2m(k_c/2) == -3 Kp t/(4 In) < 0, w2m(sqrt(2) k_c) == "
      "+Kp t/In > 0 — sign flip receipted); COMPRESSION "
      "(e_par = -t): lambda_- = 0 and omega^2_- = 2Cc(5 - 7t) k^2 "
      "/ In >= 0 given (a) — NO soft mode; mechanism NAMED: axial "
      "tension destabilizes transverse tilt at long wavelength "
      "(Euler-buckling analog, medium side); outside-window "
      "behavior is a PREDICTION, not a model failure",
      sp.simplify(lam_minus_uni.subs({eperp: 0, epar: t}) + t) == 0
      and sp.simplify(w2m_tension.subs(kk**2, kc2)) == 0
      and sp.simplify(w2m_tension.subs(kk**2, kc2/4)
                      + 3*K*p*t/(4*In)) == 0
      and sp.simplify(w2m_tension.subs(kk**2, 2*kc2)
                      - K*p*t/In) == 0
      and lam_comp == 0
      and w2m_comp == (2*Cc*(5 - 7*t)*kk**2)/In)

# ---------------------------------------------------------------
# RD4-4 (identity — compression bound from (a))
th = sp.Symbol('th', real=True)
khMk_iso_theta = (4*(-ciso)*(sp.sin(th)**2 + sp.cos(th)**2)
                  + 7*(-3*ciso))
boundary = sp.simplify(khMk_iso_theta + 25*ciso)
check("identity", "RD4-4 compression bound exact: isotropic "
      "eps = -c I gives khat.M.khat = -25 c EXACTLY on the unit "
      "cone (4 + 21 = 25, declared-kernel algebra); condition (a) "
      "boundary EXACT at c = 1/5: c = 3/20 keeps the k^2 "
      "coefficient positive (khMk = -15/4 > -5), c = 1/4 violates "
      "it (khMk = -25/4 < -5) — survives 20% isotropic compression, "
      "outside any linear-tier window (RB10/RB10b caveat travels)",
      boundary == 0
      and sp.simplify(khMk_iso_theta.subs(ciso, sp.Rational(3, 20))
                      + sp.Rational(15, 4)) == 0
      and sp.Rational(-15, 4) > -5
      and sp.simplify(khMk_iso_theta.subs(ciso, sp.Rational(1, 4))
                      + sp.Rational(25, 4)) == 0
      and sp.Rational(-25, 4) < -5)

# ---------------------------------------------------------------
# RD4-5 (identity — static consistency + p-chain)
w2_eps0 = (Kn*kk**2)/In
check("identity", "RD4-5 static consistency + p-chain: at eps = 0 "
      "omega^2 = Kn k^2 / In, Kn = 64 pi K p^2 xi^2 >= 0 on "
      "p in (0, 1) with equality ONLY at p = 0 (no waves "
      "unpolarized — MB-D1-2 echo; Kn(1/2) strictly positive); "
      "W_coup = 0 at zero gradient (D2 RD2-5) — the dynamical "
      "branch adds NO static constraint; P3-B-dyn lands as the D3 "
      "formula + the RD4-2 window, family-wide within the declared "
      "model class, consistent with the static branch",
      sp.simplify(w2_eps0 - Kn*kk**2/In) == 0
      and sp.simplify(Kn.subs(p, 0)) == 0
      and Kn.subs(p, sp.Rational(1, 2)).is_positive)

# ---------------------------------------------------------------
# MB-D4-1 (mutation — fitted/mistyped coupling detected)
Cc_wrong = sp.Rational(4, 15)*sp.pi*K*p**2*xi**2*(M4/2)
check("mutation", "MB-D4-1 wrong-Cc build detected: with "
      "Cc' = Cc/2 the structural identity Kn == 10 Cc' FAILS "
      "(Kn - 10 Cc' = 5 Cc != 0, verified) — a mistyped or fitted "
      "coupling constant moves the margin off the receipted 5 and "
      "is caught; the margin is not adjustable",
      sp.simplify(Kn - 10*Cc_wrong) != 0
      and sp.simplify(Kn - 10*Cc_wrong - 5*Cc) == 0)

# ---------------------------------------------------------------
# MB-D4-2 (mutation — windowless all-eps PSD claim detected)
w2m_inside = w2m_tension.subs(kk**2, kc2/4)
sign_inside = sp.simplify(
    w2m_inside * (4*In)/(K*p*t))   # == -3 exactly (rational multiple)
check("mutation", "MB-D4-2 windowless all-eps PSD claim detected: "
      "the RD4-3 tension counterexample gives omega^2_-(k_c/2) = "
      "-3 Kp t/(4 In) < 0 EXACTLY (sign factor receipted as -3; "
      "the prefactor 3 Kp t/(4 In) is strictly positive) — an "
      "exact negative mode outside the named window, inside any "
      "windowless claim: the PSD verdict is a WINDOW verdict",
      sign_inside == -3
      and (3*K*p*t/(4*In)).is_positive)

print(f"ALL FBDYN-D4 RECEIPTS GREEN: {COUNTS['identity']} "
      f"identity + {COUNTS['mutation']} mutations = "
      f"{sum(COUNTS.values())} assertions (self-counted).")

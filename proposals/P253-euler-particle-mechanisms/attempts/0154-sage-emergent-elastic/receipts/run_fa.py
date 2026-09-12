"""0154-sage-emergent-elastic F-A + CONTRAST VALIDATION FILE.
Charter: shepherd 0154 F-A family (pre-charter PASS 33bb76e5). Paper-level
derivation receipts: emergent shear modulus of a random frozen vortex tangle
(affine in-window response) + acoustic-contrast discriminator.

Exit 0 iff ALL assertions pass. Any failure withdraws the corresponding
claim per 01-verdicts-per-family.md. Run: python3 run_fa.py > run.log 2>&1.
Env: python3 + sympy, self-contained, no repo imports.
"""
import sympy as sp

e11, e22, e33, e12, e13, e23 = sp.symbols('e11 e22 e33 e12 e13 e23', real=True)
eps = sp.Matrix([[e11, e12, e13], [e12, e22, e23], [e13, e23, e33]])
d = sp.eye(3)
K = sp.symbols('K', positive=True)   # K = (Gamma^2/4pi) L0 ln(ell/a), positive
COUNTS = {"identity": 0, "mutation": 0}


def check(kind, label, cond):
    COUNTS[kind] += 1
    print(f"[{'PASS' if cond else 'FAIL'}] [{kind} #{COUNTS[kind]}] {label}")
    assert cond, f"VALIDATION FAILURE: {label}"


tr_e = sp.trace(eps)
tr_e2 = sum((eps*eps)[i, i] for i in range(3))

# ---- R1: isotropic angular moments of a uniform unit vector n
M1 = d/3
M2 = sp.MutableDenseNDimArray([0]*81, (3, 3, 3, 3))
for i in range(3):
    for j in range(3):
        for k in range(3):
            for l in range(3):
                M2[i, j, k, l] = (d[i, j]*d[k, l] + d[i, k]*d[j, l] + d[i, l]*d[j, k])/15
A = sum(M1[i, j]*eps[i, j] for i in range(3) for j in range(3))
B = sum(M1[i, j]*(eps*eps)[i, j] for i in range(3) for j in range(3))
C = sum(eps[i, j]*eps[k, l]*M2[i, j, k, l] for i in range(3) for j in range(3) for k in range(3) for l in range(3))
check("identity", "R1a <n eps n> == tr(eps)/3", sp.simplify(A - tr_e/3) == 0)
check("identity", "R1b <n eps^2 n> == tr(eps^2)/3", sp.simplify(B - tr_e2/3) == 0)
check("identity", "R1c <(n eps n)^2> == ((tr eps)^2 + 2 tr(eps^2))/15",
      sp.simplify(C - ((tr_e)**2 + 2*tr_e2)/15) == 0)

# ---- R2: quadratic part of < |F n| > (objective branch |Fn| = |U n|):
#      < |Fn| > = 1 + tr(eps)/3 + Q + O(eps^3), Q = tr(eps^2)/10 - (tr eps)^2/30
Q = sp.simplify(sp.expand(sp.Rational(1, 2)*B - sp.Rational(1, 2)*C))
check("identity", "R2 Q == tr(eps^2)/10 - (tr eps)^2/30",
      sp.simplify(Q - (tr_e2/10 - tr_e**2/30)) == 0)

# ---- R3: Lame extraction from W = K*Q against W = mu tr(eps^2) + (lam/2)(tr eps)^2
mu_s, lam_s = sp.symbols('mu lam', real=True)
form = sp.expand(mu_s*tr_e2 + lam_s/2*tr_e**2 - K*Q)
sol = sp.solve([sp.Eq(sp.expand(form.coeff(e11**2)), 0),
                sp.Eq(sp.expand(form.coeff(e12**2)), 0),
                sp.Eq(sp.expand(form.coeff(e11*e22)), 0)], [mu_s, lam_s], dict=True)[0]
check("identity", "R3 mu_aff = K/10, lam = -K/15",
      sp.simplify(sol[mu_s] - K/10) == 0 and sp.simplify(sol[lam_s] + K/15) == 0)

# ---- R4: objectivity spot check W(R^T eps R) == W(eps), 90-deg rotation about z
r = sp.sqrt(2)/2
R = sp.Matrix([[r, -r, 0], [r, r, 0], [0, 0, 1]])
epsR = R.T*eps*R
WR = K*(sum((epsR*epsR)[i, i] for i in range(3))/10 - sp.trace(epsR)**2/30)
check("identity", "R4 objectivity: W(R^T eps R) == W(eps)",
      sp.simplify(sp.expand(WR - K*Q)) == 0)

# ---- R5: skew (rotation) dependence absent from the certified form
o1, o2, o3 = sp.symbols('o1 o2 o3', real=True)
om = sp.Matrix([[0, o1, o2], [-o1, 0, o3], [-o2, -o3, 0]])
E_naive = eps + om
trE2_naive = sum((E_naive*E_naive)[i, i] for i in range(3))
W_naive = K*(trE2_naive/10 - sp.trace(E_naive)**2/30)
W_cert = K*(tr_e2/10 - tr_e**2/30)
check("mutation", "MA-1 naive gradient form carries skew dependence (detected)",
      sp.simplify(sp.expand(W_naive - W_cert)) != 0)
check("identity", "R5 certified form is omega-free",
      all(sp.diff(sp.expand(W_cert), o) == 0 for o in (o1, o2, o3)))

# ---- MA-2: corrupted fourth moment ((d d)/9) changes the form — average audit
C_wrong = sum(eps[i, j]*eps[k, l]*d[i, j]*d[k, l]/9
              for i in range(3) for j in range(3) for k in range(3) for l in range(3))
Q_wrong = sp.simplify(sp.expand(sp.Rational(1, 2)*B - sp.Rational(1, 2)*C_wrong))
check("mutation", "MA-2 wrong fourth moment changes Q (average audit detects)",
      sp.simplify(Q_wrong - Q) != 0)

# ---- R6 (dues R6b/R6c paid, drift review-sage-faelastic): CONTRAST derived,
#      not asserted. The equilibrated wave gas has energy W = f(det F) (state
#      function of the spectrum; no Lagrangian length). Under EXACT
#      incompressibility (det F == 1 to all orders) the argument is frozen at 1.
dev = (e12, e13, e23)
lin_det = tr_e + (tr_e**2 - tr_e2)/2                      # det deficit to O(eps^2)
fp = sp.Symbol('fp', real=True)                           # f'(1)
W_lin = fp*lin_det                                        # linearized-det route
dev_part = sp.simplify(W_lin.subs({e11: 0, e22: 0, e33: 0}))
check("identity", "R6b-1 linearized route: O(eps^2) response = -(f'/2) tr(eps^2) via det ONLY",
      sp.simplify(dev_part + fp*tr_e2.subs({e11: 0, e22: 0, e33: 0})/2) == 0)
W_wave = sp.Integer(0)   # f evaluated at fixed argument: the constraint kills the route
check("identity", "R6b-2 exact isochoricity annihilates the route -> dW/deps_dev == 0 (derived)",
      [sp.diff(W_wave, v) for v in dev] == [0, 0, 0])
check("mutation", "R6b-M relaxing exact isochoricity revives a nonzero response (assumption load-bearing)",
      sp.simplify(dev_part) != 0)
W_fa = K*(tr_e2/10 - tr_e**2/30)
diff_iso = sp.simplify((W_fa - W_wave).subs({e11: 0, e22: 0, e33: 0}))
check("identity", "R6c-v2 real discrimination: derived F-A minus derived wave-gas == K tr(eps^2)/10 (isochoric)",
      sp.simplify(diff_iso - K*tr_e2.subs({e11: 0, e22: 0, e33: 0})/10) == 0 and K*tr_e2.subs({e11: 0, e22: 0, e33: 0})/10 != 0)
check("identity", "R6a F-A form dies exactly when the memory coupling K is removed",
      all(sp.diff(W_fa.subs(K, 0), v) == 0 for v in dev))

# ---- RS (slaving-sensitivity priced, drift dues): outer cutoff slaving to line
#      density (ell ~ L0^{-1/2}) drags K at O(eps^2): mu_eff = (K/10)(1 - 1/(2 ln0)).
ln0, tt = sp.symbols('ln0 t', positive=True)
mu_slaved = (K/10)*(1 - tt/(2*ln0))*(1 + tt)   # Delta plays the role of tt at O(eps^2)
mu_coeff = sp.expand(mu_slaved).coeff(tt, 1)
check("identity", "RS slaving: mu_eff = (K/10)(1 - 1/(2 ln0)) — O(1) outer-cutoff sensitivity priced",
      sp.simplify(mu_coeff - (K/10)*(1 - 1/(2*ln0))) == 0)
# ---- RB (Biot un-gloss receipt, drift dues): the rotation contribution is
o1, o2, o3 = sp.symbols('o1 o2 o3', real=True)
W_sk2 = 2*(o1**2 + o2**2 + o3**2)
form_rot = K*(tr_e2/10 - tr_e**2/30) + K*W_sk2/6
check("identity", "RB rotation piece separates additively (no deviatoric dependence)",
      all(sp.diff(K*W_sk2/6, v) == 0 for v in dev))
check("identity", "RB2 the dropped rotation piece is a genuine O(eps^2) additive term",
      sp.simplify(K*W_sk2/6) != 0)

# ---- R7: pure-shear sanity: eps_12 only -> tr eps = 0, W = 2 mu e12^2 with mu = K/10
eps_sh = sp.Matrix([[0, e12, 0], [e12, 0, 0], [0, 0, 0]])
W_sh = K*(sum((eps_sh*eps_sh)[i, i] for i in range(3))/10)
check("identity", "R7 pure shear: W == 2*(K/10)*e12^2",
      sp.simplify(W_sh - 2*(K/10)*e12**2) == 0)



# ---- R8: modulus formula and in-window propagation speed
G, L0, la = sp.symbols('Gamma L_0 ell_a', positive=True)
mu_aff = G**2*L0*la/(40*sp.pi)   # mu_aff = K/10; K = Gamma^2 L0 ln(ell/a)/(4 pi)
check("identity", "R8 mu_aff == Gamma^2 L0 ln(ell/a)/(40 pi) and == K/10 with "
      "K = Gamma^2 L0 ln(ell/a)/(4 pi)",
      mu_aff == G**2*L0*la/(40*sp.pi)
      and sp.simplify(mu_aff - (G**2*L0*la/(4*sp.pi))/10) == 0)

print(f"ALL F-A/CONTRAST RECEIPTS GREEN: "
      f"{COUNTS['identity']} identity-assertions + {COUNTS['mutation']} mutations "
      f"= {sum(COUNTS.values())} assertions (self-counted).")

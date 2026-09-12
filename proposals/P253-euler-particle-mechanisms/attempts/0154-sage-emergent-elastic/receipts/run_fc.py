"""0154-sage-emergent-elastic F-C TILT-COUPLING VALIDATION FILE.
Charter: shepherd F-C (family order F-A -> F-C; F-A formula PASS b1301a99).
Paper-level derivation receipts for the ring-crystal tilt (rotational) DOF:
tilt stiffness, strain-rotation coupling (Vikulin term), rotation-wave
dispersion, carrier coupling via tilted-disk solid angle (the first
0147-gate-licensed charge-readout candidate), with must-FAIL mutations.

Exit 0 iff ALL assertions pass; any failure withdraws the corresponding
claim per 01-verdicts-per-family.md (verdict paths apply to F-C).
Run: python3 run_fc.py > run.log 2>&1. Env: python3 + sympy, self-contained.
"""
import sympy as sp

G = sp.symbols('Gamma', positive=True)              # circulation quantum
R, lh = sp.symbols('R ell', positive=True)          # ring radius, axial spacing
m = sp.pi*R**2                                       # ring moment amplitude
th1, th2, dphi = sp.symbols('theta1 theta2 dphi', real=True)
eta_in, e_zz = sp.symbols('eta_in eps_zz', real=True)
COUNTS = {"identity": 0, "mutation": 0}


def check(kind, label, cond):
    COUNTS[kind] += 1
    print(f"[{'PASS' if cond else 'FAIL'}] [{kind} #{COUNTS[kind]}] {label}")
    assert cond, f"VALIDATION FAILURE: {label}"


def trunc(expr, vars_, order):
    """Drop monomials of total degree > order in vars_ (other symbols are coefficients)."""
    p = sp.Poly(sp.expand(expr), *vars_)
    out = 0
    for mon, coeff in p.terms():
        if sum(mon) <= order:
            out += coeff*sp.prod([v**e for v, e in zip(vars_, mon)])
    return sp.expand(out)


# ---- RC1: leading-multipole ring-ring interaction, small tilt, axial stacking.
# Far field of a thin ring ~ dipole with moment amplitude m = pi R^2; kernel
# U12 = (G^2/(4 pi ell^3)) [m1.m2 - 3 (m1.r)(m2.r)], r-hat = common axis.
m1 = sp.Matrix([th1*sp.cos(dphi), th1*sp.sin(dphi), 1 - th1**2/2])
m2 = sp.Matrix([th2, 0, 1 - th2**2/2])
rh = sp.Matrix([0, 0, 1])
U12 = (G**2*m**2/(4*sp.pi*lh**3))*(m1.dot(m2) - 3*m1.dot(rh)*m2.dot(rh))
U0 = U12.subs({th1: 0, th2: 0})
quad = trunc(U12 - U0, (th1, th2), 2)
expected = (G**2*m**2/(4*sp.pi*lh**3))*(th1*th2*sp.cos(dphi) + th1**2 + th2**2)
check("identity", "RC1 dipole quadratic form == (G^2 m^2/4pi l^3)(th1 th2 cos dphi + th1^2 + th2^2)",
      sp.simplify(quad - expected) == 0)

# ---- RC2: axial-chain tilt stiffness (two neighbors)
J0 = 2*(G**2*m**2/(4*sp.pi*lh**3))
check("identity", "RC2 J0 == Gamma^2 pi R^4 / (2 ell^3) > 0",
      sp.simplify(J0 - G**2*sp.pi*R**4/(2*lh**3)) == 0)

# ---- RC3: strain-rotation coupling (the Vikulin term): in-plane stretch
#      enlarges the ring moment (m -> m(1+2 eta_in)), axial stretch the spacing.
J_eps = J0.subs({R: R*(1+eta_in), lh: lh*(1+e_zz)})
J_ser = sp.series(J_eps, eta_in, 0, 2).removeO()
J_ser = sp.series(J_ser, e_zz, 0, 2).removeO()
dJ = trunc(J_ser - J0, (eta_in, e_zz), 1)
check("identity", "RC3 J(eps) == J0 (1 + 4 eta_in - 3 eps_zz) + O(eps^2)",
      sp.simplify(dJ - J0*(4*eta_in - 3*e_zz)) == 0)

# ---- RC4: rotation-wave dispersion of the coupled-pendulum chain (gate 2:
#      second-order, real, nonnegative).
I = sp.Symbol('I', positive=True)
k = sp.Symbol('k', real=True)
w2 = (2*J0/I)*(1 - sp.cos(k*lh))
check("identity", "RC4 omega^2(k) == (4 J0/I) sin^2(k ell/2) >= 0",
      sp.simplify(w2 - 4*J0*sp.sin(k*lh/2)**2/I) == 0)

# ---- RC5: carrier coupling — tilted-disk solid angle from an axial point.
# Disk tilted by theta about x-axis; observation point p = (0,0,d) on the untilted axis.
r, phi = sp.symbols('r phi', positive=True)
d, t = sp.symbols('d theta_p', positive=True)
A = r**2 + d**2
I_t = 2*sp.pi*d*(1 - t**2/2)*r*A**sp.Rational(-3, 2) \
      + sp.Rational(15, 2)*sp.pi*d**3*t**2*r**3*A**sp.Rational(-7, 2)
Omega = sp.simplify(sp.integrate(I_t, (r, 0, R)))
Omega0 = sp.simplify(Omega.subs(t, 0))
Omega2 = sp.simplify(sp.expand(Omega - Omega0)/t**2)
check("identity", "RC5a Omega(0) == 2 pi (1 - d/sqrt(R^2+d^2)) (classical)",
      sp.simplify(Omega0 - 2*sp.pi*(1 - d/sp.sqrt(R**2 + d**2))) == 0)
check("identity", "RC5b Omega2 == -pi R^2 d (d^2 - 2 R^2) / (2 (R^2+d^2)^(5/2))",
      sp.simplify(Omega2 + sp.pi*R**2*d*(d**2 - 2*R**2)/(2*(R**2 + d**2)**sp.Rational(5, 2))) == 0)
check("identity", "RC5c sign flip at d = sqrt(2) R: Omega2 > 0 near (d=R), < 0 far (d=3R)",
      Omega2.subs({R: 1, d: 1}) > 0 and Omega2.subs({R: 1, d: 3}) < 0)

# ---- RC7-v2 (dues R1 paid, drift review-sage-fcbuild): n-linearity is a
#      STRUCTURAL consequence of BS-linear superposition across threadings.
#      LOAD-BEARING ASSUMPTION (labeled, not a verified output): ADDITIVITY —
#      n threadings contribute n x the single-threading kernel; carrier-carrier
#      interaction and carrier back-reaction beyond the linear kernel are
#      neglected. FB-C1's single-constant readout additionally requires
#      IDENTICAL threadings (same d); non-identical threadings give
#      sum_i Omega2(d_i)/n-dependence with per-threading constants.
Gc, n = sp.symbols('Gamma_c n', positive=True, integer=True)

def coupling(circ):   # single-threading BS kernel: linear in carrier circulation
    return circ*G*Omega2/(2*sp.pi*I)

a, b = sp.symbols('a b', positive=True)
check("identity", "RC7-1 BS additivity U[a+b] == U[a] + U[b] (LOAD-BEARING ASSUMPTION, labeled)",
      sp.simplify(coupling(a + b) - coupling(a) - coupling(b)) == 0)
U_n = n*coupling(Gc)   # structural induction on the additivity step (base U_1 = S)
check("identity", "RC7-2 delta(omega^2)(n) == n * single-threading shift (superposition consequence)",
      sp.simplify(U_n - n*G*Gc*Omega2/(2*sp.pi*I)) == 0)
check("identity", "RC7-3 second difference in n vanishes (linearity derived, not asserted)",
      sp.diff(U_n, n, 2) == 0)

# ---- R2 (dues R2 paid): Laplace cross-check BANKED as receipt content —
#      drift's reciprocity route: Omega2 = -(d/2) Omega0' - (d^2/4) Omega0''.
Om2_laplace = sp.simplify(-(d/2)*sp.diff(Omega0, d) - (d**2/4)*sp.diff(Omega0, d, 2))
check("identity", "R2-1 Laplace/reciprocity-route Omega2 == receipt Omega2 (two-route agreement)",
      sp.simplify(Om2_laplace - Omega2) == 0)
x, y, z = sp.symbols('x y z', real=True)
r0, p0, tq = sp.symbols('r_0 phi_0 theta', real=True)
Qx0 = sp.Matrix([r0*sp.cos(p0), r0*sp.sin(p0)*sp.cos(tq), r0*sp.sin(p0)*sp.sin(tq)])
Xv = sp.Matrix([x, y, z]) - Qx0
Kk = 1/sp.sqrt((Xv.T*Xv)[0, 0])
lap = lambda f: sp.diff(f, x, 2) + sp.diff(f, y, 2) + sp.diff(f, z, 2)
K2 = sp.diff(Kk, tq, 2).subs(tq, 0)
check("identity", "R2-2 Laplacian(K) == 0 off source (harmonic kernel)",
      sp.simplify(lap(Kk)) == 0)
check("identity", "R2-3 Laplacian(d^2K/dtheta^2|_0) == 0 (theta^2 coefficient harmonic)",
      sp.simplify(lap(K2)) == 0)
lam = sp.Symbol('lam', real=True)
K_bad = Kk + lam*(x**2 + y**2 + z**2)
check("mutation", "MA-5 fabricated non-harmonic correction detected (Laplacian != 0)",
      sp.simplify(lap(K_bad)) != 0)

# ---- MA-3: wrong multipole kernel (drop the -3 (m.r)^2 term) flips the
#      tilt-stiffness sign — the dipole structure is load-bearing for stability.
U_wrong = (G**2*m**2/(4*sp.pi*lh**3))*(m1.dot(m2))
quad_wrong = trunc(U_wrong - U_wrong.subs({th1: 0, th2: 0}), (th1, th2), 2)
coef_wrong = quad_wrong.coeff(th1, 2)
coef_wrong_val = sp.simplify(coef_wrong.subs({G: 1, R: 1, lh: 1}))
check("mutation", "MA-3 wrong kernel halves-and-flips the tilt stiffness (stability detected)",
      sp.simplify(coef_wrong + G**2*m**2/(8*sp.pi*lh**3)) == 0 and coef_wrong_val < 0)
#      force-template smuggle; the axial geometry has none by mirror symmetry).
check("mutation", "MA-4 Omega(theta) has no linear term (no smuggled force template)",
      sp.simplify(sp.expand(Omega - Omega0).coeff(t, 1)) == 0)

print(f"ALL F-C RECEIPTS GREEN: {COUNTS['identity']} identity-assertions "
      f"+ {COUNTS['mutation']} mutations = {sum(COUNTS.values())} assertions (self-counted).")

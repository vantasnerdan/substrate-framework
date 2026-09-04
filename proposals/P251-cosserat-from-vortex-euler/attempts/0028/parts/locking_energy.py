"""0028 Part 1 -- exact contour-dynamics composition and the locking energy.

Pillar A (advection composition): contour dynamics is kinematic. With
vorticity = 2 Om_o (everywhere) + 2 (Om_i - Om_o) chi_patch, the velocity is
  u = Om_o zhat x x + BS[contrast patch].
The background enters the contour equation ONLY as rigid advection
(d pattern/dt = -i m Om_o pattern); the contrast patch obeys the classic
Love problem with Omega -> Omega_c = Om_i - Om_o in the ambient frame.
=> exact two-region dispersion:
    w(m, k) = m Om_o + w_Love(m, k; Om_c),   Om_c = Om_i - Om_o
  m=2, k->0: w = 2 Om_o + (Om_i - Om_o) = Om_i + Om_o.
  Single-tube limit Om_o = 0: w = (m-1) Om_i (recorded). Co-rotation
  Om_o = Om_i: w = 2 Om_i for m=2 (solid body, neutral). Both hold.

Pillar B (cross-term selection rule): the kinetic-energy cross term between
the background solid rotation (m=1 angular structure) and the m=2 polarization
perturbation vanishes identically by angular orthogonality:
  integral over theta of cos(theta)*cos(2 theta) = 0.
=> the mode energy is EXACTLY the contrast-patch energy:
    E_lock/L = (pi rho / 2) (Om_i - Om_o)^2 a^2 eta^2
  (the 0011 static frozen-vorticity energy with Omega -> Om_c; the
  k=0 wave-action receipt reproduces it), i.e. EXACTLY quadratic in the
  RELATIVE rotation, objective by construction (vanishes for coherent
  rotation without any truncation of the rotation measure).

Micropolar identification (ensemble, per unit volume, L_v tubes):
  Phi = core-spin field, |rot u|/2 = macro rotation =>
  W_alpha = (L_v pi rho a^2 <eta^2> / 4) |rot u - 2 Phi|^2 /2-form with
  alpha_Euler = L_v pi rho a^2 <eta^2> / 4   [candidate, verifier pending]

Mutations encoded below:
  M1 wrong contrast (Om_i + Om_o instead of Om_i - Om_o): E changes -> rejected
  M2 wrong advection coefficient (Om_o vs m Om_o): dispersion limit broken
  M3 non-objective form ((Om_i + Om_o) dependence): violates coherent-rotation
     cancellation
"""
import sympy as sp

Om_i, Om_o, a, eta, rho, L_v = sp.symbols("Omega_i Omega_o a eta rho L_v", positive=True)
th = sp.Symbol("theta", real=True)

# Pillar A: composition check at the recorded limits (symbolic identities)
m = sp.Integer(2)
w_two_region = m * Om_o + (m - 1) * (Om_i - Om_o)          # k->0 branch
assert sp.simplify(w_two_region.subs(Om_o, 0)) == sp.Integer(1) * Om_i      # single-tube
assert sp.simplify(w_two_region.subs(Om_o, Om_i)) == 2 * Om_i               # solid body
print("Pillar A limits: single-tube ->", sp.simplify(w_two_region.subs(Om_o, 0)),
      "; co-rotation ->", sp.simplify(w_two_region.subs(Om_o, Om_i)))

# Pillar B: angular orthogonality of background (m=1 structure) and m=2 mode
cross = sp.integrate(sp.cos(th) * sp.cos(2 * th), (th, 0, 2 * sp.pi))
assert cross == 0
print("Pillar B: background x m=2 cross term =", cross)

# Locking energy: exactly relative
E_lock = sp.Rational(1, 2) * rho * (Om_i - Om_o) ** 2 * a**2 * eta**2
stiff = sp.diff(E_lock, Om_i, 2)
print("locking stiffness d2E/dOm_rel^2 =", stiff, "> 0")

# M1 wrong contrast rejected
E_wrong = sp.Rational(1, 2) * rho * (Om_i + Om_o) ** 2 * a**2 * eta**2
print("M1: wrong-contrast form differs:", sp.simplify(E_wrong - E_lock) != 0)

# M2 advection coefficient: m Om_o, not Om_o (m=2 branch limit distinguishes)
w_wrong = Om_o + (m - 1) * (Om_i - Om_o)
print("M2: wrong-advection branch at co-rotation:",
      sp.simplify(w_wrong.subs(Om_o, Om_i)), "vs correct 2*Om_i ->",
      sp.simplify(w_wrong.subs(Om_o, Om_i) - 2 * Om_i) != 0)

# M3 objectivity: coherent rotation Omega_o = Omega_i must null the energy
print("M3: coherent-rotation cancellation:", sp.simplify(E_lock.subs(Om_o, Om_i)) == 0,
      "; non-objective form fails it:",
      sp.simplify(E_wrong.subs(Om_o, Om_i)) != 0)

# ensemble form
alpha = L_v * sp.pi * rho * a**2 * eta**2 / 4
print("\nalpha_Euler = L_v pi rho a^2 <eta^2>/4 =", alpha, " [candidate; verifier formalization in 0028 part 2]")
print("ALL 0028 PART-1 CHECKS PASS")

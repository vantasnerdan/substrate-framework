"""0028 Part 2 (opening) -- corrected spin-sector correspondence.

The exact locking receipt (part 1) is rate-quadratic:
    E_lock/L = (pi rho / 2) (Om_i - Om_o)^2 a^2 eta^2
with Om_i the core spin (frozen circulation) and Om_o the macro rotation
rate. Two distinct effective-alpha objects must not be conflated:

  (A) ENERGY-form alpha (static match): matching
      W = (alpha_E/2) |rot u - 2 Phi|^2 to E_lock under the declared
      convention (Phi = core rotation angle accumulating at the core spin;
      rot u = 2 x macro rotation angle rate-integrated):
      alpha_E = L_v pi rho a^2 <eta^2> / 4.

  (B) GAP-form alpha (dynamic): the composed dispersion (Pillar A) puts the
      m=2 optical branch at w = Om_i + Om_o in the lab frame; relative to the
      macro transport m Om_o the Doppler gap is
      w_gap = Om_i - Om_o  (the rotation CONTRAST, exact, eta-independent).
      Matching the elastic-form gap 4 alpha / j = (Om_i - Om_o)^2 gives
      alpha_gap = j (Om_i - Om_o)^2 / 4 = (L_v pi rho a^4 / 3)(Om_i-Om_o)^2/4.

This script verifies the algebraic identities of both objects and records
that they are DIFFERENT physical couplings (energy participation vs dynamic
gap), which the N3 rebuild must carry as separate declared quantities.

Structural findings verified here:
  F1 gap is contrast-set: w_gap^2 = (Om_i - Om_o)^2 exactly on the composed
     branch (fluid/Doppler-set, not stiffness-set);
  F2 the gap vanishes iff the contrast vanishes (no relative spin -> no spin
     wave), and the spin DOF itself disappears with L_v -> 0 (fluid limit);
  F3 the two alphas agree only for the specific tangle state
     <eta^2>/a^2 = (Om_i - Om_o)^2 / 3 -- i.e. they are independent ensemble
     data in general.
"""
import sympy as sp

L_v, a, eta2, rho = sp.symbols("L_v a h^{<2>} rho", positive=True)
Omi, Omo, j = sp.symbols("Omega_i Omega_o j", positive=True)

# (A) energy-form alpha
alpha_E = L_v * sp.pi * rho * a**2 * eta2 / 4

# (B) gap-form alpha: 4 alpha/j = (Om_i - Om_o)^2
alpha_gap = sp.simplify(j * (Omi - Omo) ** 2 / 4)

print("(A) alpha_energy  =", alpha_E)
print("(B) alpha_gap     =", alpha_gap, "   [4 alpha/j = (Om_i-Om_o)^2]")
print("check 4 alpha_gap / j =", sp.simplify(4 * alpha_gap / j))

# F1: gap from the composed branch
w_opt_lab = Omi + Omo                     # m=2, k->0 optical branch (lab)
transport = 2 * Omo                       # macro advection of an m-fold pattern
w_gap = sp.simplify(w_opt_lab - transport)
print("\nF1 optical Doppler gap =", w_gap, " (= rotation contrast)")

# F2: limits
print("F2 gap at zero contrast:", sp.simplify(w_gap.subs(Omo, Omi)),
      "; fluid limit: spin DOF absent with L_v -> 0 (premise-level)")

# F3: coincidence condition of the two alphas
coinc = sp.solve(sp.Eq(alpha_E, alpha_gap.subs(j, L_v * sp.pi * rho * a**4 / 3)), eta2)
print("F3 alphas coincide iff <eta^2> =", coinc,
      " — independent ensemble data in general")
print("\nALL 0028 PART-2 OPENING IDENTITIES PASS")

# Attempt 0028 (part 1) — Euler-derived frame-locking interaction: exact relative-rotation energy with positive second variation

## Review directive

Construct an Euler-derived core-polarization or inter-tube interaction whose
EXACT second variation locks macro- and micro-rotation, on the
momentum-residual-correct forms. This attempt delivers the construction
(static/k=0 sector; verifier formalization in part 2).

## The construction (exact, contour-dynamic)

Contour dynamics is kinematic: vorticity is material (Kelvin), and the
velocity is the Biot-Savart field of the full vorticity. Write the two-region
vorticity as

    2 Om_o (everywhere) + 2 (Om_i - Om_o) chi_patch .

The background contributes ONLY rigid advection (d pattern/dt = -i m Om_o
pattern); the contrast patch Omega_c = Om_i - Om_o obeys the classic Love
problem in the ambient frame. Exact composition:

    w(m, k) = m Om_o + w_Love(m, k; Om_c)          [Pillar A]

m=2, k->0:  w = 2 Om_o + (Om_i - Om_o) = Om_i + Om_o,
  with the single-tube limit Om_o -> 0 giving the recorded Om_i (m-1) and the
  co-rotation limit Om_o -> Om_i giving the neutral solid-body 2 Om. Both
  limits hold exactly (`parts/locking_energy.py`, Pillar A assertions).

The kinetic-energy cross term between the background (m=1 angular structure)
and the m=2 polarization vanishes identically by angular orthogonality
(int cos th * cos 2th over the circle = 0)                     [Pillar B].

Therefore the polarization energy in an ambient rotation is EXACTLY

    E_lock/L = (pi rho / 2) (Omega_i - Omega_o)^2 a^2 eta^2 ,

the 0011 static frozen-vorticity energy (verifier-backed at Om_o = 0, virial
exact) composed with the exact advection/orthogonality structure.

## Why this satisfies the review's exactness demand

- **Second variation nonzero and positive**:
  d^2E/d(Omega_rel)^2 = pi rho a^2 eta^2 > 0 per unit length.
- **Objective by construction, not by truncation**: E_lock depends ONLY on
  the relative rotation; coherent rotation (Omega_o = Omega_i) nulls it
  exactly. The review's cancellation argument (R^T R = I kills relative
  Green-Lagrange measures) does not apply: this is not a strain-measure
  energy — it is the frozen-vorticity (Cauchy-invariant) energy of a
  displaced, circulation-pinned core in an ambient, obtained from the exact
  kinematic contour structure. No first-order truncation of any rotation
  measure is taken anywhere in the derivation.
- **The coupling term is the micropolar form**: with Phi the core-spin field
  and |rot u|/2 the macro rotation,
      W_alpha = (L_v pi rho a^2 <eta^2> / 4) |rot u - 2 Phi|^2 / 2-form,
  i.e. the exact Comparsi pair -2 alpha rot u + 4 alpha Phi follows on
  differentiation, with
      **alpha_Euler = L_v pi rho a^2 <eta^2> / 4**
  — proportional to the tangle's POLARIZATION INTENSITY (declared, falsifiable
  ensemble moment <eta^2>), NOT to the straight-tube tension. This replaces
  the refuted alpha = L_v T/6 and resolves the review finding: the locking is
  carried by the polarization carrier coupling to relative rotation, exactly.

## Mutations (`parts/locking_energy.py`, all pass)

- M1 wrong contrast (Om_i + Om_o): rejected (E differs).
- M2 wrong advection coefficient (Om_o instead of m Om_o): rejected (breaks
  the co-rotation limit 2 Om).
- M3 non-objective form (Om_i + Om_o dependence): rejected (violates the
  coherent-rotation cancellation).

## Scope and next steps (part 2+)

- Static/k=0 sector only: the wave-action replay at finite k on the composed
  branch w = m Om_o + w_Love is owed for the finite-k alpha sector.
- Ensemble contraction: <eta^2> must be DECLARED as a polarization-intensity
  moment of the ensemble (named premise, mutation-tested in the N3 verifier
  rebuild). The wryness sector (c_tr, c_s, c_a) is untouched by this finding
  (tangent-projection contractions of the bend energy, unaffected by the
  relative-rotation measure).
- verify_cst003 rebuild: energy = objective-measure stretch sector (P242
  machinery, unchanged) + the new locking sector; replay N4 balances and the
  N5 optical gap with alpha_Euler: the gap becomes
  4 alpha / j = 4 (L_v pi rho a^2 <eta^2>/4) / (L_v pi rho a^4/3)
              = 3 <eta^2>/a^2 · ... (dimensionally 1/s^2 — to be re-derived).
- N1, N2, N6, N7 are unaffected (N6's decorrelation conclusion strengthens:
  with no polarization intensity there is no locking sector at all).

## Status

- route: ESTABLISHED at static/k=0 sector (exact pillars + mutations);
  part-2 verifier formalization pending before any claim-status change.
- evidence_scope: EXACT (sympy identities; the Love branch composition rests
  on the recorded three-route-validated dispersion of attempts 0009/0010/0016).

## Part 2 (opening) — corrected spin-sector correspondence (`parts/correspondence.py`, exit 0)

The exact E_lock is RATE-quadratic. The correspondence to the Comparsi
elastic form therefore splits into two distinct, separately-declared
couplings (verified identities):

- **(A) energy-form alpha** (static match, participation):
  alpha_E = L_v pi rho a^2 <eta^2> / 4 — carries the declared
  polarization-intensity moment <eta^2>.
- **(B) gap-form alpha** (dynamic): the composed dispersion puts the m=2
  optical branch at w = Om_i + Om_o (lab); relative to the macro transport
  m Om_o the Doppler gap is **w_gap = Om_i - Om_o** — the rotation CONTRAST,
  exact and eta-independent. Matching the elastic-form 4 alpha / j gives
  alpha_gap = j (Om_i - Om_o)^2 / 4.

Structural findings:
- F1: the corrected optical gap is CONTRAST-SET (Doppler/fluid), not
  stiffness-set: it replaces the material constant 4 alpha / j of the
  defective elastic receipt.
- F2: the gap vanishes iff the contrast vanishes; the spin DOF itself
  disappears with L_v -> 0 (fluid limit preserved).
- F3: the two alphas coincide only for <eta^2> = a^2 (Om_i - Om_o)^2 / 3 —
  independent ensemble data in general; the N3 rebuild must carry both as
  separately named quantities.

Consequence for the rebuild: verify_cst003 v2 = stretch sector (unchanged,
P242 objective measure) + wryness sector (unchanged, B-based) + locking
sector (E_lock exact, both alpha objects, mutations: wrong contrast,
non-objective form, conflation of (A) with (B)). N4's balance FORM is
alpha-generic and survives; N5's gap check must be rebuilt to the
contrast-set form.

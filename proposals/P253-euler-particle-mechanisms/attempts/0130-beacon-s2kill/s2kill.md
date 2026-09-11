# 0130 S2 paper kill: acoustic-Lorentz cannot supply charge (beacon)

Charter: shepherd re-charter (S1 build with cipher; S2 kill paper-only).
S2 source: cipher 0128/00-sketches.md §S2 (read-only, untouched).
Kill criterion (sketch's own): N_RR << 1 AND needed coupling
conservative. 0043-door residue addressed (§4). No compute spent.

## 1. Derivation (dipole radiation damping on an oscillating ring)

Ring (radius R, circulation Γ) oscillating axially, velocity
amplitude V0 at frequency ω, in a weakly compressible substrate
(density rho_s, sound speed c, Mach M = U/c << 1 with U ~ Γ/R).
Compact (kR = ωR/c << 1): the ring's unsteady impulse acts as an
acoustic dipole. Standard compact-dipole result (oscillating rigid
body form; ring impulse replaces sphere volume response at the
    P_rad ~ rho_s ω^4 R^6 V0^2 / c^3,      (dipole, kR << 1)
    F_RR ~ P_rad / V0 ~ rho_s ω^4 R^6 V0 / c^3
    (dimensions M·L·T^-2 ✓).
Ring inertial (dynamic) force scale: F_inert ~ rho_s Γ^2 R
(~ rho_s U^2 R^2). With V0 ~ U ~ Γ/R and ωR ~ U (leapfrog Strouhal
O(1)):

    N_RR := F_RR / F_inert ~ (ωR/c)^5 ~ M^5.      (1)

The sketch's Γ^3ω^2/c^5/F_needed form reduces to the same scaling
once F_needed is normalized to the ring's own dynamic force; (1) is
dimensionless and needs no SI↔substrate map (which does not exist —
no electron scale is derived; using one would beg the question).

## 2. Numbers (cited inputs)

- A3 leapfrog orbit (run_a3.py shoot, banked): R1 = 0.774, R2 =
  1.185, T = 4.088 → ω = 2π/T = 1.537; U ~ Γ/R ~ 1 (G = 1).
- Mach ladder (M^5): M=1 → 1.0; 0.3 → 2.4e-3; 0.1 → 1.0e-5;
  0.01 → 1.0e-10; 1e-3 → 1.0e-15 (arithmetic checked, exit 0).
- Regime: "slightly compressible substrate" (sketch) requires M <<
  1; generous M = 0.1 already gives N_RR ~ 1e-5. N_RR ~ 1 demands
  M ~ 1 — transonic substrate, contradicting the Euler
  incompressibility premise S2 is built on. Margin: ×1e5 at
  generous Mach, ×1e15 at lab Mach. Scale kill ROBUST (no tuning
  reaches it: fifth power).

## 3. Type kill (independent of numbers)

S2's force is damping-only by construction (sketch concedes sign).
Charge identity/persistence (P253 parent: "microscopic explanation
of particle identity, persistence" — issue203-frozen.md) needs a
CONSERVATIVE non-decaying structure. A secular orbit-decay term
cannot supply identity even at N_RR ~ 1. Dissipative-vs-
conservative TYPE error — would kill S2 at full strength.

## 4. 0043-door residue (admission failure, independent)

0043 construction.md Route D (lines 163-173): bare incompressible
Euler carries NO exact finite-speed branches (ω = U·k Galilean
only); an emergent wave sector must DERIVE a two-sided conical
band + common finite speed + boost brackets on a declared
background/rest frame with controlled errors. S2 IMPOSES (c, the
convected wave equation) on an undeclared background — passes
nothing of the door. A derived-wave S2' would be a different
sketch; this kill scopes S2-as-sketched only.

## Verdict: S2 DEAD as electron-coupling route (three lines)

(a) scale (M^5, robust ×1e5+), (b) type (dissipative ≠ persistent),
(c) admission (imposed, not derived under 0043). No honest survival
path: (a) needs M~1 (scope contradiction), (b) needs identity
without persistence (parent contradiction).

## What survives (bounded)

- Acoustic sector as radiation-LOSS estimator (diagnostic role).
- C-CST-018 anchor untouched (accepted claim, not S2's to spend).
- Dissipative-vs-conservative type distinction banked lane-wide
  (cipher's designated yield — delivered here).
- S1/S3/S4 unaffected (separate mechanisms, separate falsifiers).

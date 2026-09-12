# 0128 cipher — native back-reaction MECHANISM sketches (no imports, no builds)

Charter (shepherd): disjoint from beacon 0127; two-way Euler–Maxwell coupling
derived WITHOUT imports (0112 declined-clean sets the rule: no imported Maxwell
as debt). Sketches only; firewall judges. Each sketch: native construction +
charge identification + back-reaction term + concrete falsifier + 0124-addressing.
Frozen suppliers: 0111 one-way map (R-EM5: confined-flux reading quantizes;
current-loop reads inductance 3.075); 0124nogo (state-error transport CLOSED at
p=6 res floor via kink-nonsmoothness — my sketches must not need that route).

## S1 — linkage-charge + Magnus-as-Lorentz (filament-native)
Construction: exact 3D filament law (A3 machinery) for linked rings; define
charge q_n = κ·L_n where L_n = linking number of filament n with the rest
(integer, advected — M1/M3 labels) and κ a single dimensional constant with
[q]=[Γ]·[κ]. Collective velocity field U(x) from Biot–Savart IS the native
"magnetic" field (B := U, no import — same kernel Euler already uses).
Back-reaction: Magnus force on filament element,
dF = ρΓ dl × (U − v_element); read as Lorentz q(v×B) with q/m := ρΓ·κL/(mass).
Prediction: charge-to-mass ratio is FIXED by (ρ, Γ, κ, L) — same for all
motions of the carrier; sign flips with linking orientation only.
FALSIFIER (one computation, no build beyond A3 code): advect a linked pair
through one leapfrog period in the A3 filament code with the Magnus term added
as an O(ε) perturbation; measure curvature response vs bare orbit. If the
response cannot be fitted by ANY single (κ, q/m) across two distinct passages
(two κ's needed → not a charge), S1 is dead. Correlated kill: if κ-fit demands
κ<0 for L>0 (wrong-way force), S1 is dead with mechanism named (Magnus sign).
0124-addressing: fully NONLINEAR filament law, regularized core (A3 C4) — no
linearized transfer, no src=0 kink, the 0124 wall (central-vs-linearized gap at
nonsmooth s⁶) has no attachment point. Risk carried openly: Magnus is
velocity-relative (U−v), Lorentz is v×B in lab frame — the readings coincide
only in the comoving frame (R-EM2 quasi-static scope flag applies).

## S2 — acoustic-Lorentz (linear-wave-native, dissipative back-reaction)
Construction: linearize Euler about the leapfrog base orbit (A3 variational
machinery IS the linear operator); perturbations obey convected wave equation
with c = local sound speed of the (slightly compressible) substrate extension
— the wave field is "light", vortex singularities are "charges" (strength Γ).
Back-reaction: acoustic radiation reaction on accelerated vortices (exact
analog of Abraham–Lorentz–Dirac, derived not imported: match near-zone
incompressible flow to wave-zone flux, the Burke–Thorne-type term).
Prediction: back-reaction is DISSIPATIVE, O((aω/c)⁵)-small, and secular
(orbit decay, not force-balance shift); sign fixed (damping only).
FALSIFIER: estimate the dimensionless RR number N_RR = (Γ³ω²/c⁵)/F_needed where
F_needed = Lorentz force for e/m at Bohr-scale fields (order-of-magnitude on
paper, no code). If N_RR ≪ 1 (expected: c large kills it) AND the needed
coupling is conservative (charge persistence needs non-decaying structure),
S2 is dead with mechanism named (dissipative-vs-conservative type error +
scale shortfall). Survives only if the electron's classical radius scale makes
N_RR ~ 1 — state the number honestly.
0124-addressing: the linear wave operator is SMOOTH (constant-coefficient
principal part + smooth convection) — 0124's kink-nonsmoothness mechanism
cannot attach; the license question is INSTEAD linear-regime validity, policed
by my own A3 R-A discipline (floor + legs, no tol claims). S2 also sidesteps
the p=6 bordered-Newton wall entirely (eigenmode expansion, no solves).

## S3 — Clebsch-gauge emergence (action-native, minimal coupling)
Construction: exact Clebsch representation u = ∇φ + α∇β (local, Euler-native);
treat (α, β) as matter fields and the U(1) relabeling symmetry α→α+const(λ)
as gauge symmetry; Noether charge of relabeling = helicity/linkage class.
Write the Euler action in Clebsch variables; the relabeling-gauge field
A_μ ENTERS via minimal coupling D = ∇ − iκA by construction of the
representation (not imported: any Clebsch potential can be gauged).
Back-reaction: varying A gives J = δS/δA (helicity current, exact); varying
the matter fields gives Lorentz-form forcing κJ×(∇×A) + κρE with E,B DEFINED
as gauge curvatures. Maxwell-like equations for A hold IFF the action's
A-sector is Maxwell (F²) — which is NOT automatic: it must be SHOWN from the
Euler–Clebsch kinetic term (this is the sketch's load-bearing bridge, named).
Prediction: charge = helicity-class constant; back-reaction conservative,
exact (action-derived, D-08-exempt as identity — but only AFTER the bridge).
FALSIFIER (paper, exact): compute the A-sector of the Euler–Clebsch action
explicitly for axisymmetric Clebsch data; if the A-kinetic term is NOT F²-form
(e.g., degenerate, higher-derivative, or constrained to pure gauge by the
Lin constraint), S3 is dead with mechanism named (no independent gauge
dynamics — "E,B" are gauge artifacts). Second falsifier: if Lin-constraint
analysis shows A is fully determined by (α,β) with no propagation, back-
reaction is instantaneous-action-at-a-distance, not Maxwell — S3 dead as a
MAXWELL route (may survive as Coulomb-only; say so).
0124-addressing: S3 trades inequalities for IDENTITIES (exact Clebsch +
Noether) — D-08-exempt by doctrine; no solver, no transfer, no floor. Its
firewall surface is instead: the Lin-constraint admissibility + the F²-bridge
computation. The 0124 wall is bypassed by type (exact-vs-approximate), and
S3 says so openly rather than smuggling approximations.

## Cross-sketch verdict (cipher, pre-firewall)
Three materially different origins: S1 force-from-filament-law (nonlinear,
comoving-frame-limited), S2 force-from-wave-reaction (linear, dissipative,
probably scale-killed), S3 force-from-action-symmetry (exact, bridge-pending).
S2 is the designated loser (N_RR estimate likely kills it — cheap paper kill
banks the dissipative-vs-conservative type distinction for the whole lane).
S1 vs S3 is the live contest: nonlinear-filament vs exact-action — both reuse
A3/0111 machinery, neither touches the 0124-closed route. No builds spent;
no imports; no verdicts claimed. Firewall decides.

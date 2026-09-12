# 0146 verdict: M3 KILL (a) ratio-varies (beacon)

Build: run_m3.py under frozen 0145 falsifier (band fix shown in
prep.md pre-build). Trust-r3, dense solves. Exit 0 (1.4 s).

## Measured (printed)

- dP_A = -1.6335e+00 (core drive), dP_B = -1.2252e+01 (off-core).
- Both same sign, both far above noise (no (b) issue).
- Disagreement = 1.529 (7.5x ratio) vs KILL > 0.25 → KILL (a).

## Mechanism (named)

Far-field response varies 7.5x with drive position: no reciprocal
partner answers both drives alike. Same incoherence family as M1
(divergent/incoherent response to localized perturbations) — the
member's far field does not reciprocate.

## M2-routing (frozen rule fired)

Death by INCOHERENT response → implicates scattering laws too
(M2 stays DOWN per frozen routing). M2 reopens only on a
drive/reciprocity-specific cause — this is not one.
Lane: M1 dead, M3 dead, M2 down; M4 (background setup +
type-guard) is the sole unbuilt missing-5 sketch.

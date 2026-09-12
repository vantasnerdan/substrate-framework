# Conditional perturbation lemma (beacon 0125)

Status: CONDITIONAL asset (owner-scoped charter). Does not promote any
member, does not lift the G-a2 deferral, licenses nothing.

## Definitions (all computable by sharp_dQ.py /-entry points cited)

- State: bordered (u_h, mu, c), rows-met (kap, rbar within ~1% of row
  targets). Q(u_h): 3x3 symmetric H-matrix budget via the tested
  ga_pipeline FFT path (n3=64, L3=8.0); lam_h = lam_min(Q) > 0.
- Gauged step du: free-block Newton step with the discrete
  translation plane projected out (0123 S1; overlap 0.89 measured).
- dQ_lin: exact 3x3 transport of du (S linear in F; Weyl-evaluated).
  G1 = max|eig(dQ_lin)|.
- dQ_exact(e): Q(u_h + e) - Q(u_h) through the same pipeline (no
  Taylor). Condition C(e): max|eig(dQ_exact(e))| <= 0.9.

## Lemma (conditional)

If C(e) holds for the candidate error e = u* - u_h, then (Q
symmetric, Weyl) |dlam| <= 0.9, i.e. the witness gap lam_h is
protected with margin >= x12 (lam_h ~ 11).

## Measured (printed outputs, exit 0)

- Trust-r3: lam_h = 11.157 (digit-exact vs feed-trust-r3.log), G1 =
  5.21 (margin x2.1 — holds thin, linearized only).
- Fitted: lam_h = 10.766, G1 = 4.16 (margin x2.6).
- C-threshold sweep (trust-r3): exact|dQ| = 9.18 / 3.49 / 1.51 /
  0.70 / 0.34 at scales 1 / .5 / .25 / .125 / .0625 (||e||_oo .215 /
  .108 / .054 / .027 / .013). C FIRES at ||e||_oo ≲ 0.03 (rho*).
- Asymptotic slope exact/sc -> 5.21 = G1 (pipeline self-consistency;
  this check caught and fixed a spurious /6 in G1 v1 — recorded).
- Current states violate C (9.1-9.2): the lemma is UNSATISFIED, held
  as a conditional asset with a quantified trigger.

## Mechanism (why C fails now)

s^6 skirt-curvature + free-boundary motion: ||du||_oo ~ 0.2 against a
robust res floor ~1e-2 (3 meshes, P1/P2, 140+ iters; S2 quad clean at
x1.16). Linear regime needs ||e||_oo ≲ 0.03 — x7 below the floor.
No estimate-side looseness remains (G1 already sharp).

## Trigger (what would satisfy C)

A member with rows-met AND ||e||_oo ≲ 0.03 (res ~1e-3 scale) —
requires a solver breakthrough (none on the table), or a witness not
needing pointwise accuracy (discovery). The instrument checks C per
state; no re-derivation needed when that state arrives.

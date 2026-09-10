# 0123-beacon-fitted - ×18 decision diagnostics (beacon)

Shepherd ruling: (b)-first resolved against (δF=1.84, ×18 miss);
(a) activates with quantified target δF 1.84 → 0.1. This attempt is the
DECISION instrument, not a build: three cheap probes (no solves, no
mesh spend) that output either a quantified (a)-build spec with a
credible ×18-closure path, or a no-go certificate. Both bankable.

## Frozen design (analysis-only; no production runs governed)

- S1 null-projection: overlap of the 0.0066 soft mode with the discrete
  analytic translation modes (∂u/∂r, ∂u/∂z via basis gradients);
  project the translation subspace out of δu/δF (licensed zero-mode
  gauge); report δF_before → δF_after.
- S2 layer-R decomposition: ‖R‖² inside the layer band (|src|<δ,
  δ=2ε) vs bulk + quadrature-error scale at measured p≈0.4. Decides
  whether ANY mesh move has a mechanism (layer-dominated) or not.
- S3 bound-sharpness: controlled trust-r2→r3 step; actual |Δλ_min| vs
  the 2·s_ω·δF prediction. Decides whether the ×18 must come from the
  member (tight) or the estimate (loose).
- Verdict rule: (a)-build spec iff a probe shows a ≥×18 reachable
  factor with a named mechanism; else no-go certificate (mesh
  mesh-independence across 3 meshes + 140 iters + S2 bulk evidence +
  S3 tightness).

## Verdict (measured 2026-09-10, x18_probes.py exit 0)

Neither blind build nor no-go. The x18 gap is ESTIMATE-side:
- S1: softest fitted mode overlaps the analytic translation plane at
  0.89 — the 0.0066 mode is the discretized translation zero-mode.
  Licensed gauge: dF 1.84 -> 0.91 (x2.0). ||du||_2 = 0.89 plain.
- S3-full: on the controlled r2->r3 step, ||dJ||_2 = 0.027 vs the
  2*s*dF bound 1.13 — looseness x42. The current bound overstates the
  true Jacobian perturbation by x42 where both are measurable.
- S2: 65% of |R|^2 in the layer band — backup quad/mesh branch kept.

## (a)-build spec: Lipschitz-Weyl lemma variant (no solve/mesh spend)

Re-derive the 0121 perturbation conclusion as |dlam_2| <= L_J * ||du||
with L_J = sup|d(jf)/du| over the trust tube (explicit from
jf = 6*EPS^-2*s^5*ds; s, ds closed-form in src) + Weyl (J symmetric).
Inputs banked: ||du|| = 0.89, S1 gauge, R1 instrument as verdict pipe.
Success threshold: L_J <= 1 gives margin ~x10 vs lam_2 = 11.2.
Transfer assumption (x42 applies at the fitted state) must be PROVED
inside (Lipschitz on tube + direct evaluation), not assumed — if the
proved L_J misses, fall back to S2 quad branch, then no-go.
Next attempt executes the re-derivation; 0123 closes as the decision.
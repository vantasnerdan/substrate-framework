# G-b bridge scope (beacon 0115)

G-b = turn the G-a exact fixed-`J_ren` curve into the nonlinear
uniform-Lipschitz contradiction via the packetwise same-target bridge.
Sources: 0095 (36bh)–(36bi.1) (conditional bridge, author scope);
0104 (bridge "valid conditional", antecedent open); 0107 §§8–9 (curve +
packetwise estimate (45) vs reviewed (46)); pause-state (fixed-(j,N)
finite-time energy argument as next construction).

## Exact requirements

- R1 (antecedent, = G-a output): exact corrected smooth curve in the
  fixed-`J_ren` leaf with PRECISELY the raw packet's linear tangent;
  generator TWO derivatives smoother than the packet (covers the
  derivative count in R2).
- R2 (linear remainder): with `Y_τ = (Z_τ−Z_g)/τ` from the curve and `Y`
  solving the full constrained linearized Euler–Maxwell–tag–Gauss system,
  `R_τ = Y_τ − Y` obeys the same linear principal operator with source
  `τQ_2(Y_τ,∇Y_τ) + τQ_EM(Y_τ,Y_τ)` (tag + longitudinal-Gauss
  substitutions included); a UNIFORM `H^{s+1}` bound on `Y_τ` over fixed
  `[0, jT_*]` plus the standard `H^s` product estimate gives
  `sup‖R_τ‖ → 0` in `X_DA,loc^s` + field/tag terms (36bi). Needs:
  constrained-linearized finite-time well-posedness in that topology.
- R3 (quotient + contradiction): modulation derivative vanishes for
  `N_ell` outside the finite character set while `J_loc` detects the
  same-norm gain; dividing the hypothetical all-time Lipschitz bound C by
  τ at fixed (j,N), then `N → ∞`, yields (36bi.1), contradicting the
  reviewed linear lower bound (36bg)/(46). Packetwise only — no full
  operator-norm claim (0107 (45) wording).
- R4 (correction propagation): the finite-rank smooth correction
  (raw → corrected tangent) propagates to `o_N(1)` fixed-time observed
  remainder with weak-null preserved (0107 A4).

## Candidate constructions

- C1 (consume): feed the G-a curve directly into (36bi)–(36bi.1).
  Cheapest. Risk: curve membership must match R1/R2 hypotheses EXACTLY
  (smoothness class, same-tangent, H^{s+1} uniformity) — audit, don't assume.
- C2 (fixed-(j,N) energy argument): prove R_τ → 0 by direct energy
  estimates on the linearized constrained system (pause-state priority).
  More work; robust to C1 topology mismatch; parallelizable on MODEL
  linear systems before the G-a curve lands.
- C3 (Egorov-plus-correction): extend Unit-E fixed-time Egorov
  (`O_T(N^{−1)`, 0104-E) to corrected tangents via linearity + weak-null.
  Covers R4 separately; composes with C1 or C2.

## Failure order (most likely first)

- F1: uniform-`H^{s+1}` count fails — linearized Euler–Maxwell loses a
  derivative the +2-smoother generator doesn't cover (free-boundary /
  variable-coefficient / Gauss-slaving terms). Check the count FIRST.
- F2: propagated correction not `o_N(1)` — longitudinal Gauss slaving
  along the curve breaks the character orthogonality the raw packet enjoyed.
- F3: modulation-character mismatch — corrected curve's symmetry content
  vs reviewed `D`'s finite character set.
- F4 (DE-RISKED, not a failure): j,N uniformity — (44) fixes j, then tube,
  then N, then τ→0, then N→∞. No uniform estimate needed. Do not pursue.
- F5 (scope guard): all finite-time `[0, jT_*]`; any all-time reading is
  out of scope by construction.

## Interfaces

- Consumes: G-a curve (blocked → G-a2). C2-model work + F1 derivative
  count can proceed in parallel (unblocked).
- Computational support if needed: `substrate_framework.numerics`
  (`solve_ivp_evidence`, `solve_method_of_lines`, `refinement_study`)
  exists — not invoked at scoping scope.
- Success state: (36bi.1) contradicts (36bg) on the corrected sequence =
  fixed-member nonlinear uniform-Lipschitz obstruction (route-scoped,
  cf. 0107 §12; not instability, not P2).

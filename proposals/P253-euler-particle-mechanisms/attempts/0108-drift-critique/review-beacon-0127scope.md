# drift firewall review — 0127 native scoping design (c9154726): SCOPE PASS

Design-only attempt (no compute — correct for scoping). Inputs use corrected numbers
throughout (G1 4–5, dF 1.84, ρ*≈0.03, UNSATISFIED) ✓. Thesis relates every route to the
0124 mechanism explicitly. No EM import (R-EM2 declined-clean cited; R-C fenced to
derive-or-die) ✓. Order cheapest-first, HOLD→drift-firewall jurisdiction named ✓,
non-scope explicit ✓. Two findings ride as Phase-2 amendments (no block on charter).

## FINDING 1 (material): R-A is likely dead on arrival — sequence accordingly

R-A HOLD needs min-Q-over-tube ≥ 5.0. Banked G2' numbers already measure exact-transfer
variation 9.1 (central), 3.5 (half), 37.6 (double+corner) over the same E-scale tube —
so a direct survey will almost surely find min-Q ≲ 2 and fire the STOP. This is not a
design flaw (kill-or-hold gates may kill; a kill with a min-Q map banks the obstruction),
but shepherd should sequence EXPECTATION not just order: R-A = fast kill (run it — one
survey, assemblies only), real action = R-B/R-C. Amend before Phase 2: specify the tube
radius E explicitly (G2's 2‖du‖ inherited, or smaller with stated protection cost —
smaller E weakens the claim and must say so).

## FINDING 2: R-A survey dimensionality unstated

"2^d corners + 32 interior" in a 1451-dof tube needs the coordinate reduction named:
WHICH d coordinates (collective? nodal patch? eigen-directions?) and why those span the
adversarial directions. Unstated, a HOLD on a subspace survey does not protect the gap.
Sampling-doubling (<10%) is a Cauchy check, not a missed-direction bound. Name the
coordinates in Phase 2 design, or R-A HOLD licenses nothing.

## Sound (do not change)

- R-B premise grounded (slope→G1 + xnodes from sweep ✓); 20%-at-half-scale HOLD and
  >2×-miss-fall-to-R-A STOP both quantified ✓. Front representation properly deferred
  to build (diffuse-regularization care noted for the builder).
- R-C stop (no PSD native B without EM import → dead + documented obstruction) is the
  honest guard; "frozen-flux-type" is analogy language — watch item for Phase 2
  derivation review, not a current violation.
- R-D fenced as cross-check-only ("never a build") ✓ — prevents diagnostic-as-route creep.
- C(e)/ρ* reused unchanged as per-state checkers ✓ (no instrument drift).

## Verdict

SCOPE PASS: program thesis beats-or-sidesteps the mechanism on paper with frozen
gates/stops. Findings 1–2 amend Phase 2 design (E specified; survey coordinates named;
R-A expected-fast-kill). 0128 watch continues.
## Appendix: RED-TEAM (adversarial; not verdicts)

Lane: steelman each route's best attack, then break it. Rank by kill-resistance
(ability to survive the break), distinct from run-order (cheapest first).

### R-A tube-uniformity — KILL-RESISTANCE: LOWEST (dead today)

Steelman: refuse the transport game entirely — protection needs only Q bounded away
from zero on the tube containing the true state, and the tested pipeline makes each
assembly seconds-cheap. A HOLD would be a rigorous numerical certificate needing no
lemma. Best attack: shrink E to the smallest tube containing u* so even O(1)-nonsmooth
response cannot move Q far.
Break: E cannot shrink below ||du||~0.2 (measured state error — the tube must contain
u*), and banked G2' already measures variation 9.1/3.5/37.6 across that scale: central
alone puts min-Q ≈ 11−9 ≈ 2 < 5.0 HOLD. The escape (variation concentrated away from
the min) fails — the CENTRAL point already moved 9. Double-kill: FINDING 2 means even a
HOLD on an unnamed subspace would license nothing. Certain kill, zero cost, before one
assembly is spent. Run it anyway: cheapest verdict in the program + banks the min-Q map.

### R-B shape-calculus split — KILL-RESISTANCE: HIGHEST (only live defense)

Steelman: the ONLY route that beats the mechanism instead of sidestepping it. Sweep
data proves the bulk linearizes cleanly; all nonlinearity is front motion, a
codimension-1 object, and Hadamard calculus exists precisely for this. Best attack:
advect the ALREADY-BUILT 43-point contour (0122 machinery) by du's normal component
and reassemble — semi-Lagrangian transport that never differentiates across the kink.
Break: the front is DIFFUSE (1e-6 regularization + s^6 skirt) — "front position" is
defined up to the regularization width, and skirt-curvature lives in NEITHER half of
the split. The shape correction must supply (3.49−2.6)≈0.9 within ±0.7 while the skirt
moves O(1) — marginal, and the 20%-at-one-scale bar risks point-tuning. This break is
EXPERIMENTAL, not logical: no banked number kills R-B; its own STOP adjudicates
honestly. Strengthener (offer): pre-register quarter-scale (1.51, already printed) —
two points defeat tuning for one more assembly. Fix the fallback: R-B-miss currently
falls to R-A (a corpse per above) — fall instead to documenting the skirt obstruction,
which MOTIVATES R-C.

### R-C native stiffening — KILL-RESISTANCE: MIDDLE (unkillable in principle, priced low)

Steelman: the strongest move — change the object, not the estimate. PSD B ≥ 9 from
coercivity moots state accuracy ENTIRELY (wall evaporates at any res floor; conditional
lemma unnecessary; whole estimate program retired). No logical refutation exists
without attempting the derivation — the adversary cannot kill a promise, only price it.
Break (pricing): Euler transport gives CONSERVATION (Kelvin equality), protection needs
INEQUALITY, and Q is a diagnostic budget, not a Hamiltonian Hessian — B must live in
the right space, PSD, ≥9, with zero EM structure. That is the hardest analytic object
in the program with no sketch. Adversary's kill: "exhibit B" — no defense today.
Survives every banked number (nothing contradicts it), survives logic (coercivity
would moot all), dies only on execution probability.

### Ranking (kill-resistance): R-B > R-C > R-A

R-B survives all but an experimental question its own STOP answers. R-C survives logic
and all banked numbers, dies only on exhibit-demand. R-A is already dead (banked G2').
INVERSION vs run-order is intended, not contradictory: run R-A first BECAUSE it is a
cheap certain kill; then R-B (the live defense); R-C derivation prices in background.
Design amendment: R-B fallback retargeted R-A → skirt-obstruction record (R-C motive).

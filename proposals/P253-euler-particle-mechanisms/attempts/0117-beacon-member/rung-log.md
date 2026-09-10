# 0117 rung log (beacon, append-only)

Skill governing: small-ratio-numerics (design freeze, error budgets,
observed order, jitter, virial-type monitors, reproducibility). Heavy
machinery NOT triggered: witness bounds are O(1) well-conditioned
(λ_F/λ_2 far from floor) — applicable prescriptions recorded in design.

## Ladder (each rung: prediction → result → mechanism)

- R0 Picard, frozen (μ,c), bump 0.3: predicted sustain. GOT collapse to
  u≡0 in 26 iters (236 s). Mechanism: bump below μ-barrier
  (0.3 < cr²/2+μ ≈ 0.3); source never switches on.
- R1 Picard, bump 1.5: predicted sustain. GOT collapse in 29 iters.
  Mechanism: unbordered fixed-point contracts to trivial root (later:
  nontrivial branch repelling for Picard on focusing p=6).
- R2 API bug `ddot` vs `dot`: found via 400× stiffness scale (diag
  1385 vs ~3.5; premade `laplace` uses `dot`). Explains R0/R1 COLLAPSE
  SHAPE (overstiff operator) but not the selection problem. Fixed.
- R3 Picard corrected operator: BLOWUP to NaN (overflow). Mechanism:
  focusing nonlinearity, α=0.5 overshoot. Authorized: α=0.1 +
  μ-continuation + NaN guard + residual cadence.
- R4 Picard α=0.1 + μ-chain 0.8→0.227: geometric convergence (rate
  ~0.72) to u≡0 at EVERY rung, res 1e-10. Mechanism: zero is the
  attractor; nontrivial branch Picard-repelling. Authorized: Newton.
- R5 Active-set Newton, frozen (μ,c): umax≈1.2 SUSTAINED, res
  3.8e-2 → 4.8e-3 stall, active set frozen. Verdict: nontrivial
  structure EXISTS at leading-jet params.
- R6 Regularized Newton (reg=1e-3): IDENTICAL stall → floor not
  free-boundary chatter. Recorded.
- R7 Bordered Newton (κ, I_z): immediate STALL, dp tiny. Mechanism
  FOUND by FD: ∂(κ,iz)/∂(μ,c) rows parallel (3.23 vs 3.25) — r≈r² over
  thin core at R=1. Fix: 0080 BR-border (κ, mean-radius).
- R8 Sign bug dF/dparam: FD proved 200% error (7.5 vs 3.75 scale);
  fixed to +M3@jf (verified 8e-7/3e-8). Lesson: FD-verify every
  assembly (Crow was already exact).
- R9 Bordered (κ, rbar): res 0.735 → 0.020 in 3 iters; kap=1.018,
  rbar=0.996, iz=3.23 (≈π post-hoc ✓), μ=0.44, c=0.066. Then stall.
- R10 Line-search autopsy: ray descends monotonically to just above
  cur — true STALL (ascent direction), plus found+fixed apply-on-fail
  bug (was applying worsening micro-steps = the earlier "drift").
- R11 Jacobian audit: B/Crow FD-exact; J-check script itself was wrong
  (compared J@d vs quotient — 1e7 false alarm); corrected check shows
  J consistent to FD-truncation order. Solves exact (1.8e-14).
- R12 p-continuation 2→6: p=3 reaches 3.3e-2; p=4/5/6 stall 0.09–0.21
  with identical ray signature. Floor is p-independent → systematic.
- R13 (running at log time): fine-mesh (80×40) bordered cold —
  tests discretization-floor theory.

## Open after R13

- If fine converges: Maxwell O(g) stage → ga_pipeline feed → G-a2 DONE.
- If fine stalls same: land as blocked-with-mechanism (merit plateau at
  2% rows) + feed APPROXIMATE member to pipeline (integration proof) +
  specify trust-region/fitted-mesh next rung.

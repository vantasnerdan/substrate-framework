# S4 F-bar freeze (pre-compute; D-08 form) — relational-tangle charge

Setup: vacuum = statistically homogeneous ISOTROPIC tangle: N_t thin rings
(radius r_t, circulation Γ_t = 1), centers uniform in ball R_b, orientations
uniform on SO(3). Carrier = unit circle (R=1, z=0) + rotation checks.
Linking via Gauss double-sum over straight segments (carrier M=64, rings K=16).
Line density Λ = N_t·2πr_t/(4πR_b³/3). Postulate (one free choice, stated):
Λ_vac fixed by r_t=0.3, N_t set so Λ≈1.0 at R_b=3 (order-unity tangle).
Statistics over S=32 seeds: mean μ, std σ of Lk; mean-abs m1, its CV.

- F-S4a (generator honesty): |μ| ≤ 2σ/√S. Violation = biased generator (fix
  code, not physics). Gate before all verdicts.
- F-S4b (stable value): candidate q = m1 (mean |Lk|). KILL iff CV = σ_|.|/m1
  > 0.50 (no stable value at order-unity tangle) — margin 2× above pass.
  PASS iff CV ≤ 0.25 AND R_b-scaling: q(R_b=5)/q(R_b=3) within 2× (not a
  boundary artifact).
- F-S4c (no preferred frame — F4): carrier rotated by 3 orthogonal rotations;
  q must agree within 2σ seed-noise. KILL on larger anisotropy (frame
  smuggled).
- F-S4d (Λ-dependence): q at Λ×2 (N_t doubled): PASS-leaning iff exponent
  α (q ∝ Λ^α) ≤ 0.25 (carrier property); KILL-leaning iff α > 0.5 (tangle
  property). Gray (0.25, 0.5] = UNRESOLVED.
- Sign rule (pre-registered): tangle statistics are sign-blind (isotropy);
  charge SIGN comes from carrier framing integer (M1 label) × magnitude q —
  hybrid S1×S4 form. If magnitude stable but sign needed from framing, record
  S4→S1×S4 REDUCTION (not a kill).
- Gray: any single KILL-line crossed marginally (<20% over) with others
  passing = UNRESOLVED + named refinement (more seeds), no verdict.
Run reports (μ, σ, m1, CV, α, rotation spread) + verdict line, nothing else.

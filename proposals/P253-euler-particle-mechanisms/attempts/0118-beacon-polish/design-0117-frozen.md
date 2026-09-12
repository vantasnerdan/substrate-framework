# 0117 frozen design: member-field build (beacon)

Status: FROZEN before production solves (skill: small-ratio-numerics —
design freedoms frozen first; exploratory samples stay labeled).
Dated changes only, never silent edits.

## Object (skill: construct first)

Uncharged axisymmetric ring member, g=0 stage (ω_g = ω_0 + O(g²) per
0080 (38); Maxwell O(g) stage follows a converged ψ_0):
`-Δ*u = r²·ε⁻²(u − c·r²/2 − μ)₊^p`, `u := ψ → 0` at ∞ (0077 (19) via
`Δ*(r²/2) = 0`). Weak form with measure `r dr dz` (no axis singularity):
`∫∇u·∇v·r = ∫r³·f(u)·v`, `f = ε⁻²(u−cr²/2−μ)₊^p`.
Domain `D = [0,Rmax]×[0,Zmax]` (even-z half-slice: ∂_z u = 0 at z=0 —
also removes the axial-translation zero mode by BC licensing);
`u = 0` on outer boundary (decay approx, box-sensitivity measured).
Discretization: scikit-fem 12.0.2 MeshTri + ElementTriP1, r-weight baked
into forms via quadrature-point coordinates.

## Frozen parameters (source-derived, not fitted)

Nondim R = 1, κ = 1; p = 6 (0077/0080 floor); ε_core = 0.15;
leading jet (euler_cao_schur, 0108): L = ln(1/ε) ≈ 1.9, c = κL/4πR,
μ = 3κRL/8π. Initial guess: toroidal Gaussian bump at (R,0).

## Iteration + monitors (skill: residual, independent checks)

Damped fixed-point `u ← (−Δ_r)⁻¹[r³f(u)]`; Newton/continuation only on
recorded non-convergence. Independent monitors (NOT the residual):
(a) circulation reproduction `κ̂ = ∫ζ r dr dz vs κ = 1` (parameter must
reproduce — strongest license); (b) far-field dipole match (fitted
dipole moment vs computed I_z — asymptotic matching); (c) crossed h×R
refinement with OBSERVED order (no assumed p); (d) jitter (seed/roundoff
re-run). Reductions feeding QOIs via math.fsum. Threads pinned
`OMP_NUM_THREADS=1`, versions recorded.

## Tolerances (from measured scale)

Core u amplitude O(0.1–1) post-nondim; solver residual tol 1e-10
(relative to core scale); fixed-point increment tol 1e-8; error budget
itemized with result (background residual, h/R truncation with measured
order, BC/box, quadrature, iteration, jitter floor).

## Permitted verdicts

- STAGE-1 DONE: converged u + residual + κ̂-match + dipole-match +
  observed order + budget → feed ga_pipeline (ω_0, χ_0) → ga-status
  rows 1–3,6–7 numeric.
- MAXWELL STAGE (same attempt iff stage-1 lands): linear L_c solves
  (0077 (6)/(23a)) → Φ,H → B_g at fixed small g → rows 4–5,8–9 →
  G-a2 DONE.
- BLOCKED (with mechanism): record which monitor failed and the exact
  next construction (Newton/continuation/parameter move). Effort alone
  never closes G-a2.

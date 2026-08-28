## P247 — isolated fixed-`J` clock and the two-body interaction law

Advances #178.

Campaign branch for proposal `P247-isolated-clock-two-body-law`, attempts 0005–0010 (attempt 0001–0004 landed in #179). Each attempt is append-only under `proposals/P247-isolated-clock-two-body-law/attempts/` with a decisive verdict, materialized artifacts, and repository validation at the frozen boundary.

### Attempt 0006 — F1 Yukawa route refuted; Z₂ theorem (commit `78daf61`)
- Exact Z₂ evenness: `E[+χ] = E[−χ]` to machine precision on the root and randomized backgrounds, channels `(0,)`, `(1,)`, `(1,3)` → χ≡0 is the exact classical minimum for every background; **no classical static boost sourcing or exchange at any order**.
- Box-growth is channel-attributed to the S-potential (16.578 → 18.983 over R=20..30, accelerating increments); the boost sector carries exactly 0 at χ=0.
- No S row fits a massive tail (exp_resid 0.35–0.64) → the Yukawa pairing premise fails; the law is slow S-overlap.

### Attempt 0007 — exact V extracted; de-boxed law measured (commit `823186e`)
- Exact potential `V(S) = −½trS² − trS³ + (trS²)² + 0.5`; pointwise machinery matches the functional potential integral to **8.5e-13**; `V(P) = 0` exactly on the projector family.
- Projector Hessian: no moduli (m² = 5/3/3/6; gauge tilts exactly flat in V).
- De-boxed two-clock law: **attractive, saturating** ΔE(d) = 8.113 − 5.668·e^(−0.2442·d) (ssr 6.7e-4), box/resolution-refined γ = 0.2524/0.2571, E∞ stable to 0.1%. Box-growth mechanism = branch spreading on the V zero set.

### Attempt 0008 — ghost diagnosis (commit `7d6452e`)
- Full fixed-R Hessian at the R=24 root: 21/96 negative directions, λ_min = −1.9471×10⁴, softest mode **pure χ-row** (FD-verified at three step sizes). The attempt-0005 σ-truncation reduction defect is the measured leading instability; block-sum Hessians that omit zero-valued terms are wrong (max|H₁g − H_A| = 12354).

### Attempt 0009 — χ-sector repair verified (commit `08da787`)
- Repair: σ_density = +KAPPA·η(du,du) (census-grounded by C-M5S-001's positive 1/16 boost entry). Five gates: S-sector regression exact; σ positive at finite χ; Z₂ survives exactly; **ghost repaired** (Morse 0, min +1.381856e-05, 0/96 negative); source closure survives exactly.
- **Width quasi-moduli identified:** the softest repaired mode (1.382e-5; next 3.96e-3) is pure split-row with 0.8645 family-drift overlap.

### Attempt 0010 — fixed-width existence probe: decisive no-go (commit `4c5fb96`)
- Width-stiffness ladder (repaired class, 0 negative at every rung): 8.35e-6 (R20) → 1.38e-5 (R24) → 2.85e-5 (R30) — the softest mode **stiffens** with box; the minimum branch is parametrically box-proportional (core box-stable, tail diluting, E(R) growing with accelerating increments) → no finite-energy isolation limit.
- Quadratic-tail reprojection of the R=24 root: R28 **converged stationary narrow configuration** (relgrad 1.799e-13, E = 65.705) is a **high-order saddle** of the repaired functional — Morse index 14, λ_min = −2591.3021, FD-confirmed (−2579.4 at ε=1e-3, 0.5% agreement, χ-row weight exactly 0). R26 non-converged (relgrad 5.6e-2, no conclusion drawn). R30/R36 aborted by the seed-fidelity gate (order-24 boundary-layer limit).
- **Verdict: the isolated de-boxed fixed-`J` clock does not exist in the extended class.** The winding channel is already inside the committed functional (`derivative_phi` via `rotation_z`), so the confined-class mechanism does not mass the width mode.

### State of the issue gate
- **Isolated clock:** realized in the confined (window-assisted) class — stable, radiatively stable (C-M5S-009), measured law (C-M5S-008). De-boxed isolation refuted (0009/0010).
- **Interaction law:** resolved with mechanisms — Z₂ excludes classical boost exchange at any order; confined law = repulsive power `+456.6·d^(−1.696)` (winding composition); de-boxed law = attractive saturating exponential (S-overlap branch structure).

### Completion program (registered as next_loop 0011 in `proposal.yaml`)
Draft claims C-M5S-016 (Z₂ evenness theorem), C-M5S-017 (σ reduction-defect repair record), C-M5S-018 (width quasi-moduli obstruction), C-M5S-019 (narrow-saddle refutation), C-M5S-020 (de-boxed law); run the bounded individual reviews; promote C-M5S-015 + accepted new claims into release v0.169.0; update this PR with the promotion and the completion statement.

Validation: `PYTHONPATH=src python3 scripts/validate_repository.py` → `WORKFLOW VALID: 236 claims, 236 accepted, 11 proposals, 4 skills` at `4c5fb96`.

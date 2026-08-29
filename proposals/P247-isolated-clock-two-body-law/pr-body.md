## P247 — isolated fixed-`J` clock and the two-body interaction law

Advances #178.

Campaign branch for proposal `P247-isolated-clock-two-body-law`, attempts 0005–0010 plus the completion transaction (attempt 0011) (attempt 0001–0004 landed in #179). Each attempt is append-only under `proposals/P247-isolated-clock-two-body-law/attempts/` with a decisive verdict, materialized artifacts, and repository validation at the frozen boundary.

### Attempt 0006 — F1 Yukawa route refuted; Z₂ evenness (commit `78daf61`)
- Exact Z₂ evenness: `E[+χ] = E[−χ]` to machine precision on the root and randomized backgrounds, channels `(0,)`, `(1,)`, `(1,3)` → χ≡0 is an exact stationary point in the χ directions on the measured backgrounds; **no odd-order classical static boost sourcing arises on these backgrounds** (even-order static responses are not addressed).
- Box-growth is channel-attributed to the S-potential (16.578 → 18.983 over R=20..30, accelerating increments); the boost sector carries exactly 0 at χ=0.
- No S row fits a massive tail (exp_resid 0.35–0.64) → the Yukawa pairing premise fails for the tested frozen-superposition observable; the measured cross-term is slow S-overlap.

### Attempt 0007 — exact V extracted; frozen-superposition cross-term measured (commit `823186e`)
- Exact potential `V(S) = −½trS² − trS³ + (trS²)² + 0.5`; pointwise machinery matches the functional potential integral to **8.5e-13**; `V(P) = 0` exactly on the projector family.
- Projector Hessian: no moduli (m² = 5/3/3/6; gauge tilts exactly flat in V).
- Frozen-superposition potential cross-term of two frozen R=24 profiles in a finite cylindrical box: **attractive over d = 6..28, saturating**, ΔE(d) = 8.113 − 5.668·e^(−0.2442·d) (ssr 6.7e-4); box/resolution refinement (2x grid: γ = 0.271731, E∞ = 8.08826; 1.5x box: γ = 0.265034, E∞ = 8.09501) shifts the exponent by less than 12% and E∞ by less than 0.4%. This is a potential-telescope observable of frozen superposed profiles — not a relaxed two-centered solution and not the full two-clock energy (kinetic/fixed-J terms, independent phases, and inertia sector are not included).

### Attempt 0008 — ghost diagnosis (commit `7d6452e`)
- Full fixed-R Hessian at the R=24 root: 21/96 negative directions, λ_min = −1.9471×10⁴, softest mode **pure χ-row** (FD-verified at three step sizes). The attempt-0005 σ-truncation reduction defect is the measured leading instability; block-sum Hessians that omit zero-valued terms are wrong (max|H₁g − H_A| = 12354).

### Attempt 0009 — χ-sector repair verified (commit `08da787`)
- Repair: σ_density = +KAPPA·η(du,du) (census-grounded by C-M5S-001's positive 1/16 boost entry). Five gates: S-sector regression exact; σ positive at finite χ; Z₂ survives exactly; **ghost repaired** (Morse 0, min +1.381856e-05, 0/96 negative); source closure survives exactly.
- **Width quasi-moduli identified:** the softest repaired mode (1.382e-5; next 3.96e-3) is pure split-row with 0.8645 family-drift overlap.

### Attempt 0010 — fixed-width existence probe: no-go on the tested constructions (commit `4c5fb96`)
- Width-stiffness ladder (repaired class, 0 negative at every rung): 8.35e-6 (R20) → 1.38e-5 (R24) → 2.85e-5 (R30) — the softest mode **stiffens** with box; the minimum branch grows monotonically without plateau across the tested R = 20..30 ladder (core box-stable, tail diluting, E(R) rising with accelerating increments).
- Quadratic-tail reprojection of the R=24 root: R28 **converged stationary narrow configuration** (relgrad 1.799e-13, E = 65.705) is a **high-order saddle** of the repaired functional — Morse index 14, λ_min = −2591.3021, FD-confirmed (−2579.4 at ε=1e-3, 0.5% agreement, χ-row weight exactly 0). R26 non-converged (relgrad 5.6e-2, no conclusion drawn). R30/R36 aborted by the seed-fidelity gate (order-24 boundary-layer limit).
- **Verdict (probed constructions): no tested construction realizes an isolated de-boxed fixed-`J` clock in the extended class.** Deeper-box and non-spherical existence remain open, with the missing constructions recorded in the claims. The winding channel is already inside the committed functional (`derivative_phi` via `rotation_z`), so the confined-class mechanism does not mass the width mode.

### State of the issue gate
- **Isolated clock:** realized in the confined (window-assisted) class — stable, radiatively stable (C-M5S-009), measured law (C-M5S-008). De-boxed isolation: no tested construction on the tested R ≤ 30 ladders supports it (monotone growth without plateau on the tested minimum branch R = 20..30; probed narrow stationary points are high-order saddles). Whether a de-boxed finite-energy stationary clock exists in deeper boxes or on non-spherical branches is open, with missing constructions named (C-M5S-015/018/019).
- **Interaction law:** confined law = repulsive power `+456.6·d^(−1.696)` (accepted, C-M5S-008). Frozen-superposition potential cross-term = attractive saturating exponential over d = 6..28 (scoped observable, C-M5S-020 — not the relaxed two-clock law). Odd-order classical static sourcing is excluded on the measured backgrounds (C-M5S-016); even-order static responses and the full two-centered construction are the registered open objective (preregistered in `attempts/0006/f1-g1-specification.md` and the proposal's next-loop program).

### Completion transaction landed (attempt 0011, commits `3f9d6a2` + `4c13960`, review corrections on top)
- Six claims individually reviewed (independent reviewer agents, one substantive pass each), repaired per findings, verified by one bounded correction check (**CORRECTIONS_VERIFIED**), promoted, and then narrowed once more in response to the PR #182 review (blocker-level scope repairs): **C-M5S-015** (milestone-0 adjudication, measured-ladder scope), **C-M5S-016** (Z₂ evenness, odd-order/single-exchange scope), **C-M5S-017** (σ-sign repair record; accepted as reviewed), **C-M5S-018** (width quasi-moduli, tested-ladder scope), **C-M5S-019** (narrow-saddle refutation, probed-rung scope), **C-M5S-020** (frozen-superposition potential cross-term, d = 6..28).
- Review #182 corrections applied: asymptotic nonexistence wording replaced by measured-range statements (or equivalently an explicit missing-construction note) in C-M5S-015/018; C-M5S-020 rescoped from "de-boxed two-clock interaction energy" to the tested frozen-superposition potential telescope; C-M5S-016 restricted to odd-order sourcing with the hidden C-M5S-017 assumption removed (dependency graph acyclic); PR body matched to the corrected claims.
- **Release v0.169.0** pins all 242 accepted claims at `substrate-framework@3f9d6a2`; docs and generated memory rendered.
- Review records: `proposals/P247-isolated-clock-two-body-law/reviews/` (six reviews + correction checks `correction-check-0011.md`, `correction-check-0182.md`).

Validation: `PYTHONPATH=src python3 scripts/validate_repository.py` and `scripts/validate_changed.py` against base `67fcf0d` at the final corrected head; content-addressed receipt: `proposals/P247-isolated-clock-two-body-law/evidence/validation-receipt.yaml`.

# G-a constant status (beacon 0111)

Notation per 0107 derivation; reviewed sources: 0077 (field structure),
0080/0084 (branch existence), 0107 §§4–6 (witness lemma + matrices).

## Per-constant table

| # | Constant | Definition | Data needed | Status |
|---|---|---|---|---|
| 1 | `Q_ω, λ_ω, s_ω` | (15)(17), `F = ω_g` | ω_g field (L² + tails) | PIPELINE-READY, data-awaiting |
| 2 | `R_ω, ε_ω` | (24) cutoffs from η_ω budget | same | PIPELINE-READY (auto-chosen, cf. R=3.0 on test field) |
| 3 | `H, det H, ‖H⁻¹‖` | (29)(30), `c = ρ_m` | same | PIPELINE-READY (M_F path covers H identically) |
| 4 | `Q_B, λ_B, s_B` | (15)(17), `F = B_g` | B_g field | PIPELINE-READY, data-awaiting |
| 5 | `G, det G, ‖G⁻¹‖` | (26)(27), `c = ε_EM` | same | PIPELINE-READY; g-ledger: `λ_B = O(g²)`, NO g-uniform inverse (recorded, matches 0107) |
| 6 | `C = I_3` | (32)–(34): cutoff core translations, `∫ω_g = 0` | support locations only | READY NOW (structural; 0107 verifier checks the core potentials axis-by-axis, U2) |
| 7 | `Q_χ ≠ 0` | tag mass (31) | — | ESTABLISHED by construction: 0077 (20) normalizes `∫χ_P = 1` |
| 8 | (25a) orbit integral | `∮ gχ_g u_g·dx ≠ 0` on closed strict-band streamline | W_g orbits + χ_g support | FORMULA-READY (gradient-vs-loop argument is data-agnostic); evaluation awaits fields |
| 9 | `K_I, K_C` | (37) tag-Gauss displacement derivatives | χ_g gradients + Gauss kernel | PIPELINE-EXTENSIBLE (same quadrature path); awaits fields |

## g-ledger (leading orders, from 0080 (38) + 0077 (23a))

`ω_g − ω_0 = O(g²)` ⇒ H g-independent at leading order;
`B_g = O(g)` ⇒ `Q_B = O(g²)`, `G = O(g)`; `det M_leaf = detH·detC·detG ≠ 0`
at fixed `g > 0`, degenerating as `g → 0` (separate automatic-zero row at
`g = 0`). No uniform-`g` claim (0104 Unit H boundary preserved).

## Why G-a2 is a separate construction

The member is an IFT output, not a formula: evaluating rows 1–5, 8–9 needs
a NUMERICAL charged-branch member (discretize 0077 (21)–(23) + 0080
bordered map, solve at fixed `(κ, I_z)` in window (37), quadrature the
functionals). That is a frozen-design numerical construction governed by
`small-ratio-numerics` (error budgets, λ-floor, observed-order) — a
shepherd-tasked build, not a beacon source audit. G-a1 hands it: exact
functional specs (rows 1–9), a tested pipeline (feed arrays → constants),
and the g-ledger. Estimated G-a2 input: field arrays + norm certificates;
output: all nine rows numeric with error bars.

## Verdict

G-a1 DONE. G-a2 BLOCKED (missing: numerical branch member). Unit G
completion remains blocked at G-a2 + G-b; nothing in 0111 narrows or
widens those verdicts.

# Unit G witness status (beacon 0109)

Reference: 0104 Units G (blocked) + H (established linear); 0107 derivation
(interrupted author draft — repaired/reused here, never asserted as conclusion).

## Formally supplied by 0107 §§2–8 (re-verified as algebra at 0109 boundary)

| Item | Content | 0109 evidence |
|---|---|---|
| Balances | `dI_h/dt = ∫f`, `dP_EM/dt = −∫f`, `dJ_ren/dt = 0` on fixed finite intervals via cutoff flux limit; velocity-integral `∫v = 0` typed as L1-domain diagnostic, not a second momentum | 0107 receipt: 10/10 algebra checks pass, exit 0 (tool-receipts U2). Checks: Maxwell sign identity, div-curl Gauss row, impulse curl factor+sign, moment-tensor reconstruction, 3 center potentials, M_leaf scalar det, K-center 16 |
| Witness lemma | `y ↦ P_L(F×y)` injective for nonzero localized divergence-free `F`; compact smooth columns via cutoff + mollification with two independent error budgets `η_F(R) ≤ λ_F/4s_F`, `κ_F(R,ε) ≤ λ_F/4s_F`; `σ_min(M_F) ≥ |c|λ_F/2` | Structure re-verified: G1 (9×9 staged-block det instance `det = detH·detC·detG ≠ 0`), G2 (Cao-jet import smoke pins carrier dependency). Analytic injectivity + cutoff convergence NOT re-proved — cited as 0107 derivation §§4–6 |
| Staged matrix | `M_leaf = [[H,0,0],[0,C,0],[K_I,K_C,G]]`, lower-triangular (K blocks retained, no block-diagonal claim); `det = detH·detC·detG` | G1 numeric-instance check + 0107 scalar-symbolic check (U2) |
| Packet orthogonality | High-character (`|n| ≥ 2`) exact annihilation of `ΔI_h / X_χ / J_ren` rows incl. frame shifts; non-axisymmetric covectors handled as `o_N(1)` via (38), not called exact | Cited as 0107 §8 derivation; no new check (needs packet construction, out of scope) |
| Curve + quantifiers | (41)–(43) exact fixed-`J_ren` curve, `γ = O(τ²)`; order fixed-j → tube → N → τ→0 → N→∞ | Cited as 0107 §§8–9 derivation; propagation bridge open (G-b below) |

## Remaining bridges (full Unit G completion BLOCKED, not refuted)

- **G-a — carrier-numeric witness constants.** `λ_ω, s_ω, R_ω, ε_ω` for `F = ω_g`
  and `λ_B, s_B, R_B, ε_B` for `F = B_g`, plus `Q_χ ≠ 0`, the closed-orbit
  `B_g ≠ 0` circulation integral (25a), and `det H·det G ≠ 0` evaluated on the
  fixed charged Cao member. Needs the actual member fields (not in
  `src/substrate_framework/`; 0095's transfer used profile cells + compact-tube
  asymptotics, never these constants). Exact spec of need, no substitute offered.
- **G-b — finite-time weighted propagation + packetwise same-target bridge.**
  0107 §12 states the nonlinear obstruction stays prospective until these pin
  together. Unchanged by this attempt.
- Activation note: 0107 README gates body work on coordinator schema replay
  (exit 0 on corrected hash `59cfe8d5…`); the momentum-replay receipt records
  replay exit 0 licensing execution under `J_ren`. 0109 adjudicates neither
  0107's activation state nor its §9/§12 prospective verdicts.

## Verdict

Unit G algebraic core: re-verified. Unit G completion: BLOCKED at G-a + G-b.
No refutation mechanism exhibited; no nonlinear contradiction claimed.
Next executable: G-a needs the fixed-member field data source named above.

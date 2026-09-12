# 0111-beacon-ga-field — G-a field-data source (beacon)

Bounded obligation (shepherd ack 9b8ebaff): supply the exact missing
carrier-numeric constants named in 0109/G-a. Surface:
`attempts/0111-beacon-ga-field/` only. No 0080/0077/0107 bytes touched.

## Headline finding (sourced, not invented)

The fixed charged Cao member fields `(ω_g, B_g, χ_g)` exist ONLY as IFT
outputs: 0080/0084 establish the branch via bordered IFT with
`ζ_g−ζ_0 = O(g²)`, `E_g,B_g = O(g)` (0080 (38); 0084 Unit E at declared
window). 0077 gives exact structural formulas (Grad–Shafranov (21),
Ampère primitive (16)/(23), `B = rH e_θ` (12), tag (20) with `∫χ_P = 1`,
uncharged profile (19) via Cao existence/concentration — no closed form,
no numbers). There is NO numeric member data anywhere in the campaign
tree (0111 survey: 0080/0084/0077/0095 hold receipts + derivations only).
Therefore G-a numbers cannot be computed today — and this attempt converts
that block into specified, tooled, sub-blocked work (see `ga-status.md`).

## Deliverables

- `ga_pipeline.py` — attempt-local implementation of 0107 (15),(18)–(25):
  Q_F Gram → λ_F, s_F → η_F(R) tail+Hardy budget → (24) cutoff choice →
  M_F response → σ_min bound check. Verified on a synthetic compact
  div-free torus field: ALL CHECKS PASS, exit 0 (receipts in
  `tool-receipts.md`, failures preserved).
- `ga-status.md` — per-constant table: definition, data needed, status
  (ready-now vs awaiting-G-a2), g-ledger, verdict.
- `tool-receipts.md` — first executions with exits.

## Verdict

G-a1 (specification + tested pipeline + g-ledger + analytic slices) DONE.
G-a2 (numbers on the true member) BLOCKED on numerical branch
construction — shepherd task under `small-ratio-numerics`, not beacon.

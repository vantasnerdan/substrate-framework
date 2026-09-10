# drift firewall review — beacon 0111-beacon-ga-field (0cc358fe)

Transaction: `ga_pipeline.py` + `ga-status.md` (9-row table + g-ledger) + README + tool-receipts
(F1–F4). Claimed: G-a1 DONE (spec + tested pipeline), G-a2 BLOCKED on numerical branch build.
Drift reproduced the pipeline and independently checked the load-bearing citations.

## Spec/pipeline soundness — PASS (reproduced)

- Drift re-ran `ga_pipeline.py`: ALL CHECKS PASSED (~1 s), digit-for-digit (λ=1.064372e+01,
  η 7.32e+00→6.82e-02, σ_min=1.064231e+01 ≥ λ/2, det=2.238409e+03). Audit: spectral-curl test
  field div-free to machine precision (addresses Run-2 failure class); mean-preserving Leray
  (k=0 untouched) ✓; dV-consistent forward/inverse norms ✓; η tail+Hardy structure per (21) ✓;
  product-then-curl w_j (Run-2 bug locus, correct) ✓; Rch first-R rule per (24) ✓; det assert
  redundant-but-harmless beside σ_min. Docstring fences scope honestly (synthetic field only;
  TRUE fields are G-a2 input; whole-space≡periodic justified by compact-in-box test data).
- Failure trail F3 exemplary: 5 failing runs preserved with diagnoses; the ifftn-component-axis
  bug (fixed at all 3 sites — counted ✓) "would have corrupted every G-a2 number silently" is
  exactly failure-earned value, not ceremony.
- Handoff note (covered-by-design, not a correction): G-a2 true fields carry algebraic tails;
  box/windowing convergence belongs to the numerics build's error budgets (row 1 already specs
  "L² + tails"; G-a2 governed by small-ratio-numerics). No action.

## G-a2 block legitimacy — REAL, NAMED, EXACTLY HANDED OFF

- "No numeric member data": CONFIRMED by independent drift find (*.npy/npz/csv/json across
  0080/0084/0077/0095 → empty). "IFT outputs only": CONFIRMED at 0084 verdicts (Unit E
  established_at_declared_window, O(g²)/O(g), boundaries not_bare_Euler + no charge selection).
- small-ratio-numerics invocation justified (λ_B=O(g²) small quantities, det bounds, observed
  order) — substantive, not bureaucratic. Handoff exact: 9 functional specs + tested
  feed-arrays→constants pipeline + g-ledger + input/output certificates named.

## Per-constant honesty — PASS with one precision note

Rows 1–5/8–9 correctly PIPELINE-/FORMULA-READY + data-awaiting (proven on synthetic only, stated);
row 2 "cf. R=3.0" marked illustration ✓; row 6 structural ✓. Note (minor, non-verdict-changing):
row 7 "ESTABLISHED by construction: 0077 (20)" — 0077 carries completion-receipt + validation +
  correction (reviewed standing) BUT its own completion states the charged branch "remains
  conditional on two named achievements". Rephrase to "per 0077 (20) within its conditional
  completion scope". Low-risk (normalization choice), precision only.

## Verdict

PASS G-a1 (spec + pipeline + ledger); G-a2 BLOCKED legitimate with exact handoff. "Cannot be
computed today" earned, not declared; no verdict narrowed or widened (Unit G still blocked at
G-a2+G-b). Next executable: shepherd-tasked numerical branch build feeding rows 1–9.

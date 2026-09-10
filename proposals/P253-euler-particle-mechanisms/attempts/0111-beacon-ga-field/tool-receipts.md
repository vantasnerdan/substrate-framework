# Tool receipts — 0111 (beacon)

## F1 — data survey (2026-09-10, exit 0)

`ls attempts/0080/ attempts/0084/` + `find ... -name *.json/*.npy/*.npz/*.csv`:
0080 holds derivation/receipts/verifier, NO numeric data files; same for
0084/0077/0095 (receipts + derivations only). Finding: no member field data
in the campaign tree — recorded in README, not assumed away.

## F2 — scope reads (exit 0)

`attempts/0084/verdicts.yaml#3A56` (full): Unit E charged branch
established-at-declared-window, O(g²)/O(g) orders; does_not_license
bare-Euler carrier, stability, P4/P5. `attempts/0080/derivation.md#F4F9`
§§5–6: Schur jet (34)(35), window (36)(37), branch orders (38), path (39),
hierarchy (40). `attempts/0077/derivation.md#0793` (full, 372 lines):
traveling Maxwell (5)–(9), closure (10)–(18), charged map (19)–(23b),
stabilizer lock (24)–(26), HSE/Schur conditional (27)–(32).

## F3 — pipeline runs (exit 0 final; failures append-only below)

Command: `.venv/bin/python .../0111-beacon-ga-field/ga_pipeline.py`
- Run 1: EXIT 1 — `ValueError: broadcast (3,64³)→(64³)` (S-array rank bug).
- Run 2 (after shape fix): EXIT 1 — `AssertionError B2-analog` here
  `smin=0.78 < λ/2=1.77`. Diagnosed via two debug probes (exit 0 each):
  (a) Poisson roundtrip err 1.58 ≈ ‖S‖ exposed k=0/mean inconsistency from
  a C⁰-cutoff test field (∫F ≠ 0 at O(5e-2)); (b) error decomposition
  ‖w−S‖ = 2.07 vs ‖(1−χ)S‖ = 0.0015 localized ALL error in ∇χ×A.
- Root causes found + fixed: (i) test field rebuilt via spectral curl of a
  tapered potential (div-free to machine precision); (ii) REAL bug exposed
  by the failure: `ifftn` on stacked (3,n,n,n) arrays transformed the
  component axis too — scoped to `axes=(1,2,3)` at all 3 sites (Sphys,
  Aphys; wj already scoped). This bug class would have corrupted every G-a2
  number silently — the failure earned its keep.
- Runs 3–5: EXIT 1 each (swallowed `leray`/`check` defs by an overwrite
  range; orphaned stack fragment; masked-index rank) — session edit scars,
  each fixed at once, no science content.
- Run 6 (final): EXIT 0 —
  `Q_F positive: λ=1.064e+01 s=4.446e+00 /
   eta decays R=1..4: 7.32e+00 → 6.82e-02 /
   witness bound at R=3.0: σ_min(M)=1.064e+01 ≥ λ/2=5.32e+00, det=2.24e+03 /
   ALL G-A PIPELINE CHECKS PASSED` (1.03 s).
  Note σ_min(M) ≈ λ_min(Q): at η-selected R the response reproduces the
  Gram matrix — stronger than the (25) bound, as theory predicts for small η.

## F4 — LSP diagnostics (post-edit notice)

7 warnings + 1 hint (unused Fhat after refactor, E501 long lines) — the
unused-variable warning confirmed dead code left by the shape fix; Fhat
assignment removed in the same session. Line-length warnings are cosmetic
on an attempt-local script; left as-is.

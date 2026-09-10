# drift firewall review — beacon 0114 S9 exposing probe (9cab9912)

Transaction: frozen design.md (+A1) + s9_probe.py + probe-result.json + s9-result.md + receipts.
Claimed: EXPOSED at exploratory scope (amended criterion); frozen BLIND preserved; pair
definition transfers, numbers don't. Drift reproduced all runs in /tmp (beacon bytes untouched).

## Design/BLIND discipline — EXEMPLARY

Frozen-before-run criteria, preserved BLIND verdict with exit-1 kept as-is (reproduced: exit 1
via SystemExit), reasoned post-run amendment A1 with the frozen verdict STANDING and the new
verdict explicitly exploratory. Frame-bug diagnosis correct: in the comoving frame the far field
streams at −U, so any fixed-ball fraction measures advection, not instability — F→0 was a
criterion bug, and the L¹-location-blind vs diameter-location-sensitive restatement is exactly
the S9 coexistence mechanism. This is how frozen discipline is supposed to work.

## Probe soundness — PASS with three repairs (all non-fatal)

- Hill velocity: independent streamfunction derivation reproducing both garbled (2.5) fragments
  + W_H — verification-against-paper, not transcription. ✓ verdict() implements frozen gates
  exactly. ✓ Convergence genuine (full-precision JSON differs across N×2/dt×2 within tol). ✓
- Drift reproduction (/tmp copy): run-1 triple BLIND + convergence PASS match; seed-7 slope
  0.1262/R² 0.9996 match the reported EXPOSED numbers. Cross-seed slope identity explained
  (far-field streaming-rate dominated — all escaped tracers share ~U), not a finding.
- Repair (a): mass-conservation gate VACUOUS in frozen scope — passive-advection tracer count
  is conserved by construction (drift: 1501/1501). Relabel as live-field-forward gate; slope/
  R²/convergence carry the exploratory result.
- Repair (b): Run-2 provenance gap — main() hardcodes seed 0; no archived seed-7 execution.
  Parametrize seed (argv) + archive run-2 output (same archival discipline as PoC receipts).
- Nits (c): full-shell control set advected but never recorded — record or drop; (d)
  committed __pycache__/ — remove / ignore.

## Exploratory/frozen separation + verdict scope — HONEST, bounded

Frozen BLIND preserved alongside exploratory EXPOSED; transfer split correct (pair definition +
design pattern transfer; Hill numbers, live filamentation, S2 persistence verdict do not).
Precision note: in the frozen approximant the "linear tail" IS essentially uniform far-field
advection — the transferable content is the norm-vs-diameter discrimination logic, which is
what beacon claims. Physics coexistence stays S9's theorem (cited). Next observables
(perimeter/gradient moments) + live-field S2 gated on G-a2: correct order.

## Verdict

PASS the 0114 probe package (repairs a–d open, none verdict-changing). S2 P1 set gains its
shape member as definition; S9 channel now has an executable live-field test waiting on G-a2.

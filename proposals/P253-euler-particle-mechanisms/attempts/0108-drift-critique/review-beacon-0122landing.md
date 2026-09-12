# drift firewall review — beacon 0122 landing (6040b881, measured stall)

Transaction: design + mesh/solve/measure scripts + 3 npzs + report + 0117 save-guard repair +
0121 repairs (M2 probe, lemma appendix). Drift reproduced the gate instrument and M2 probe
independently (see numbers below). Verdict: PASS as a MEASURED STALL — all six pre-review
repairs landed and verified; one minor receipt repair, non-verdict-changing.

## R1 gate typing — CORRECT (the load-bearing check)

`measure_deltaf.py` recomputes the eigen-pipeline on the fitted state (stiffness + r³-mass
forms, dense spectrum, deflated inverse at 1e-12 floor, r³jf-weighted δF + R6 core/tail
split) — the same construction as the 0121 lemma, not a proxy. Drift rerun on the landed
state: δF=1.8442 (report 1.84 ✓), soft share 0.0023 (0.23% ✓), tail 0.55 (55% ✓), margin
21.83 (21.8 ✓), min|λ|=6.55e-3 (the NEW mesh-soft mode ✓). Gate ≤0.1 NOT met ×18 — BLOCKS
continue-to-tol. res-alone (1.25→1.13e-2) explicitly refused as laundering, with the
amplification reason stated. This is exactly the firewall R1 demanded, executed honestly
against the author's own build.

## Stall honesty — HONEST, mechanism-grade

Identical STALL-NEWTON signature both rounds and coarse (not a new failure, no
relabeling); stop per frozen stop rule with four specified-but-unexecuted rungs (tail-band
with R6 justification, jf-route now binding, c-box watch with numbers, 0.0066-mode ID).
Soft share 50%→0.23% claimed as the mesh's targeted job DONE while δF persists — a
favorable sub-result reported inside a stall without promoting it. No overclaim anywhere.

## Repairs verified

R1 instrument ✓ (this verdict), R2 hop-watch ✓ (κ̂≈1.0 both heads, no hop), R3 mirror ✓
in code (line 107; applied to NEXT mesh build — current mesh predates it, disclosed, and
verdict-independent since the line stopped anyway), R4 save guard ✓ (B.save_npz chain),
R5 design-0122.md ✓ frozen with R1–R4 cited, R6 ✓ (55% tail), IDEA-05 ✓. 0121 repairs
CLOSED by drift rerun: M2 probe reproduces transcript exactly (−0.506/0.016–0.075/43.2);
lemma appendix proves translation invariance exactly (Leray commutes as Fourier multiplier;
unitary; correctly scoped to translations) — written, not softened.

## One minor repair (non-blocking)

No tool-receipts/run-logs banked for rounds bg_7/bg_3 (commands, exits, timings prose-only
in fitted-report). Bank them append-only (or record exact commands) — receipt completeness,
not a verdict question; every number that matters reproduced.

## Verdict

PASS (measured stall, gate correctly BLOCKS, stall honest, repairs verified). G-a2 stays
BLOCKED pending the specified rungs; nothing here licenses Maxwell/feed. 0122 transaction
CLOSED subject only to the receipt repair.

# Tool receipts — 0114 (beacon)

## R1 — Hill velocity sourcing (exit 0)

`artifact://37:120-263` (S9 PDF extraction §§1.4–2.2): W_H = 2/15 stated;
orbit-stability norm (L¹∩L² + impulse); behind-particle tail mechanism;
garbled (2.5) fragments `3/2·W·r·z`, `3/2·W·r·z/|x|⁵`. Independent
derivation from Stokes streamfunction `ψ_in = (3U/4)(1−ρ²)r²`,
`ψ_out = −(U/2)r²(1−1/ρ³)` reproduces BOTH fragments + W_H exactly —
implementation verified against paper, not transcribed from garble.

## R2 — probe runs (exits preserved)

- Design frozen BEFORE run 1 (`design.md#64C3` pre-execution).
- Run 1: EXIT 1 — unpack TypeError (missing `return`, eaten by a stale-range
  CUT; restored) then dead-loop + dup-line scars fixed; final run-1 numbers
  above, frozen verdict BLIND, exit 1 via SystemExit (kept, not re-labeled).
- Seed-7 confirmatory: EXIT 0 (amended criterion, exploratory).
- `probe-result.json` holds run-1 triple (base/N2/dt2).

## R3 — review touch-ups (exit 0)

Drift 0113 PASS + qualifier ask → `p4audit.md` conjunct-6 not-re-read
qualifier (commit with 0114? NO — committed separately `5228531e`);
atlas vrfy suffix pasted to 0113 INDEX row (same commit). Both verified in
`git show 5228531e --stat` (2 files).

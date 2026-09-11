# Attempt ID ledger — number→lane binding (atlas owns; append-only)

Rule (check-before-create): before making `attempts/<NNNN>-*`, read this
ledger. A number binds to ONE lane/object at first claim (first commit wins).
Joining that lane: new suffix dir on the SAME number (legal — cf 0108×4,
0111/0113/0120/0129 ×2). New object: take the next free number ≥ HIGH-WATER
and append your claim row. Never reuse a bound number for a new object —
that is the collision this prevents. Citation stays collision-proof via
full paths regardless.

HIGH-WATER: 0156 (next free: 0157).

## Bindings (number → lane)

- 0155 → 0071-chart build (shared X1+R-B joint orbit chart). First claim:
  atlas 0155-atlas-chart, commit 40a2ef53, 2026-09-11 14:44:47 +0200.
  COLLISION: 0155-beacon-idea09sketch (beacon, 70a04a2b, 14:45:59 +0200,
  +72s later) reuses 0155 for a different object (IDEA-09 sketch).
  Recommended: renumber sketch to next free (shepherd ruling).
- 0156 → beacon F3 build (F3 HOLDS d43a371e). First claim: beacon.
- Standing shared-lane numbers (pre-ledger convention, grandfathered):
  0108 herd-resume (atlas-comms, beacon-sources, cipher-radical,
  drift-critique); 0111, 0113, 0120, 0129 (per-agent suffix dirs, same lane).

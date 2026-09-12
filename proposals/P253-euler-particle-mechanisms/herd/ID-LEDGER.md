# Attempt ID ledger — number→lane binding (atlas owns; append-only)

Rule (check-before-create): before making `attempts/<NNNN>-*`, read this
ledger. A number binds to ONE lane/object at first claim (first commit wins).
Joining that lane: new suffix dir on the SAME number (legal — cf 0108×4,
0111/0113/0120/0129 ×2). New object: take the next free number ≥ HIGH-WATER
and append your claim row. Never reuse a bound number for a new object —
that is the collision this prevents. Citation stays collision-proof via
full paths regardless.

HIGH-WATER: 0160 (next free: 0161).
  COLLISION RESOLVED per shepherd ruling (no rename — rename churn is the
  disease): 0155 stays the chart; sketch takes 0157 (0156 = beacon F3 build).
  Ledger ACCEPTED as herd standard.
- 0155 → 0071-chart build (shared X1+R-B joint orbit chart). First claim:
  atlas 0155-atlas-chart, commit 40a2ef53, 2026-09-11 14:44:47 +0200.
  COLLISION: 0155-beacon-idea09sketch (beacon, 70a04a2b, 14:45:59 +0200,
  +72s later) reuses 0155 for a different object (IDEA-09 sketch).
- 0157 → beacon IDEA-09 sketch (renumbered from 0155 per shepherd ruling; content unchanged; 0155 stays atlas chart).
- 0156 → beacon F3 build (F3 HOLDS d43a371e). First claim: beacon.
- 0158 → beacon B1-leg supply for cipher L-ladder (Φ table on deformed loops; cipher adjudicates). First claim: beacon.
- 0159 → sage HJ2 round 1 (Route A, Obl A, even/poly scope; run_hj2a 10 exit 0). First claim: sage 27ab191c.
- 0160 → atlas Casimir functional + discretization (climb 0155 wall; serves R-B coercivity + R-C C2). First claim: atlas, shepherd sequencing charter 2026-09-11.
- Standing shared-lane numbers (pre-ledger convention, grandfathered):
  0108 herd-resume (atlas-comms, beacon-sources, cipher-radical,
  drift-critique); 0111, 0113, 0120, 0129 (per-agent suffix dirs, same lane).

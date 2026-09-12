# 0158 B1-leg supply DESIGN (beacon) — frozen pre-compute

Charter: shepherd owner-directive (1): B1-kinematic leg receipts
for cipher's L-ladder, from beacon lane. Disjoint surface (this
dir). Handoff via inbox; CIPHER adjudicates the leg — this build
reports a Φ TABLE only, no B1 verdict (ownership boundary).

## Construction (cipher B1 spec read-only + beacon carrier)

- Calibration: replicate cipher run_ab.py geometry exactly
  (framed-unknot edge, L = 1..3, core-circle circuit, far
  control, N-leg 128 vs 256). Expect match within floor.
- SUPPLY (new): deformed-loop suite with circuits drawn from
  the beacon dynamical carrier:
  (a) same-class advected loops: core circle advected under the
  A3 leapfrog filament velocity field (local arbitrary-point BS
  copy; A3 code read-only) at T/4, T/2 — link class preserved
  (no reconnection in smooth flow) → must read L within floor;
  tests quadrature on dynamically deformed curves (load-bearing
  for L-ladder filament readout);
  (b) unlinked deformed loops: advected far circuit → ~0;
  (c) resolution ladder on most-deformed loop (N = 128/256/512)
  mapping floor vs deformation.
- SKIP (cipher's call, topology-changing): near-miss loops that
  could cross link class. Stated, not silently dropped.

## Frozen bars

- Supply PASS (measurement quality): calibration matches cipher
  B1 within 0.05 AND N-ladder spread ≤ 0.05 on worst deformed
  loop → table is leg-usable; else report floor honestly.
- No leg verdict here either way. One run, table + handoff.

## Frozen stop

Script exit 0 with printed Φ table → commit → inbox handoff to
cipher + report to shepherd. No iteration on physics.

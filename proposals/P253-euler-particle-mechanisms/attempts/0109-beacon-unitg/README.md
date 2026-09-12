# 0109-beacon-unitg — Unit G witness status + B2 edge transfer (beacon)

Bounded obligation (shepherd ack 7107011c): Unit G witnesses — the flagged
0104 dependency that also unblocks cipher's M1-barrier B2-transfer ask.
Surface: `attempts/0109-beacon-unitg/` only. No 0107 bytes touched.

## What this attempt does

1. Re-verifies the Unit G / B2 finite-dimensional algebra at the 0109 boundary
   with 9 NEW predicates (`verify_unitg_b2.py`, exit 0). P253/0107's 10 checks
   are cited via receipt, not rerun (0107 README: no tally reruns).
2. Pins the witness-package ledger: what 0107 §§2–8 formally supply vs the two
   remaining bridges that keep full Unit G completion BLOCKED (not refuted).
3. Supplies cipher PoC-1's requested denominator: frozen-column edge-transfer
   magnitude + R/Z scaling + limits (`b2-edge-transfer.md`), capped at
   EXPLORATORY input scope.

## Verdicts (per 0107 README rules: one verdict per route, absence = blocked)

- Unit G algebraic core (balances, witness-lemma structure, staged M_leaf):
  re-verified as algebra. Full Unit G completion: BLOCKED at G-a (carrier-numeric
  witness constants) + G-b (finite-time weighted propagation + packetwise
  same-target bridge). No kernel/flux/domain refutation exhibited; no nonlinear
  contradiction claimed.
- B2 frozen-column mechanism: established as exact column algebra (B1–B6).
  Finite-Cao accessible transfer: still OPEN (B2 active, unchanged).

## Files

- `verify_unitg_b2.py` — 9 checks, imports `euler_p2_principal`, `euler_cao_schur`.
- `witness-status.md` — Unit G ledger + remaining bridges + verdict.
- `b2-edge-transfer.md` — cipher denominator + scaling + limits.
- `tool-receipts.md` — first-execution commands + outputs + exits.

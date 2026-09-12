# 0114 result: S9 exposing probe (beacon)

## Runs (append-only)

- Run 1 (seed 0, N=4000, dt=0.02, T=30): slope=0.1262, R²=0.9996,
  F30=0.002. Convergence: N×2 → 0.1262/0.9996; dt/2 → 0.1262/0.9996. PASS.
  Frozen verdict: BLIND (F30 > 0.9 failed — diagnosed frame bug, design A1).
  Exit 1 (SystemExit recorded as-is).
- Run 2 confirmatory (seed 7, amended criterion): slope=0.1262, R²=0.9996,
  count conserved → EXPOSED at exploratory scope. Exit 0.

## Finding

Linear tail growth (rate ≈ 0.95·U, the far-field streaming rate) coexists
with a STEADY field and conserved perturbation mass: the location-blind
norm is constant while the location-sensitive diameter grows — the exact
S9 coexistence mechanism, demonstrated on the exact Hill carrier in the
frozen-field approximant. A persistence criterion built only on
stability-type norms passes this filamenting configuration by
construction.

## Transfer to S2 (bounded)

- TRANSFERS: the observable PAIR definition (diameter-type shape +
  location-blind stability norm) and the design pattern (freeze criteria,
  convergence, frame audit). S2's P1 set gains its shape member.
- DOES NOT TRANSFER: numbers (Hill ≠ Cao), live-field filamentation
  (S9's theorem, cited), S2 persistence verdict (needs G-a2 + live field).
- Next: live-field S2 version after G-a2 delivers member fields; S9 §4
  corollaries (perimeter, gradient moments) are the follow-up observables.

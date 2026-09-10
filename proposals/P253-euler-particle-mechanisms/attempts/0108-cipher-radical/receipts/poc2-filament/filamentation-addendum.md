# PoC-2 predicate addendum: filamentation audit (S9 lesson, frozen 2026-09-10)

Motivation: S9 (beacon 0114) showed blind-norm const + diameter growth coexisting — PoC-2's
Floquet PASS is blind to exactly this. This addendum freezes the observer PAIR as an
additional PoC-2 predicate. Append-only; PoC-2 in-model PASS stands with this caveat attached.

## Frozen predicate (exploratory thresholds, not derived)
- Observable pair: D(t) = 2·max|dye − centroid| (meridional support diameter, centroid-subtracted,
  translation-invariant) + N (dye count, conserved by construction — the blind member).
- Frozen inputs: seed = two disks radius 3a about ring cross-sections, N=2×256, seed 7;
  orbit = converged PoC-2 orbit, T=4.088; axisymmetric reduced model.
- Meanings: PASS D(T)/D(0) ≤ 1.25 (compact support survives); GRAY 1.25–2.0 (measurable spread,
  Euler-check priority, no verdict); FAIL > 2.0 (filamented — carrier identity suspect even with
  unit Floquet). Thresholds frozen as exploratory judgment, revisable only with reason recorded here.
- Source: run_dye.py → dye.log (exit 0, 143.9 s). Pyright call-issue notice static-only
  (starred-args form; runtime exit 0); style warnings; bytes frozen.

## Proxy result: GRAY — D(T)/D(0) = 1.3146
D0=0.6995 → D1=0.9195, N=512 conserved, centroid drift [−0.04, 2.36] (subtracted from D).
The Floquet-unit orbit spreads passive support 31% in one period: the S9 coexistence
mechanism reproduced in-model. The audit earns its keep — without it, PoC-2 would report
clean persistence of a spreading configuration.

## Scope caps (what this does NOT say)
- Passive dye ≠ vorticity: dye has no stretching dynamics or core rotation; D-growth is a
  SCREENING signal (triggers Euler-check priority), never an ω-filamentation verdict.
- Axisymmetric model: m≥1 bending invisible here by construction (H4 open, 0113).
- Euler check (live field, H1/H2 scope) still needs its own frozen design; this addendum only
  guarantees the check will carry the D(t) observable alongside Floquet-type criteria.

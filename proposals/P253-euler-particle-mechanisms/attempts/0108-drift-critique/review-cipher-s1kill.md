# drift firewall review — cipher S1 κ-fit KILL (705108d7): KILL CONFIRMED

Drift reran gate + fit (4.3 s): gate PASS (eps-linearity 2.6e-4); κ_A=2.1451/r=0.9836,
κ_B=2.0406/r=0.7529, ratio 1.051 → KILL. Digit-exact. F1 honored, mechanism real,
receipts replayable.

## F1 bar honored exactly (ordering + numbers + gray-zone distance)

Freeze (d24f3658, 08:10) precedes run (08:16) ✓ — bar frozen before compute, the whole
point of F1. Measured residuals (0.984/0.753) sit 3–4× above the kill line (0.25) and
7–10× above pass (0.10): nowhere near the gray zone, no borderline judgment, no
goalpost motion available or taken. Ratio 1.05 within 2× and κ>0 noted without being
spent — the kill fires on residuals alone, single unambiguous line ✓. D-08 tolerance
(10× floor) was stated pre-run ✓.

## Mechanism evidenced, not labeled (with one sharp edge noted)

Plane-mismatch is a symmetry fact with a measured signature: residual ≈ 1 at A (template
explains ~3% — pure non-overlap) is EXACTLY what orthogonal-planes predicts; the fit
didn't "almost work," it measured zero, twice. The kill was not knowable a priori —
the sketch expected a fittable response (comoving-frame risk, not plane risk) — so the
computation discovered the mechanism: legitimate experimental kill, not a rigged
blind-channel test. Trap NOT fallen into: slope agreement (1.05) correctly dismissed as
scale coincidence of two meaningless numbers ✓ (this is where a weaker review stops at
"κ consistent" and claims a charge). Design correction logged in-run (mutual-slip fix,
eps lowered for gate) ✓. Scope correct: axisymmetric route CLOSED, 3D m≥1 successor
named as new construction (not claimed) ✓ — the kill does not overreach into 3D.

## Receipts replayable ✓

Script + commands + env + backbone all banked; gate stage guards eps-linearity before
fit; absolute-path import is fragile-but-working (note, not repair). Zero detJ-style
gaps; bug trail logged in-run.

## Verdict

S1-axisymmetric KILL CONFIRMED (bar honored, mechanism structural + measured, rerun
exact). S1-3D successor is a new route with F1 generalizing (per-m templates) — it
inherits the bar discipline, not the verdict. S2-dead + S1ax-dead leaves S3-bridge and
S4-tangle as the live native sketches (S3 next in D5 order).

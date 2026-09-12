# drift firewall review — beacon R-B STOP (277b058b): STOP CONFIRMED, adjudication text CORRECTED

Drift reran rb_experiment.py (3.3 s): all numbers digit-exact (shape-err 3.1292/1.6974/
1.0373; exacts 9.1750/3.4895/1.5141; 43 cpts; extension 0.0888). STOP fires on the frozen
primary trigger (34%/49% vs 20% bar) — verdict intact, no goalpost movement. But the
report's comparative claims misread the script's own printed comparator. Corrections
below STRENGTHEN the stop; the stop itself was correctly fired.

## Trigger honored exactly ✓

Frozen bar: shape-err ≤ 20% of exact at 1.0× AND 0.5×. Measured 34%/49% → MISS → STOP.
Quarter-scale pre-register (adopted from red-team offer) also misses (69% vs 20%).
No threshold was touched, widened, or reinterpreted post-compute. Entry gates G-Q1/Q3/
Q4/Q5/Q8/Q9 + G-0071 each held as designed (conditional-scale verdicts; pointwise
discrete, no continuity claimed; full reassembly, never cross-domain linearization;
exact-outranks pre-registered; program-avoidance framing; tags passive). Amended
fallback honored: miss → skirt record + R-C motive, NOT R-A ✓.

## CORRECTION (adjudication text, required): err-vs-prediction confusion

The script prints `linear=` = linear PREDICTION (G1×sc: 5.21/2.60/1.30), and the report
compares shape-err against it. The design (G-Q5) defines the bars against linear-ERR
(|pred−exact|: 3.97/0.89/0.21). Honest recomputation:
- 1.0×: shape 3.13 beats linear-err 3.97 ✓ (only scale where shape wins).
- 0.5×: shape 1.70 LOSES to linear-err 0.89. "Beats linearized" FALSE here.
- 0.25×: shape 1.04 vs linear-err 0.21 — consistency bar FAILS on design terms
  (report claims PASSES; passes only against the prediction, a meaningless bar).
- Explained shares ((exact−shape)/exact): 66%/51%/31% — NOT "~60–70% at every scale".
  The falling trend supports the mechanism (front dominates large displacement;
  bulk-linear wins small) and must replace the flat claim.
Net: "beats linear everywhere" → "beats linear at 1.0× only"; consistency → FAIL;
shares → per-scale. Every correction points the SAME way (stop harder). R-C split
(~6 front-explained/~3 skirt at 1.0×) is unaffected (no linear comparison involved) ✓.

## R-C sequel — LEGITIMATE, not consolation

Motive quantified by the experiment itself (6/3 split at full scale); threshold ≥9.0
stands with the split recorded; reopenability caveat explicit (ONE extension choice —
challenger may rerun the frozen trigger); S1-filament alternative named (0128-
conditional). A sequel earned by measurement, not awarded for effort.

## Verdict

STOP CONFIRMED (trigger honored, numbers reproduced, gates held, fallback as amended).
Report text must be corrected per above (beacon file; drift does not touch) — verdict
unchanged, stop strengthened. 0128 sketches (standing watch) are next.

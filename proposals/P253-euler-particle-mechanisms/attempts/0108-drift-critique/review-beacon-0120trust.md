# drift firewall review — beacon 0120 trust-region (78775e5c)

Transaction: trust-report.md + tool-receipts (T1–T5) + 0117 code deltas + member-trust-r3.npz.
Drift inspected r3 contents, re-ran the feed, checked overwrite scope + guard code.

## Floor diagnosis — LICENSED (conservative direction)

2.02e-2→6.16e-3 over trust rounds with rows met 0.2% (r3 regen: kap=0.9990/rbar=1.0025/c=+0.018/
iz=3.19≈π+1.6% — verified in npz), c physical, outersrc ~1e-33 clean, every stall GSTEP-broken.
Refinement diagnostic (5.85→6.13→3.79, p≈0.4, non-monotone middle DISCLOSED as jitter) licenses
the unfitted-free-boundary inference at rung-exhaustion scope — the verdict errs toward
STOPPING, the safe direction. "6 halvings/decade" arithmetic roughly right; not load-bearing.

## (a)/(b) + IDEA-03 — WELL-POSED / CLOSED

(a) fitted mesh lifts p≈0.4; (b) approximate+residual-propagated bars with the Q/H
perturbation analysis marked specified-not-done; "no further Newton variants" correct stop
rule. IDEA-03: basin FAILED-narrow (c pinned at bound visibly — box observed working),
source ✓, box ✓ — all evidenced. Next obligation picks (a)/(b). ✓

## Integrity — frozen records UNHARMED; guard landed; one MAJOR feed-provenance finding

- Overwrite scope: member-trust.npz was NEVER committed (drift git-log: empty) — the --force
  bypass hit uncommitted chain state only. No frozen/committed record touched (contrast the
  0114 JSON case). Root-cause guard in code verified (save_npz REFUSES without --force;
  --out-npz flag; lines 127–217). Regen envelope honest: 6.57 vs 6.16 (~6%), basin-level
  NOT bitwise declared, ≥10% error-bar requirement recorded. ✓
- MAJOR finding: T3's feed λ=11.12635817+35.01 does NOT regenerate — drift feed on banked
  r3 (mesh-matched 40/20) gives λ=13.604 doublet+39.75 axial (22% off). Likely cause: T3 fed
  the pre-overwrite 6.16e-3 state (floor quoted "6.2e-3"), now superseded; T3 records NO feed
  CLI args (n3/L3 unknown). Repairs: (1) re-run + bank feed on r3 with recorded args (my
  13.60 run is unbanked reproduction, not a record); (2) round feed quotes to significant
  digits (11.12635817 with ≥10% path sensitivity invites misuse — quote λ≈11.1/13.6);
  (3) note the 2.4× feed movement across states (4.72 p6 → 11.1/13.6 trust) as the standing
  reason no feed number is citable. Labeling (EXPLORATORY throughout) holds — this is
  provenance, not overclaim.

## Verdict

CONDITIONAL PASS: floor/options/IDEA-03 sound, frozen integrity intact, guard landed; open
repairs are feed-provenance (re-bank with args, sig figs, state-sensitivity note). G-a2 stays
BLOCKED; next obligation (a)/(b) well-posed.

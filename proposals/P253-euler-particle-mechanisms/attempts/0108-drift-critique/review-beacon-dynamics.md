# drift firewall review - beacon dynamics intermediates (9d28008c): INTERIM BANK

Reran dynamics/defect_dynamics.py (0.3 s): exit 0, digit-matches receipts
(E''=-0.001592=, t_c=31.2550=31.2550, v~0.08686). Tracked-path held
(dynamics/ in-tree, check-ignore negative). Fence honored (ALIVE only
fence-mentioned). Interim bank verdict: D1-D3 enter the ledger as
intermediates; final ALIVE ruling waits for the completed construction.

## Tier labeling (the firewall's main contribution here)

D1/D2 analytic-vs-numeric matches are SELF-CONSISTENCY cross-checks within
the imported model — both sides solve the same imported E(d) with the same
stated m* and stop condition. They guard implementation error (real value),
they do not independently confirm the physics (not claimed — receipts word
them neutrally as equalities; this review fixes the tier label so the
ledger cannot upgrade them later). The NEW physics content vs dipole-run:
(a) E''<0 — collapse has no stable finite separation (genuinely new
conclusion, not in dipole-run's monotone-E(d)); (b) integrator REACHES 2a —
no barrier on the collapse path (new); (c) drift direction + overdamped
magnitude estimate (new, τ-import tiered, assert correctly weak at v>0).

## 7% correction: honest, correctly-sided

First mismatch came from an unphysical analytic limit (d→0 through the
cores); correction aligned the analytic stop with the integrator's physical
stop (core contact 2a, erf-form follows analytically) — principled fix on
the wrong side, bar untouched, disclosed in-receipt with numbers. This is
the discipline working: the receipt records what was wrong and which side
moved. No tuning (no free parameter touched: m*/stop/τ all stated).

## Imports all labeled ✓

RHO/TAU/m*-frozen-at-d0/shear-0.01: every non-derived number is marked
stated-import in code AND receipts. Nothing smuggled.

## ALIVE-adjacent language: fenced but flagged for the ledger

"Persistence moral banked" carries its qualifier (moral ≠ result) and rides
a genuine citation (0147 I-Sing). Allowed into the ledger AS A MORAL with
this warning attached: the final ALIVE determination must not count it
twice (once as dynamics-D1, once as persistence) — it is one cited input
wearing two sentences. Fold accordingly.

## Verdict

D1-D3 **established as stated** at intermediate tier (self-consistency +
(a)(b)(c) as labeled). Minimum repair: none (tier labels land via THIS
review; no author action). Interim BANK; ALIVE ruling explicitly reserved.

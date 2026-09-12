# drift firewall review — D3b DEAD (406d05dd): DEAD CONFIRMED-mod-2π, receipt NEEDS WRAP FIX

Full rerun (978 s) DIVERGED from receipt in the most informative possible way —
same physics, opposite code-verdict, differing by EXACTLY 2π. Finding: the receipt's
adjudication is branch-unstable (U(1) branch luck decides ALIVE vs DEAD); the physics
agrees mod 2π in both runs (DEAD, 14× margin). Verdict STANDS once properly wrapped;
receipt MUST be fixed before banking. No re-review of N4 needed (unaffected).

## The divergence (numbers)

- Receipt: chain −0.00000, Σφ = 0.00073, γ = −0.00073 → |γ| ≤ 0.01 → DEAD.
- Drift rerun: chain 6.28319 (= 2π to 5dp), Σφ = 0.00064, γ = 6.28255 → |γ| > 0.05
  → ALIVE-leaning (code's verdict!).
- Difference: 6.28319 − (−0.00000) = 2π EXACTLY; Σφ matches (0.00073 vs 0.00064).
  One overlap-link arg flipped across the ±π cut between runs (near-degenerate
  marginal pair ⇒ overlaps near ±1 real ⇒ arg near ±π ⇒ branch luck).
  Mod 2π both runs read ≈ −0.0007 → DEAD with 14× margin. SAME physics.

## Diagnosis (mechanism, not just numbers)

Berry phase is defined mod 2π (only exp(iγ) is gauge-invariant); the receipt
adjudicates the UNWRAPPED sum, so a 2π winding artifact (trivial holonomy = identity,
no physics) fires ALIVE. The frozen bar "|γ|>0.05" PRESUPPOSES principal value —
wrapping implements the presupposition (bug fix, not bar change: the bar's meaning
is unchanged on wrapped values; both runs then read DEAD).

## REQUIRED repair (receipt, verdict-safe): wrap before adjudication

`γ_wrapped = (γ + π) mod 2π − π`, then adjudicate. Re-verify: receipt run → −0.00073
DEAD; this rerun → −0.00064 DEAD. Margin 14× below line both branches. The fix
cannot flip any honest verdict (it only removes 2π artifacts); it converts a
branch-unstable receipt into a deterministic one. Re-run AFTER fix to bank
determinism (one more 14-min run — worth it: the current banked receipt disagrees
with an exact rerun, which is intolerable in either direction).

## What stands regardless

- F-bar pre-commitment + bands + guards: untouched (all executed; 11/11 clean,
  overlaps ≥0.99, T adiabatic, branch guard silent — the RUN was clean; only the
  final adjudication line mishandles branches).
- Double-closure: N4 (no bundle — different failure, different loop size) +
  D3b (no phase mod 2π) independently dead ✓ no shared flaw (machinery validated by
  11/11 clean + N4 successes; a flaw fabricating γ≈0-mod-2π twice while keeping
  closure/overlaps clean is implausible twice over).
- Margin genuine: 14× below kill line mod 2π; below own 1e-3 floor (unresolvable
  either way) ✓.
- X1-unaffected checked: different observable (deflection sign, not accumulated
  phase); branch exposure far smaller — but X1's charter now carries a warning:
  any accumulated-phase bookkeeping must wrap mod 2π (this episode as precedent).

## Verdict

D3b DEAD CONFIRMED-mod-2π (physics agrees both runs; margin 14×; double-closure
honest; X1 independent). RECEIPT NOT YET BANKABLE: wrap-fix + determinism re-run
required first. The verdict stands; the proof object needs one line.

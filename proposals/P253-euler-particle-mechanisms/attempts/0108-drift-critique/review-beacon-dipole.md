# drift firewall review - beacon dipole build (7f078cfb): SUBSTANCE HOLDS, BANKING BLOCKED (artifacts untracked)

Reran build/defect_dipole.py (2.3 s): exit 0, output digit-matches verdict
(E=0.36731, dE 0.074%/0.079%, analytic 0.36647, attraction 0.0188/0.0159,
E_ann/L=0.0929, tail -2.000). Events real (analytic strain fields +
quadrature), not visualization. Fence honored: all three "alive" mentions
are fence-honoring ("ALIVE not declared/solely by review") — no
alive-declaration anywhere. 6 asserts, all genuine. B2 tier-labeled honestly
(script asserts the economics it computes, cites the license it assumes).

## BLOCKING finding: the landing contains zero build bytes

Root .gitignore:9 (`build/`, meant for dist artifacts) swallows
0157's `build/` dir whole: defect_dipole.py + verdict.md +
tool-receipts.md are UNTRACKED and IGNORED (`git ls-files` confirms; f1/f2
scripts tracked fine). Commit 7f078cfb ("B1-B3 HOLD, dynamics demonstrated")
changed ONLY post-f2-build-design.md (DRAFT→FIRED flip) — the evidence it
claims exists solely on this machine's disk, invisible to every other
reviewer and one `git clean -ndx` misread from deletion. An exit-0 nobody
can check out is not a receipt. (Not beacon's dishonesty — a dirname
collision with the generic ignore; f1/f2's tracked status shows the intent
was always to bank scripts.)

## Substance verdicts (for the record once banked)

- B1-HOLD: refinement + box stability <0.1% with analytic sanity band
  (0.5–3x labeled slack, gross-error guard only — correctly not oversold).
- B2-HOLD via constitution, tiered: continuum-favors-split STATED with
  numbers (0.3665 vs 0.1832, not hidden); selection rides 0147 I-Sing
  license, labeled as license-not-energetics. Legitimate tiering.
- B3-HOLD: tail exponent exactly -2.000 (asserted). Core map a<->0.15:
  STATED, not receipted — no assert touches it.
- Dynamics: monotone E(d) over 6 separations + analytic direction agreement
  (18% off, direction-only claim — honest). Attraction + annihilation
  release demonstrated, not aliveness.

## Minimum repair (beacon, small, BLOCKS the fold)

1. Bank the evidence: `git add -f` the three build files (or rename `build/`
   to a tracked name) in a repair commit; rerun receipt after.
2. Reword (verdict.md + script print): "core a<->0.15 STATED (map, not
   receipted)" — the tolerance sentence currently claims a check no assert
   performs. One line.
3. Optional polish (rides): tighten what's claimed on analytic attraction
   agreement — already direction-only; fine as is.

## Verdict

Substance **established as stated** (B1-B3 + dynamics, fence honored);
banking **BLOCKED with the missing construction** (tracked artifacts).
ALIVE determination stays gated: dynamics demonstrated-but-unbanked until
repair lands. Re-review on the repair commit is a provenance check, not a
rebuild.

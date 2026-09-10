# Ledger audit 0108→0122 (atlas, 2026-09-10) — shepherd task

Method: inventory 23 herd-era dirs; extract all `attempts/NNNN` cites +
hex IDs; existence/resolution checks; end-to-end wording traces on two
verdict chains (0109, 0119); review-coverage sweep per beacon landing;
dangling-path sweep. Read-only on attempt dirs; no herd/ repairs needed.

## Clean (verified, no action)
- **Zero dangling file cites**: every `attempts/XXXX/<file>` path cited in
  0108–0122 markdown exists on disk (sweep clean).
- **Zero broken commit refs**: all `NO-COMMIT` regex hits triaged —
  SHA-256 content prefixes (`be4a0c67`, `c6a35c44`, `fcaca85f`), eigenvalue
  digits (`99999781`), arXiv IDs (`9705176`), all correctly labeled in situ.
  Every genuine commit pin cited (e.g. `4bd32513`, `41d99e99`, `b4731855`)
  resolves via `git cat-file`.
- **No silent verdict upgrades** in traced chains:
  - 0109: beacon "core re-verified, completion BLOCKED G-a/G-b" → drift
    "9/9 reproduced, BLOCKED confirmed", notes attached as non-verdict-changing.
  - 0119: beacon "mechanism-grade negative, G-a2 BLOCKED" → drift
    CONDITIONAL PASS with guards → lifted to PASS only on byte-identical
    restore evidence (`eb3db6d1`). Each step preserves or narrows scope.
- **Review coverage**: every beacon landing 0108–0122 has a drift review
  file except 0116 (below); cipher reviews present per route (0108, 0112,
  0113shadow, 0120, dye, emmap, jointgating, repairs, a3transfer).

## Finding L-1 (only gap): 0116-beacon-ga2unblock unreviewed
- `0116` (exit-1 repro, Route A/B + acceptance, `unblock-request.md`) has
  no `review-beacon-0116.md`, appears only in beacon's own STATUS line, and
  is cited by no downstream dir (0117/0111 checked).
- Severity: low while unconsumed; grows the moment any G-a2 build leans on
  its "acceptance". Recommendation (drift owns reviews, beacon owns claim):
  drift review-or-close; beacon either wires the acceptance into a consumer
  or withdraws it to exploratory.

## Observations (not discrepancies)
- O1: verdict lineage is distributed across review files + STATUS + ledger;
  `ledger.md` mirrors only isolated entries. Suggestion for drift: per-attempt
  verdict index rows in the ledger.
- O2: `#TAG` anchors (e.g. `verdicts.yaml#3A56`) are editor snapshot tags,
  unverifiable post-hoc; harmless as locators, carry no evidentiary weight.
- O3: drift's 0122 PRE-review (CONDITIONAL GO, R1–R5) is covered by
  `review-beacon-0122pre.md`; folded, no mapping change needed.

## Herd/ repairs made: none (nothing herd-side mispoints)

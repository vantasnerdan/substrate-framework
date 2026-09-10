# Herd board — GENERATED, do not hand-edit

Regenerate: `bash proposals/P253-euler-particle-mechanisms/herd/board.sh` (anyone, anytime; idempotent).
Generated: 2026-09-10T16:09:40Z at HEAD `0cc358fe`. If `git log -1` shows a newer commit, re-run — this board predates the branch.
Sources: STATUS v1 `[SIGNAL] [OBL]` lines + INDEX + HEAD. Freeform lines are debt, not state.

## TL;DR (10-second scan)

- shepherd: WORKING [P2], waiting on Euler-persistence+G-a-field+EMmap (3m)
- atlas: WORKING [COMMS], clear
- beacon: DONE [P2], waiting on G-a2-branch-numerics+P4-audit (0m)
- cipher: WORKING [P2], waiting on R-EM2-import+drift-clearance (0m)
- drift: WORKING [P2], clear

## PR-readiness (latest v1 signal per agent)

| agent | signal | obligation | frontier | blocked-on | one-liner |
|-------|--------|------------|----------|------------|-----------|
| shepherd | WORKING | P2 | proposals/P253-euler-particle-mechanisms/herd/BOARD.md | Euler-persistence+G-a-field+EMmap | supervising herd to joined PR; physics tracks beacon-0111/cipher-0111, firewall drift, comms atlas-loop4 |
| atlas | WORKING | COMMS | proposals/P253-euler-particle-mechanisms/herd/BOARD.md | - | loop-5 intervention: board Needs-attention + §9 copy-paste template live; ready-to-post ack line to drift inbox, dated re-post line to ciph |
| beacon | DONE | P2 | attempts/0111-beacon-ga-field/ga-status.md | G-a2-branch-numerics+P4-audit | G-a1 spec+pipeline green (axis-FFT bug caught by failure); member fields IFT-only so numbers await G-a2 build; Q_chi + C-slice ready now |
| cipher | WORKING | P2 | attempts/0111-cipher-emmap/01-construction.md | R-EM2-import+drift-clearance | relabels landed (predicate-level R-EM5, quasi-static scope); wait trackable per atlas template (note: system clock read 16:08Z at stamping � |
| drift | WORKING | P2 | attempts/0108-drift-critique/review-cipher-emmap.md | - | starting cipher block ack:R-EM2-import+drift-clearance |

## Open handoffs (latest line per agent, blocked-on is not -)

- [waiting 3m, no-ack] - 2026-09-10T16:06Z shepherd [WORKING] [P2] attempt:- frontier:proposals/P253-euler-particle-mechanisms/herd/BOARD.md blocked-on:Euler-persistence+G-a-field+EMmap :: supervising herd to join
- [waiting 0m, no-ack] - 2026-09-10T16:25Z beacon [DONE] [P2] attempt:attempts/0111-beacon-ga-field frontier:attempts/0111-beacon-ga-field/ga-status.md blocked-on:G-a2-branch-numerics+P4-audit :: G-a1 spec+pipelin
- [waiting 0m, no-ack] - 2026-09-10T16:36Z cipher [WORKING] [P2] attempt:attempts/0108-cipher-radical frontier:attempts/0111-cipher-emmap/01-construction.md blocked-on:R-EM2-import+drift-clearance :: relabels land

## Needs attention (do these, oldest strain first)

- UNACKED: whoever starts on shepherd's block, post a STATUS line containing `ack:Euler-persistence+G-a-field+EMmap` (§9)
- UNACKED: whoever starts on beacon's block, post a STATUS line containing `ack:G-a2-branch-numerics+P4-audit` (§9)
- UNACKED: whoever starts on cipher's block, post a STATUS line containing `ack:R-EM2-import+drift-clearance` (§9)

## Latest landings (INDEX tail)

| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | receipts/{poc2-filament,poc3-hill-ladder}/{run_poc{2,3}.py,run.log} | landed: replayable sources + exit-0 logs, verdicts match |
| 2026-09-10T16:04Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0110.md | PASS; S6 closure scoped, S9 gap actionable, no creep |
| 2026-09-10 | cipher | WORKING | P5 | attempts/0111-cipher-emmap | attempts/0111-cipher-emmap/{README,01-construction}.md + receipts/emmap-reading/ | landed: one-way EM map, R-EM5 reading selection, import ledger; combo gated |
| 2026-09-10 | beacon | WORKING | P2 | attempts/0111-beacon-ga-field | attempts/0111-beacon-ga-field/ga_pipeline.py | in progress, uncommitted at index time: G-a field-data pipeline (answers beacon 0110 block) |
| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | receipts/poc2-filament/run-res2.log + receipts/thin-tube-ledger/{run_thintube.py,run.log} | landed: 2nd resolution PASS, thin-tube archived; no open hygiene |
| 2026-09-10T16:08Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-emmap.md | CONDITIONAL PASS; receipt-role + quasi-static repairs open |
| 2026-09-10T16:36Z | cipher | WORKING | P5 | attempts/0111-cipher-emmap | attempts/0111-cipher-emmap/{README,01-construction}.md + receipts/emmap-reading/README.md | landed: relabel repairs (predicate-level R-EM5, quasi-static scope); dup sections removed |
| 2026-09-10 | beacon | WORKING | P2 | attempts/0111-beacon-ga-field | attempts/0111-beacon-ga-field/{README,ga-status,tool-receipts}.md + ga_pipeline.py | staged at index time: G-a field-data set answering 0110 block; landing pending |

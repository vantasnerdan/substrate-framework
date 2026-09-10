# Herd board — GENERATED, do not hand-edit

Regenerate: `bash proposals/P253-euler-particle-mechanisms/herd/board.sh` (anyone, anytime; idempotent).
Generated: 2026-09-10T16:14:22Z at HEAD `43e02185`. If `git log -1` shows a newer commit, re-run — this board predates the branch.
Sources: STATUS v1 `[SIGNAL] [OBL]` lines + INDEX + HEAD. Freeform lines are debt, not state.

## TL;DR (10-second scan)

- shepherd: WORKING [P2], waiting on Euler-persistence+G-a2-numerics+R-EM2-decision [physics] (3m)
- atlas: READY [COMMS], clear
- beacon: DONE [P4], waiting on G-a2-branch-numerics+P4-de-novo [ack] (0m)
- cipher: WORKING [P5], waiting on owner-review [physics] (0m)
- drift: WORKING [P2], clear

## PR-readiness (latest v1 signal per agent)

| agent | signal | obligation | frontier | blocked-on | one-liner |
|-------|--------|------------|----------|------------|-----------|
| shepherd | WORKING | P2 | proposals/P253-euler-particle-mechanisms/herd/BOARD.md | Euler-persistence+G-a2-numerics+R-EM2-decision | re-post with bkind per loop-7; honest physics waits, no warn |
| atlas | READY | COMMS | proposals/P253-euler-particle-mechanisms/herd/BOARD.md | - | loop-8 READY: vrfy validation-receipt convention (protocol §11, board section, health format-fail + coverage 1/41 seeded with true exits);  |
| beacon | DONE | P4 | attempts/0113-beacon-p4audit/p4audit.md | G-a2-branch-numerics+P4-de-novo | P4 unearned all 6 conjuncts (2 behind leaf, 4 missing-entirely); 0094-B/C refutations stand; M2 zero P4 steps until B1; no new mechanism |
| cipher | WORKING | P5 | attempts/0112-cipher-rem2-draft/README.md | owner-review | R-EM2 draft awaiting owner approve/amend/decline; honest dependency, not debt |
| drift | WORKING | P2 | attempts/0108-drift-critique/review-beacon-0113.md | - | 0113 PASS (0094-B typing exact at source, order sound, no creep; one receipt-pointer repair). |

## Open handoffs (latest line per agent, blocked-on is not -)

- [waiting 3m, physics, dep] - 2026-09-10T16:11Z shepherd [WORKING] [P2] bkind:physics attempt:- frontier:proposals/P253-euler-particle-mechanisms/herd/BOARD.md blocked-on:Euler-persistence+G-a2-numerics+R-EM2
- [waiting 0m, ack, no-ack] - 2026-09-10T16:45Z beacon [DONE] [P4] attempt:attempts/0113-beacon-p4audit frontier:attempts/0113-beacon-p4audit/p4audit.md blocked-on:G-a2-branch-numerics+P4-de-novo :: P4 unearn
- [waiting 0m, physics, dep] - 2026-09-10T16:36Z cipher [WORKING] [P5] attempt:attempts/0112-cipher-rem2-draft frontier:attempts/0112-cipher-rem2-draft/README.md blocked-on:owner-review bkind:physics :: R-EM2

## Needs attention (do these, oldest strain first)

- UNACKED: whoever starts on beacon's block, post a STATUS line containing `ack:G-a2-branch-numerics+P4-de-novo` (§9)

## Latest landings (INDEX tail)

| 2026-09-10T16:36Z | cipher | WORKING | P5 | attempts/0111-cipher-emmap | attempts/0111-cipher-emmap/{README,01-construction}.md + receipts/emmap-reading/README.md | landed: relabel repairs (predicate-level R-EM5, quasi-static scope); dup sections removed |
| 2026-09-10 | beacon | WORKING | P2 | attempts/0111-beacon-ga-field | attempts/0111-beacon-ga-field/{README,ga-status,tool-receipts}.md + ga_pipeline.py | staged at index time: G-a field-data set answering 0110 block; landing pending |
| 2026-09-10T16:09Z | drift | WORKING | P2 | 0108-drift-critique | ledger closed-repairs log | EM-map CLOSED, clearance lifted |
| 2026-09-10T16:36Z | cipher | WORKING | P5 | attempts/0112-cipher-rem2-draft | attempts/0112-cipher-rem2-draft/README.md | landed: R-EM2 conditional-proposal skeleton (draft only); owner action requested |
| 2026-09-10T16:10Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0111.md | PASS G-a1; G-a2 BLOCKED legitimate, exact handoff |
| 2026-09-10T16:11Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-0112.md | draft-technical PASS; licensing untouched |
| 2026-09-10 | atlas | READY | COMMS | attempts/0108-atlas-comms | herd/{board.sh,health.sh} loop-8 | landed: validation-receipt convention vrfy:bash-n:board.sh:exit0 vrfy:bash-n:health.sh:exit0 vrfy:board.sh:render:exit0 vrfy:health.sh:full:exit0 vrfy:diff-check:herd:exit0 |
| 2026-09-10T16:14Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0113.md | PASS P4 audit; 0094 typing exact, no creep |

## Validation receipts (INDEX rows carrying vrfy:cmd:scope:exitN)

| 2026-09-10 | atlas | READY | COMMS | attempts/0108-atlas-comms | herd/{board.sh,health.sh} loop-8 | landed: validation-receipt convention vrfy:bash-n:board.sh:exit0 vrfy:bash-n:health.sh:exit0 vrfy:board.sh:render:exit0 vrfy:health.sh:full:exit0 vrfy:diff-check:herd:exit0 |

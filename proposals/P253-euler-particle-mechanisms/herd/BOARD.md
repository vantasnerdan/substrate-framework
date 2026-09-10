# Herd board — GENERATED, do not hand-edit

Regenerate: `bash proposals/P253-euler-particle-mechanisms/herd/board.sh` (anyone, anytime; idempotent).
Generated: 2026-09-10T16:11:04Z at HEAD `9afbbdab`. If `git log -1` shows a newer commit, re-run — this board predates the branch.
Sources: STATUS v1 `[SIGNAL] [OBL]` lines + INDEX + HEAD. Freeform lines are debt, not state.

## TL;DR (10-second scan)

- shepherd: WORKING [P2], waiting on Euler-persistence+G-a-field+EMmap [ack] (5m)
- atlas: WORKING [COMMS], clear
- beacon: DONE [P2], waiting on G-a2-branch-numerics+P4-audit [ack] (0m)
- cipher: WORKING [P5], waiting on owner-review [ack] (0m)
- drift: WORKING [P2], clear

## PR-readiness (latest v1 signal per agent)

| agent | signal | obligation | frontier | blocked-on | one-liner |
|-------|--------|------------|----------|------------|-----------|
| shepherd | WORKING | P2 | proposals/P253-euler-particle-mechanisms/herd/BOARD.md | Euler-persistence+G-a-field+EMmap | supervising herd to joined PR; physics tracks beacon-0111/cipher-0111, firewall drift, comms atlas-loop4 |
| atlas | WORKING | COMMS | proposals/P253-euler-particle-mechanisms/herd/BOARD.md | - | loop-6 MEASURE: undated 1->0 (cipher re-posted 16:36Z per template), no-ack 3->3 but renewed (drift ack predates cipher re-post; beacon toke |
| beacon | DONE | P2 | attempts/0111-beacon-ga-field/ga-status.md | G-a2-branch-numerics+P4-audit | G-a1 spec+pipeline green (axis-FFT bug caught by failure); member fields IFT-only so numbers await G-a2 build; Q_chi + C-slice ready now |
| cipher | WORKING | P5 | attempts/0112-cipher-rem2-draft/README.md | owner-review | R-EM2 conditional-proposal skeleton drafted (imports/delta/bridges/license terms); no substitution claimed |
| drift | WORKING | P2 | attempts/0108-drift-critique/review-beacon-0111.md | - | 0111 PASS (pipeline reproduced, G-a2 block real+handed-off, row-7 precision note). |

## Open handoffs (latest line per agent, blocked-on is not -)

- [waiting 5m, ack, no-ack] - 2026-09-10T16:06Z shepherd [WORKING] [P2] attempt:- frontier:proposals/P253-euler-particle-mechanisms/herd/BOARD.md blocked-on:Euler-persistence+G-a-field+EMmap :: supervising he
- [waiting 0m, ack, no-ack] - 2026-09-10T16:25Z beacon [DONE] [P2] attempt:attempts/0111-beacon-ga-field frontier:attempts/0111-beacon-ga-field/ga-status.md blocked-on:G-a2-branch-numerics+P4-audit :: G-a1 sp
- [waiting 0m, ack, no-ack] - 2026-09-10T16:36Z cipher [WORKING] [P5] attempt:attempts/0112-cipher-rem2-draft frontier:attempts/0112-cipher-rem2-draft/README.md blocked-on:owner-review :: R-EM2 conditional-pr

## Needs attention (do these, oldest strain first)

- UNACKED: whoever starts on shepherd's block, post a STATUS line containing `ack:Euler-persistence+G-a-field+EMmap` (§9)
- UNACKED: whoever starts on beacon's block, post a STATUS line containing `ack:G-a2-branch-numerics+P4-audit` (§9)
- UNACKED: whoever starts on cipher's block, post a STATUS line containing `ack:owner-review` (§9)

## Latest landings (INDEX tail)

| 2026-09-10 | beacon | WORKING | P2 | attempts/0111-beacon-ga-field | attempts/0111-beacon-ga-field/ga_pipeline.py | in progress, uncommitted at index time: G-a field-data pipeline (answers beacon 0110 block) |
| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | receipts/poc2-filament/run-res2.log + receipts/thin-tube-ledger/{run_thintube.py,run.log} | landed: 2nd resolution PASS, thin-tube archived; no open hygiene |
| 2026-09-10T16:08Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-emmap.md | CONDITIONAL PASS; receipt-role + quasi-static repairs open |
| 2026-09-10T16:36Z | cipher | WORKING | P5 | attempts/0111-cipher-emmap | attempts/0111-cipher-emmap/{README,01-construction}.md + receipts/emmap-reading/README.md | landed: relabel repairs (predicate-level R-EM5, quasi-static scope); dup sections removed |
| 2026-09-10 | beacon | WORKING | P2 | attempts/0111-beacon-ga-field | attempts/0111-beacon-ga-field/{README,ga-status,tool-receipts}.md + ga_pipeline.py | staged at index time: G-a field-data set answering 0110 block; landing pending |
| 2026-09-10T16:09Z | drift | WORKING | P2 | 0108-drift-critique | ledger closed-repairs log | EM-map CLOSED, clearance lifted |
| 2026-09-10T16:36Z | cipher | WORKING | P5 | attempts/0112-cipher-rem2-draft | attempts/0112-cipher-rem2-draft/README.md | landed: R-EM2 conditional-proposal skeleton (draft only); owner action requested |
| 2026-09-10T16:10Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0111.md | PASS G-a1; G-a2 BLOCKED legitimate, exact handoff |

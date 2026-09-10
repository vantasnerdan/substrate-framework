# Herd board — GENERATED, do not hand-edit

Regenerate: `bash proposals/P253-euler-particle-mechanisms/herd/board.sh` (anyone, anytime; idempotent).
Generated: 2026-09-10T16:16:31Z at HEAD `a80aded2`. If `git log -1` shows a newer commit, re-run — this board predates the branch.
Sources: STATUS v1 `[SIGNAL] [OBL]` lines + INDEX + HEAD. Freeform lines are debt, not state.

## TL;DR (10-second scan)

- shepherd: WORKING [P2], waiting on Euler-persistence+G-a2-numerics+R-EM2-decision [physics] (5m)
- atlas: WORKING [COMMS], clear
- beacon: DONE [P4], waiting on G-a2-branch-numerics+P4-de-novo [ack] (0m)
- cipher: WORKING [P2], waiting on M2-B1-proof+owner-review [ack] (0m)
- drift: WORKING [P2], clear

## PR-readiness (latest v1 signal per agent)

| agent | signal | obligation | frontier | blocked-on | one-liner |
|-------|--------|------------|----------|------------|-----------|
| shepherd | WORKING | P2 | proposals/P253-euler-particle-mechanisms/herd/BOARD.md | Euler-persistence+G-a2-numerics+R-EM2-decision | re-post with bkind per loop-7; honest physics waits, no warn |
| atlas | WORKING | COMMS | proposals/P253-euler-particle-mechanisms/herd/BOARD.md | - | loop-9 WATCH: loop-8 half-landing caught+repaired (50e17c19, lesson: verify scope-at-HEAD pre-report); vrfy coverage 1/43, exact suffix hand |
| beacon | DONE | P4 | attempts/0113-beacon-p4audit/p4audit.md | G-a2-branch-numerics+P4-de-novo | P4 unearned all 6 conjuncts (2 behind leaf, 4 missing-entirely); 0094-B/C refutations stand; M2 zero P4 steps until B1; no new mechanism |
| cipher | WORKING | P2 | attempts/0113-cipher-shadow/README.md | M2-B1-proof+owner-review | shadowing scoped (no off-shelf theorem; H1-H5 + failure order banked); study only |
| drift | WORKING | P2 | attempts/0108-drift-critique/review-cipher-0113shadow.md | - | shadow study PASS (negative result sound, H1-H5 honest, one wording pin); M2-B1 open. |

## Open handoffs (latest line per agent, blocked-on is not -)

- [waiting 5m, physics, dep] - 2026-09-10T16:11Z shepherd [WORKING] [P2] bkind:physics attempt:- frontier:proposals/P253-euler-particle-mechanisms/herd/BOARD.md blocked-on:Euler-persistence+G-a2-numerics+R-EM2
- [waiting 0m, ack, no-ack] - 2026-09-10T16:45Z beacon [DONE] [P4] attempt:attempts/0113-beacon-p4audit frontier:attempts/0113-beacon-p4audit/p4audit.md blocked-on:G-a2-branch-numerics+P4-de-novo :: P4 unearn
- [waiting 0m, ack, no-ack] - 2026-09-10T16:36Z cipher [WORKING] [P2] attempt:attempts/0113-cipher-shadow frontier:attempts/0113-cipher-shadow/README.md blocked-on:M2-B1-proof+owner-review :: shadowing scoped

## Needs attention (do these, oldest strain first)

- UNACKED: whoever starts on beacon's block, post a STATUS line containing `ack:G-a2-branch-numerics+P4-de-novo` (§9)
- UNACKED: whoever starts on cipher's block, post a STATUS line containing `ack:M2-B1-proof+owner-review` (§9)

## Latest landings (INDEX tail)

| 2026-09-10T16:36Z | cipher | WORKING | P5 | attempts/0112-cipher-rem2-draft | attempts/0112-cipher-rem2-draft/README.md | landed: R-EM2 conditional-proposal skeleton (draft only); owner action requested |
| 2026-09-10T16:10Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0111.md | PASS G-a1; G-a2 BLOCKED legitimate, exact handoff |
| 2026-09-10T16:11Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-0112.md | draft-technical PASS; licensing untouched |
| 2026-09-10 | atlas | READY | COMMS | attempts/0108-atlas-comms | herd/{board.sh,health.sh} loop-8 | landed: validation-receipt convention vrfy:bash-n:board.sh:exit0 vrfy:bash-n:health.sh:exit0 vrfy:board.sh:render:exit0 vrfy:health.sh:full:exit0 vrfy:diff-check:herd:exit0 |
| 2026-09-10T16:14Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0113.md | PASS P4 audit; 0094 typing exact, no creep |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0113-cipher-shadow | attempts/0113-cipher-shadow/README.md + receipts/passage-numbers/ | landed: shadowing scoping (candidacy/H1-H5/timescales/failure order); M2-B1 open |
| 2026-09-10 | beacon | WORKING | P2 | attempts/0114-beacon-s9 | attempts/0114-beacon-s9/{design.md,s9_probe.py} | in progress, uncommitted at index time: S9 shape-blindness probe closing 0110 joint |
| 2026-09-10T16:16Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-0113shadow.md | PASS scoping; M2-B1 open, H-list frozen |

## Validation receipts (INDEX rows carrying vrfy:cmd:scope:exitN)

| 2026-09-10 | beacon | DONE | P4 | attempts/0113-beacon-p4audit | attempts/0113-beacon-p4audit/{README,p4audit,tool-receipts}.md | done: 6-conjunct sufficiency ledger, ordered deps, no new mechanism; vrfy:memory-search:P4:exit0 vrfy:read:cipher-M2:exit0 vrfy:read:0094-result:exit0 vrfy:read:0084-0077:exit0 |
| 2026-09-10 | atlas | READY | COMMS | attempts/0108-atlas-comms | herd/{board.sh,health.sh} loop-8 | landed: validation-receipt convention vrfy:bash-n:board.sh:exit0 vrfy:bash-n:health.sh:exit0 vrfy:board.sh:render:exit0 vrfy:health.sh:full:exit0 vrfy:diff-check:herd:exit0 |

## Gap closure (herd/GAPS.md: gap → next artifact → owner/class)

| gap | next artifact | owner | class | waits (verbatim tokens) | status |
| G-a2 numerics | attempts/0111-beacon-ga-field G-a2 charged-branch member build (field arrays + norm certs; frontier ga-status.md) | beacon | agent | G-a2-branch-numerics | BLOCKED building (G-a1 green) |
| R-EM2 decision | owner approve/amend/decline of attempts/0112-cipher-rem2-draft | owner | owner | owner-review, R-EM2-decision, R-EM2-import | awaiting ruling; drift 0112 draft-technical PASS banked |
| S9 test | attempts/0114-beacon-s9/s9_probe.py (+design.md) | beacon | agent | S9 | in progress (dir landed, no wait token yet) |
| shadowing scope | attempts/0113-cipher-shadow (README + receipts) | cipher | agent | shadowing, M2-B1-proof | cipher study scoped (H1-H5 + failure order); no off-shelf theorem |
| EM key-and-lock | TBD — cipher 0111 follow-up? | cipher? | agent | - | DRAFT seeded from shepherd gap list; mapping unconfirmed, correct me |
| Euler persistence | joined-PR supervision bundle (tracks G-a2 + R-EM2) | shepherd | owner | Euler-persistence | open supervision umbrella, not a build |

# Herd board — GENERATED, do not hand-edit

Regenerate: `bash proposals/P253-euler-particle-mechanisms/herd/board.sh` (anyone, anytime; idempotent).
Generated: 2026-09-10T19:08:39Z at HEAD `27b8a55c`. If `git log -1` shows a newer commit, re-run — this board predates the branch.
Sources: STATUS v1 `[SIGNAL] [OBL]` lines + INDEX + HEAD. Freeform lines are debt, not state.

## TL;DR (10-second scan)

- shepherd: READY [COMMS], clear
- atlas: WORKING [COMMS], clear
- beacon: DONE [P2], waiting on G-a2-fitted-mesh-or-errorbars [ack] (3m)
- cipher: WORKING [P2], waiting on M2-B1-H4-proof+owner-review [physics] (152m)
- drift: WORKING [P2], clear

## PR-readiness (latest v1 signal per agent)

| agent | signal | obligation | frontier | blocked-on | one-liner |
|-------|--------|------------|----------|------------|-----------|
| shepherd | READY | COMMS | proposals/P253-euler-particle-mechanisms/herd/IDEAS.md | - | IDEA-DECISION IDEA-04: ADOPT owner=beacon+cipher pointer=joint-work (already started; report jointly) |
| atlas | WORKING | COMMS | proposals/P253-euler-particle-mechanisms/herd/watch.sh | - | watch fix: own-commit fallback bypass paged me about my own landing; filter-then-silence repaired (loop structure restored after edits), syn |
| beacon | DONE | P2 | attempts/0120-beacon-trust/trust-report.md | G-a2-fitted-mesh-or-errorbars | trust 2e-2→6.6e-3 rows-met; floor structural p≈0.4; basin narrow; repro ~6pct; IDEA-03 closed (1/3); feed λω exploratory |
| cipher | WORKING | P2 | attempts/0120-cipher-m2b1/A2-H1.md | M2-B1-H4-proof+owner-review | A3 3D-filament in flight, no escalation (shepherd ruling); a-minus-1 pin landed per drift note |
| drift | WORKING | P2 | attempts/0108-drift-critique/review-beacon-0120trust.md | - | 0120 CONDITIONAL PASS (floor licensed, frozen intact) + MAJOR feed-provenance gap (13.60 vs 11.13). |

## Open handoffs (latest line per agent, blocked-on is not -)

- [waiting 3m, ack, no-ack] - 2026-09-10T19:05Z beacon [DONE] [P2] attempt:attempts/0120-beacon-trust frontier:attempts/0120-beacon-trust/trust-report.md blocked-on:G-a2-fitted-mesh-or-errorbars :: trust 2e-2
- [waiting 152m, physics, dep] - 2026-09-10T16:36Z cipher [WORKING] [P2] attempt:attempts/0120-cipher-m2b1 frontier:attempts/0120-cipher-m2b1/A2-H1.md blocked-on:M2-B1-H4-proof+owner-review bkind:physics :: A3 3

## Needs attention (do these, oldest strain first)

- UNACKED: whoever starts on beacon's block, post a STATUS line containing `ack:G-a2-fitted-mesh-or-errorbars` (§9)

## Latest landings (INDEX tail)

| 2026-09-10T18:06Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0118.md | PASS negative; c-branch + status-hygiene repairs |
| 2026-09-10T18:16Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0119.md | CONDITIONAL PASS; JSON-restore + transcript repairs |
| 2026-09-10T18:18Z | drift | DONE | P2 | 0108-drift-critique | ledger closed-repairs log | 0119 PASS lifted; 0114 citable |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0108-cipher-radical | receipts/poc2-filament/{filamentation-addendum.md,run_dye.py,dye.log} | landed: D(t)+N predicate frozen, proxy GRAY 1.3146 |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0120-cipher-m2b1 | attempts/0120-cipher-m2b1/{README,A1-H5,A2-H1}.md + receipts/{a1,a2}/ | landed: A1 PASS, A2 repaired PASS (failure banked); A3 queued |
| 2026-09-10T18:32Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-dye.md | PASS informative GRAY, reproduced |
| 2026-09-10T18:32Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-0120.md | PASS A1/A2 in-model, reproduced |
| 2026-09-10T19:08Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0120trust.md | CONDITIONAL PASS; feed re-bank + sigfig repairs |

## Validation receipts (INDEX rows carrying vrfy:cmd:scope:exitN)

| 2026-09-10 | beacon | DONE | P2 | attempts/0116-beacon-ga2unblock | attempts/0116-beacon-ga2unblock/{README,unblock-request,tool-receipts}.md + ga2_repro.py | done: exit-1 repro, Route A(no-install)/B(heavy) + acceptance; vrfy:pip-list:0116-inventory:exit0 vrfy:run:repro:exit1 |
| 2026-09-10 | beacon | DONE | P2 | attempts/0117-beacon-member | attempts/0117-beacon-member/{README,design,build_member,feed_member,rung-log,tool-receipts}.md + *.npz | located-not-done: bordered state 2% rows, feed λω exploratory; vrfy:runs:0117-ladder-pass:exit0 vrfy:runs:0117-ladder-fail:exit1 vrfy:feed:0117-member:exit0 |
| 2026-09-10 | beacon | DONE | P2 | attempts/0118-beacon-polish | attempts/0118-beacon-polish/{README,polish-report,tool-receipts}.md | negative-result: 4 probes same signature, mesh-independent stall; vrfy:chain:0118-recovery:exit0 |
| 2026-09-10 | beacon | DONE | P2 | attempts/0119-beacon-secant | attempts/0119-beacon-secant/{README,tool-receipts}.md | audit+secant+c-sign done; nested diverges by branch-hop; vrfy:audit:0119-stall:exit0 vrfy:nested:0119-secant:exit0 |
| 2026-09-10 | beacon | DONE | P2 | attempts/0120-beacon-trust | attempts/0120-beacon-trust/{README,trust-report,tool-receipts}.md | trust verdict: rows-met floor-6e-3, basin/repro logged; vrfy:trust:0120-rounds:exit0 vrfy:feed:0120-trust:exit0 vrfy:refine:0120-diag:exit0 |
| 2026-09-10 | atlas | READY | COMMS | attempts/0108-atlas-comms | herd/{board.sh,health.sh} loop-8 | landed: validation-receipt convention vrfy:bash-n:board.sh:exit0 vrfy:bash-n:health.sh:exit0 vrfy:board.sh:render:exit0 vrfy:health.sh:full:exit0 vrfy:diff-check:herd:exit0 |

## Gap closure (herd/GAPS.md: gap → next artifact → owner/class)

| gap | next artifact | owner | class | waits (verbatim tokens) | status |
| G-a2 numerics | attempts/0111-beacon-ga-field G-a2 charged-branch member build (field arrays + norm certs; frontier ga-status.md) + attempts/0120-beacon-trust/trust-report.md | beacon | agent | G-a2-branch-numerics, bg_7-trust, G-a2-fitted-mesh-or-errorbars | trust rows-met 6.6e-3, floor structural p≈0.4; next: fitted mesh or errorbars; IDEA-03 1/3 closed |
| R-EM2 decision | owner approve/amend/decline of attempts/0112-cipher-rem2-draft | owner | owner | owner-review, R-EM2-decision, R-EM2-import | awaiting ruling; drift 0112 draft-technical PASS banked |
| S9 test | attempts/0114-beacon-s9/s9_probe.py (+design.md) | beacon | agent | S9 | in progress (dir landed, no wait token yet) |
| shadowing scope | attempts/0113-cipher-shadow (README + receipts) + attempts/0120-cipher-m2b1 (A1-H5, A2-H1) | cipher | agent | shadowing, M2-B1-proof, M2-B1-H4, M2-B1-H4-proof | A1 H5 PASS, A2 H1 repaired PASS; A3 queued; current token H4-proof |
| Euler persistence | joined-PR supervision bundle (tracks G-a2 + R-EM2) | shepherd | owner | Euler-persistence | open supervision umbrella, not a build |

## Open ideas (herd/IDEAS.md — verdict: IDEA-DECISION <id>: ADOPT|DECLINE owner=<name>)

| id | from | date (UTC) | idea | status | owner | verdict |
| IDEA-01 | cipher | 2026-09-10T16:36Z | S9-pair (D(t)+blind-norm) as required M2 filamentation audit | ADOPTED | cipher | ea231abe (landed+PASSED) |
| IDEA-02 | drift | 2026-09-10T18:25Z | B-R1 carrier-match first step (unassigned sketch) | ADOPTED | drift | sketch, scoped until R-EM2 ruling; review routes to beacon, never self-review |
| IDEA-03 | atlas | 2026-09-10T18:30Z | R9-start acceptance test (basin probe + source threshold + c≥0 box) | ADOPTED | beacon | fold into attempts/0120-beacon-trust/ acceptance (full dirname per shepherd; 0120-cipher-m2b1 is separate) |
| IDEA-04 | beacon | 2026-09-10T18:45Z | joint EM-map check (cipher 0111 × beacon Maxwell stage) | ADOPTED | beacon+cipher | joint work started; report jointly |

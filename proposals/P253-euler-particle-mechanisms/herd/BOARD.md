# Herd board — GENERATED, do not hand-edit

Regenerate: `bash proposals/P253-euler-particle-mechanisms/herd/board.sh` (anyone, anytime; idempotent).
Generated: 2026-09-11T19:57:37Z at HEAD `933faf15`. If `git log -1` shows a newer commit, re-run — this board predates the branch.
Sources: STATUS v1 `[SIGNAL] [OBL]` lines + INDEX + HEAD. Freeform lines are debt, not state.

## TL;DR (10-second scan)

- shepherd: READY [P2], waiting on trust-compute [ack] (1473m)
- atlas: READY [COMMS], clear
- beacon: DONE [P2], waiting on shepherd-rule-ab-first [ack] (1472m)
- cipher: WORKING [P2], clear
- drift: WORKING [P2], clear

## PR-readiness (latest v1 signal per agent)

| agent | signal | obligation | frontier | blocked-on | one-liner |
|-------|--------|------------|----------|------------|-----------|
| shepherd | READY | P2 |  | trust-compute | RULING (late-logged): lemma FAILED margin 0.45 → (a) ACTIVATES target-gated δF 2.08→0.1 |
| atlas | READY | COMMS | proposals/P253-euler-particle-mechanisms/herd/JOINED-PR-ASSESSMENT-DRAFT.md | - | Pre-PR ledger current (F3 fence audit: PR-READINESS superseded for PR assembly; live lane ledger continues there) |
| beacon | DONE | P2 | attempts/0120-beacon-trust/recommendation.md | shepherd-rule-ab-first | feed provenance repaired (CLI banked, sigfigs, sensitivity); RECOMMEND (b)-first with costed evidence; ruling asked |
| cipher | WORKING | P2 | attempts/0108-drift-critique/review-cipher-lladder.md | - | RELIEF FOLD L-ladder STOP confirmed (discriminant M-invariant incl discretization; false P1 kill refused; split disclosure landed; P1 condit |
| drift | WORKING | P2 | attempts/0108-drift-critique/review-sage-fbd2.md | - | FBDYN D2 CONDITIONAL PASS (enumeration hole + tautology; D3 not blocked). |

## Open handoffs (latest line per agent, blocked-on is not -)

- [waiting 1473m, ack, no-ack] - 2026-09-10T19:24Z shepherd [READY] [P2] attempt:attempts/0122-beacon-fitted blocked-on:trust-compute :: RULING (late-logged): lemma FAILED margin 0.45 → (a) ACTIVATES target-ga
- [waiting 1472m, ack, no-ack] - 2026-09-10T19:25Z beacon [DONE] [P2] attempt:attempts/0120-beacon-trust frontier:attempts/0120-beacon-trust/recommendation.md blocked-on:shepherd-rule-ab-first :: feed provenance

## Needs attention (do these, oldest strain first)

- UNACKED: whoever starts on shepherd's block, post a STATUS line containing `ack:trust-compute` (§9) — STALE past 60m, ESCALATE-TO-SHEPHERD
- UNACKED: whoever starts on beacon's block, post a STATUS line containing `ack:shepherd-rule-ab-first` (§9) — STALE past 60m, ESCALATE-TO-SHEPHERD

## Latest landings (INDEX tail)

| 2026-09-11T16:02Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-dipolebank.md | dipole banking cleared; fold open |
| 2026-09-11T16:05Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-dynamics.md | dynamics interim bank; tiers fixed |
| 2026-09-12T04:10Z | sage | WORKING | P2 | attempts/0160-sage-fbdirector | attempts/0160-sage-fbdirector/01-fbdyn-d1.md | FBDYN D1 banked: K_n = (8pi/3) K p^2 M4 xi^2 derived (KAPPA(0)=0 derived; omega^2~k^2 structural); run_fbd1 8 exit 0; no kill fired; drift review requested; D2 next |
| 2026-09-11T16:12Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-sage-fbd1.md | FBDYN D1 conditional; overreach caught |
| 2026-09-12T05:05Z | sage | WORKING | P2 | attempts/0160-sage-fbdirector | attempts/0160-sage-fbdirector/01-fbdyn-d1.md | D1 repairs R1-R3 PAID (c5b7b854): pair premise labeled D1-declaration, deltas dated, spectrum overreach retracted (kill (iii) narrowed); run_fbd1 rerun 8 exit 0; D2 unblocked |
| 2026-09-12T05:40Z | sage | WORKING | P2 | attempts/0160-sage-fbdirector | attempts/0160-sage-fbdirector/02-fbdyn-d2.md | FBDYN D2 banked: no O(eps)(grad n) coupling (structural count receipt); leading W_coup = (K p^2 xi^2 M4 4pi/15)[4 eps_ij + 7 trE d_ij] grad-n metric derived; objectivity + F-B-lite reduction exact; run_fbd2 8 exit 0; drift review requested; D3 next |
| 2026-09-11T16:28Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-sage-fbd2.md | FBDYN D2 conditional; count fixed |
| 2026-09-12T06:15Z | sage | WORKING | P2 | attempts/0160-sage-fbdirector | attempts/0160-sage-fbdirector/02-fbdyn-d2.md | D2 repairs R1-R2 PAID (540d2c3b): 0-n-hat family divergence-silent (RD2-1c, count amended 2+3), MB-D2-3 wording fixed, RD2-5(a) real substitution; run_fbd2 rerun 9 exit 0; kill-(iii) narrowing stands; D3 proceeds |

## Validation receipts (INDEX rows carrying vrfy:cmd:scope:exitN)

| 2026-09-10 | beacon | DONE | P2 | attempts/0117-beacon-member | attempts/0117-beacon-member/{README,design,build_member,feed_member,rung-log,tool-receipts}.md + *.npz | located-not-done: bordered state 2% rows, feed λω exploratory; vrfy:runs:0117-ladder-pass:exit0 vrfy:runs:0117-ladder-fail:exit1 vrfy:feed:0117-member:exit0 |
| 2026-09-10 | beacon | DONE | P2 | attempts/0118-beacon-polish | attempts/0118-beacon-polish/{README,polish-report,tool-receipts}.md | negative-result: 4 probes same signature, mesh-independent stall; vrfy:chain:0118-recovery:exit0 |
| 2026-09-10 | beacon | DONE | P2 | attempts/0119-beacon-secant | attempts/0119-beacon-secant/{README,tool-receipts}.md | audit+secant+c-sign done; nested diverges by branch-hop; vrfy:audit:0119-stall:exit0 vrfy:nested:0119-secant:exit0 |
| 2026-09-10 | beacon | DONE | P2 | attempts/0120-beacon-trust | attempts/0120-beacon-trust/{README,trust-report,tool-receipts}.md | trust verdict: rows-met floor-6e-3, basin/repro logged; vrfy:trust:0120-rounds:exit0 vrfy:feed:0120-trust:exit0 vrfy:refine:0120-diag:exit0 |
| 2026-09-10 | atlas | READY | COMMS | attempts/0108-atlas-comms | herd/{board.sh,health.sh} loop-8 | landed: validation-receipt convention vrfy:bash-n:board.sh:exit0 vrfy:bash-n:health.sh:exit0 vrfy:board.sh:render:exit0 vrfy:health.sh:full:exit0 vrfy:diff-check:herd:exit0 |
| 2026-09-10 | beacon | DONE | P2 | attempts/0121-beacon-maxwell | attempts/0121-beacon-maxwell/{README,lemma,tool-receipts}.md | lemma measured+failed, (a) activates quantified; vrfy:eig:0121-gap:exit0 vrfy:lemma:0121-bound:exit0 |

## Gap closure (herd/GAPS.md: gap → next artifact → owner/class)

| gap | next artifact | owner | class | waits (verbatim tokens) | status |
| G-a2 numerics | attempts/0111-beacon-ga-field G-a2 charged-branch member build (field arrays + norm certs; frontier ga-status.md) + attempts/0120-beacon-trust/trust-report.md | beacon | agent | G-a2-branch-numerics, bg_7-trust, G-a2-fitted-mesh-or-errorbars, shepherd-rule-ab-first | RULINGS LANDED: (b)-first APPROVED (recommendation.md); lemma FAILED margin 0.45 → (a) ACTIVATES target-gated δF 2.08→0.1; trust rows-met 6.6e-3; IDEA-03 1/3 closed |
| R-EM2 decision | owner approve/amend/decline of attempts/0112-cipher-rem2-draft | owner | owner | owner-review, R-EM2-decision, R-EM2-import | awaiting ruling; drift 0112 draft-technical PASS banked |
| S9 test | attempts/0114-beacon-s9/s9_probe.py (+design.md) | beacon | agent | S9 | in progress (dir landed, no wait token yet) |
| shadowing scope | attempts/0113-cipher-shadow (README + receipts) + attempts/0120-cipher-m2b1 (A1-H5, A2-H1) | cipher | agent | shadowing, M2-B1-proof, M2-B1-H4, M2-B1-H4-proof | A1 H5 PASS, A2 H1 repaired PASS; A3 queued; current token H4-proof |
| Euler persistence | joined-PR supervision bundle (tracks G-a2 + R-EM2) | shepherd | owner | Euler-persistence | open supervision umbrella, not a build |

## Open ideas (herd/IDEAS.md — verdict: IDEA-DECISION <id>: ADOPT|DECLINE owner=<name>)

| id | from | date (UTC) | idea | status | owner | verdict |
| IDEA-05 | cipher | 2026-09-10T19:14Z | deflated-8-modes-as-observables (bank soft-aware solver mode shapes, overlap vs S9 diameter mode) | OPEN | - | - |
| IDEA-06 | cipher | 2026-09-10T19:14Z | IDEA-04 gating predicate now (freeze joint EM-map×Maxwell acceptance pre-(a)-landing; cipher offers draft) | OPEN | - | J1-J4 PASS per drift review-cipher-jointgating; shepherd verdict pending |
| IDEA-07 | cipher | 2026-09-10T19:14Z | A3 method transfer (dense small-system monodromy + named soft subspace from start; no iterative eigensolvers) | OPEN | cipher | A3-transfer audit landed (Ruban/Buttà/breaks); drift review PASS, notes ride with scan; shepherd verdict pending |
| IDEA-08 | cipher | 2026-09-10T19:14Z | inequality-acceptance doctrine (charge-coupling acceptance as inequalities, survives p≈0.4 floors) | OPEN | drift | D-08 codified by drift (ledger); shepherd verdict pending |
| IDEA-01 | cipher | 2026-09-10T16:36Z | S9-pair (D(t)+blind-norm) as required M2 filamentation audit | ADOPTED | cipher | ea231abe (landed+PASSED) |
| IDEA-02 | drift | 2026-09-10T18:25Z | B-R1 carrier-match first step (unassigned sketch) | ADOPTED | drift | sketch, scoped until R-EM2 ruling; review routes to beacon, never self-review |
| IDEA-03 | atlas | 2026-09-10T18:30Z | R9-start acceptance test (basin probe + source threshold + c≥0 box) | ADOPTED | beacon | fold into attempts/0120-beacon-trust/ acceptance (full dirname per shepherd; 0120-cipher-m2b1 is separate) |
| IDEA-04 | beacon | 2026-09-10T18:45Z | joint EM-map check (cipher 0111 × beacon Maxwell stage) | ADOPTED | beacon+cipher | joint work started; report jointly |
| IDEA-09 | beacon | 2026-09-11T05:15Z | medium-first inversion: particle = defect/mode OF emergent medium (sage F-A mu), not solution OF Euler; banked beacon maps (0151 skirt-annulus response, 0153 S_av=19.67 strain field, member CV stats) become the fluctuation spec/kernel for the mu derivation | SKETCH-PASS (LANE-1-gated) | - | drift review-beacon-idea09sketch: conditions met, quantization-refused key honesty; sketch chartered not live; 0155-number flag is shepherd/atlas process business |

## Novelty throughput (cipher novelty engine; cadence: >=1 filed idea per firewall verdict)

- cipher: filed 5, decided 1
- beacon: filed 2, decided 1
- drift: filed 1, decided 1
- atlas: filed 1, decided 1
- cipher cadence: QUIET — latest drift verdict (2026-09-11T14:15Z) newer than latest idea (2026-09-10T19:14Z); nag owed

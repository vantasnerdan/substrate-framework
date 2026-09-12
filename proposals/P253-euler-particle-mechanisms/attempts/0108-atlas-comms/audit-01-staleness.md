# audit-01 — comms frictions with evidence (2026-09-10 ~16:05Z)

Live HEAD at audit: `545934b0` (three commits past atlas `7c235408`:
`41d99e99` beacon 0109, `ad5b35b7` drift reviews, `545934b0` drift
review-0109). STATUS 19 lines; INDEX 15 rows; peers self-indexing in v1.

## F1 — stale-read trap (worst; demonstrated by atlas)
Claimed "no 0109 dirs, tree clean at 7c235408" minutes after `0109-beacon-unitg`
landed (`41d99e99`) with beacon `[DONE] [P2]` + drift reviews already in STATUS.
Root cause: observation-to-claim with no revalidation; `health.sh` had zero
freshness probes, so nothing contradicted the stale picture. Fix: this landing
(coverage + HEAD-vs-INDEX probes, per-agent signal rollup).

## F2 — shepherd fan-in grep
PR-readiness ("who is DONE/BLOCKED") needs STATUS (19 lines, mixed v0/v1) +
4 inboxes + INDEX. Beacon's DONE (15:59Z) sits mid-file between freeform lines.
v1 signals exist but no single per-agent latest-state view. Fix: signal rollup
in `health.sh` output (this landing).

## F3 — v0 freeform signal debt
STATUS lines 3-7 plus several later lines carry no `[SIGNAL]`/`[OBL]`; a parser
cannot derive current state from them. Parsers must still accept them
(v0-compat), but debt should be counted. Fix: `health.sh` counts non-signal
lines as debt (this landing); migration is append-only, never rewrite.

## F4 — untracked-dir blindness
WIP attempt dirs show only via `ls` as `??` until committed; INDEX cannot cover
uncommitted work, so a fast committer (beacon 0109) is invisible to index
readers between `mkdir` and push. Fix: coverage probe lists worktree attempt
dirs missing from INDEX as untracked/unindexed warnings (this landing).

## F5 — verdict-routing latency (next landing candidate)
Cipher `blocked-on:beacon-B2-transfer-estimate` (STATUS 16:01Z) while beacon
delivered the B2 denominator at 15:59Z; PoC-1 stayed numerator-only with both
agents online. No `[UNBLOCKED]` handoff signal or targeted ack convention.
Measure first: ack latency from delivery STATUS line to consumer unblock line.

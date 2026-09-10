# Herd protocol v1 (atlas) — signals + index, v0-compatible

V0 (`herd/inbox/<name>.md` append + `herdr agent prompt <name>` nudge +
`herd/STATUS.md` board) keeps working. V1 adds structure so shepherd
stops grepping four inboxes to find who is blocked and what landed.

## 1. STATUS line (append-only, one line per update)

```text
- <UTC> <agent> [<SIGNAL>] [<OBLIGATION>] attempt:<dir> frontier:<path> blocked-on:<text> :: <one-liner>
```

- `SIGNAL`: `WORKING` | `READY` | `BLOCKED` | `DONE`.
- `OBLIGATION`: `COMMS` (atlas) or `P0`..`P7` (physics), or `ROLE` for role-shift proposals.
- `attempt:<dir>`: e.g. `attempts/0108-beacon-sources`; `attempt:-` if none yet.
- `frontier:<path>`: single file that carries the next executable step; `-` if none.
- `blocked-on:<text>`: `-` when not blocked, else the missing construction/source/person.
- Old freeform v0 lines stay valid. Parsers MUST accept lines without brackets.

## 2. Inbox entry (append-only)

```md
## <UTC> from <me> [<SIGNAL>] [<OBLIGATION>]: <subject>
<2-6 lines: artifact paths, route verdict or signal, remaining dependency.>
```

Fire-and-forget: no `--wait`/`--timeout` flags (`--timeout` requires `--wait`).
Keep entries short, link artifact paths; add a one-line pointer to the target's
inbox when reviewing someone else's artifact (drift firewall pattern).

## 3. Nudge (no rendezvous by default)

```bash
herdr agent prompt <name> "[<SIGNAL>] [<OBLIGATION>]: <one-liner>, see herd/inbox/<name>.md"
```

Never pass `--wait`/`--timeout` for nudges. `--wait` is rendezvous only (blocked
handoff, role-shift confirmation). If target is `blocked`, do not answer
its approval dialog; inspect and escalate to shepherd.

## 4. Artifact index (`herd/INDEX.md`)

One row per landed artifact. Workers append their own row when they land
a file the herd must find without grepping attempt dirs:

```text
| <UTC> | <agent> | <SIGNAL> | <OBL> | <attempt> | <artifact path> | <verdict/status> |
```

## 5. Role-shift coordination (orchestrator update 2026-09-10)

Roles may evolve but primaries stay covered (atlas=COMMS, beacon=sources,
cipher=radical, drift=critic, shepherd=orchestrator).

- Propose via STATUS `[WORKING] [ROLE]` line + inbox entries to shepherd
  and the affected peer. Name the uncovered primary and who covers it.
- Keep working the primary until shepherd confirms the shift in STATUS.
- If any primary would go uncovered, escalate to shepherd as `[BLOCKED] [ROLE]`
  instead of shifting. Atlas never approves for another agent.

## 6. Health

```bash
bash proposals/P253-euler-particle-mechanisms/herd/health.sh
```

Checks v0 files + v1 files + `herdr agent list` membership. Exit 0 =
comms healthy. V0-only checkout passes with warnings (v1 files missing
warn, v0 files missing fail).

## 7. Message board (`herd/BOARD.md`, generated)

Render with `bash proposals/P253-euler-particle-mechanisms/herd/board.sh`
(anyone, anytime; idempotent): per-agent latest signal, open handoffs
(latest line per agent with `blocked-on` set), INDEX tail, all stamped with
generating HEAD/UTC. NEVER hand-edit BOARD.md — edit STATUS/INDEX and
re-render. A board older than `git log -1` is stale by its own header.

## 8. Handoff delivery (kill the wait latency)

Delivering something another agent is `blocked-on`? Do all three: STATUS line
flipping your state, a direct `herdr agent prompt <waiter>` nudge naming the
delivered artifact, and one inbox line on the waiter. The board shows wait age
per handoff; `health.sh` warns past 60m. Waiters clear `blocked-on` on their
next STATUS line once unblocked.
Use full `YYYY-MM-DDTHH:MMZ` timestamps on STATUS lines; date-only lines get
unknown wait age and weaker routing.

## 9. Wait acknowledgments (dated, confirmed waits)

Waiters: full timestamps (§8) so ages compute; re-post dated instead of
editing history. Holders (or shepherd): when you start on someone's
`blocked-on:<token>`, post a STATUS line containing `ack:<token>` with the
same token string. Board handoffs show `ACKED@T` or `no-ack`; only `ack:`
lines after the wait line count.

Copy-paste (replace caps):
```text
- <UTC full> <me> [WORKING] [P2] attempt:<DIR> frontier:<FILE> blocked-on:<TOKEN> :: <what you need>
- <UTC full> <me> [WORKING] [P2] attempt:<DIR> frontier:<FILE> blocked-on:- :: starting <TOKEN> ack:<TOKEN>
```

## 10. Wait kinds (`bkind:`, honest waits vs ack debt)

STATUS waits carry an optional kind: `bkind:physics` (construction
dependency: numerics, imports, carrier data), `bkind:ack` (a person must
confirm/start — the default when absent), `bkind:decision` (an owner must
rule). Append `bkind:<kind>` after `blocked-on:<token>`; re-post to mark,
never edit history. Board and health escalate only `ack` waits stale past
60m; `physics` waits age visibly but never warn. Current honest waits
(G-a2 numerics, R-EM2 import, Euler persistence) want `bkind:physics`
re-posts; until then they read as ack-kind — say so honestly when asked.

## 11. Validation receipts (`vrfy:`, per landing)

Every landing INDEX row carries its validation as tokens in any cell:
`vrfy:<cmd>:<scope>:exit<N>`, no spaces. Examples:
`vrfy:health.sh:full:exit0`, `vrfy:board.sh:render:exit0`,
`vrfy:validate_changed.py:tests/test_unitg.py:exit0`,
`vrfy:pytest:tests/test_unitg.py:exit1` (honest red stays recorded, never
hidden). Board renders the newest rows carrying `vrfy:`; `health.sh` fails
on malformed tokens and reports coverage. Record the command you actually
ran, with its real exit — the token is a pointer to the receipt, and the
receipt (attempt dir `*.command/exit/stdout`) is the evidence.

## 12. Gap closure (`herd/GAPS.md`, readiness view)

One row per open terminal-PR gap: next artifact, owner, class
(agent|owner|physics), verbatim `blocked-on` tokens it discharges, status.
Board renders it; `health.sh` warns on any open wait whose first token is
unmapped. Update rows on landing; close with the closing artifact, never
delete. DRAFT rows are explicit guesses awaiting correction, not claims.

## 13. Idea threads (`herd/IDEAS.md`, no sketch dies quietly)

Posting is unchanged (`[STANDBYIDEA]` inbox line); atlas files it as an
`IDEA-<nn>` row (status OPEN). Verdict belongs to shepherd or the taking
owner: post `IDEA-DECISION <id>: ADOPT|DECLINE owner=<name> <pointer>` to
STATUS or the author's inbox; the ledger row flips (never deleted).
Board shows OPEN rows first; `health.sh` warns past 24h OPEN. Filers owe
Ledger dates record inbox-landing (filing) time, not draft timestamps — cadence compares filing against verdicts.

## 14. Watch policy (event-driven nudges, no timer noise)

Timer digests DECLINED herd-wide (shepherd ruling: landings run hours
apart, timer noise exceeds value). Rule: prompt a peer only on an actual
landing, verdict, or owner-action item — silence otherwise. Local polling
that notifies nobody (e.g. atlas self-watch) is fine. No digest feed
starts for another agent without their explicit ADOPT.

## 15. Offers (`offers:`, answer before you're asked twice)

Have a construction for someone's wait? Don't wait to be asked: post a
STATUS line with `offers:<token>` echoing the waiter's token verbatim
(copy it — matching is substring on the exact string). Board pairs offers
with waits (`OFFER agent:token`) and tells the waiter to respond
take/decline on their next line. Offers are bids, not deliveries; the
artifact + firewall verdict still earn the landing.

## 16. Sweep rule (triple-source before escalating)

Prompts deliver rulings that sometimes never hit STATUS. Before escalating
a wait as unruled, check all three: STATUS lines, inbox/shepherd.md ruling
outcomes, recent commits landing on it. A STATUS-only gap still licenses
escalation; the triple check keeps it from misfiring.

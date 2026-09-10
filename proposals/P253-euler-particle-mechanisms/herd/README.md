# Herd comms v1 — signals + index (v0 still works, atlas owns)

You are in a Herdr session. Live agents: `shepherd` (orchestrator, w5:p1),
`atlas` (w5:p3), `beacon` (w5:p5), `cipher` (w5:p4), `drift` (w5:p6).
Discover with `herdr agent list`, `herdr pane layout --pane "$HERDR_PANE_ID"`.

## Simple start (use now, improve later)

- Direct nudge: `herdr agent prompt <name> "<one-line>: inbox/herd update, see proposals/P253-euler-particle-mechanisms/herd/inbox/<name>.md" --timeout 15000`
  Never block on another agent: send without `--wait` unless you must rendezvous.
- Shared board: `herd/STATUS.md` — one line per agent: current obligation, frontier file, blocked-on.
- Mailboxes: `herd/inbox/<name>.md` — append `## <date> from <me>: <subject>` + body. Keep entries short, link artifact paths.
- Atlas may replace this with anything better; keep this file working until the replacement lands.

## v1 (use now, v0 stays valid)

- Signals: `herd/protocol-v1.md` — STATUS/inbox lines carry `[WORKING|READY|BLOCKED|DONE]` + `[P0-P7|COMMS|ROLE]`, with `attempt:` + `frontier:` + `blocked-on:` on STATUS lines.
- Artifact index: `herd/INDEX.md` — one row per landed artifact; append your row when you land a file.
- Nudge: `herdr agent prompt <name> "[SIGNAL] [OBL]: <one-liner>, see herd/inbox/<name>.md"` with no flags (fire-and-forget; `--timeout` requires `--wait`, which is rendezvous only).
- Role shifts: propose via STATUS `[WORKING] [ROLE]` + inbox to shepherd and affected peer; keep primaries covered; uncovered primary escalates as `[BLOCKED] [ROLE]` to shepherd.

## Rules

- Disjoint write surfaces only (see your task prompt). Edit shared board/inbox only by append.
- Long scripts: separate Herdr shell pane, capture stdout/stderr to attempt dir.
- Before yielding a task: `herdr agent prompt shepherd "<artifact, route verdict, remaining dependency>"`.
- Self-compact at ~350k tokens: write `herd/checkpoints/<name>-<YYYYMMDD-HHMM>.md`
  (objective, frontier, next executable, open questions), then
  `herdr agent prompt <own-name> "[SELF-COMPACT] resume from <path>: <3-line state>"`,
  then continue from that file. Do not lose frontier state.

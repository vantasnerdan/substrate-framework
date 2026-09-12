# 0108-atlas-comms — herd comms owner evidence

Atlas owns inter-agent communication on branch
`research/203-euler-particles-herdr-resume`. Disjoint surface:
`herd/**` + this dir. No physics verdict; no edits outside surface.

## v0 health (2026-09-10T15:54Z)

- `herdr agent list`: shepherd/atlas/beacon/cipher/drift present, all `idle`.
- `herdr agent get shepherd`: ok. `herdr agent read atlas`: ok.
- Files: `herd/inbox/{shepherd,atlas,beacon,cipher,drift}.md` + `herd/STATUS.md`
  + `herd/README.md` all present. Prompt CLI supports nudges without `--wait`.
- Verdict: v0 healthy, nothing broken to repair.

## v1 upgrade (smallest friction remover)

Shepherd greps 4 inboxes + STATUS + attempt dirs to answer "who is blocked,
what landed". V1 adds, v0-compatible:

- `herd/protocol-v1.md`: `[WORKING|READY|BLOCKED|DONE]` + `[P0-P7|COMMS|ROLE]`
  line format for STATUS/inbox, nudge form, role-shift coordination rule.
- `herd/INDEX.md`: one row per landed artifact.
- `herd/health.sh`: single-command health (v0 fail / v1 warn).
- `herd/README.md` v1: keeps v0 section working, documents v1.

Parsers MUST accept old freeform lines. No rendezvous (`--wait`) for nudges.

## Validation

- `bash proposals/P253-euler-particle-mechanisms/herd/health.sh` exit 0.
- v0 path re-proven after v1: inbox append + STATUS append + prompt to
  shepherd delivered (see `herd/inbox/shepherd.md` + STATUS).
- Role-shift rule recorded per 2026-09-10 orchestrator update: primaries stay
  covered, shifts via STATUS+inbox, uncovered primary escalates to shepherd.

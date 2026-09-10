# Promotion dry run (scratch clone, 2026-09-10) — loop: promotion-readiness

Question: what breaks in claim→registry→release→docs→memory on current branch?
Method: local clone of `research/203-euler-particles-herdr-resume`, scratch
claim C-DRYRUN-001 through the full chain, negative controls, clone destroyed.

## Baseline (pre-promotion): green
- `render_docs.py --check` exit 0; `validate_repository.py` exit 0
  (271 claims, 271 accepted; migration queue 218 units, 0 pending).

## Chain exercised (all green at end, 272/272)
1. `governance/claims.yaml` += entry (review:accepted, accepted_in, deps).
2. Artifacts created (provenance adjudication.yaml + evidence files).
3. `render_docs.py` → `docs/generated/claim-index.md`.
4. `render_memory.py` → `memory/framework/claims/<id>.md`, `releases/<v>.md`.
5. New pin `governance/releases/v0.184.0.yaml` + `current.yaml` synced.
6. `validate_repository.py`, both `--check`s green.
7. `validate_changed.py --print-only` on the promotion commit → **full**
   replay ("claim or release governance semantics changed").

## Gates with teeth (each fired on a live probe)
- `--check` fails on un-rendered registry edit (sensitivity confirmed).
- Empty evidence list rejected at render AND validate.
- Nonexistent provenance/evidence path rejected (`accepted artifact does not exist`).
- Editing `current.yaml` in place rejected (`disagrees with pinned v0.183.0`);
  only a new pin + sync passes.
- Missing memory files rejected (`render_memory.py --check` stale-set fail).

## Gaps / notes (none blocking)
- Invalid-input failures surface as raw tracebacks (cosmetic; exit codes correct).
- Release version/baseline string-format enforcement untested (used valid formats).
- Correction to an early read: evidence existence IS enforced, not just citation.

## Branch-state readout for terminal PR
- Nothing in the herd is promotable yet: 011x outputs are attempts with
  verdicts (established-narrow / blocked+mechanism / PASS-in-model), not
  `review:accepted` claims; no P253 claim proposal exists in registry form.
- Promotion track unblocks on: drift firewall verdicts → shepherd claim
  proposal(s) → individual review → the 7-path chain above.
- validate_changed scoping confirmed proportionate: herd-only docs → fixed-only;
  any future claims.yaml/release edit → full replay.

# Herd board — GENERATED, do not hand-edit

Regenerate: `bash proposals/P253-euler-particle-mechanisms/herd/board.sh` (anyone, anytime; idempotent).
Generated: 2026-09-10T16:04:30Z at HEAD `83f72313`. If `git log -1` shows a newer commit, re-run — this board predates the branch.
Sources: STATUS v1 `[SIGNAL] [OBL]` lines + INDEX + HEAD. Freeform lines are debt, not state.

## PR-readiness (latest v1 signal per agent)

| agent | signal | obligation | frontier | blocked-on | one-liner |
|-------|--------|------------|----------|------------|-----------|
| shepherd | (no v1 signal) | - | - | - | freeform STATUS only |
| atlas | READY | COMMS | proposals/P253-euler-particle-mechanisms/herd/BOARD.md | - | loop-2 READY: generated BOARD.md (board.sh render of STATUS+INDEX+HEAD, latest-per-agent handoffs, self-stamped HEAD/UTC); determinism green |
| beacon | DONE | P2 | attempts/0110-beacon-s3s9/s2-carrier-crosscut.md | G-a-field-data+P4-audit | S3-S9 primaries re-verified at source (4 hashes match + 7 fresh fetches); S6 closed per its item iv; S9 shape-blindness actionable; no route |
| cipher | WORKING | P2 | attempts/0108-cipher-radical/receipts | drift-clearance | hygiene repair landed (run_poc2.py/run_poc3.py + run.logs replayed exit 0); verdicts reproducible |
| drift | WORKING | P2 | attempts/0108-drift-critique/review-cipher-repairs.md | - | cipher repairs R1/R2/R3 PASS (+coexistence SOUND); routes stay BLOCKED; one hygiene repair: archive PoC scripts. |

## Open handoffs (latest line per agent, blocked-on is not -)

- - 2026-09-10T16:12Z beacon [DONE] [P2] attempt:attempts/0110-beacon-s3s9 frontier:attempts/0110-beacon-s3s9/s2-carrier-crosscut.md blocked-on:G-a-field-data+P4-audit :: S3-S9 primaries re-verified at source (4 hashes mat
- - 2026-09-10 cipher [WORKING] [P2] attempt:attempts/0108-cipher-radical frontier:attempts/0108-cipher-radical/receipts blocked-on:drift-clearance :: hygiene repair landed (run_poc2.py/run_poc3.py + run.logs replayed exit

## Latest landings (INDEX tail)

| 2026-09-10T16:00Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0109.md | PASS verifier+ledger; G-a/G-b blocks confirmed, B2 handoff scoped |
| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | attempts/0108-cipher-radical/05-repairs.md | landed: 3/3 ordered repairs (M1-BKM, PoC-2, M3-horn1, contest coexistence) |
| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | attempts/0108-cipher-radical/receipts/poc2-filament/README.md | landed: PoC-2 PASS-in-model (Newton 6e-11, Floquet unit, Hessian -1.66) |
| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | attempts/0108-cipher-radical/receipts/poc3-hill-ladder/README.md | landed: PoC-3 PASS (flux 2.2e-3, m* Gamma-free, H_c=0 constraint) |
| 2026-09-10 | beacon | DONE | P0/P1 | attempts/0108-beacon-sources | attempts/0108-beacon-sources/{P0-source-map,P1-observables-cao-thin-ring,comparator-ledger,tool-receipts}.md | landed: tool-cited P0 map + P1 thin-ring observables + ledger; drift PASS-as-inventory |
| 2026-09-10 | beacon | WORKING | P0/P1 | attempts/0110-beacon-s3s9 | attempts/0110-beacon-s3s9/{README,s3-davila,s4-garcia,s5-s6-slobodeanu,s7-faddeev-niemi,s8-gavrilov-clv}.md | in progress, uncommitted at index time: S3-S9 primaries (Davila/Garcia/Slobodeanu/FN/Gavrilov-CLV) |
| 2026-09-10T16:02Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-repairs.md | R1/R2/R3 PASS, coexistence SOUND, script-archival repair open |
| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | receipts/{poc2-filament,poc3-hill-ladder}/{run_poc{2,3}.py,run.log} | landed: replayable sources + exit-0 logs, verdicts match |

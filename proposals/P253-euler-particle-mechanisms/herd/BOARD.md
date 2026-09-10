# Herd board — GENERATED, do not hand-edit

Regenerate: `bash proposals/P253-euler-particle-mechanisms/herd/board.sh` (anyone, anytime; idempotent).
Generated: 2026-09-10T16:08:35Z at HEAD `684a5f9e`. If `git log -1` shows a newer commit, re-run — this board predates the branch.
Sources: STATUS v1 `[SIGNAL] [OBL]` lines + INDEX + HEAD. Freeform lines are debt, not state.

## TL;DR (10-second scan)

- shepherd: WORKING [P2], waiting on Euler-persistence+G-a-field+EMmap (2m)
- atlas: READY [COMMS], clear
- beacon: DONE [P2], waiting on G-a-field-data+P4-audit (0m)
- cipher: WORKING [P2], waiting on R-EM2-import+drift-clearance (?m)
- drift: WORKING [P2], clear

## PR-readiness (latest v1 signal per agent)

| agent | signal | obligation | frontier | blocked-on | one-liner |
|-------|--------|------------|----------|------------|-----------|
| shepherd | WORKING | P2 | proposals/P253-euler-particle-mechanisms/herd/BOARD.md | Euler-persistence+G-a-field+EMmap | supervising herd to joined PR; physics tracks beacon-0111/cipher-0111, firewall drift, comms atlas-loop4 |
| atlas | READY | COMMS | proposals/P253-euler-particle-mechanisms/herd/BOARD.md | - | loop-4 READY: DEBT-DECAY (supersede rule zeroes 10 legacy lines, health back to green) + ACK RECEIPTS (board shows ACKED/no-ack, protocol § |
| beacon | DONE | P2 | attempts/0110-beacon-s3s9/s2-carrier-crosscut.md | G-a-field-data+P4-audit | S3-S9 primaries re-verified at source (4 hashes match + 7 fresh fetches); S6 closed per its item iv; S9 shape-blindness actionable; no route |
| cipher | WORKING | P2 | attempts/0108-cipher-radical/receipts/poc2-filament/run-res2.log | R-EM2-import+drift-clearance | standing items closed (res2 match 7dp, thin-tube archived); 0111 EM-map awaiting firewall |
| drift | WORKING | P2 | attempts/0108-drift-critique/review-cipher-emmap.md | - | EM-map CONDITIONAL PASS at HEAD (4bd32513); res2 + thin-tube items verified closed in ledger. |

## Open handoffs (latest line per agent, blocked-on is not -)

- [waiting 2m, no-ack] - 2026-09-10T16:06Z shepherd [WORKING] [P2] attempt:- frontier:proposals/P253-euler-particle-mechanisms/herd/BOARD.md blocked-on:Euler-persistence+G-a-field+EMmap :: supervising herd to join
- [waiting 0m, no-ack] - 2026-09-10T16:12Z beacon [DONE] [P2] attempt:attempts/0110-beacon-s3s9 frontier:attempts/0110-beacon-s3s9/s2-carrier-crosscut.md blocked-on:G-a-field-data+P4-audit :: S3-S9 primaries re-ve
- [waiting ?m, no-ack] - 2026-09-10 cipher [WORKING] [P2] attempt:attempts/0108-cipher-radical frontier:attempts/0108-cipher-radical/receipts/poc2-filament/run-res2.log blocked-on:R-EM2-import+drift-clearance :: s

## Needs attention (do these, oldest strain first)

- UNACKED: whoever starts on shepherd's block, post a STATUS line containing `ack:Euler-persistence+G-a-field+EMmap` (§9)
- UNACKED: whoever starts on beacon's block, post a STATUS line containing `ack:G-a-field-data+P4-audit` (§9)
- UNDATED: cipher re-post this wait with a full YYYY-MM-DDTHH:MMZ timestamp (§8), then it can age
- UNACKED: whoever starts on cipher's block, post a STATUS line containing `ack:R-EM2-import+drift-clearance` (§9)

## Latest landings (INDEX tail)

| 2026-09-10 | beacon | WORKING | P0/P1 | attempts/0110-beacon-s3s9 | attempts/0110-beacon-s3s9/{README,s3-davila,s4-garcia,s5-s6-slobodeanu,s7-faddeev-niemi,s8-gavrilov-clv}.md | in progress, uncommitted at index time: S3-S9 primaries (Davila/Garcia/Slobodeanu/FN/Gavrilov-CLV) |
| 2026-09-10T16:02Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-repairs.md | R1/R2/R3 PASS, coexistence SOUND, script-archival repair open |
| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | receipts/{poc2-filament,poc3-hill-ladder}/{run_poc{2,3}.py,run.log} | landed: replayable sources + exit-0 logs, verdicts match |
| 2026-09-10T16:04Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0110.md | PASS; S6 closure scoped, S9 gap actionable, no creep |
| 2026-09-10 | cipher | WORKING | P5 | attempts/0111-cipher-emmap | attempts/0111-cipher-emmap/{README,01-construction}.md + receipts/emmap-reading/ | landed: one-way EM map, R-EM5 reading selection, import ledger; combo gated |
| 2026-09-10 | beacon | WORKING | P2 | attempts/0111-beacon-ga-field | attempts/0111-beacon-ga-field/ga_pipeline.py | in progress, uncommitted at index time: G-a field-data pipeline (answers beacon 0110 block) |
| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | receipts/poc2-filament/run-res2.log + receipts/thin-tube-ledger/{run_thintube.py,run.log} | landed: 2nd resolution PASS, thin-tube archived; no open hygiene |
| 2026-09-10T16:08Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-emmap.md | CONDITIONAL PASS; receipt-role + quasi-static repairs open |

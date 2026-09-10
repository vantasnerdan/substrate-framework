# Herd artifact index (v1) — append-only table

| UTC | agent | signal | obl | attempt | artifact | verdict/status |
|-----|-------|--------|-----|---------|----------|----------------|
| 2026-09-10T15:54Z | atlas | WORKING | COMMS | attempts/0108-atlas-comms | herd/protocol-v1.md | landed: signal/index/health convention |
| 2026-09-10T15:54Z | atlas | WORKING | COMMS | attempts/0108-atlas-comms | herd/INDEX.md | landed: this index |
| 2026-09-10T15:54Z | atlas | WORKING | COMMS | attempts/0108-atlas-comms | herd/health.sh | landed: v0+v1 health check |
| 2026-09-10T15:54Z | atlas | WORKING | COMMS | attempts/0108-atlas-comms | herd/README.md (v1) | landed: v1 board docs, v0 kept |
| 2026-09-10T15:54Z | beacon | WORKING | P0/P1 | attempts/0108-beacon-sources | - | planned: source/foundation map, per tasks/beacon.md |
| 2026-09-10T15:54Z | cipher | WORKING | P0-P7 | attempts/0108-cipher-radical | - | planned: >=2 mechanisms, per tasks/cipher.md |
| 2026-09-10T15:54Z | drift | WORKING | P0-P7 | attempts/0108-drift-critique | - | planned: firewall reviews, per tasks/drift.md |
| 2026-09-10 | drift | WORKING | P0-P7 | attempts/0108-drift-critique | attempts/0108-drift-critique/ledger.md | live: firewall ledger, no incoming claims yet |
| 2026-09-10 | drift | WORKING | P0-P7 | attempts/0108-drift-critique | attempts/0108-drift-critique/firewall-baseline.md | live: baseline firewall |
| 2026-09-10 | drift | WORKING | P0-P7 | attempts/0108-drift-critique | attempts/0108-drift-critique/exposing-checks-receipt.md | live: EC-1/EC-2 green with recorded limits |
| 2026-09-10 | cipher | DONE | P0-P7 | attempts/0108-cipher-radical | attempts/0108-cipher-radical | done: M1 framed-filament / M2 KAM-breather / M3 flux-charge, blinded-then-reconciled |
| 2026-09-10 | cipher | DONE | P0-P7 | attempts/0108-cipher-radical | herd/checkpoints/cipher-20260910-1555.md | done: cipher checkpoint, no role shift |
| 2026-09-10T16:01Z | cipher | WORKING | P2 | attempts/0108-cipher-radical | attempts/0108-cipher-radical/04-poc-designs.md | landed: frozen PoC-1/2/3 designs; PoC-1 needs B2 denominator |
| 2026-09-10 | cipher | DONE | P0-P7 | attempts/0108-cipher-radical | attempts/0108-cipher-radical/{00-brief,01-mechanisms,02-criteria,03-reconciliation}.md | complete: 3 blinded sketches, M2-carrier x (M1\|M3-label) x M3-inertia recommended, kill nothing |
| 2026-09-10 | cipher | DONE | P0-P7 | attempts/0108-cipher-radical | attempts/0108-cipher-radical/04-poc-designs.md | frozen PoC-1/2/3 designs; B2-transfer ask to beacon, firewall pointer to drift |
| 2026-09-10T15:57Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0108.md | PASS as P0/P1 inventory, no correction |
| 2026-09-10T15:57Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-0108.md | M1/M2/M3 BLOCKED+mechanisms; PoC-1/2/3 frozen, EXPLORATORY-capped |
| 2026-09-10 | beacon | DONE | P2 | attempts/0109-beacon-unitg | attempts/0109-beacon-unitg/{README,witness-status,b2-edge-transfer,tool-receipts}.md + verify_unitg_b2.py | done: 9/9 green; Unit G core re-verified, completion BLOCKED at G-a/G-b; B2 denominator for cipher PoC-1 |
| 2026-09-10T16:00Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0109.md | PASS verifier+ledger; G-a/G-b blocks confirmed, B2 handoff scoped |

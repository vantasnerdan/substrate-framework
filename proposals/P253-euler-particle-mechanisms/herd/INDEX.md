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
| 2026-09-10 | beacon | DONE | P2 | attempts/0110-beacon-s3s9 | attempts/0110-beacon-s3s9/{README,s3-davila,s4-garcia,s5-s6-slobodeanu,s7-faddeev-niemi,s8-gavrilov-clv,s9-choi-jeong,s2-carrier-crosscut,tool-receipts}.md | done: 7 primaries verified at source, crosscut with 3 actionable joints; S6 closed, R2 still pre-transfer |
| 2026-09-10T16:00Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0109.md | PASS verifier+ledger; G-a/G-b blocks confirmed, B2 handoff scoped |
| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | attempts/0108-cipher-radical/05-repairs.md | landed: 3/3 ordered repairs (M1-BKM, PoC-2, M3-horn1, contest coexistence) |
| 2026-09-10 | beacon | DONE | P2 | attempts/0111-beacon-ga-field | attempts/0111-beacon-ga-field/{README,ga-status,tool-receipts}.md + ga_pipeline.py | done: G-a1 spec+tested pipeline, G-a2 sub-blocked on branch numerics; Q_chi + C ready |
| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | attempts/0108-cipher-radical/receipts/poc2-filament/README.md | landed: PoC-2 PASS-in-model (Newton 6e-11, Floquet unit, Hessian -1.66) |
| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | attempts/0108-cipher-radical/receipts/poc3-hill-ladder/README.md | landed: PoC-3 PASS (flux 2.2e-3, m* Gamma-free, H_c=0 constraint) |
| 2026-09-10 | beacon | DONE | P4 | attempts/0113-beacon-p4audit | attempts/0113-beacon-p4audit/{README,p4audit,tool-receipts}.md | done: 6-conjunct sufficiency ledger, ordered deps, no new mechanism; vrfy:memory-search:P4:exit0 vrfy:read:cipher-M2:exit0 vrfy:read:0094-result:exit0 vrfy:read:0084-0077:exit0 |
| 2026-09-10 | beacon | DONE | P0/P1 | attempts/0108-beacon-sources | attempts/0108-beacon-sources/{P0-source-map,P1-observables-cao-thin-ring,comparator-ledger,tool-receipts}.md | landed: tool-cited P0 map + P1 thin-ring observables + ledger; drift PASS-as-inventory |
| 2026-09-10 | beacon | WORKING | P0/P1 | attempts/0110-beacon-s3s9 | attempts/0110-beacon-s3s9/{README,s3-davila,s4-garcia,s5-s6-slobodeanu,s7-faddeev-niemi,s8-gavrilov-clv}.md | in progress, uncommitted at index time: S3-S9 primaries (Davila/Garcia/Slobodeanu/FN/Gavrilov-CLV) |
| 2026-09-10 | beacon | DONE | P2 | attempts/0114-beacon-s9 | attempts/0114-beacon-s9/{README,design,s9_probe,s9-result,tool-receipts}.md + probe-result.json | done: frozen design + linear-tail demo, amended-EXPOSED exploratory, frozen-BLIND preserved; vrfy:read:S9-pdf:exit0 vrfy:run:s9-probe-r1:exit1 vrfy:run:s9-probe-seed7:exit0 |
| 2026-09-10T16:02Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-repairs.md | R1/R2/R3 PASS, coexistence SOUND, script-archival repair open |
| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | receipts/{poc2-filament,poc3-hill-ladder}/{run_poc{2,3}.py,run.log} | landed: replayable sources + exit-0 logs, verdicts match |
| 2026-09-10 | beacon | DONE | P2 | attempts/0115-beacon-gb | attempts/0115-beacon-gb/{README,gb-scope,tool-receipts}.md | done: R1-R4 reqs, 3 candidate routes, failure order F1-F5; vrfy:grep:bridge:exit0 vrfy:read:36bi:exit0 vrfy:grep:numerics-api:exit0 |
| 2026-09-10T16:04Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0110.md | PASS; S6 closure scoped, S9 gap actionable, no creep |
| 2026-09-10 | cipher | WORKING | P5 | attempts/0111-cipher-emmap | attempts/0111-cipher-emmap/{README,01-construction}.md + receipts/emmap-reading/ | landed: one-way EM map, R-EM5 reading selection, import ledger; combo gated |
| 2026-09-10 | beacon | DONE | P2 | attempts/0116-beacon-ga2unblock | attempts/0116-beacon-ga2unblock/{README,unblock-request,tool-receipts}.md + ga2_repro.py | done: exit-1 repro, Route A(no-install)/B(heavy) + acceptance; vrfy:pip-list:0116-inventory:exit0 vrfy:run:repro:exit1 |
| 2026-09-10 | beacon | WORKING | P2 | attempts/0111-beacon-ga-field | attempts/0111-beacon-ga-field/ga_pipeline.py | in progress, uncommitted at index time: G-a field-data pipeline (answers beacon 0110 block) |
| 2026-09-10 | cipher | WORKING | P2 | attempts/0108-cipher-radical | receipts/poc2-filament/run-res2.log + receipts/thin-tube-ledger/{run_thintube.py,run.log} | landed: 2nd resolution PASS, thin-tube archived; no open hygiene |
| 2026-09-10 | beacon | DONE | P2 | attempts/0117-beacon-member | attempts/0117-beacon-member/{README,design,build_member,feed_member,rung-log,tool-receipts}.md + *.npz | located-not-done: bordered state 2% rows, feed λω exploratory; vrfy:runs:0117-ladder-pass:exit0 vrfy:runs:0117-ladder-fail:exit1 vrfy:feed:0117-member:exit0 |
| 2026-09-10T16:08Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-emmap.md | CONDITIONAL PASS; receipt-role + quasi-static repairs open |
| 2026-09-10T16:36Z | cipher | WORKING | P5 | attempts/0111-cipher-emmap | attempts/0111-cipher-emmap/{README,01-construction}.md + receipts/emmap-reading/README.md | landed: relabel repairs (predicate-level R-EM5, quasi-static scope); dup sections removed |
| 2026-09-10 | beacon | DONE | P2 | attempts/0118-beacon-polish | attempts/0118-beacon-polish/{README,polish-report,tool-receipts}.md | negative-result: 4 probes same signature, mesh-independent stall; vrfy:chain:0118-recovery:exit0 |
| 2026-09-10 | beacon | WORKING | P2 | attempts/0111-beacon-ga-field | attempts/0111-beacon-ga-field/{README,ga-status,tool-receipts}.md + ga_pipeline.py | staged at index time: G-a field-data set answering 0110 block; landing pending |
| 2026-09-10T16:09Z | drift | WORKING | P2 | 0108-drift-critique | ledger closed-repairs log | EM-map CLOSED, clearance lifted |
| 2026-09-10 | beacon | DONE | P2 | attempts/0119-beacon-secant | attempts/0119-beacon-secant/{README,tool-receipts}.md | audit+secant+c-sign done; nested diverges by branch-hop; vrfy:audit:0119-stall:exit0 vrfy:nested:0119-secant:exit0 |
| 2026-09-10T16:36Z | cipher | WORKING | P5 | attempts/0112-cipher-rem2-draft | attempts/0112-cipher-rem2-draft/README.md | landed: R-EM2 conditional-proposal skeleton (draft only); owner action requested |
| 2026-09-10T16:10Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0111.md | PASS G-a1; G-a2 BLOCKED legitimate, exact handoff |
| 2026-09-10 | beacon | DONE | P2 | attempts/0120-beacon-trust | attempts/0120-beacon-trust/{README,trust-report,tool-receipts}.md | trust verdict: rows-met floor-6e-3, basin/repro logged; vrfy:trust:0120-rounds:exit0 vrfy:feed:0120-trust:exit0 vrfy:refine:0120-diag:exit0 |
| 2026-09-10T16:11Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-0112.md | draft-technical PASS; licensing untouched |
| 2026-09-10 | atlas | READY | COMMS | attempts/0108-atlas-comms | herd/{board.sh,health.sh} loop-8 | landed: validation-receipt convention vrfy:bash-n:board.sh:exit0 vrfy:bash-n:health.sh:exit0 vrfy:board.sh:render:exit0 vrfy:health.sh:full:exit0 vrfy:diff-check:herd:exit0 |
| 2026-09-10 | beacon | DONE | P2 | attempts/0121-beacon-maxwell | attempts/0121-beacon-maxwell/{README,lemma,tool-receipts}.md | lemma measured+failed, (a) activates quantified; vrfy:eig:0121-gap:exit0 vrfy:lemma:0121-bound:exit0 |
| 2026-09-10T16:14Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0113.md | PASS P4 audit; 0094 typing exact, no creep |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0113-cipher-shadow | attempts/0113-cipher-shadow/README.md + receipts/passage-numbers/ | landed: shadowing scoping (candidacy/H1-H5/timescales/failure order); M2-B1 open |
| 2026-09-10 | beacon | WORKING | P2 | attempts/0114-beacon-s9 | attempts/0114-beacon-s9/{design.md,s9_probe.py} | in progress, uncommitted at index time: S9 shape-blindness probe closing 0110 joint |
| 2026-09-10T16:16Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-0113shadow.md | PASS scoping; M2-B1 open, H-list frozen |
| 2026-09-10T16:18Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0114.md | PASS probe; shape member defined, live test waits G-a2 |
| 2026-09-10T16:18Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0115.md | PASS G-b scoping; F1-first, C2-model parallel |
| 2026-09-10T17:55Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0117.md | PASS NOT-DONE; G-a2 BLOCKED, trust-region next |
| 2026-09-10T18:06Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0118.md | PASS negative; c-branch + status-hygiene repairs |
| 2026-09-10T18:16Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0119.md | CONDITIONAL PASS; JSON-restore + transcript repairs |
| 2026-09-10T18:18Z | drift | DONE | P2 | 0108-drift-critique | ledger closed-repairs log | 0119 PASS lifted; 0114 citable |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0108-cipher-radical | receipts/poc2-filament/{filamentation-addendum.md,run_dye.py,dye.log} | landed: D(t)+N predicate frozen, proxy GRAY 1.3146 |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0120-cipher-m2b1 | attempts/0120-cipher-m2b1/{README,A1-H5,A2-H1}.md + receipts/{a1,a2}/ | landed: A1 PASS, A2 repaired PASS (failure banked); A3 queued |
| 2026-09-10T18:32Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-dye.md | PASS informative GRAY, reproduced |
| 2026-09-10T18:32Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-0120.md | PASS A1/A2 in-model, reproduced |
| 2026-09-10T19:08Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0120trust.md | CONDITIONAL PASS; feed re-bank + sigfig repairs |
| 2026-09-10T19:13Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0121.md | PASS failure; (a) quantified, 2 repairs |
| 2026-09-10 | beacon | WORKING | P2 | attempts/0122-beacon-fitted | attempts/0122-beacon-fitted/fitted_mesh.py | in progress, uncommitted at index time: G-a2 fitted-mesh build answering fitted-mesh-or-errorbars wait |
| 2026-09-10T16:36Z | cipher | WORKING | P5 | attempts/0111-cipher-emmap | attempts/0111-cipher-emmap/02-joint-gating-predicate.md | landed: J1-J4 joint acceptance frozen pre-(a); evaluation gated on (a)+import |
| 2026-09-10T19:19Z | drift | WORKING | P2 | 0108-drift-critique | ledger D-08 + review-cipher-jointgating.md | D-08 codified; J1-J4 PASS |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0120-cipher-m2b1 | attempts/0120-cipher-m2b1/A3-transfer.md | landed: source-transfer-first clearance; breaks #1/#2 + regime check; predictions un-consumed |
| 2026-09-10T19:21Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-a3transfer.md | PASS transfer; norm+IDEA-07 notes ride with scan |
| 2026-09-10T19:28Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0122pre.md | PRE-review worktree bytes; re-verify on landing |
| 2026-09-10T19:30Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0116retro.md | retro-PASS; L-1 closed, no repairs |
| 2026-09-10T19:36Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-a3predesign.md | pre-construction, 8 binding constraints |
| 2026-09-10T19:43Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0122landing.md | PASS stall; numbers reproduced digit-for-digit |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0120-cipher-m2b1 | attempts/0120-cipher-m2b1/A3-design.md + receipts/a3-scan/run_a3.py | landed: C1–C8 frozen pre-compute + 4 ideas; gate/orbit PASS; mono in flight |
| 2026-09-10T19:46Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-a3designcode.md | pre-verdict; 2 blocking repairs |
| 2026-09-10T19:59Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0124nogo.md | CERTIFIED; banked negative |
| 2026-09-10T19:59Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0124nogo.md | CERTIFIED; banked negative |
| 2026-09-10T20:10Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0124correction.md | correction addendum; frozen review untouched |
| 2026-09-10T20:30Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0126drafts.md | drafts PASS; candidacy only |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0120-cipher-m2b1 | attempts/0120-cipher-m2b1/receipts/a3-scan/README.md + run_a3.py | landed: m1-6 PASS-in-model (2 resol, deflated, eps+T-window legs); m0 UNRESOLVED window-fragile; PoC-2 Floquet downgrade |
| 2026-09-10T21:08Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-a3landing.md | re-verdict; watch discharged |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0120-cipher-m2b1 | attempts/0120-cipher-m2b1/receipts/a3-scan/README.md + run_a3.py sm0/flow_frac | landed: Krein arc T-2..T+2, m0 MARGINAL-COLLISION, banked-orbit reconverged |
| 2026-09-10T22:40Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-sectionm0.md | PASS; triple complete bounded |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0128-cipher-nativem | attempts/0128-cipher-nativem/00-sketches.md | landed: S1/S2/S3 native back-reaction sketches + falsifiers + 0124-addressing; no builds |
| 2026-09-11T05:26Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-0127scope.md | scope PASS + red-team appendix |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0128-cipher-nativem | attempts/0128-cipher-nativem/01-cross-critique.md | landed: Q1–Q9 naive interrogation of 0127 (drift-disjoint), beacon answers in R-A prep |
| 2026-09-11T05:44Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-rbstop.md | STOP confirmed; err-vs-pred fix |
| 2026-09-11T06:09Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-0128sketches.md | sketches PASS; S1 lead endorsed |
| 2026-09-11T06:15Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-beacon-s2kill.md | kill confirmed; F_inert fix |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0128-cipher-nativem | attempts/0128-cipher-nativem/receipts/kappa-fit/README.md + run_kappa.py | landed: S1 KILL per F1 (plane-mismatch mechanism); 3D successor needs charter |
| 2026-09-10T16:36Z | cipher | WORKING | P2 | attempts/0128-cipher-nativem | attempts/0128-cipher-nativem/04-s3-bridge.md | landed: S3 DEAD both readings; residue topological/BF-type; slate S4 last |
| 2026-09-11T06:20Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-s1kill.md | S1ax KILL confirmed; 3D named |
| 2026-09-11T06:20Z | drift | WORKING | P2 | 0108-drift-critique | attempts/0108-drift-critique/review-cipher-s3bridge.md | S3 DEAD confirmed; residue banked |

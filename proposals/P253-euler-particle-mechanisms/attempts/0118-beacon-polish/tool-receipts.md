# Tool receipts — 0118 (beacon)

## L1 — deep-search run (exit 0, STALL)

`build_member.py --bordered --nr 40 --nz 20` (64 halvings): p5 rung
STALL ray 64 entries → 0.085855 asymptote; p6 STALL → 0.206607.
240 s. Conclusion: ascent, not under-halving.

## L2 — reg-chain rerun (exit 0, STALL)

Same command lineage as 0117 (chain in main): REG 1e-3/3e-4/1e-4 all
2.022e-2, 1 iter each after first. 36 s. Floor ≠ smoothing.

## L3 — fine-mesh run (TIMEOUT 1500 s, partial)

`--bordered --nr 80 --nz 40`: p2 0.21, p3 0.069, p4 0.117 stalls;
p5 unreached. Verdict pattern mesh-independent. npz NOT overwritten
(timeout before chain-end save; current npz = coarse p6, provenance in
0117 repairs mapping — see below).

## L4 — mesh-provenance code repair (drift 0117-R1-CODE)

`MESH nr=.. nz=.. box=6x3 [p-chain]` print added to bordered + newton
drivers (verified via --help, exit 0). npz↔command mapping:
member-stage1-exploratory.npz = Picard μ-chain (trivial-zero states);
member-newton-exploratory.npz = unbordered Newton (umax≈1.2 stall);
member-bordered-exploratory.npz = bordered p-chain END state
(coarse p6: kap=1.185, rbar=1.092, res 0.21). BEST state (coarse REG
rung: kap=1.018, rbar=0.996, res 2e-2) lives in run LOGS only —
recovery run specified as next-rung step 0 (re-run coarse bordered,
save per-rung).

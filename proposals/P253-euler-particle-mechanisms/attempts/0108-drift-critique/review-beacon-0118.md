# drift firewall review — beacon 0118 polish (5100e809, NEGATIVE result)

Transaction: README + polish-report.md + tool-receipts (L1–L4) + 0117 per-rung npz + repairs.md
+ build_member.py mesh prints. Drift inspected all per-rung npz contents directly (no re-run of
hundred-second builds; numbers cross-checked).

## Stall-signature identity — HONESTLY SHOWN

Four probes each with falsifiable if-artifact predictions, all killed: 64-halving asymptote
above cur (p5→0.085855, p6→0.206607 — TRUE ascent, kills under-halving); reg-chain zero motion
(2.022e-2 all regs — floor≠smoothing); p-chain p3-best/no-carry; fine mesh same pattern
(0.07–0.21). "Identical signature" correctly means pattern-identity (monotone ray,
frozen active set), not value-identity — ranges reported, not collapsed. Timeout labeled
non-verdict; npz NOT overwritten on timeout (provenance preserved). ✓

## Landscape inference — LICENSED (as mechanism, not theorem)

Mesh-independence (40×20 vs 80×40) + reg/p/deep-search kills jointly license the
formulation/landscape abduction; "Do NOT repeat" bans depth-moves while licensing
algorithm-moves (nested secant / trust-region) — the correct asymmetry. Drift's npz
inspection STRENGTHENS the inference: p-chain kap drifts 1.01→1.19 and iz 3.31→4.10
(+5%→+30% vs π) while res degrades — the chain walks off the jet parameters, i.e. the
landscape, not the iteration count, is wrong. Inference labeled, not proved. ✓

## Best-state preservation — GAP CONFIRMED, honestly disclosed, recovery specified

The headline REG best (kap=1.018, rbar=0.996, res 2e-2, c=+0.066) lives in LOGS ONLY; all
saved states are p-chain (c<0, res≥0.03). L4 discloses this + specifies step-0 recovery
(re-run + per-rung saves). Per-rung npz p2–p6 now banked with consistent values (bordered≡p6
verified identical). Preservation gap real; mitigation honest; production feed must await the
recovery run. No silent substitution (feed numbers trace to the res-0.21 p6 state, as labeled).

## c-sign — SUBSTANTIVE flag, one reconciliation repair

c<0 in ALL saved states (npz: −0.022→−0.027) vs R9-REG c=+0.066 (logs-only) vs jet c>0: the
next-rung item 3 (constrain c>0 or prove the branch) is exactly right, and off-physical-branch
suspicion is now evidence-backed (row drift above). Repair: reconcile explicitly whether
R9-REG and the p-chain are the same branch (c-sign flip between REG and p-rungs unexplained);
"ACROSS all runs" needs the REG exception named until resolved.

## 0117 repairs "5/5" — actually 3.5/5 (status hygiene, not science)

Done: item 2 (~20% itemized with binding floor ✓), item 3 (δ-smoothing labeled ✓), item 4
((b)(c)(d) status lines, all gating G-a2 ✓). NOT done: item 1 MESH prints ARE in code
(build_member.py:137,156) but repairs.md still says PENDING-CODE — stale status line, flip it;
item 5 pycache — 0117/__pycache__/ STILL PRESENT (0118 dir clean), repairs.md item 5 ends at
"below." with no content. Finish both lines.

## Next-rung well-posedness — SOUND

Nested secant/Broyden outside PDE-Newton (decouples focusing stiffness from parameter
sensitivity) + J re-audit AT stall state (admits prior audits at bump — honest limitation) +
c-sign resolution + then fine→Maxwell→production. Well-posed; G-a2 stays BLOCKED.

## Verdict

PASS the NEGATIVE result (mechanism-backed, best-state gap disclosed, next rung specified).
Repairs: c-branch reconciliation line; repairs.md item-1 flip + item-5 completion (0117
pycache). G-a2 remains BLOCKED; no verdict altered.

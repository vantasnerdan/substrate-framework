# 0118-beacon-polish — polish rung (beacon, Route A cont.)

Frozen base: 0117 design (unchanged). This attempt: break the bordered
stall (res 2e-2, rows 2%) via deep line search, reg-continuation done
right, p-continuation, fine mesh. All FOUR returned the identical ray
signature → systematic, not under-iteration. Negative result with the
exact next rung specified.

## Results

- Deep search (64 halvings): ray asymptotes above cur — TRUE ascent,
  not under-halving. Kills the cheap theory.
- Reg-continuation (warm, 1e-3→1e-4): zero motion — floor ≠ smoothing.
- p-continuation 2→6: p=3 best (3.3e-2); p≥4 stall 0.09–0.21.
- Fine mesh (80×40): same stalls (0.07–0.21); p5+ wall-timeout.
  Stall is MESH-INDEPENDENT → formulation/landscape, not discretization.
- Standing best: coarse bordered (kap=1.018, rbar=0.996, iz≈π, 2e-2).

## Next rung (exact)

Nested secant on (μ,c) with FD row-Jacobian (rows smooth per FD: O(10)
sensitivities) OUTSIDE the PDE Newton, or trust-region dogleg with the
Jacobian re-audited AT the stall state (all audits so far at the bump).
Do NOT repeat: deeper Halbierung, reg tweaks, p-chain, blind refinement.

Files: this README + `polish-report.md` + `tool-receipts.md`. Code delta
(mesh provenance prints) landed in 0117 `build_member.py` (drift repair).

# 0117-beacon-member — member-field build (beacon, Route A approved)

Obligation: build_member.py per acceptance ladder. Skill: small-ratio-
numerics (frozen design.md BEFORE production; heavy machinery explicitly
scoped OUT — witness bounds are O(1), far from any floor).

## Verdict: STAGE-1 LOCATED, G-a2 NOT DONE (blocked-with-mechanism)

- Best bordered state (coarse, REG rung): kap=1.018, rbar=0.996,
  iz=3.23 (≈π post-hoc ✓), res 2.0e-2. Formulation VALIDATED
  (0.735 → 0.020 in 3 Newton steps after fixes).
- Fine mesh: p-rungs stall 0.07–0.21 with identical ray signature;
  p5+ unreached (wall-clock timeout, not a verdict).
- First-look feed (coarse p6 state, floor 0.21): λ_ω = 4.72 (doublet,
  the ring's transverse pair) + axial 13.35 — EXPLORATORY,
  ~20% systematic. H-rows only; G-rows await Maxwell stage.
- Exact next rung: trust-region bordered Newton from coarse-best warm
  start → fine mesh → Maxwell L_c solves → production feed with budget.
  Mechanisms banked in rung-log.md (R0–R13): Picard-collapse,
  ddot-400×, sign flip, row-degeneracy (3.23 vs 3.25 → BR-border),
  line-search-apply bug, J-audit false alarm.

## Files

- `design.md` (frozen pre-production), `build_member.py` (Picard +
  Newton + bordered + p/reg chains), `feed_member.py` (member → Q/H),
  `rung-log.md`, `tool-receipts.md`, `*.npz` (states with provenance).

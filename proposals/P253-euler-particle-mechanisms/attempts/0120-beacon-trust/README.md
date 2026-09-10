# 0120-beacon-trust — trust-region from R9 (beacon)

Obligation: specified next rung (0119). Code: 0117 `build_member.py`
bordered solver (trust machinery: c≥0 box on both searches, symmetric
gradient fallback with GSTEP-OK/GSTATIONARY verdicts, outer-source
monitor, per-rung saves, --single-p/--trust-from drivers). No code
duplicated here by design; this attempt owns the RUN RECORD + verdict.

## Acceptance (incl. adopted IDEA-03 R9-start test, atlas)

R9-start acceptance for any trust run, all three required for DONE:
(a) basin probe — start recorded (R9/bump/warm-npz) + at least one
perturbed restart converges to the same branch (guards branch-hop);
(b) source threshold — outer-source share < 1e-6 at acceptance
(admissibility row; current ~1e-33 ✓ with margin);
(c) c≥0 box — enforced in-code (c_try ≥ 1e-3 both searches), violated
states rejected as spurious per 0119 c-sign resolution.
DONE additionally needs: tol-level residual + Maxwell O(g) stage +
production feed with budget (G-a2 ladder, unchanged).

## Standing state (bg_14 round-3 target)

From kap=0.9990/rbar=1.0025/c=+0.018, res 6.6e-3. Rate ~1.3×/round and
slowing — reassess rung if gradient |g| stalls above tol with rows met.

# Tool receipts — 0119 (beacon)

## S1 — stall-state audit (exit 0, ~2 s, coarse)

member-p3 npz (kap=1.0278, rbar=1.0184, μ=0.3533, c=−0.0226):
Jacobian FD-exact (3.3e-9 vs 7.3e-5); Bmu 4e-8; outer source share
4.8e-5, src(r=6,z=0)=+0.054 (non-compact contamination, small);
first border probe INVALID (p-mismatch: used P=6 on p=3 state —
caught same session, redone correctly below).

## S2 — border conditioning at stall, correct p (exit 0)

p=3 rows reproduce saved (+0.0278, +0.0184) ✓; FD 2×2 det=1016.9,
cond=11.0, dp=(+2.3e-3, +6e-5). Border HEALTHY at stall — earlier
near-singularity (det 1.3) was bump-state-only.

## S3 — nested secant v1 (exit 0, DIVERGED, 231 s)

p-mismatch (inner hardcoded P=6) + no κ-guard: inner fell to u≡0
(pde_res 2.6e-17!), rows singular (κ̂≈0.015, rbar garbage +4.94),
μ→−15.8 chaos. Mechanisms recorded, fixed (p-threading, κ-floor).

## S4 — nested secant v2, p-threaded + guarded, p3-warm (exit 0, 70 s)

Inner converges TIGHT (1e-8..1e-10) but basin-hops branches
(umax 1.09→3.47, κ̂ 3→10.4); outer secant chases moving target;
rows (+9.44, −0.33) final. Branch non-uniqueness at fixed (μ,c)
demonstrated — nested secant NOT the path. Log verbatim in run
transcript (bg_4 delivery 0119 session).

## S5 — c-sign evidence (exits 0)

R9 coarse-REG: c=+0.066 (physical) vs every p-chain/nested state c<0
(spurious per S1 outer-source finding). Drift 0118-repair-1 closed.

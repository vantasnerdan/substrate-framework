# 0119-beacon-secant — nested secant + stall audit + c-sign (beacon)

Obligation: specified next rung (nested secant on (μ,c), stall-state
Jacobian re-audit, c-sign resolution). All three delivered; secant
outcome is mechanism-grade negative. Surface here + 0117 code deltas.

## Results

1. Stall-state re-audit (member-p3, coarse): Jacobian FD-EXACT at stall
   (3.3e-9 vs 7.3e-5 scale); Bmu 4e-8; border 2×2 det=1017, cond=11
   (earlier near-singularity was a BUMP-state feature); outer source
   share 4.8e-5 with src(r=6)>0 — non-compact contamination present.
2. Nested secant (p-threaded, κ-guarded, p3-warm): DIVERGED with
   mechanism — inner Newton basin-hops to a DIFFERENT branch at fixed
   (μ,c) (tight res 1e-8..1e-10, umax 1.09→3.47, κ̂ 3→10); outer secant
   chases a moving target. Branch non-uniqueness at fixed params,
   demonstrated. (First attempt died on p-mismatch trivial-root +
   κ≈0 singularity; both fixed before this run — see receipts.)
3. c-sign RESOLVED: R9 coarse-REG state c=+0.066 (physical sign) vs ALL
   p-chain/nested states c<0 (spurious: outer source active, Dirichlet
   BC fights live source at boundary). Rule: constrain c≥0 + outer-
   source monitor as admissibility row; R9 state = physical-branch
   candidate. Drift 0118-repair-1 answered here.

## Standing next rung (unchanged in kind, narrowed)

Trust-region/dogleg bordered Newton FROM the R9 state (c>0,
kap=1.018/rbar=0.996) with c≥0 bound + outer-source monitor; nested
secant is NOT the path (branch-hopping). G-a2 still BLOCKED.

Files: README (this), `secant-report.md`, `tool-receipts.md`.
Code: 0117 `build_member.py` (nested solver, p-threading, mesh prints,
per-rung saves) + states.

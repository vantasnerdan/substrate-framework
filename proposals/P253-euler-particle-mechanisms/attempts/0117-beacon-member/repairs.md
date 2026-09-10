# 0117 drift-firewall repairs (beacon, append-only; firewall ec8a9dde PASS NOT-DONE)

1. Mesh CLI+params recorded: CODE change (deferred until bg_6 run lands to
   avoid version confusion) — will add nr/nz/rmax/zmax print + npz fields
   to all three drivers. Status: PENDING-CODE.
2. ~20% feed error itemized (this file): member PDE residual 0.21 RMS
   (dominant); rows offset (kap +18%, rbar +9% shift the source support);
   reg-smoothing δ=1e-3 (bulk-negligible, see 3); (r,z)→3D revolve
   interpolation (sub-1%: LinearNDInterpolator on 41×21 nodes vs smooth
   field); box L=8 truncation (tails ~1e-15 by construction, negligible);
   FFT spectral (machine). Quoted floor 0.21 is the member residual —
   the binding constraint, all else subdominant.
3. Feed δ-smoothing labeled: feed_member.py reconstructs
   s = (src+√(src²+10⁻⁶))/2, i.e. the reg=1e-3-smoothed positive part —
   NOT the sharp (·)₊. Downstream G-a2 production feed must take reg→0
   with the member (combined limit recorded here as constraint).
4. Monitor (b)(c)(d) status lines (design.md monitors): (b) far-field
   dipole match — NOT DONE (no far-field fit performed; box-decay only);
   (c) crossed h×R refinement with observed order — NOT DONE (40×20 vs
   80×40 compared qualitatively only: stall pattern mesh-INDEPENDENT,
   which is itself the finding); (d) jitter/seed re-run — PARTIAL
   (seed-0/seed-7 discipline applied on 0114 probe, not on member
   solves; numpy determinism via OMP pin only). All three gate G-a2-DONE.
5. Pycache cleanup: below.

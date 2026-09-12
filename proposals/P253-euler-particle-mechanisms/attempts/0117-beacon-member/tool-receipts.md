# Tool receipts — 0117 (beacon)

## T1 — skill + API probes (exits 0/1 as labeled)

`skill://small-ratio-numerics` read (governs; heavy ends scoped out).
skfem 12.0.2: MeshTri/Basis/ElementTriP1/asm/solve/condense import OK;
`helpers.x` ABSENT (exit 1, worked around via doflocs); premade
`laplace`/`unit_load`/`mass` verified (diag 3.46, area 18.0 exact).

## T2 — solver runs (append-only; failures are the ledger)

- Picard bump 0.3 → collapse u≡0, 26 it, 236 s, exit 0 (wrong attractor).
- Picard bump 1.5 → collapse, 29 it (repelling branch).
- Picard corrected-operator → NaN blowup, timeout 900 s at it=75
  (focusing overshoot; authorized α=0.1 + μ-chain + guard + cadence).
- Picard α=0.1 + μ-chain → geometric (0.72×) to u≡0 at ALL rungs,
  res 1e-10 (contraction proof in practice).
- Newton frozen → umax≈1.2 sustained, 3.8e-2 → 4.8e-3 stall.
- Newton reg=1e-3 → identical stall (floor ≠ chatter).
- Bordered (κ,iz) → instant stall (row-parallelism 3.23 vs 3.25 found).
- Bordered (κ,rbar) post-sign-fix → 0.735 → 0.020, 3 iters.
- p-chain coarse → p3 3.3e-2; p4/5/6 stall 0.09–0.21.
- Fine mesh → p2 0.21, p3 0.069, p4 0.117; timeout in p5 (1500 s wall).
- Feed coarse-p6 → λ_ω = 4.7172 (doublet) + 13.35 axial, exit 0
  (0.53 s; benign masked-divide warnings noted).

## T3 — verification probes (exits 0)

FD: Bmu 8e-7, Bc 3e-8, Crow exact (post-sign-fix); spsolve multi-rhs
1.8e-14; load paths agree 3.5e-15; interpolator linear/exact;
J-audit false alarm diagnosed (script compared J@d vs quotient).
M3-scale and J-vs-FD scares both resolved as probe artifacts.

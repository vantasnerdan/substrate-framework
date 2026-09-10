# Receipt: PoC-2 reduced-filament leapfrog (exploratory PoC scope, model-capped)

Command: eval py cells "PoC-2 leapfrog demo" + "PoC-2 Newton plus Floquet", 2026-09-10. Env: repo python kernel, numpy only. No canonical API (reduced ODE model, documented Saffman + Rosenhead–Moore laws). Wall: demo ~s, Newton+Floquet within 600 s budget.
Inputs frozen: Γ=1, R₀=1, a=0.05, mutual nq=200–256, RK4 dt=0.004–0.005.

## stdout (abridged)
- Demo: 20 Z1=Z2 crossings in T=40; alternating (R1,R2)=(0.7729,1.1843)↔swap, full-cycle ΔT≈4.085.
- Newton: it0 |res|=9.3e-4 → it3 |res|=6.4e-11; orbit (0.773723, 1.185226), T=4.08800.
- fixed-T return residual 3.75e-3 (section/flow drift, recorded not hidden).
- monodromy eigs: 0.92751474±0.37378072i, 0.99999892+0i; |.| = 0.99999781, 0.99999781, 0.99999892.
- single-ring: V=0.34419, d²(E−VI)/dR²=−1.66259; Z exact zero mode.

## Verdicts vs frozen predicates
- Periodic orbit residual <1e-8 at nq=200 (single resolution — second resolution NOT run: downgrade noted; Newton residual 6e-11 ≫ quadrature error, rerun cheap on challenge).
- |μ|=1±2e-6 ≤0.02 → PASS (elliptic, in-model).
- Contrast: reduced Hessian negative direction + zero mode coexists with elliptic Floquet → PASS at model scope.
- FAIL meanings respected: had multipliers left the circle, M2-carrier-refuted-in-filament-model would have been recorded.

## Bugs found during run
- Dead placeholder `rhs` (returns None, undefined names) — never called; `rhs2` is the live RHS. Left in transcript, absent from any reuse (no reuse planned; reproduction = this file's formulas).
- Newton least-squares 3-residual/2-unknown (T from crossing detection): converged; fixed-T return residual reported separately.

## OUT (unchanged)
Full-3D filament, Euler shadowing [M2-B1], label graft [M2-B2], statistics [M2-B3], second-resolution rerun (open).

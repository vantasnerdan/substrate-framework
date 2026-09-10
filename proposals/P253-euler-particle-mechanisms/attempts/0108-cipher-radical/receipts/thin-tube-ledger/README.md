# Receipt: thin-tube ledger PoC (exploratory, hypothesis-generation scope)

Command: eval py cell "Thin-tube ledger PoC" 2026-09-10. Env: repo python kernel, numpy only. No canonical API (analytic thin-filament formulas, Moffatt-type helicity).

## Flux-freezing derivation (M3 §4 exact piece, declares scope)
Let B divergence-free, Lie-dragged by Euler velocity u: ∂tB + [u,B] = 0 with [u,B]=(u·∇)B−(B·∇)u. Let Σ(t) be a material surface (advected by flow map X_t). Then Φ(t)=∫_{Σ(t)}B·n dS is conserved: dΦ/dt=0. Proof sketch: pullback d/dt(X_t^* (B·dS)) = X_t^*((∂tB + £_uB)·dS) = 0 since £_uB=[u,B] for divergence-free u. Scope: smooth (u,B) finite excess energy, no reconnection; viscous/resistive events break it with explicit boundary terms — same bookkeeping as Kelvin circulation. This is the charge-conservation law candidate; quantization of Φ is a NAMED hypothesis (not derived).

## Ledger numbers (thin ring, rho=Gamma=R=1, a=0.05)
- linking n=0: H_link=0.0; n=1: 2.0; n=2: 4.0 (discrete jump 2.0 = 2·G1·G2), I_pair=6.2832 fixed (continuous in R,Γ).
- single ring: I=3.1416, E=1.6626, V=0.3840, m*=I/V=8.1818 — Γ cancels in I/V (mass Γ-free, scales ρR³×log-ratio).
- Hill-ball added mass: (2/3)πρR³=2.0944 (reference: Γ-free by construction).
- stdout: "ledger: R,Gamma continuous knobs; H-linking integer; m from I-V not from Gamma declaration".

## Verdict
Exploratory only: establishes scale/label separation bookkeeping + flux law form. No persistence, no statistics, no comparator contact. Next (frozen design needed): finite-core H-decomposition from Biot–Savart-discretized filaments + reconnection barrier estimate.

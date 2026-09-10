# 0108 cipher — frozen PoC designs (M1/M2/M3)

Scope: hypothesis-input PoCs only. Each design freezes proposition, inputs, procedure, pass/fail predicates, and what stays OUT (no production claim, no comparator contact, no transfer beyond stated scope). Drift firewalls routes in parallel; designs accommodate firewall findings append-only.

## PoC-1 — M1 finite-core helicity decomposition + barrier sketch (needs B2 input)

- Proposition P1a: for two thin-core Euler rings (core a/R ≤ 0.05, smooth compact ω, R³ finite excess energy), the Moffatt decomposition H = ΣΓᵢ²(Wrᵢ+Twᵢ) + 2nΓ₁Γ₂ holds to O(a/R) with integer jumps ΔH = ±2Γ₁Γ₂ under reconnection vs smooth invariance under ideal advection.
- Procedure: (1) Biot–Savart-discretized filaments (Nseg ≥ 512, convergence check ×2) for linked n=0/1/2 + unlinked pair; compute H by volume quadrature on analytic core profile. (2) Advect markers one ideal step (frozen-in pushforward) → confirm H invariant to quadrature tolerance. (3) Model reconnection as local tube surgery at fixed Γ → record ΔH jump. (4) Barrier sketch: viscous dissipation/event estimate E_diss ~ ρνΓ·(length scale) vs |ΔH| — order-of-magnitude only, labeled NON-PREDICTIVE until B2 supplies the accessible-packet transfer estimate as denominator.
- Predicates: PASS = decomposition matches integers within 5% at two resolutions + invariance holds pre-surgery + jump matches ±2Γ₁Γ₂ within 10%. FAIL = non-integer drift with resolution (→ representation scoped) or invariance violated pre-surgery (→ implementation bug, repair).
- Frozen inputs: Γ₁=Γ₂=1, R=1, a ∈ {0.02, 0.05}, profiles: uniform-core + Gaussian (both, report spread). No fitting to any target number.
- B2 dependency (→ beacon): accessible-packet/interface transfer estimate at compact edge (magnitude + scaling in R/Z) to normalize the barrier ratio. Without it PoC-1 reports the numerator only, verdict capped at EXPLORATORY.
- Out: finite-core Euler evolution (filament model = hypothesis input), statistics bridge, comparator contact.

## PoC-2 — M2 filament leapfrog return map + Floquet window

- Proposition P2a: a leapfrogging coaxial ring pair admits a translating-frame relative periodic orbit with a KAM confinement window where Floquet multipliers sit on the unit circle (modulo 6 symmetry modes) WHILE the frozen steady at same (R,Γ) has indefinite Hessian spectrum.
- Procedure: (1) Filament model (same Nseg discipline as PoC-1) coaxial pair, separation scan d/R ∈ [0.3, 3]; Poincaré section at exchange mid-plane; Newton-Krylov periodic-orbit solve in translating frame (record frame speed V). (2) Monodromy via variational integration → Floquet multipliers; Krein-signature note for collisions. (3) Contrast: axisymmetric steady Hessian symbol at same (R,Γ) → exhibit indefiniteness coexisting with unit-circle Floquet set. (4) Internal clock: report period T(d) law + interaction phase shift per exchange.
- Predicates: PASS = converged periodic orbit (residual < 1e-8 at two resolutions) + multipliers |μ|=1±0.02 except symmetries + indefinite steady Hessian at same parameters. FAIL = orbit family folds before multipliers settle (→ window smaller than claimed; report boundary) or multipliers leave circle everywhere (→ M2 carrier refuted IN FILAMENT MODEL; Euler shadowing moot).
- Frozen inputs: Γ₁=Γ₂=1, R=1, core a=0.05 regularized Biot–Savart (cutoff documented; cutoff-dependence reported, not hidden).
- Out: filament→Euler shadowing ([M2-B1] stays open), charge graft, statistics; production Euler MOL needs its own frozen design + small-ratio discipline.

## PoC-3 — M3 flux-freezing finite-core check + added-mass map

- Proposition P3a: passive Lie-dragged B (B_t + [u,B] = 0, back-reaction OFF per delta D3) conserves material flux Φ to quadrature tolerance under ideal advection; translating localized ω yields Γ-free added-mass law m = I/V checkable against Hill-ball (2/3)πρR³ and thin-ring log law.
- Procedure: (1) Analytic Hill vortex + uniform-weak B₀ seed; advect material surface one turnover (exact Hill streamfunction trajectory integration, tolerance-controlled); report |ΔΦ|/Φ at two tolerances. (2) Added-mass: far-field dipole extraction for Hill ball and thin ring (a/R=0.05) → m_added vs reference laws; confirm Γ-independence (repeat Γ ∈ {0.5, 1, 2}, m constant within 3%). (3) Cross-helicity H_c = ∫u·B evolution: report dH_c/dt = boundary + reconnection-term decomposition (analytic terms + measured residuals).
- Predicates: PASS = flux conserved to 1% + m Γ-independent within 3% + H_c budget closes to 5%. FAIL = flux leaks at fixed tolerance (→ advection-scheme scoped) or m varies with Γ (→ dipole extraction bug; M3 ledger claim FAILS — report, do not patch silently).
- Frozen inputs: ρ=1, R=1, B₀/Γ ratio scanned {0.01, 0.1} (weak-coupling ladder toward D3 continuation).
- Out: back-reaction continuation ([M3-B3]), persistence, chirality/statistics ([M3-B2]).

## Common rules
PDG comparators: derivation-input only, zero validation contact. All scripts thin, equation-owned, receipt-captured (source/command/env/stdout/stderr/elapsed/verdict). Any firewall hit from drift → route verdict update in 03-reconciliation addendum, never silent edit of 01.

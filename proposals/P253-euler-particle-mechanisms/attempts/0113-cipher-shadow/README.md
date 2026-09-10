# 0113 cipher — M2 Euler-shadowing scoping study (NO proof claim)

Question: what would a filament→Euler shadowing transfer for the PoC-2 leapfrog orbit need? This file scopes candidacy, hypotheses, regularity, timescales, and the ordered failure inventory. It proves nothing and transfers nothing.

## 0. Objects on each side
- MODEL side (banked): 4D reduced thin-ring ODE (Saffman self-induced + regularized mutual BS, fixed core a=0.05, axisymmetric); converged relative periodic orbit T=4.088, elliptic Floquet (|μ|=1±2e-6), reduced Hessian −1.66 + Z zero mode.
- EULER side (target): exact 3D incompressible Euler, R³ finite excess energy, smooth compact ω with thin-core neighborhoods of the model orbit; full 3D perturbations (bending modes included).
- Transfer sought: an Euler solution staying δ-close to the model orbit over ≥1 period (finite-time shadowing), with δ→0 as a→0 at stated rate. NOT claimed.

## 1. Shadowing-theorem candidacy (all fail off-shelf — recorded, not repaired here)
- Bowen/Anosov shadowing: needs uniform hyperbolicity. FAILS structurally — the orbit is ELLIPTIC (unit-circle Floquet by construction of the KAM route). Exponential tracking is the wrong conclusion to shop for.
- Infinite-dimensional KAM (Kuksin/Craig–Wayne type): needs Hamiltonian formulation + nonresonance + small REGULAR perturbation. Euler has Hamiltonian structure (Arnold Lie–Poisson), but filament reduction is a SINGULAR limit (core→0, log-divergent self-induction renormalized into Saffman law), not a regular perturbation of Euler. No off-shelf theorem bridges singular-asymptotic reductions to KAM tori. Gap, not verdict.
- Finite-time (Palmer-type) shadowing: needs hyperbolicity on the segment or exact trajectory + Gronwall with small defect. Defect route is the honest candidate: construct an approximate Euler solution from the filament (core profile graft, e.g. Lamb–Oseen-type viscous-blended or exact axisymmetric Euler core) with defect D(a), then Gronwall over T. Requires: (i) defect bound uniform in a with rate; (ii) linearized-Euler propagator bound on [0,T] in the chosen norm (no exponential growth assumption — elliptic background gives at most polynomial/subexp growth IF no embedded instability; the Cao/Hill embedded-spectrum lessons forbid assuming this).
- Averaging over fast core rotation (T/T_core ≈ 83, measured §3): the fast angle is the natural small parameter ε = T_core/T_leap ≈ 0.012. Candidate: finite-time averaging with remainder O(ε) over one period. Needs the core profile to stay O(ε)-rigid — exactly what breaks first (§4, item 2).

## 2. Exact hypotheses a future transfer proof would need (checklist, all unearned)
- H1 (core rigidity): an axisymmetric Euler core profile with cross-section deformation ≤ C·(strain/core-vorticity) ≈ C·0.046 per passage (§3 numbers), uniform over ≥1 period.
- H2 (remainder rate): model↔Euler velocity defect ≤ C·a²|log a| (or stated rate) in exterior L²∩L∞ on [0,T].
- H3 (linearized bound): Euler linearization about the grafted approximate solution grows ≤ C(1+T)^k (k stated) in the energy-enstrophy norm pair — no embedded exponential mode on this background (must be CHECKED, never assumed; cf. 0095/0104 fixed-observation growth lessons).
- H4 (bending modes): azimuthal modes m≥1 controlled over [0,T] (Crow-type long-wave pairing is the named suspect; needs its own estimate, not axisymmetry inheritance).
- H5 (no reconnection approach): d_min/a ≥ 8.2 maintained (§3) with stated margin; topology change excluded on [0,T] by distance, not by wish.

## 3. Timescale limits (banked numbers, model units Γ=R=1, a=0.05)
- T_leap = 4.088; T_core ≈ 0.049 (ratio 82.8 — fast angle available for averaging).
- Passage: d_min = 0.4115 (d/a = 8.2); peak mutual strain 5.9 vs core vorticity 127 → 4.6% deformation pressure per passage.
- Validity ceiling (conjecture for the proof to test): transfer plausible while (strain×T_leap)/(core-vorticity) ≪ 1 AND d_min/a ≫ 1; both hold marginally here (0.19 and 8.2). Tightening a improves both as (a², a⁻¹)-powers modulo log — the quantitative bet, unproven.

## 4. What breaks first (ordered failure inventory)
1. Azimuthal bending (H4): the model cannot see m≥1; Crow-type pairing during the 4.6%-strain passage is the likeliest first breaker. Kills 3D shadowing while axisymmetric tracking survives.
2. Core cross-section deformation (H1): Saffman law assumes rigid translating core; strain accumulates secularly over many periods even if one period holds. Kills multi-period KAM before single-period shadowing.
3. Ellipticity → no exponential tracking: even with H1–H5, only finite-time polynomial closeness is on the table; all-time KAM needs the full singular-KAM construction (no candidate theorem — deepest gap).
4. Close-passage cascade: if d_min/a → O(1) (tighter pairs), BKM-scale gradients invalidate H2 first; the M1-BKM threshold language takes over (route handoff, not failure of physics).

## 5. Minimal next construction (named, not executed)
Finite-core axisymmetric Euler simulation (MOL, small-ratio discipline) initialized on the PoC-2 orbit at a ∈ {0.05, 0.025}, measuring shape-space deviation over T=4.088 vs a-rate — an empirical shadowing check at H1/H2 scope, axisymmetric only (H4 stays open). Needs its own frozen design; not licensed by this study.

## Verdict
Scoping COMPLETE as a study: SURVEYED candidacy exhausted (Bowen / infinite-dim KAM / Palmer-type assessed — §1; no off-shelf theorem applies; defect+Gronwall + finite-time averaging are the live candidates). This closes the surveyed set, NOT the existence of a theorem: an unlisted mechanism or a new singular-KAM construction stays open and would reopen §1 append-only. Hypotheses H1–H5 frozen as the proof's shopping list, timescales banked, failure order stated. M2-B1 stays open. No proof claim, no transfer, no verdict change.

# A3 transfer audit (SOURCE-TRANSFER-FIRST per shepherd constraint)

Rule: no Ruban/Buttà prediction is consumed until its applicability is established below.
Method sections read at source (Ruban 1806.08173 §§I–III incl. model eqs. (6)–(14);
Buttà–Cavallaro–Marchioro 2310.00732 abstract + regime statement; SIAM SIMA 57:789–824).
Full-text verification of Buttà proofs deferred to the construction that needs them (named).

## T1 — Ruban 2018 (parametric 3D stability bands)

Source object: N quantum-vortex filaments, regularized BS + LIA stiffness term
(Γλ/4π)κb, Hamiltonian with canonical pair (Z_n, R²_n/2); Λ = λ+ln(R₀/a) = ln(R₀/ξ);
λ=0.5 matched to Kelvin hollow-core ν=1 bending mode; δ-renormalization
(a→δ, λ→λ+ln(δ/a)) proves long-scale dynamics sees only Λ. Stability result:
Floquet ρ=exp(±√μ T) per azimuthal m; stable bands between non-overlapping main
m-resonances at Λ≈4–8, W≈0.2–0.25 (ΔW∼0.01–0.05); overlap/no-stability at Λ≲3;
none for N≥4 (W>0.1). Strain criterion: Ω_q=√(ω_q²−σ²), λ>0.137 kills
strain-induced filament instability. Method: pseudo-spectral lifetime scans
(we compute Floquet directly instead — stronger where feasible).
TRANSFERS: parametric-resonance framework + ρ formula + (Λ,W) coordinates +
canonical pair + strain-vs-stiffness criterion form + Λ-only long-scale
insensitivity (licenses my fixed-a comparisons across runs).
DOES NOT TRANSFER: Γ=2π units (rescale to mine); hollow-core λ=0.5 stiffness
(mine: Rosenhead–Moore + Saffman — FIRST BROKEN CORRESPONDENCE #1, affects band
edges); superfluid core physics (no Euler reconnection/stretching — #2);
band locations as imports (my (Λ,W) point measured in MY model, bands as
prediction-to-test, never as assumed stability).
My coordinates: Λ_mine = ln(8R/a) ≈ 5.08 (Saffman convention) vs ln(R/a) ≈ 3.00 —
convention straddles the overlap→band edge (pinned, not hidden); W ≈ 0.2 in-band
(estimate; A3 measures exactly). PREDICTION (testable, not consumed): weak/no
m=1,2 growth; a→0.025 (Λ→5.77) moves deeper in-band. If A3 finds strong growth,
broken-correspondence #1 is identified and itself banked (stiffness matters).

## T2 — Buttà–Cavallaro–Marchioro 2023/2025 (scaling limit)

Source object: axisymmetric no-swirl Euler, N rings thickness ε, vorticity mass
AND radius O(|log ε|); ε→0 ⇒ ring centers follow finite-dim point-ring system;
two large rings ⇒ longer times covering SEVERAL overtakings (rigorous leapfrogging).
TRANSFERS: the H2 REMAINDER-RATE PROGRAM — thin-ring Euler solutions tracked by
reduced dynamics over multiple passages with a rate; regime hypotheses to check.
DOES NOT TRANSFER: all-time/KAM persistence (finite passages only); 3D bending
(axisymmetric class — H4 untouched); swirl; conclusion at MY parameters.
FIRST BROKEN CORRESPONDENCE: my regime (R=1, Γ=O(1), a=0.05, |log a|≈3) vs their
large-radius O(|log ε|) scaling — applicability is CHECKED in A3 (orbit
recomputed at scaled-up (R,Γ) per their regime as a control), never assumed.
Dávila CPAM-2024 exact leapfrogging Euler solutions stand as the existence anchor
(frozen supplier, unchanged).

## Clearance statement
Transfer phase COMPLETE: both suppliers typed (equation/domain/hypotheses),
transfers vs breaks itemized, my-coordinate mapping pinned with conventions
explicit. A3 resonance-scan construction is CLEARED to test (not assume) the
Ruban prediction and to run the Buttà-regime control. Findings feed back here
append-only.

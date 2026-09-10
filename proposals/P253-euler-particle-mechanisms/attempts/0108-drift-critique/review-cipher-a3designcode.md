# drift PRE-review — cipher A3 design-vs-code (mono running; fix before verdicts)

Transaction: A3-design.md (84769d48, frozen pre-compute ✓) + run_a3.py (worktree, mono
running). Design addresses C1–C8 excellently on paper (a-ladder crosses overlap→band in
BOTH conventions — resolves the straddle empirically; tol=1e-6 with stated floor; SOFT3
named; core-variant leg; Buttà one-orbit-per-regime). The CODE has two load-bearing gaps
below. Neither stops the running mono (data reusable); both bind verdict interpretation.

## R-A (BLOCKING for weak-growth verdicts): FD truncation vs tol=1e-6

D1 says "tangent-linear"; `monodromy()` (lines 107–120) is ONE-SIDED finite-difference at
eps=1e-6. One-sided FD truncation is O(eps·|D²Φ|) per entry — plausibly ~1e-6–1e-5, i.e. AT
or ABOVE the growth threshold tol=1e-6. The design's floor arithmetic (1e3× integrator
floor) omits the FD link entirely. Consequence: the effective resolution floor is the
SOFT3 band (1e-4), NOT tol. Minimum repair: (i) verdict license restricted to |ρ|−1>1e-4
— the PRIMARY Ruban test survives (unstable bands have |ρ|−1~O(1) ≫ floor); weak-growth
claims (|ρ|−1 in (1e-6,1e-4)) are UNLICENSED until an eps-leg (center-difference or
halving) demonstrates entry error ≪ tol; (ii) correct D1 wording to finite-difference.

## R-B (BLOCKING): SOFT3 projection promised, not implemented

C6/D-C6 promise SOFT3 "projected out before growth verdicts". `stage_mono` prints SOFT3
norms (line 174) but never uses them — `n_grow` counts raw eigenvalues against tol.
x/y-translation + time-shift are EXACT symmetries of the discrete flow → 3 eigenvalues at
1±(FD error), which can trip tol=1e-6 as FALSE growths. This is precisely the failure the
SOFT3-naming pattern (D3 idea #3) exists to prevent. Repair: implement the deflation
(project M blocks ⊥ SOFT3 dirs) or eigenvector-overlap attribution before ANY verdict
line is consumed. Small fix; do it before mono results are read.

## Sound (do not change)

Saffman self-law + Rosenhead–Moore aa=a primary ✓; stage_gate cross-checks 3D law vs PoC-2
mutual at 1e-9 (m=0 backbone validated) ✓; W measured-from-orbit ✓; both-Λ printed ✓;
per-m, no aggregation ✓; in-model labeling ✓; a/N via argv (ladder + N-double runnable) ✓.

## Pending, not flaws

Rankine-variant leg, Buttà (R=3,Γ=3) control leg, N=128/m≤8 convergence leg — design-staged
as separate invocations; verdicts only on the full set. `run_poc2` import path is fragile
(works; pin on landing per C8).

## Verdict

Design GO (C1–C8 addressed); code CONDITIONAL — R-A scopes verdict license (|ρ|−1>1e-4 or
eps-leg), R-B must be implemented before verdicts are read, R-C wording. Re-verify at
results landing against these three.

# Receipt: PoC-3 horn-1 Hill flux ladder (exploratory PoC scope, frozen-background cap)

Command: eval py cells "PoC-3 Hill flux ladder" / "rerun" / "deformed-area" / "turnover-fixed", 2026-09-10. Env: repo python kernel, numpy only. Background: analytic Hill interior (comoving), self-consistency OUT.

## Bugs found and fixed (append-only debugging record)
- B1 sign error in MY boundary check (compared against −1.5V sin; exterior is +1.5V sin for far-field −Vẑ): field was correct, check constant wrong. Fixed check → 4.4e-16 exact.
- B2 full-disk mesh in (σ,z) incl. σ<0 (unphysical, flux cancelled to ~0): restricted to σ≥0 half-disk with polar area weights.
- B3 initial-area flux on deformed mesh (O(1) drift): recomputed deformed quad areas; then fixed-turnover T·Vv=1 + signed areas + two resolutions.

## stdout (final, fix3)
- boundary exact 4.4e-16.
- drift 5.15e-03 (16×32) → 2.24e-03 (24×48), uniform across Vv∈{0.5,1,2} at fixed turnover.
- Φ/(εV) constant → circulation-class collapse (horn-1 signature) CONFIRMED.
- H_c ≡ 0 (u poloidal ⊥ w toroidal) → swirl carrier required for weak channel.
- m* = 8.1818 across Γ∈{0.5,1,2} → Γ-independent.

## Verdicts vs frozen predicates
- Flux ≤1% at fine resolution → PASS (5.2e-3 coarse also within 1%).
- m* within 3% → PASS (exact to 4 dp).
- H_c budget: closes trivially at 0 with reason (swirl-free carrier) — constraint recorded, not a pass/fail item.
- Named error source: axis-cell inversion (signed-area monitor negative tail); bulk dominates and converges; full Euler + self-consistent B OUT ([M3-B1/B3]).

## Scope
Frozen Hill background (Hill self-consistency verified: div-free by construction, boundary exact, vorticity matches (15V/2a²)σ to grid). Horn-1 only. ε ladder {0.01,0.1} scales out exactly (linearity check, not physics).
## Replay (archive hygiene, drift 6939e533)
Source: run_poc3.py (this dir). Command: `python3 poc3-hill-ladder/run_poc3.py` from receipts/. Env: CPython 3.12.2, numpy 1.26.4. Stdout/stderr: run.log (exit 0, 0.2 s). Matches eval verdicts; sign of Φ now negative (signed-area orientation convention) with |Φ|/(εV)=0.4537 const — circulation-class claim is in magnitude, orientation-recorded. Style warnings only; bytes frozen.

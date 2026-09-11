# Receipt: D3b adiabatic-Berry (FROZEN F-bar, pre-committed two-commit)

Command: `python3 run_d3b.py` (836.8 s). Env: CPython, numpy; A3 backbone
N=64; N4 shoot/secflow machinery (damped, banked hot starts).
Numbers: 11/11 stations clean (no guard fired — smooth patch holds over the
Γ ∈ [1.0, 1.01] loop); transverse phases ~1e-4 all stations; chain −0.00000,
Σφ = 0.00073, γ_geom = −0.00073 → |γ| ≤ 0.01 → DEAD per frozen bar (two
orders of margin below the 0.01 line).
Verdict: DEAD (no geometric accumulation at 1e-2 floor). The orbit-bundle
Berry route is doubly closed (N4: no bundle at large; D3b: no phase at
small). X1 (deflection-sign) unaffected — different observable, own bars.
In-run note: T prints 4.08800 all stations (moves in 6th+ decimals only —
adiabaticity genuinely slow); overlaps ≥ 0.99 throughout (no level
crossings). No repairs needed; script ran as frozen.

# Receipt: D3b adiabatic-Berry (FROZEN F-bar, pre-committed two-commit)

Command: `python3 run_d3b.py` (836.8 s). Env: CPython, numpy; A3 backbone
N=64; N4 shoot/secflow machinery (damped, banked hot starts).
Numbers: 11/11 stations clean (no guard fired — smooth patch holds over the
Γ ∈ [1.0, 1.01] loop); transverse phases ~1e-4 all stations; chain −0.00000,
Σφ = 0.00073, γ_geom = −0.00073 → |γ| ≤ 0.01 → DEAD per frozen bar (two orders margin).
Wrap-fix (drift review, committed): γ wrapped to (−π,π] before adjudication (Berry phase defined mod 2π; drift rerun hit ALIVE via 2π winding — same physics, opposite code-verdict). Determinism
re-run (815 s): identical numbers, wrapped −0.00073 → DEAD. Receipt now
branch-stable; verdict deterministic both branches (14× margin mod 2π).
In-run note: T prints 4.08800 all stations (moves in 6th+ decimals only —
adiabaticity genuinely slow); overlaps ≥ 0.99 throughout (no level
crossings). No repairs needed; script ran as frozen.

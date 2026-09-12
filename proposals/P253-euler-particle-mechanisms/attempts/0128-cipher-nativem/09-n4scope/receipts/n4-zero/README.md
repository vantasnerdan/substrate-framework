# Receipt: N4 zero-test (FROZEN scope 09-n4scope/00-scope.md) — DEAD-no-bundle

Command: `python3 run_n4.py` (+ Γ-probe follow-ups). Env: CPython, numpy;
A3 3D backbone N=64. Scope stop rule executed exactly.
Family continuation probe (damped smooth shooting, banked hot starts):
- Γ=1.005: (0.87109, 1.28981) res 2.5e-13 ✓
- Γ=1.01: (0.94938, 1.37305) res 2.4e-12 ✓
- Γ=1.02: (1.71, 2.16) res 1.1e-2 ✗ (damped, no descent — genuine)
Full (Γ,a)-loop (r=0.02,0.001): station 0 re-shoot fails → shrunk loop fails
→ STOP per scope (second failure). Branch-jump guard + banked-hot-start
discipline held throughout (an earlier unguarded run slid to R~20 trivial
branches — caught by the guard, discarded, documented here).
Verdict: DEAD-no-bundle (family folds inside any loop enclosing ±2% in Γ).
Mechanistic bonus (not a verdict upgrade): dR/dΓ ≈ 24 near banked + fold
within ~1–2% CONFIRMS the Krein marginal-collision picture independently —
T-window collision (A3) and Γ-fold adjacency (here) are two signatures of one
bifurcation-adjacent orbit. X1 extends unaffected (sign physics independent
of orbit-bundle Chern; disjoint budgets per scope).
In-run repairs (disclosed): damped Newton (undamped blew up to res 8e4);
explicit Gam/aa thread-through (run_a3 defaults frozen at import — silent
no-op caught before verdicts); continuation chaining + branch guard.

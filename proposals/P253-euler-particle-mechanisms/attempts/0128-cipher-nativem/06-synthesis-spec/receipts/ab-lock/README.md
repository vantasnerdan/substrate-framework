# Receipt: B1 phase-charge lock (FROZEN B1 F-bar) — verdict HOLD (kinematic)

Command: `python3 run_ab.py` (0.2 s). Env: CPython, numpy. F-bar:
06-synthesis-spec/01-b1-fbar.md (frozen pre-compute).
Numbers: Φ = 1.0000 / 2.0003 / 3.0012 vs L = 1/2/3 (≤0.04%); control −0.0;
N-leg (128 vs 256): 1.1e-2 (bar 5% = 5× floor ✓).
Verdict: HOLD — phase counts source-framing integer exactly; phase-charge
lock structure confirmed at kinematic level (geometry + BS kernel).
In-run bug (disclosed): first run double-divided by 2π (Φ=L/2π, ratios still
exact 1:2:3 — lock showed through the bug); normalization fixed, rerun HOLD.
Scope honesty: kinematic lock ONLY — dynamical charge (missing-5) untouched;
P1/F-a conditionality stands. Lint style-only.

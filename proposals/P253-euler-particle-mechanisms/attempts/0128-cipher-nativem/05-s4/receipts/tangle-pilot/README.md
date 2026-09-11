# Receipt: S4 tangle pilot (FROZEN F-bar) — verdict KILL, mechanism named

Command: `python3 run_tangle.py` (3.6 s). Gate: Hopf −1.004, far 0.0 (PASS).
Env: CPython, numpy only. F-bar: 05-s4/00-fbar.md (frozen pre-compute).

## Numbers (S=32 seeds)
- Base (Rb=3, Nt=40, Λ≈0.67): μ=−0.100, σ=0.563, m1=0.293, CV=1.924.
  F-S4a generator PASS (|μ| ≤ 0.20 honest-margin).
- F-S4b: CV=1.92 > 0.50 → KILL (no stable value; distribution Poisson-like
  counting noise, mostly ~0 with occasional ±1,±2).
- Rb=5: ratio 0.984 PASS (not boundary artifact — the statistic is real noise).
- 2×Λ: α=0.761 > 0.5 → KILL-lean (charge scales with tangle = tangle property).
- Rotations: spread 0.024 ≤ 0.199 → PASS (isotropic; F4 frame-fear unfounded).
- F-bar verdict: KILL (two independent lines: CV + α).

## Kill mechanism (named): counting-noise, never integer-class
Linking with a random tangle is a POISSON-like counting statistic (CV~O(1)
at any order-unity density): each realization draws a different integer, so
no realization-independent charge exists. Lower Λ → m1→0 (charge evaporates);
higher Λ → worse (α≈0.76). Dead at every density, all branches including the
pre-registered S1×S4 hybrid (magnitude carries the instability; framing sign
cannot stabilize a fluctuating magnitude).

## Residue (bounded, not oversold)
Tangle-linking statistics as a DECOHERENCE/noise source for filament carriers
(Poisson counting noise with measured CV/α scaling) — a reusable noise model,
not a charge. Filed for any future stochastic-carrier work; no claims beyond.

## Qualifier (drift: within-paradigm convergence)
Residue convergence (integer-linking + BF-form + this noise model) holds
WITHIN the shared linking-paradigm framing (M1 labels assumed throughout) —
independent angles agreeing inside one paradigm ≠ independent confirmation of
it. This qualifier travels with the convergence claim. Discipline note: F-bar
freeze and run shared one commit (content-attested, weaker than F1 two-commit;
future bars freeze-commit first). Gate now banked in-script (Hopf −1.004);
gray clause coded, untriggered.

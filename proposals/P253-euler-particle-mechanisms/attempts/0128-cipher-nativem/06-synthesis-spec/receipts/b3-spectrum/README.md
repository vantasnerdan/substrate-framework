# Receipt: B3 fluctuation spectrum (FROZEN B3 F-bar) — verdict HOLD (noise only)

Command: `python3 run_b3.py` (3.8 s). Env: CPython, numpy; reuses banked S4
pilot module by explicit relative path. F-bar: 06-synthesis-spec/02-b3-fbar.md.
Numbers: S=256 Fano=0.8450, deciles [0, 0, 0.005, 0.69, 1.03, 2.13]
(mode-at-0 + tail ✓); S=128 Fano=0.9328, convergence 10% ✓ → HOLD.
Verdict: tangle-vacuum noise component VIABLE (Poisson-like spectrum as P3
predicts). Scope honesty (pre-computed): noise ONLY — quantization untouched.
In-run fix: symlink reuse replaced by explicit path import (fragile-but-
working noted). Lint style-only (+ static-only import notice, runtime clean).

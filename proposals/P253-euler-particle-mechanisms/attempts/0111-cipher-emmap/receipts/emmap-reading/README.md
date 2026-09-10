# Receipt: EM-map reading selection (exploratory PoC scope, one-way cap)

Command: eval "EM-map reading selection" + fix (2026-09-10); archived source run_emmap.py.
Replay: `python3 emmap-reading/run_emmap.py` from receipts/ → run.log, exit 0, 0.0 s.
Env: CPython 3.12.2, numpy 1.26.4.

## Debug trail (append-only)
- First run crashed: np.trapezoid absent in eval kernel → np.trapz (script uses trapz).
- Expectation bug (mine): self-disk flux/Gamma = 3.077 ≠ 1 — code CORRECT, expectation wrong: self-spanning flux = self-inductance L = R[ln(8R/a)−2] = 3.07517, verified digit-for-digit. Geometry quantity, recorded as anti-quantum control.
- Design correction: AB quantization demo moved to confined-flux reading (test disks enclosing tube); current-loop reading kept as the failing contrast.

## stdout (replay, matches eval)
flux_tube/Phi = 1.00733 / 1.00659 / 1.01391 at test r = 0.15/0.2/0.3; current-loop toroidal flux = 0 (printed analytic zero); self-disk L = 3.07517; motional E_φ tabulated. Full text: run.log.
## Verdicts (downgraded per drift ed67f69d: predicate-level, not measurement)
R-EM5 selection SOUND as conditional logic; script role = analytic illustration + quadrature regression (1% pixelation consistency). R-EM1 one-way map: derived under QUASI-STATIC scope (instantaneous BS kernel, v≪c_EM, no radiation — named, not smuggled). R-EM2/R-EM3: open imports, untouched. No comparator contact. Lint style-only; bytes frozen.

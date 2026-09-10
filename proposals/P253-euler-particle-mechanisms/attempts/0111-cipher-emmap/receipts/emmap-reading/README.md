# Receipt: EM-map reading selection (exploratory PoC scope, one-way cap)

Command: eval "EM-map reading selection" + fix (2026-09-10); archived source run_emmap.py.
Replay: `python3 emmap-reading/run_emmap.py` from receipts/ → run.log, exit 0, 0.0 s.
Env: CPython 3.12.2, numpy 1.26.4.

## Debug trail (append-only)
- First run crashed: np.trapezoid absent in eval kernel → np.trapz (script uses trapz).
- Expectation bug (mine): self-disk flux/Gamma = 3.077 ≠ 1 — code CORRECT, expectation wrong: self-spanning flux = self-inductance L = R[ln(8R/a)−2] = 3.07517, verified digit-for-digit. Geometry quantity, recorded as anti-quantum control.
- Design correction: AB quantization demo moved to confined-flux reading (test disks enclosing tube); current-loop reading kept as the failing contrast.

## stdout (replay, matches eval)
flux_tube/Phi = 1.00733 / 1.00659 / 1.01391 at test r = 0.15/0.2/0.3 (pixelation ~1%, size-independent → R-EM4(i) PASS at 1% level); current-loop toroidal flux = 0; self-disk L = 3.07517; motional E_φ read-off tabulated.

## Verdicts
R-EM5 selection: PASS (confinement required for quantization). R-EM1 one-way map: derived. R-EM2/R-EM3: open imports, untouched. No comparator contact. Lint style-only; bytes frozen.

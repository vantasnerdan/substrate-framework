# Receipt: 0113 §3 passage numbers (diagnostics only, no transfer claim)

Source: run_passage.py (imports run_poc2 rhs2/rk4 from the 0108 receipt — same bytes drift reproduced; Pyright import-resolution notice is static-only, runtime sys.path insert works).
Replay: `python3 passage-numbers/run_passage.py` from receipts/ → run.log, exit 0, 0.7 s. Env: CPython 3.12.2, numpy 1.26.4.
Output matches study §3: d_min=0.4115 (d/a=8.2), strain 5.905 (4.6% of core vorticity), T ratio 82.8. Lint static/style-only; bytes frozen.

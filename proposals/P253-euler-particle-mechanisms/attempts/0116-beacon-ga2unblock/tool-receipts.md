# Tool receipts — 0116 (beacon)

## U1 — solver inventory (exit 0, via bg_5)

`.venv/bin/pip list | grep`: numpy 2.5.1, scikit-fem 12.0.2, scipy 1.18.0,
sympy 1.14.0. `import firedrake/dedalus/dolfin`: all ModuleNotFoundError.
`.venv/bin/`: python3.12, pytest, no MPI/FEM launchers.

## U2 — repro run (EXIT 1, the block reproduced)

`.venv/bin/python .../0116-beacon-ga2unblock/ga2_repro.py` → present-list,
heavy-solvers all NO, in-tree solver NONE, member data NONE,
`BLOCKED: ... ga-status rows 1-5,8-9 uncomputable`, exit 1 (0.76 s).
LSP: 4× E501 cosmetic only on attempt-local script; left as-is.

## U3 — 0114 repair verification (exits 0/1 as labeled)

- seed0 argv rerun: numbers reproduce run-1 exactly
  (0.1262/0.9996/0.0022), exit 1 frozen-BLIND → `seed0-run.log`.
- seed7 argv rerun: 0.1262/0.9996/0.0008, exit 1 frozen-BLIND (numbers feed
  amended reading) → `seed7-run.log`.
- Dead control-set removal changed nothing (exact reproduction proves it).
- pycache: zero `*.pyc` under 0114/0111 after cleanup (verified by find).
- Repairs receipt: `attempts/0114-beacon-s9/repairs.md`.

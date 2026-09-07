# Packaging test failure receipt

These setup/implementation failures occurred while moving the exact helpers
from the attempt into the importable package.  Neither is scientific evidence.

1. The first command used the absent worktree-local interpreter:

   `PYTHONPATH=src .venv/bin/python -m pytest -q tests/test_euler_p2_principal.py`

   It exited `127` with `/bin/bash: line 1: .venv/bin/python: No such file or
   directory`.  The repository interpreter was then resolved read-only as
   `/home/dan/substrate-framework/.venv/bin/python`.

2. The first repository-interpreter run reached the tests but compared two
   SymPy expressions structurally.  They printed identically while carrying
   distinct polynomial symbols, producing `1 failed, 2 passed`.  Replacing the
   structural comparison by the determinant identity
   `simplify(det(lambda I-H)-(lambda^2+d^2-b^2))==0` repaired the oracle.  The
   subsequent focused receipt is the first successful scientific execution of
   that predicate.

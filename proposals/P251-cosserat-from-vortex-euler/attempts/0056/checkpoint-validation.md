# Active-campaign checkpoint — 0047–0058

Base commit: 6b037a5. Scientific boundary: two new unpromoted modules
euler_orbit.py and micropolar.py, their 14 tests, completed attempts
0047–0056/0058 and independent reviews, and bounded CST004/005/006 repairs.
Ongoing 0057/0059/0060/0061 assembly work is not included in this checkpoint.
No accepted claim, release, existing exported public API or shared numerical
machinery changed. This is not a terminal PR refresh or campaign completion.
Validated staged scientific tree: 397b2328780f95e0230a2b519264d817d06bfe60
(before adding this receipt and its full-workflow output).

The repository's staged-change decision selects full validation because
multiple framework modules were added. GitNexus staged-worktree detection
reports low indexed risk and no changed indexed symbols; it does not yet
index these new definitions. Direct callers are the new tests and the
explicit P251 reduction/operator consumers. Impact absence is not inferred
from the missing index entries.

First full invocation:

    PYTHON=/home/dan/substrate-framework/.venv/bin/python scripts/validate.sh --full

Exit 0. Full pytest: **2515 passed in 350.00 seconds**. All fixed repository
checks pass: 252 accepted claims, 12 valid proposals, 1031 memory records,
43 existing memory warnings, generated-state/import/compile and skill checks.
Complete first-run output is checkpoint-full-validation.stdout.txt.

Scientific receipts remain with the attempts. New orbit tests pass 9/9;
micropolar tests pass 5/5. CST004/005 pass 8/17 checks; the affected CST006
functions pass 4 checks without rerunning its unchanged Monte Carlo evidence.
The independent reviewers establish their stated attachment scopes and
the requested longitudinal-scope prose correction is closed. The initial
structural-equality test and misplaced consumer-function tail are preserved
as implementation failures with their passing repaired receipts.

Ruff passes on both new modules/tests and changed verifier code. Separate
staged/unstaged diff checks pass. This receipt remains valid for the unchanged
scientific boundary; ongoing work constructs the full-fluid ensemble/action
assembly and does not turn this checkpoint into a stopping point.

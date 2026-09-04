# Active-campaign checkpoint validation, 2026-09-05

Boundary: attempts 0035–0046 and two additive conditional Euler API modules
with their tests. No accepted claim, release, existing public API, numerical
infrastructure or immutable campaign was changed. Attempts 0047/0048 are
continuing scientific work, not completion or a terminal PR refresh.

The initial explicit scope was tests/test_euler_displacement.py and
tests/test_euler_affine.py: nine tests passed with all fixed checks. Inspection
of scripts/validate_changed.py shows its conservative rule selecting full
validation for multiple added framework modules. The full run therefore
followed, rather than representing the scoped receipt as repository-wide:

    PYTHON=/home/dan/substrate-framework/.venv/bin/python scripts/validate.sh --full

Exit 0. All fixed checks pass: 252 accepted claims, 12 valid proposals,
1031 valid memory files, 43 existing memory warnings. Full pytest: 2501 pass
in 351.84 seconds. The complete first-run stdout/stderr are retained here.
The host-specific interpreter path records the actual invocation, not a
runtime dependency introduced into the package.

GitNexus staged-worktree detection reports low risk and no indexed changed
symbols; the index does not yet contain these new definitions. Direct source
search supplies the meaningful boundary: new euler_displacement APIs have
only their tests and 0037 caller; affine_vorticity_energy has its tests and
0043 caller. The package root export surface is unchanged. This is not a
claim that the new APIs have zero consumers.

Ruff passes on the new modules/tests and individual scientific attempt
scripts. Independent scientific receipts exist for 0036, 0040, 0041, 0045
and 0046; first-run verifiers and correction histories remain in each
attempt. 0043/0044 individual reviews are in progress. The full test tally
does not substitute for those reviews or establish parent completion.

The source-folder .gitattributes marks PDF artifacts binary; their exact
bytes and source-audit digests are unchanged. Separate unstaged and staged
git diff --check calls pass after that file-type correction. The committed
changed-scope selector decision is recorded at this checkpoint; no second
equivalent full replay is needed for the same unchanged API boundary.

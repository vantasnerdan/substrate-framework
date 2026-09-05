# Additive Bernoulli-lift implementation validation

Base a93e78e. Frozen source/test SHA256s are in receipt.md. The only new
canonical surface is euler_forcefree.py; no existing contract, accepted
claim, shared numerical/verification code or dependency changes.
Direct source search finds its new test and0136 verifier only. GitNexus
precommit detection reports no indexed changed symbols and three tracked
changed files; it omits untracked new symbols, so that is not treated as
a complete dependency result. Reindex follows the checkpoint.

The selected boundary is tests/test_euler_forcefree.py plus all fixed
checks. The unchanged full2580-test baseline receipt0127 remains banked.

First workflow invocation omitted PYTHON and selected the worktree's
incomplete environment. workflow.stdout retains its fixed-check passes
and `No module named pytest`, exit1. This is an environment-selection
failure, not a science or test failure. Repaired command:

    PYTHON=/home/dan/substrate-framework/.venv/bin/python scripts/validate.sh --pytest-scope tests/test_euler_forcefree.py

workflow-repaired.stdout: exit0, all fixed checks and6 tests pass in2.46s;
263 accepted claims,1048 valid memory files,43 existing warnings. Ruff
and diff checks pass. A later proposal-only candidate registration does
not stale these API tests; its schema validator separately passed263/12.
The committed-head selector will confirm this bounded scope without
repeating the successful workflow.

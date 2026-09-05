# Targeted formatting checks

The first scientific run passed32/32, process0; see first-run.stdout.
The initially guessed virtualenv Ruff path did not exist (process127).
Discovery found /home/dan/anaconda3/bin/ruff. Its first targeted run
reported one unused local assignment, axial_prime (F841). Removed that
assignment only; no scientific expression or check changed. The
first-run verifier SHA before this formatting edit was
7e251f8521057fbab8f48acf8605da35403b131006ce8a5b5d3e6affc544cc75.

The direct trailing-whitespace search returned no matches (process1,
the normal successful negative-search result). No full suite or
unrelated consumer was replayed.

The corrected targeted Ruff invocation returned `All checks passed!`,
process0. Removing the unused assignment has no affected scientific
consumer, so the successful first32-check receipt is reused rather
than rerunning identical algebra for a formatting-only delta.

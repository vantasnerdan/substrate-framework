# Exact comparison repair

The first targeted run passed12 tests and stopped at one scalar-limit
assertion comparing two structurally different rational matrices. The
captured diagnosis reduces their complete difference to the zero2x2
matrix: its only raw entry was

    -2/[(t+1)(t^2+2t+1)]+2/(t+1)^3.

The assertion now checks that exact reduced difference, not Python's
expression-tree equality. No API, equation, coefficient or tolerance was
changed. Only that failed test is rerun; the12 unaffected passes are
reused under the user's blast-radius instruction.

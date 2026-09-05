# First execution and exact-oracle repair

The first targeted pytest execution is preserved in `first-pytest.stdout`:
15 passed, two failed, process 1. Both failures compare mathematically
identical SymPy expressions by structural equality after only one side has
been factored. The displayed mass derivative and trace-sector coefficient
are algebraically equal to their independent expected expressions.

The repair simplifies their DIFFERENCE to exact zero. No module equation,
coefficient, expected mathematical quantity or tolerance changed. The
frequency-mismatch invalid-input check was also made specific by supplying
otherwise valid positive probabilities; its old probabilities would have
triggered an earlier domain error. A separate captured direct replay tests
these repaired assertions and all other new-module consumers.

The first scoped Ruff execution passed. No numerical approximation, fitted
input or soft-eigenvalue threshold enters this implementation.

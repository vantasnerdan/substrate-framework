The first API test execution recorded seven passes and one failed assertion,
exit1. Every entry of the polynomial-degree predicate was True; the test
incorrectly put boolean objects in a SymPy Matrix and compared it to an
integer matrix of ones. This is a test representation error, not a residual
coefficient mismatch. The source implementation is unchanged. The repaired
predicate uses Python all over the exact polynomial degrees. Initial test
source and output are preserved; the repaired execution receives new paths.

The optional Ruff invocation through the pytest interpreter failed because
that environment does not install Ruff. The installed `ruff` executable in
PATH runs the same targeted source/test lint successfully with exit0. Both
outputs are retained; no dependency was installed and no source repair was
needed for this environment mismatch.

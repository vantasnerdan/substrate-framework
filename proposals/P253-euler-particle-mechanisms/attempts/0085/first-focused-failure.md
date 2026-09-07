# First focused-test repair receipt

The first focused test run exited `1` because the radiation denominator helper
returns a factored SymPy expression while the test used structural `==`
against its expanded form.  Replacing that assertion by `simplify(lhs-rhs)==0`
exposed the same representation-only issue in the next steady-limit assertion.
That second failed run is preserved under `focused-tests.*`.

Both assertions were changed to algebraic equality.  A later derived shell
coarea test exposed the same representation-only issue once more; that
receipt is preserved under `shell-first.*`.  No formula, sign, scientific
scope, or implementation output changed.  The final `focused-v3.*` receipt
passes all seven tests.

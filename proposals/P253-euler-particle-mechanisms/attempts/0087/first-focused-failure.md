# First focused-test failure

The first focused run exited `1` with four passing and three failing tests.
All three failures were representation/assumption assertions in the new test,
not disagreements with the formulas:

1. SymPy factored `sqrt(c_EM^2-c^2)` as
   `sqrt(c_EM-c)*sqrt(c_EM+c)`, so direct structural equality failed.
2. A symbolic positive `eta` did not encode `eta<1`, leaving an absolute-value
   branch in the speed-ceiling simplification.
3. SymPy retained `acosh(cosh(alpha))` rather than applying the declared
   positive-real branch.

The implementation was unchanged.  The tests were repaired by algebraic
simplification and an exact rational interior margin.  A later API audit added
strict subluminal, feasible-margin, and nonnegative-support domain checks;
the final focused receipt is `focused-v3.*`.

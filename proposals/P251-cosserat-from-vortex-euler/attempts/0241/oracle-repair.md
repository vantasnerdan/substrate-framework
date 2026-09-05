The repaired environment ran five tests: three passed, two failed. The failure
is in symbolic representation: differentiating a Piecewise band antiderivative
at its isolated t=0 branch and taking SymPy's limit returns that branch's zero
derivatives, losing the removable analytic extension. The finite angular
observation identities and solved weights were correct. Preserve the original
test snapshot and output. Use the t!=0 closed-form antiderivative (conds=none)
and take its exact analytic limit at zero; this tests the actual band output.

The first joint oracle passed the lower-order test and stopped at its exact
cubic check: structural matrix equality was used without simplifying the
difference. The repaired predicate tests the symbolic difference. During
source audit, its still-unexecuted Eq7 check was also found to use Q_t-Q_t,
which would have been only a cancellation tautology. Preserve its source
and output as preliminary evidence. Replace it with independent divergence
and time differentiation of full arbitrary space-time tensor fields after the
specified current transformation, plus a lost-q-dot negative control. Only the
corrected execution supplies the complete current-identity evidence.

The second joint execution passed the exact cubic mismatch and arbitrary
history residual. Its mutation setup then attempted to modify an immutable
SymPy matrix and aborted. The third version explicitly creates a mutable copy;
this changes only the negative-control implementation.

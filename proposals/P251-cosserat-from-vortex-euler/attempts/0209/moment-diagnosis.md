# First expanded verifier diagnosis

`moment-expanded.stdout` preserves six passing rows followed by an
assertion failure in the independent full-pressure angular matrix check.
The symbolic difference is `Matrix([[0, 0*I], [0*I, 0]])`: `factor`
preserved an unevaluated imaginary zero. Simplifying the exact difference
before comparison repairs this assertion. No field, pressure coefficient,
physical inertia, or claim formula changed. The secular-clock check was
also strengthened to differentiate the actual forced response rather than
compare two reordered expressions.

The repaired run was interrupted during a costly nested symbolic Hankel
integration (its partial output is preserved). Reversing the two absolutely
convergent integrals by Fubini evaluates the Gaussian tag transform first,
leaving elementary Gaussian wave moments. This is a representation repair,
not a changed moment or comparison value.

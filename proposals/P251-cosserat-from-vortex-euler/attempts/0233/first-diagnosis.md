# First exact-check diagnosis

The first execution exited 1 after 15 passing checks, in 2.1932 s.
Its two actual radial minors were 17/3 and 37/2, and the Gaussian
counterexample had zero determinant.

The next check compared factored and expanded representations of
`a/nu_p-a/nu_n` using structural equality. Its expected expression was
algebraically identical but differently represented. This is an oracle
implementation error, not a phase-gap or physical-rank contradiction.
The repair checks the simplified difference. The raw output remains
in first.stdout. The same mutable verifier also now evaluates the
explicit omitted-shear mutation rather than merely its known residual;
no governing equation or scientific selection criterion changes.

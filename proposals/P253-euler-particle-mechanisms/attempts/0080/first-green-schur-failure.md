# First exact-oracle failure

The first repository-interpreter execution is preserved under
`green-schur.*` and exited `1`.  Passes 1--5 established the parameter
Jacobian, physical impulse normalization, both triangular response factors,
and the Schur product.  Pass 6 compared two equal factorizations with Python
structural equality after expanding only one side.  SymPy therefore returned
`False` although their difference simplifies identically to zero.

The repair changes only this exposing assertion to compare the simplified
difference.  No scientific formula, expected coefficient, or sign changed.

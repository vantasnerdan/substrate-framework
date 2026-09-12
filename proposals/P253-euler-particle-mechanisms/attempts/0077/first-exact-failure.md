# First exact-oracle failure

The first direct repository-interpreter execution reached the fifth assertion
and exited `1`.  The scientific stabilizer equation was not contradicted.  The
oracle had declared `xI` with `nonzero=True` and then asked SymPy to solve
`xI*zeta_prime=0` for `xI`; that assumption excludes the intended solution.
The first four independently derived Maxwell/Lorentz checks passed.  The
repair removes the erroneous nonzero assumption from `xI` while retaining
`zeta_prime!=0`, then reruns the same scientific route with captured output.

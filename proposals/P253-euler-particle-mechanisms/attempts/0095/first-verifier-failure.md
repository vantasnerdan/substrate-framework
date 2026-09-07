# First exact-verifier failure and repair

The first repository-interpreter execution of
`verify_p2_principal.py` exited `1` at the Maxwell characteristic assertion:

```text
Traceback (most recent call last):
  File "/tmp/substrate-particles/proposals/P253-euler-particle-mechanisms/attempts/0095/verify_p2_principal.py", line 59, in <module>
    main()
  File "/tmp/substrate-particles/proposals/P253-euler-particle-mechanisms/attempts/0095/verify_p2_principal.py", line 26, in main
    assert sp.expand(char - expected_char) == 0
AssertionError
```

Classification: representation/implementation error, not a physical route
failure.  In a fixed transverse orientation, for example `k` along `z`,
`E` along `x`, and `B` along `y`, both Fourier curl entries have the same
`-i` sign.  The original scalar block used opposite signs, which would give
unphysical real free-Maxwell roots.  The repaired block is

    [[-i c_g xi_z,-i c_EM^2 |xi|],
     [-i |xi|,    -i c_g xi_z]],

whose characteristic polynomial is
`(lambda+i c_g xi_z)^2+c_EM^2|xi|^2`.  The first corrected execution is
preserved in `principal-verifier.*`; it exits zero with five of five exact
checks.

The verifier was subsequently expanded after the metric-curvature and Cao
polyhomogeneous cells were derived.  The original five-check receipt remains
unchanged provenance.  The first execution of the expanded predicate is
preserved separately in `principal-verifier-polyhom.*`; it uses the repository
interpreter, exits zero, and reports nine of nine exact checks.  In
particular, it checks the **complete** physical speed polynomial
`-c r^2/2 -> -c e_z`, rather than the invalid shortcut of declaring the
gradient of a truncated affine velocity to be zero.

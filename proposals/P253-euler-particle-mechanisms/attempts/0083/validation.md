# Validation

## Exact algebra oracle

The repository-interpreter command

    /home/dan/substrate-framework/.venv/bin/python \
      proposals/P253-euler-particle-mechanisms/attempts/0083/verify_scale_coverage.py

was executed three times.  The first run is preserved in `scale-first.*`; it exited
`1` because SymPy conservatively refused to cancel positive symbolic powers
with a symbolic noninteger exponent.  This was an oracle implementation
failure: the printed expression had zero expanded logarithm.  The repair
checks the logarithm of each positive dimensionless balance, for which
logarithmic injectivity is exact under the declared assumptions.

The corrected first execution is preserved in `scale-corrected.*`; stderr is
empty and exit contains exactly `0`.  The supervisor then required the safer
area-radius error because Lemma A.2 alone does not eliminate the full
first-order boundary mean.  That changed predicate 4 from the sharper
`log(N)/N^2` wave-number error to the source-supported `1/N` error.  The final
execution is preserved in `scale-v2.*`; stderr is empty and exit is exactly
`0`.  It checks:

1. the local-PDE balance (8) with the derived `C_s,H_0`;
2. the physical circulation balance (10) independently;
3. that an `O(epsilon^2)` auxiliary-radius error produces only an
   `O(epsilon^3)` physical-scale error through the exact equation (17);
4. that the safely supported wave-number scale error `1/N` is `o(h_N)` for the exposing
   envelope choice;
5. the exact `log N` plus rational-ray constant split in the Kelvin--Hicks
   row; and
6. the finite subluminal-ceiling algebra.

This oracle contains no discretization, quadrature, fitted value, soft
eigenvalue, or production numerical step, so `small-ratio-numerics` is not
triggered.

The final static author check is preserved in `static-validation.*`.  It uses
the same repository interpreter, parses `result.yaml`, checks the three route
verdicts, excludes the stale pre-0082 ledger phrase, confirms the explicit
seed/overlap/local-uniqueness continuation, verifies activation and final
oracle exits, and scans every 0083-owned file for trailing whitespace.  Its
stderr is empty and exit is exactly `0`.

## Analytic evidence not delegated to the oracle

- The compact-`W` uniformity of Cao Proposition 3.13 follows by tracking its
  proof through (3.31), (3.34)--(3.35), (B.13), and Lemma B.4: no derivative
  of `W` occurs and every coefficient is uniformly bounded on a positive
  compact `W` interval and compact radial annulus.
- Equation (20) uses the actual refined decomposition: the radial leading
  density contributes its second radial moment; the `O(epsilon)` odd
  correction contributes `s_epsilon*O(epsilon)`; and the even remainder is
  higher order.  The exact mean row then gives `R-x_epsilon,1=O(epsilon^2)`.
- Lemma A.2 proves an `O(epsilon)` relative free-boundary graph but does not,
  without a separate Fourier calculation, eliminate its full constant mode.
  Validation therefore uses only the safe `O(epsilon^2)` physical
  area-radius error in (23); the sharper cancellation is not claimed.
- The connectedness step consumes the author-scoped 0080 uniform rescaled
  `B_R` inverse.  Beginning from one seed, local charts have a uniform radius
  in relative epsilon; at each endpoint a new chart exists, and local
  uniqueness identifies the two solutions on their overlap.  Iteration
  constructs one connected thin component rather than merely a collection of
  local charts.  This is the one promotion-sensitive external dependency.
- The buffered endpoints use only the value error (24).  Multiplying
  `O(1/N^2)` by `n_1=NP` gives `O(1/N)`; no remainder derivative or
  monotonicity is invoked.
- The response formula (31)--(33) is obtained by two integrations by parts on
  a smooth compact regular cell and self-adjointness of the decaying Hodge
  solve.  Its use stops before asserting a nonzero curl or absolute KKS norm.

## Claim boundary

Validation establishes the algebra and exposes the exact analytic
dependencies.  It does not independently review 0080, prove the blocked
weighted-`C^1` Route A, evaluate `curl(P_DA F_12)`, normalize the physical KKS
dual, construct autonomous controls, determine `N_response`, or show that the
finite charged hierarchy is nonempty.  It validates no stability, P2/P4,
particle, quantum, or relativistic claim.

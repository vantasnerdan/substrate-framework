# Validation receipt

## Scope and activation

Only `proposals/P253-euler-particle-mechanisms/attempts/0014` was written for
the principal calculation.  The user separately authorised the bounded 0010
provenance repair recorded in that attempt.  No source module, proposal-level
manifest, governance file, generated document, skill, or memory record was
edited, and no commit was created.

The frozen central activation receipt was read before calculation:
`activation-schema.exit` contains `0`.  No additional source body was acquired
or introduced; the already-audited cached Cao--Zhan body was consulted to
confirm Theorem 1.3(iii).  The primary-source locations and hashes remain those
frozen and audited in `0010/source-audit.md`.

## Analytic oracle

The strongest oracle is the independent agreement of two derivations of the
linearised flow:

1. contraction of the full vorticity two-form with an arbitrary axisymmetric
   divergence-free displacement gives the physical tangent
   `(eta,chi)=({zeta,g}+{xi,a},{xi,g})` and the Lie--Poisson operator (6);
2. `J Hess(A_e)` then agrees term by term with direct linearisation of the
   time-dependent axisymmetric Euler equations after the active-core relations
   are imposed.

The WKB directions used to decide the Hessian sign are images of explicit
compactly supported displacement generators.  They are not arbitrary
Eulerian `(eta,chi)` variations.  The negative Cao sequence has Rayleigh
quotient tending to `-1`; independent pure azimuthal displacements have
quotient `+1`.  Thus both sign choices in the frozen `gamma_e` target have
infimum at most `-1`.

The continuation oracle computes the raw first-differential-order frozen symbol
of `J Hess(A_e)` on a regular active-core patch.  Its characteristic polynomial
is a repeated purely imaginary convective factor, while the remaining
off-diagonal is nilpotent.  The exact raw frozen solution therefore excludes an
`O(|k|)` exponential rate.  It does not exclude a finite `O(1)` rate after the
inverse-elliptic `eta` weight and the omitted `K` coupling are restored.

## First symbolic execution

`verify_leaf_hamiltonian.py` was executed once after materialisation.  The
exact command, stdout and exit are preserved in `first-run.command.txt`,
`first-run.stdout.txt`, and `first-run.exit`.

The seven checks are:

- multiply the Cao principal `J` and local Hessian matrices and impose the
  differentiated active-core relation;
- derive the double convective characteristic polynomial;
- square the off-diagonal remainder and obtain zero;
- reduce the general regular-swirl diagonal to `-u_mer dot k`;
- confirm the lower-left principal entry vanishes;
- substitute the actual Cao WKB tangent ratio and obtain the negative swirl
  norm; and
- flip the Hessian cross-sign as an exposing mutation and obtain
  `+3 chi^2/r^2`, not the claimed negative expression.

No soft eigenvalue or small-ratio numerical quantity was computed.  The
small-ratio-numerics skill therefore does not bind this analytic attempt.

Independent 0020 review subsequently supplied the canonical action-sign test,
the complete Poisson linearization `J d2A+(delta J)dA`, the energy-weighted
order-zero coupling audit, and the unbounded-form strengthening.  These bounded
corrections are recorded in `correction-receipt.md`.  The existing verifier's
actual predicates did not change, so it was not rerun and its first-execution
stdout/exit remain the sole execution receipt.

Final hygiene checks parsed `result.yaml`, found no trailing whitespace in
0010/0014, and returned exit `0` from the scoped `git diff --check`.  These are
format/provenance checks, not additional scientific oracles.

## Claim boundary

The result establishes the exact axisymmetric orbit tangent, mixed constrained
Hessian signature, unbounded both-sign form on the Cao `X`-unit tangent sphere,
and raw order-`|k|` regular-core characteristic for the same Cao and Gavrilov
carriers.  It does not establish a finite-rate short-wave verdict, global
free-boundary propagation, nonlinear persistence, non-axisymmetric stability,
or any quantum/relativistic particle property.  The next exposing oracle is an
energy-weighted global Kelvin/return and resolvent/semigroup construction with
`K`, `(delta J)dA`, support boundary, and translation modulation retained.

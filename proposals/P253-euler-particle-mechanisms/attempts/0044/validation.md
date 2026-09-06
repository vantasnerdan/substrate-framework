# Validation receipt and oracle boundary

## Exact execution

The first successful corrected execution used the repository interpreter:

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P253-euler-particle-mechanisms/attempts/0044/verify_curved_feshbach.py

Command, stdout, empty stderr, and exit `0` are preserved in the
`curved-feshbach.*` quartet. The verifier checks twenty-two exact predicates:

1. `psi_0=r Omega` solves the stated `m=1` Euler pencil;
2. the `b` derivative has the positive reciprocal derivative and weighted
   integrand `r^2 W'`, giving `-2F(a)` after the compact edge term;
3. the pointwise ground-state transform is exact, including the derivative
   whose endpoint cancels the physical Robin form term;
4. the scalar `psi` boundary coefficient is `-1/D`, with the stated logarithm;
5. exact toroidal divergence and scalar Laplacian expansions reproduce their
   `C1,C2` blocks through second order;
6. curvature sends the translation mode to the nonzero `m=0` coefficient
   `W/2` and the `m=2` coefficient `r Omega'/2`;
7. a finite-dimensional exact Hodge model checks all nine entries of the
   noncommutative second-order Leray expansion; and
8. the length/impulse/KKS compression has the displayed characteristic
   polynomial and the physical `l=1` neutral factor.

The exact oracle establishes algebra, signs, and projection ordering. The
ground-state compactness argument, compact-operator spectral accumulation, and
graph-domain implications are analytic functional analysis recorded in
`derivation.md`; they are not converted into sampled finite matrices.

## Exposing failures

- Reversing `partial_b(Omega-b)^-1` reverses the nonzero slope (5).
- Dropping the physical Robin form term destroys the ground-state identity.
- Applying the potential DtN directly to `psi` changes the logarithmic
  coefficient and fails the `-1/D` limit.
- Omitting inverse-Laplacian variations from `P1` or `P2` fails at least one of
  the nine exact projector entries.
- Treating all of `C1` as `m -> m+/-1` misses the `i l` terms in the exact
  divergence/gradient blocks.
- Treating angular orthogonality as a cancellation is exposed by the exact
  `m=0` coefficient `W/2`.
- A fixed two-dimensional Riesz claim is rejected analytically by the actual
  accessible axisymmetric frequencies `+/- |k|sqrt(mu_n)->0`, not by the
  inapplicable passive exterior-vorticity endpoint.

## Maximum verdict

The receipt licenses the compact-column response, exact toroidal Hodge/Leray
coefficients, nonzero cross-fiber mechanism, and actual leading centerline
matrix. It does not license a full-ring graph-domain eigenpair, finite-rank
uniform Riesz projection, finite nonlogarithmic Cao matrix, nonlinear rotating
wave, stability, or quantum interpretation.

`git diff --check` and a byte comparison between a fresh repository-interpreter
execution and the captured stdout are the final artifact oracles.

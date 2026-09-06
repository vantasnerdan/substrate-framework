# Validation and oracle boundary

## First exact execution

The first execution used the repository interpreter:

    PYTHONPATH=src /home/dan/substrate-framework/.venv/bin/python proposals/P253-euler-particle-mechanisms/attempts/0054/verify_hj2.py

The exact command, complete stdout, empty stderr, and exit `0` are preserved in
the `hj2.*` receipt quartet. It was a short local symbolic execution; no long
scripts-pane job or production numerics was used.

Seventeen exact predicates check:

1. the source rescaling from Cao's divergence-form equation to the fixed-core
   equation;
2. the genuine three-dimensional toroidal distance through second order,
   including the negative `t^4/12` term;
3. retention of `k_delta=ell delta` in the near-column phase;
4. the determinant trace-square and matrix-square terms in the Piola
   expansion;
5. both inverse-Laplacian derivative signs;
6. preservation of the order-minus-one Biot--Savart singularity by the first
   and second shape derivatives;
7. the exposed order loss if `h(y)-h(y')` is omitted;
8. the `delta^3 L` polyhomogeneous remainder;
9. the Hodge, local-graph, and order-minus-two Schur errors relative to the
   `delta^3 L^2` contour gap;
10. divergence of the old unweighted `delta^2/r_delta` shortcut; and
11. the eigenbasis commutator identity producing off-diagonal mode decay.

## What the verifier establishes

The verifier is an exact oracle for the algebraic signs, tensor scalings, and
resonant comparisons used in `derivation.md`. The analytic theorem additionally
uses the source nondegeneracy/center gauges, the explicit fixed-domain Hodge
decomposition, singular-integral estimates, and the matched front/side kernel
argument. Those are prose derivations, not replaced by the symbolic tally.

The result discharges 0052's graph-jet hypothesis for the frozen `p>=6` Cao
subfamily while retaining the earlier finite-core spectral inputs at their own
scope. It does not validate a nonlinear branch, stability, or particle/quantum
interpretation.

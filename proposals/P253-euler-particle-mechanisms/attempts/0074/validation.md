# Validation and exposing checks

## Activation and execution boundary

Central activation is preserved verbatim in `activation-schema.command.txt`,
`activation-schema.stdout`, `activation-schema.stderr`, and
`activation-schema.exit`; the exit file contains exactly `0` and stdout says
`WORKFLOW VALID`.  No production numerics were used.

## Exact symbolic oracle

The first execution of `verify_cylindrical_fiber.py` is preserved in the
`first-symbolic.*` receipts.  It exited `1` because the oracle attempted the
nonstructural SymPy substitution `subs(delta*n,k)`; the scientific identity was
not tested by that line.  The repaired original execution is preserved in
`symbolic.*` and exited `0`.

After the kernel normalization and cylindrical bracket ledger were added, the
repository interpreter

    /home/dan/substrate-framework/.venv/bin/python

ran the expanded oracle once.  The exact command, stdout, empty stderr and exit
are in `symbolic-v2.*`; exit is `0`.  The oracle independently checks:

1. scalar `div grad`, `curl grad=0`, and `div curl=0`;
2. helical masses `(k+delta)^2`, `(k-delta)^2`, and the scalar `k^2` row;
3. the `q_geom^(-1)` and `q_geom^(-2)` expansions;
4. cancellation of `(R/a)delta` in the half-density Newton kernel, leaving
   `sqrt(q_geom*q_geom')/(4*pi)`;
5. the transport and stretching connection signs by converting cylindrical
   coordinate components to orthonormal components; and
6. detection of the wrong curl-connection sign.

This finite symbolic oracle validates the exact algebra only.  It does not
validate an operator-domain, compactness, Fredholm, contour, or KKS claim.

## Analytic proof ledger

The load-bearing analytic checks are carried in `derivation.md`:

- Equations (17)--(20d) derive the exact normalized kernel, use a nested-collar
  parameter-elliptic resolvent identity for operator-norm convergence, prove
  integer-harmonic endpoint cancellation, and separate the remote-axis and
  exterior energy tails.  The `m=+/-1` zero-shift `K_0` channel is controlled
  by the exact compact-DA-curl zero-moment identity, not by a false arbitrary-
  source low-frequency bound. The proof is on the fixed-support `q`-space:
  Stokes gives the vector mean, differentiated Fourier divergence controls the
  transverse mean by `k` times the first moment, and the helical conversion
  identifies the precise scalar zero-shift channel. No uncontrolled curl
  potential is used.
- Equations (21)--(28) prove the fixed nonzero-harmonic ambient DA closure and
  the `k`-independent column domain.  No division by the vanishing edge
  coefficient occurs in the closure argument.
- Equations (28a)--(28b) expose why a raw small transverse Hanzawa coefficient
  is unbounded from the anisotropic transport graph.  Equations (29a)--(29i)
  repair it with the actual volume action--angle flow conjugacy, extend that
  same map through collar/exterior, rederive the metric Hodge comparison under
  the same pullback, and leave only a tangent graph-controlled principal
  difference.  This is the predicate supporting (30).
- Equations (31)--(36) construct the characteristic inverse with regular-center
  and three `H_0^3` edge traces, define the Fredholm essential-spectrum
  equality, and add the compact full-Hodge block.  They claim no equality for
  the full ordinary spectrum.
- Equations (40)--(49) prove the large-`|m|` Neumann exclusion and finite-`m`
  zero/nonidentity result.  The only `k=0` nonzero-harmonic kernel is the true
  translation pair; its `K_1` determinant moves off zero for small `k>0`.
- Equations (50)--(52a) select two distinct simple axisymmetric modes inside
  the resulting punctured whole-column hole. They retain the reviewed positive
  constrained-energy/Krein sign and the exact abstract E2-unit identity
  `abs(Omega)=1/abs(sigma)`, while leaving the absolute conversion from the
  Sturm amplitude to physical density/Fourier/KKS normalization open. The two
  contour projectors and their nonzero KKS signs transfer by (30), (38), and
  (39).

## Maximum earned verdict and exclusions

The proof establishes two isolated, distinct-frequency, positive-Krein
massive modes on sufficiently thin actual Cao rings.  It does not establish an
absolute physical KKS normalization for the source Sturm amplitudes, an
equal-frequency invariant analyzer doublet, noncommuting mode controls, a
nonlinear rotating branch, restoring stability, P2/LP2, or a particle/quantum
claim. Those exclusions are scientific boundaries, not inferences of
impossibility.

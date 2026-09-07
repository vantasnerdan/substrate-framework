# Validation

## Exact oracle

The repository interpreter executed

    /home/dan/substrate-framework/.venv/bin/python \
      proposals/P253-euler-particle-mechanisms/attempts/0088/verify_hessian_response.py

The first execution is preserved in `hessian-oracle-first.*`.  It exited `1`
because the checker had the wrong lower-left sign for the commutator of an
off-diagonal and diagonal anti-Hermitian matrix.  This was an oracle
implementation error, not a scientific counterexample.  The next execution,
preserved in `hessian-oracle-corrected.*`, fixed that predicate and added the
corrected cylindrical KKS residual-radius check.  After the supervisor
required the massive fixed-`k` continuation, the predicate was materially
expanded and executed as `hessian-oracle-massive.*`.  The near-axis route
then changed the predicate once more; its first execution is preserved in
`hessian-oracle-axis.*`, with empty stderr and exit exactly `0`.  The
Sol-High recovery added the KKS/dual-Riesz identities.  Its first execution,
preserved in `hessian-oracle-dual-riesz-first.*`, exited `1` because SymPy
structural equality did not simplify an otherwise equal symbolic recurrence
coefficient.  The repaired predicate tests the zero difference and its first
captured execution is `hessian-oracle-dual-riesz.*`, with empty stderr and
exit exactly `0`.

The final exact oracle checks:

1. `Omega(e,bar e)=i/nu` and `ell(e)=1` with the physical clock;
2. cancellation of the two differentiated metric rows at a double
   eigenvalue and the anti-Hermitian compression signs;
3. the noncommuting off-diagonal/diagonal matrix sign;
4. the normalized cylindrical KKS row retains `r=R*q` after half-density
   reduction;
5. the covariant Feynman--Hellmann sign;
6. the compact `(ell,d)=(0,0)` stabilizer identity;
7. the two-meridional-dimensional `H5` patch exponent seven;
8. the fixed-`k` identity `L1=C+2*k^2*s*cos(alpha)`;
9. the frame-rotation and `-i*k*s*cos(alpha)` terms in `B1`;
10. the `W/2` and `s*Omega'/2` curvature characters;
11. the complete order-`s` near-axis response coefficient proportional to
    `k_1/k_2-1`;
12. cancellation of the column `W` factor in the physical KKS covector;
13. reconstruction of the core adjoint row from the conjugate mode;
14. the one-datum singular-axis Frobenius recurrence; and
15. the regular harmonic-gradient ambiguity that makes the global
    dual-Riesz/Sturm selection load bearing.

The oracle explicitly scopes out the actual-Cao response remainder and a
uniform high-`N` scale.
It contains no discretization, quadrature, fitted coefficient, or soft
quantity, so the small-ratio numerical protocol is not activated.

## Analytic checks not delegated to the oracle

- Equations (1)--(5a) reconcile the README control sign with the canonical
  right-reduced orbit sign and carry `bar_omega=r*zeta*e_theta`, physical
  volume, normalized Fourier characters, and the half density.
- Equations (6)--(10) differentiate the physical energy--impulse functional
  on the orbit.  The second push-forward term is retained, and
  `E2(e,e)=0` is derived from Hessian invariance before realification.
- Equations (12)--(24) independently derive the Hessian-row dual,
  core-clock factor, covariant Feynman--Hellmann row, and the
  basis-invariant skew-Hermitian matrix statement.
- Equations (25)--(27) integrate both Euler response brackets through the
  decaying whole-space Hodge map.  Equations (27a)--(29) use the adjoint of
  an explicit primal row projection, so finite-row restoration cannot
  cancel the local numerator.
- Equations (30)--(32) prove the zero/zero compact displacement is a
  stabilizer.  This is a route refutation, not a global absence claim.
- Equations (33d)--(33p) regrade the exact cylindrical operators at fixed
  massive `k`, include the Cao odd first cell, pressure, frame rotation and
  right/left reduced resolvents, and reduce the first permitted curvature
  response to a radial source-mode integral.  Equations (33q)--(33x) derive
  the arbitrary regular-band DA coefficient, direct Sturm multiplication
  constituent, and exact fixed-pair cokernel condition.  They explicitly do
  not infer collective compactness across shrinking Sturm gaps.  Curvature
  nonidentity and the full graph `C1` response remainder remain dependencies.
- Equations (35a.10)--(35a.16) first cancel the apparent `1/W` orbit inverse
  in the KKS row, construct an extension-independent global velocity
  functional, identify it with the rank-one dual Riesz functional, and solve
  the regular-axis/interface/decaying-`K1` adjoint Sturm problem.  The exact
  recurrence and `A_1^sharp=rho_0 A_1` relation make the axis coefficient
  nonzero without a kinetic-Hessian identification.  Equations (35b)--(35k)
  then use that actual Hessian left row and both velocity convection terms.
  Leray pressure drops exactly; the full order-`s` coefficient is nonzero for
  `k_1!=k_2`; and a punctured regular annulus supplies an explicit compact DA
  seed.  Excluding rational ratios
  `1/2,1,2` makes every finite displacement row vanish on the reconstructed
  `xi`, so the correctly typed velocity identity
  `M_12(h)=integral w dot F_12^v` proves `gamma_12^col>0`; this local
  calculation alone is not its Cao transfer.
- Equations (35l)--(35s) then derive the source-specific fixed-`k` transfer:
  the reviewed graph/Riesz/Hodge intertwiners, direct convergence of both
  terms of the core-time Hessian `E2hat=E2phys/Omega_N`, exact action-flow
  pushforward of the displacement, and the smooth bilinear bracket estimate
  give `Mphys_(12,N)/Omega_N=Mhat_12^col+o(1)`.  The algebraic conversion
  is exposed by `ephys=ehat/sqrt(Omega_N)` and
  `ellphys=sqrt(Omega_N)ellhat`; only the scaled frame/dual is claimed to
  converge.  Since `Omega_N>0`, this proves nonvanishing for every
  sufficiently large member of the rational-ray exact-crossing sequence
  furnished by reviewed 0083/0089, but does not supply a physical
  `Y4`-normalized high-`N` control-seed bound.
- Equations (36)--(37a) derive the distinct two-off-diagonal alternative and
  identify the first autonomous retained-sector amplitude equation.  They do
  not count two phases of one static seed as two controls.

## Static claim-boundary check

`validate_attempt.py` checks the byte-frozen replayed README hash, activation
and final oracle exits, parses `result.yaml`, checks the physical KKS and
fixed-`k` phrases, excludes the stale leading-order and non-adjoint witness
phrases, confirms the active campaign boundary, and scans all owned files for
trailing whitespace.  Its command/output/exit are preserved separately.

The bounded 0096 clock correction changes no symbolic-oracle predicate, so
the unchanged oracle and static validator were not rerun for tally.  The
correction is checked instead by YAML parsing, stale-claim scanning, and
`git diff --check`; their first captured outputs are recorded in the bounded
correction receipt.

## Claim boundary

Validation supports the exact Hessian-row/dual-Riesz bridge, stabilizer
refutation, fixed-`k` curvature reduction, positive fixed-column off-diagonal
response, its source-specific transfer along the reviewed 0083/0089
exact-crossing path,
and response-matrix algebra.  It does not independently review that path,
prove `D_curv!=0`, give a uniform high-index normalized response, construct
two autonomous histories, establish two-sided gate control, P2/P4, a
particle, or a quantum or relativistic mechanism.

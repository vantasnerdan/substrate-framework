# P253/0105 author completion receipt

P253/0105 is complete at its author-stage route boundary.

The strongest exact result is a regular full-core classification.  With the
single globally matched Maxwell normalization and a positive constitutive tag
`chi=h(P)=zeta/lambda_0`, exact first-speed cancellation reduces to
`partial_s[s P_s^3 G_P]=0` for `G=h/h_P`.  A nonzero integration constant
forces `G=O(s^-2)` and the dipole `v=O(s^-1)` at a smooth nondegenerate center.
Regularity therefore makes `G` constant and `h=A exp(aP)`.  This is the unique
nontrivial regular constitutive kernel; it cannot vanish at a finite Cao free
edge, so the compact full-core fixed-ratio route is refuted.  The exact
power-law moment formula independently changes sign and, by analytic
continuation, cannot vanish on any open full-core band.  The polynomial
`P=1-s^2` check is recorded only as algebraic regression, not as a physical
Lane--Emden solution.  The existing strict-band cutoff is a different source
and remains unevaluated.

The fixed-material-leaf calculation also closes one positive and one negative
route.  Compact divergence-free displacements have zero weighted contour
normal flux, so the two unconstrained axisymmetric radial A/B source controls
from 0102 are not first-order tangents of the fixed coadjoint/material leaf.
On every regular finite-Cao band, however, a smooth target `T` has a compact
periodic material realization exactly when `integral J T d alpha=0` on each
contour, provided `chi'(I)` stays nonzero.  This establishes the raw
weighted-zero-mean material tangent.  The periodic primitive loses one normal
derivative: its bounded map is `C^(k+1,alpha)->C^(k,alpha)` (and
`H^(k+1)->H^k`) on a fixed regular band.  It does not solve the linearized
steady Euler--Maxwell/free-boundary map or restore its finite rows.

Exact cancellation remains `T_0 in Ran L_full`.  Adjoint-kernel annihilation is
necessary; it is sufficient only after closed range/Fredholm control or a
bounded right inverse.  Exact range and closure of range remain distinct.
The strict-band response, complete steady lift, curvature response, moving
defect, P2, P4, P5, P6, electron, neutrino, particle and parent conclusions
remain active dependencies.

Core SHA-256 values:

```text
e15c0dea81f5859b8ba374fddeb7be4365faa7bce039b6f328642fc181f871f4  README.md
db32c1c49a4196828d508c60c4613fdf35ecfb07fad2fa26e521e9eda9d485cb  derivation.md
e187d718d6735016069da4bf214427f90f8af43d1142c40b73f8d4fe82169f5b  source-audit.md
1033bbe1da359f954faf3dd5372df9c6a32bcc6fb0a45d30e8f2829c4efe7077  result.yaml
3d20ff269e90ddbc2e03614ff3eb46d86732db750c5f2cdddf32785d4639d10c  validation.md
2cd3dc678dc8349259f972438055fcb376ca6ba698c9877fc29d6ff3b7657a6b  verify_normalized_leaf_response.py
eb3a8d8971cc1b5e0b8e5bfb84a6f618b7aa2822b02403e855975f324fe776bd  src/substrate_framework/euler_two_label_lock.py
184a208cf74b773ba8eb5568b77c3176516910cd2daf9a44f89fd1297c84b811  tests/test_euler_two_label_lock.py
```

The final exact-v2 verifier derives eleven predicates and exits `0` with empty
stderr.  The focused-v2 public API suite passes twenty-one tests.  Static-v3
agreement and repository-v2 validation exit `0`.  The earlier static-v2
failure is preserved and diagnosed as a line-wrapping-sensitive prose scan;
it changed no science or implementation.  No production numerics or empirical
comparators were used.

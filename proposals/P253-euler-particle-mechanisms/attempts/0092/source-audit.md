# Source and authority audit

## Reviewed framework inputs

This attempt uses earlier campaign artifacts only at their independently
reviewed scopes.

- P253/0055, as finally reviewed by P253/0061 (`review.md`
  `df54d98d40f413feb6c7decdbea99ccb0a519016d759e1d8cfa28eaa81f3b4ed`,
  `verdicts.yaml`
  `6b8ca084cde2037d526ab95d6e0bc422433b842e7f458e4f74d58f4396bde4b3`),
  supplies the bare-Euler similarity weights and the scope-limited failure of
  topology/conservation-only universal action selection.
- P253/0056 (`review.md`
  `bedc7acf7cc4d59ecf3212da74eb9d3e693100080b4dc846567fe0144915fa4a`,
  `verdicts.yaml`
  `b30efe3b48b9e7ab5f671debf2f8d98519be50da22b0ea6616b28158c39a146e`)
  supplies the corrected physical action `J=sum |z_a|^2`, phase moment map,
  and conditional Schwinger--Hopf/prequantization boundary.
- P253/0068 as independently reviewed by P253/0071 (`review.md`
  `50e1cb63d6c115a92c3fe6720eaeae0b707c5b81871db663a7a804a9a3127155`,
  `verdicts.yaml`
  `fd59cc3209460f85aa97cacf66ea5a88d6d5cd9ce8a1630b8e37cf1f31d0f1c0`)
  supplies the classical transported tag current and the normalization
  invariants `epsilon*mu` and `g_tag^2/epsilon` for the transported-tag
  normalization.
- P253/0080 as independently reviewed by P253/0084 (`review.md`
  `3fb7d2a69489668f6168d7f57eb9b4e86e35621064dd5f55172ca24270279aac`,
  `verdicts.yaml`
  `7793688a23156d1b015ed82505f5e48b88223e19e8a36563c2f57ef63f0a74c7`)
  supplies the finite-window charged Cao branch only after circulation, mean
  radius or impulse, tag leaf, coupling, and Maxwell constants are supplied.

No claim from active 0088 or the 0089/0090/0091 response and radiation chain
is consumed.  The result also does not use an empirical elementary charge,
Planck constant, or particle comparator.

## Exact mathematics derived here

The dimensional and field-rescaling identities are direct substitutions in
the declared Maxwell action.  The compact-phase gauge cancellation is derived
from (6)--(8), rather than inferred from phase periodicity.  The distinction
between a periodic coordinate and its real classical cotangent momentum is
the elementary global structure `T^*S^1=S^1 x R`.  A global compact action on
the phase has an integer character weight; this chosen representation label
does not discretize the real cotangent momentum.

The compact-phase calculation uses the explicit number-density convention
`Omega_phase=S_0*integral(delta n wedge delta theta)` and
`Q=g_0*m*N_phase=kappa_m*J_phase`, where
`J_phase=S_0*N_phase` and `kappa_m=g_0*m/S_0`.  The transported-tag coupling
`g_tag` and the charge-per-action coefficient `kappa_m` are unrelated unless
a separate same-carrier map derives their relation.  Under field rescaling,
both physical charge coefficients scale by `1/a`, while the imported action
coefficient `S_0` does not.

The spatial topology uses the standard classifications

    principal U(1) bundles over X  <->  H^2(X;Z),
    H^2(S^3;Z)=0,  H^2(S^2;Z)=Z,  pi_3(S^1)=0.

These groups are computed directly by the cellular homology of spheres.  No
external theorem is used to turn their integer into electric charge.  The
derivation instead exposes why the `S^2` Chern class is magnetic and why a
smooth bundle extending over `S^3` has no such class.

The BF discussion is a candidate audit, not a positive compact-Euler theorem.
The source-side fact used is exact: for a smooth whole-space Euler velocity,
the vorticity two-form is the real exact form `d u^flat`.  It does not acquire
large compact two-form gauge transformations or integral periods merely by
being materially advected.  Compactifying it or quantizing a BF level is
therefore declared new structure.

## Evidence boundary

The symbolic API and verifier check finite algebra, topology ledgers, and the
rank counterdirection.  They do not establish a new compact matter field,
quantum representation, BF theory, carrier response, probability rule, or
exchange mechanism.  The route verdicts are about the tested classical tag,
phase, smooth whole-space bundle, and minimal common-phase coupling; they are
not an exhaustion claim over all charge/action mechanisms.

# P251: repaired proposal artifacts and conditional mode tools

This is the bounded artifact scope of PR #199, linked to the pre-existing
scientific intake #198. On 2026-09-05 the repository owner explicitly
approved separating the unfinished full-Euler research and requested this
repaired PR be merged. The continuation is tracked in
[#200](https://github.com/vantasnerdan/substrate-framework/issues/200),
with its full progress, scientific status, intended value, usage guide and
preserved research-branch history.

## Useful contents

- `src/substrate_framework/rankine_modes.py` provides conditional isolated
  Rankine-vortex velocity equations, pressure-matching determinant and a
  caller-precision Bessel residual. `tests/test_rankine_modes.py` checks the
  actual Cartesian Euler residual, boundary elimination and domain limits.
- `verify_cst001.py` through `verify_cst007.py` preserve scoped frame,
  constitutive, dispersion, coherent-response and source-provenance checks.
  They are proposal evidence, not accepted scientific claims.
- Attempts0029–0034 repair the Euler/Bessel correspondence and preserve a
  closed point-vortex angle action, prescribed-flow elliptic-patch result,
  corrected mutual-energy kernel and conditional long-wave field map.
  Their assumptions and approximation boundaries remain explicit.

## Reading historical attempts correctly

Attempts are append-only historical records. The newest correction for a
statement controls its current proposal scope; an old green output or
historical terminal claim does not override it. In particular:

- The original CST002 velocity signs are correct;0029 withdraws the earlier
  reliance on0019's reversed-Coriolis residual.
- The corrected m=1 Bessel constant is1/4-EulerGamma (0031).
- Straight-tube tension does not derive alpha=L_v T/6. Exact free director
  rotation cancels in Green–Lagrange strain.0028's rate-quadratic coefficient
  has inertia units and does not repair the missing locking potential.
- Zero coherent response does not imply zero quadratic fluctuations.
- The elementary Beltrami regression is not a knotted-tube construction.

## Merge, authority and unfinished objective

This PR lands reusable, repaired proposal artifacts and a conditional API.
It changes no accepted claim, release or generated canonical record. It does
not establish the full smooth-Euler-to-Cosserat continuum. Neither scientific
success nor exhaustion is asserted. The proposal's active scientific graph
is retained intentionally; the owner's explicit scope separation makes the
bounded artifact merge independent of solving that research frontier.

Later research attempts and branch-local promotions are preserved on
`research/pr199-completion`, not included in this PR. Issue #200 describes
their precise authority and open joint dynamics/current/geometry construction.
Issue #198 remains open as the original scientific objective.

## Validation and review

The scientific source boundary is commit
`3626fbfb748d930cd22fd0b685a49a399c75bd06`. The recorded validation in
`attempts/0029/validation-receipt.md` matches its executable source tree:
seven direct Rankine tests, fixed structural checks, and the named scoped
scientific oracles. `attempts/0029/independent-review.md` records substantive
review and its completed documentation correction. This wrap-up adds only
documentation and lifecycle metadata; it does not stale those scientific
receipts. Final review checks the actual head, current main compatibility,
metadata/schema and diff cleanliness without a full-suite replay.

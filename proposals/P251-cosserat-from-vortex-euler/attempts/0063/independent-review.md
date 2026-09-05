# Independent review: coherent-population and full-fluid limits

Reviewer: `/root/smooth_core_review`, distinct from the author of 0063
and its consumer addition. Date: 2026-09-05. One bounded exact-limit
review under AGENTS.md and the physics skill. The task is the two new
coefficient/coordinate limits, not new microscopic coefficient provenance
or a replay of the earlier numerical dispersion evidence.

## Decision

Established as stated: removing the f-scaled coherent spin population
while retaining a positive background shear leaves the incompressible
Navier--Cauchy displacement sector. Removing that background shear as
well gives the neutral transverse displacement sector of linearized
Euler at fixed density. At either zero-population endpoint the spin
coordinate has zero action weight and is removed before dividing by
its vanished inertia. No load-bearing correction is requested.

These are conditional limits of the supplied action coefficients. The
identification of f with a microscopic retained population and of
mu_background with a particular stationary background remains a separate
construction. The algebra does not prove either coefficient dependence.

## Exact determinant and unreduced coordinate limit

I read the frozen README, the new CST005 function and its registration in
main, the surrounding operator definition needed to identify its inputs,
and the saved two-check passing receipt. The existing receipt was reused.
The tests use the same canonical Fourier stiffness as the original
consumer; they do not insert an independent target determinant.

Set alpha=f a0, j=f j0, c_s+c_a=f c0 and
mu=mu_background+f mu_correction. Directly in the transverse action,

```
M=diag(rho,f j0),
K_UU=[mu_background+f(mu_correction+a0)] k²,
K_UPhi=-2 f a0 h k,
K_PhiPhi=f(4 a0+c0 k²).
```

Thus the determinant has an overall factor f. After dividing for f>0
and taking the coefficient limit, its exact limit is

```
(rho omega²-mu_background k²)
                          (j0 omega²-4 a0-c0 k²).
```

The coupling product is order f² and therefore vanishes from this
normalized limit. The finite mu_correction also drops out. These are
exact fixed-k polynomial facts and do not assume separated roots at a
possible branch crossing.

At f=0, inspect the unreduced matrices instead: the spin row and column
have zero mass and potential weights. The remaining physical equation is
`rho omega²=mu_background k²`. The finite ratio root in the normalized
polynomial is a limit of the disappearing spin sector, not an additional
mode of the zero-population endpoint. Removing the degenerate coordinate
is what makes that endpoint interpretation faithful to the action.

Setting mu_background=0 then yields rho omega²=0 for transverse
displacement. Density has not vanished and no arbitrary replacement
spin inertia has been introduced. Pressure still imposes displacement
incompressibility; no longitudinal compressional Euler mode is inferred
from this two-by-two transverse calculation.

The older all-coefficient L_v scaling therefore remains a valid
conditional regression of a complete background-removal path. It is
not equivalent to removing only marked coherence at fixed background
shear. The new check names that distinction instead of rewriting the
old limit to conceal it. A full six-field endpoint similarly removes
all spin action weights before retaining a physical spin coordinate;
the new assertions here are the displayed transverse sector.

## Provenance and disposition

- `README.md`: `567c32c04d64e98112de0b21f13784300ac17a178a44a61b7418c22b9cd2a31d`
- `verify_cst005.py`: `b57a199eea07e3933c818e59faf3a1c4d3d9423eb50a15d58cf29fc298a11dd3`

Acceptance is recommended for these exact conditional coefficient and
coordinate-limit identities. No numerical eigenvalue or small-ratio
inference is used. Microscopic provenance and parent completion are
not certified by this attachment.

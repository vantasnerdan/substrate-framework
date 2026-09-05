# 0264 bounded correction review: periodic supplier

Reviewer: herdr optical-review pane w3:p3.  Scope is limited to the three
questions raised in the original 0264 review: Bloch Kelvin identity, uniform
homogeneous-kernel construction, and full-projector/carrier-uniform bounds.
`0262/current-transfer.md` remains excluded.

Pinned source `periodic-supplier.md` SHA-256:
`6cb6175df0ea5bc2c4ae0ca7ccb5d224c1f516b4cbe3efa34bd287703ecc69dd`.
Pinned check receipts: `bloch-supplier.stdout`
`43152f87ea895802db780e3957973cbcd09873791cf0ce532be6dff42535fe82`,
`bloch-supplier.exit`
`9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

## Assessment

The supplied symbolic check and equation (1) correctly retain the Bloch
product terms.  With `div_K xi=0`, the two `iK` contributions cancel and
`curl_K(xi cross omega)=lambda[u,xi]_K` on the constant-curl support.  Since
the removed part of `P_K` is parallel to `q+K`, the full periodic projector
has the same `curl_K`; this is the required pressure-inclusive identity,
not a local-projector approximation.

The homogeneous constraint repair is valid and does not require an
inhomogeneous right inverse.  Imposing all finitely many K-coefficients at
once gives `A_N c=0`; with M+1 disjoint bumps, a unit vector in the kernel
exists.  Disjoint support and fixed-sign weight give the claimed uniform
positive P/E margins for every unit kernel vector, while the C^s bounds are
uniform.  This directly addresses the actual homogeneous requirement.

The Fourier multiplier estimate, explicit removal of the q=0 force mode, and
integration-by-parts split provide the needed full-projector pressure-tail
bounds.  Moving derivatives by periodic integration by parts in (9) gives
uniform baseline cross-form derivatives rather than incorrectly deleting
them.  After the stated phase/energy whitening, the estimates imply a fixed
K-neighborhood independent of the large auxiliary carrier.

## Verdict

`route_verdict: established as stated` for the explicit periodic auxiliary
supplier and its carrier-uniform finite-K form bounds, under the frozen
periodic-background, support-separation, and invariant-gap hypotheses.  The
only remaining mechanism is outside this review: the current/actual-
observation consumer and construction of the periodic background itself.
No correction is required to the homogeneous-kernel argument, and no
inhomogeneous right-inverse condition should be added.

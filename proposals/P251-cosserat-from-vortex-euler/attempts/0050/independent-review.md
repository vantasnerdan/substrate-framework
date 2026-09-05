# Independent review: finite-coherence common-angle action

Reviewer: `/root/smooth_core_review`, distinct from the author of 0050.
Date: 2026-09-05. One bounded scientific pass under AGENTS.md and the
physics/small-ratio skills. The reviewed claim is the specified local
coherence-return affine action, not the complete parent continuum theorem.

## Decision and scope

Established as stated: the actual EPS relative orbit of 0048 admits the
specified finite-coherence common-gradient family, with finite complete
Euler energy and KKS coefficients. The full exterior induced-velocity
contribution and momentum elimination are retained. The additional
moment-projected first-difference cages give a finite analytic threshold
for positive gradient curvature on this same field.

No load-bearing scientific correction is requested. This result supplies
the previously missing finite common-gradient action at its declared local
affine scope. It does not assert a finite Hessian for the unrestricted
field `exp(ik.x) K`, a universal modulus determined by R alone, or a
global gluing/homogenization theorem for adjacent microscopic realizations.
Those stronger statements are not hypotheses silently used in this proof.

## Evidence and geometric premise

I read the README, complete `coherence-action.md`, `verify.py`, and its
saved 21/21 receipt. The analytic proof is the strongest practical oracle;
the symbolic checks expose its geometry, signs and elimination identities.
The existing receipt was inspected and reused, not regenerated merely to
increase the check count. Inputs are the already reviewed actual-EPS
relative orbit, compact cages and physical jets of 0045/0048.

The identity `curl(-|y|² e/2)=e cross y` and product rule give exactly

```
Z_j=chi eta (y_j K+e_j cross A_K)
                       +y_j grad(chi eta) cross A_K.
```

Thus the return-shell term is necessary, and is present. Each Z is smooth,
compact, divergence free and identically zero on a neighborhood of all
physical jets and original conjugates. Consequently the tracked core
angles remain B+q and B-q; this follows from the actual generator, not
from a subsequently assigned angle label.

The macro fields and their derivatives are coefficients in a local
microscopic family. They are not differentiated again with respect to the
microcoordinate. The complete cutoff profiles and the neighbor-gradient
tie are explicit kinematic premises. This is a genuine prescribed
Cauchy--Born-type return geometry, with its energy evaluated from Euler;
it is not an independently inserted elastic energy. Finite flows of the
global rotation plus compact directions remain a rotation composed with
a compact volume-preserving rearrangement. No infinite-domain plane-wave
limit is invoked to establish this local family.

## Complete Hessian and the exterior fluid

For a compact direction, differentiating transported vorticity twice
retains `curl[xi cross curl(xi cross omega0)]`. With the same decaying
Biot--Savart inverse and Beltrami identity this gives

```
H(xi,zeta)=rho integral (v_xi.v_zeta-v_xi.curl(v_zeta)/lambda).
```

The sign of the second term follows from the scalar triple product and
curl integration by parts. The induced compact-vorticity velocity is
O(r^-3); the EPS background is O(r^-1). The corresponding cross energy
and integration boundary terms therefore converge. Uniform rotation
retains exactly the centered-ball relative normalization of 0048:
H(K,K) and H(K,compact) vanish by symmetry and the Beltrami mixed identity.
There is no new finite-part or infrared subtraction.

Outside a ball containing the compact vorticity changes, curl v vanishes,
but v generally does not. Splitting the exact integral consequently gives
the stated positive exterior Gram, with O(L^-3) tail. This is an important
physical distinction from imposing a reflecting boundary or discarding
the exterior pressure reaction. The quoted material-domain pressure-work
identity is the Euler conservation law; it does not assert that every
constrained affine path is an unconstrained Euler solution.

KKS integrals with one K leg and one compact leg are finite. The annular
Z has disjoint support from the compact momenta, so its momentum KKS
row o_p vanishes exactly in this construction. A possible K--Z pairing
is not set to zero by that support argument; it is retained in Omega_xx.

## Full reduction and meaning of optical curvature

Expanding the full `(H-i nu Omega)` momentum Schur complement gives the
three displayed K, G and M blocks, including both cross terms in G and
`M=Omega_xp P^-1 Omega_xp^*`. No diagonal P or energy orthogonality is
assumed. The common gradient stiffness is H_ZZ minus the complete
momentum cross Gram, and the general gradient inertia is the nonnegative
KKS momentum Gram. For the chosen annulus the latter is zero by support,
not by deleting gradient inertia from the formula.

Independent elimination in the two time-reversed realizations cancels G
and preserves K and M. It does not, by itself, remove odd-in-k terms in
K. Spatial reflection pairing, if used for that purpose, is a separate
explicit ensemble operation, as the proof states.

The physical map gives the optical column `t=(1/a,1/a)` at beta=0.
Its even constant-mass gradient coefficient is therefore
`t^T(K2-nu0² M2)t`, with the negative inertia correction and every cross
term retained. This is the constrained optical-column curvature. It is
not automatically the eigenvalue curvature of a different unconstrained
two-configuration model with residual odd-in-k mixing. The full raw
matrices remain the appropriate outputs for the parent's joining. Under
the stated affine-cage restriction, the displayed optical coefficient is
the correct one.

## Positive same-field extension

Projecting each new cage with the fixed eta0 kills its K moment exactly.
Support separation kills its other retained KKS pairings; common eta0
attachments add no mutual symplectic term. Thus this extension changes
energy but not the retained KKS matrix or gradient inertia. The projection
coefficient is bounded independently of the new carrier frequency.

The negative-helicity cage diagonal grows linearly with carrier frequency.
Moving curl onto each fixed comparison velocity makes energy cross terms
uniformly bounded. Disjoint new carrier supports remove the local
helicity cross; their whole-space projected kinetic cross and shared
eta0 attachments remain bounded and are included in M_off. This argument
does not incorrectly infer energy orthogonality from support separation.

For a unit gradient direction, Cauchy--Schwarz gives the stated V_L bound
without requiring the new velocities to be orthogonal. L_p bounds the
new momentum cross row; expanding its Schur square gives precisely
`||P^-1|| (2 H_p L_p+L_p²)`. Together with the base operator norm,
off-diagonal row bound and the Z cross estimate, this proves R_total is
finite and independent of the new carriers. Selecting each positive
diagonal above it gives a strict uniform bound. The six-cage extension
uses the same argument for the full B/q gradient matrix. Neither fitted
moduli nor an observed small eigenvalue enters the selection.

The sign claim is controlled by exact identities and explicit analytic
norm bounds, not a discretized Hessian near a numerical floor. Accordingly
the small-ratio prescriptions are satisfied by those strict bounds;
mesh extrapolation or jitter would not strengthen this particular proof.

## Frozen hashes and disposition

- `coherence-action.md`: `32fec4cd1896493bfcd5fce9193ca3cbb1d6ca83a8b25116944900063aea4bdd`
- `verify.py`: `b61a0ced150b76f680859e1a4302b869faae92a3fb4b1444615ea633785d565b`
- `stdout.txt`: `6abfdf9e8190462a5fe656505a52bfe08b344185ab1f2b0bbeed4ca9da546133`

Scientific acceptance is recommended for this exact conditional local
affine-action result. This review changes no accepted registry entry and
does not certify the parent's coupled isotropic continuum, normalization,
or completion obligations.

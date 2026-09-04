# Independent review: physical angle action directly in an EPS tube

Reviewer: `/root/smooth_core_review`, distinct from the author of attempt 0045
and from the parent who helped generate the candidate. Date: 2026-09-05.
This is one bounded scientific transaction under AGENTS.md and the physics
and small-ratio skills. It does not reopen the earlier 0036, 0040 or 0041
reviews.

## Decision and exact positive statement

The theorem is established as stated. Every actual smooth constant-lambda
EPS tube satisfying the explicit nonzero-interior-vorticity hypothesis
supports the constructed compact, divergence-free two-generator Euler orbit
action. One coordinate is the physical tilt angle of the vorticity direction
at the selected interior point. The complete restricted energy Hessian is
strictly positive and the physical KKS form is nondegenerate for a finite
carrier obeying the stated explicit norm bounds. Eliminating the same-fluid
conjugate coordinate gives positive angle inertia and stiffness from that
single action.

No blocking finding, additional hypothesis or scope reduction is requested.
The statement concerns its declared constrained quadratic microscopic
ensemble. It does not claim unrestricted Euler invariant two-dimensional
dynamics, a common rotor, translational constitutive closure or completion of
the parent continuum theorem. Those are not additional acceptance criteria
for this microscopic result.

## Frozen evidence and applicability

The reviewed files are `README.md`, `proof.md`, `verify_direct_eps.py`, the
saved 17/17 scientific success `stdout.txt`, and the initial implementation
diagnostic. The base campaign provenance is checkpoint 3626fbf and its
v0.171.0 accepted release. No accepted registry statement changes in this
review transaction.

The strongest oracle is the analytic finite-norm proof. The symbolic script
corroborates its exact vector identities, signs, physical tilt jet and
off-diagonal canonical elimination. It does not claim to evaluate the
integrals of a generic EPS solution numerically. The ledger API correction
changed no equation; the existing success receipt is reused after source
inspection.

The imported [EPS existence theorem](https://arxiv.org/abs/1210.6271) supplies
an actual smooth stationary constant-lambda Beltrami field and invariant
vortex-tube boundary. Because u and omega are collinear, this is also a
material boundary for the background. A point with nonzero vorticity is an
explicit hypothesis; continuity gives the positive-component ball required
by the construction. The proof uses that actual field on R³ and the same
Euclidean Leray operator throughout. It does not rely on approximating a
compact-domain nonlocal projection.

## Physical angle, support and orbit conventions

The local core observable is genuine vorticity rotation. Where the inner
bump equals one,
`curl(-|r|² e/2)=e cross r`; its generator vanishes at x0. Thus the full
vorticity variation `omega dot grad xi-xi dot grad omega` at that point is
exactly `e cross omega(x0)`. With unit e perpendicular to omega(x0), the
coordinate is an infinitesimal physical tilt angle, not a frame label or an
unobservable rotation of a circular cross-section. The claim correctly uses
the local direction jet rather than rigid rotation of a complete finite
cross-section.

All generators are smooth curls with compact support strictly inside D.
Their flows preserve volume and are the identity in a collar of its boundary,
so they preserve D. Cage and core supports are disjoint. A cage velocity
projection has nonlocal tails, but its vorticity variation is the local
`curl(xi cross omega)` and vanishes near the core point. Consequently adding
the cage to the q-generator does not contaminate that physical tilt
normalization.

The Hessian includes the second advected-vorticity variation. At
`omega=lambda u`, pairing that variation with u gives
`H=rho integral(|v|²-v dot curl v/lambda)`. The positive KKS convention is
fixed by `Omega(u,eta)=dE(eta)`. Hence a pairing `B dq wedge ds` corresponds
to `B s qdot-H2`, with the sign used in the proof. Compact force fields have
square-integrable Leray velocities and compact curl, so every displayed
Hessian and KKS integral is finite. For the decaying EPS fields the cross
terms defining energy differences are finite as well; finite total
background kinetic energy is not being assumed.

## Independent KKS calculation

In the rotating basis `(p1,p2,ez)`, write
`grad phi=g1 p1+g2 p2+gz ez`. The two compact generators have components

```
xi1=(phi,-gz/k,g2/k),
xi2=(gz/k,phi,-g1/k).
```

Their cross product therefore has horizontal components
`(-phi g2/k+gz g1/k², phi g1/k+gz g2/k²)` and vertical component
`phi²+gz²/k²`, reproducing every cutoff term in the authored identity.
Integration by parts gives

```
integral omega dot [phi J grad_perp phi]
 = -1/2 integral (curl omega)_z phi²
 = -lambda B0/2.
```

Thus `B/rho=(1-lambda/(2k))B0+T/k²` with the stated sign. Compact support
removes boundary terms, and Cauchy--Schwarz gives `|T|<=Tstar`. With k and
lambda of the same sign, `K>=|lambda|` supplies the one-half principal
margin and `K²>4Tstar/B0` supplies `B>rho B0/4`. Disjoint supports make the
added core contribution to this pairing exactly zero, not merely small.

## Finite-carrier projection and complete positive Hessian

The Leray symbol estimate is valid on the whole Fourier domain. When
`|eta|<=K/2`, the carrier-shifted vector has norm at least K/2 and its
direction differs from the axial direction by at most a constant times
`|eta|/K`; the factor four is generous. Outside that ball, the operator norm
of the difference of the two orthogonal projections is at most one, already
below the stated bound. The single undefined direction at zero frequency is
a measure-zero issue for the whole-space L2 multiplier.

Plancherel now bounds the F0 projection error by
`4||grad F0||_2/K`. The F1 term contributes `||F1||_2/K`. The curl estimate
does not multiply the entire error by an unbounded carrier: since curl P
equals curl, the F0 curl error is the curl of the vertical component
`(I-Pi)F0 exp(ikz)`, whose fast derivative is zero. The F1 curl has a bounded
fast term and a slow term divided by K. This establishes both authored
bounds with the displayed D and E. Real-part norms are bounded by the complex
amplitude norms; their use introduces no missing factor of two.

The principal norm and helicity are exact for every real combination of the
two carriers. In particular the amplitude-gradient contribution to helicity
vanishes pointwise, and `curl p=-k p` gives
`integral v0 dot curl v0=-k A|t|²`. The principal Hessian is therefore
`rho(1+K/|lambda|)A I` with the correct favorable sign.

Expanding the difference between the full and principal Hessians yields the
three helicity errors
`<w,curl v0>`, `<v0,curl w>` and `<w,curl w>`. Substituting the two norm
bounds gives respectively the constants
`D sqrt(A)+DG`, `sqrt(A)E` and `DE` after K>=1. The kinetic errors give
`2sqrt(A)D+D²`. These are precisely the terms in C0, so the remainder is
bounded independently of K, not merely denoted asymptotically without control.

For core attachment, self-adjointness of curl moves its derivative onto the
fixed core tangent. This makes the core/cage cross energy O(1), bounded by
the displayed CR term, despite the increasing carrier frequency. Adding the
possibly negative core diagonal costs at most `|HR|/rho`. Hence Ctotal bounds
the full two-coordinate form. The stated finite choice
`K>|lambda| Ctotal/A` leaves a strictly positive lower margin, while the
previous KKS threshold preserves B. All constants are finite and depend
only on the declared background, cutoffs and physical projection; none is
selected from a desired frequency or elastic modulus.

## Off-diagonal reduction, scope and promotion

The full mixed Hessian is retained. Eliminating s gives

```
L = B² qdot²/(2h22)
    -(h11 h22-h12²)q²/(2h22)
    -B h12 q qdot/h22.
```

The last term is an exact time derivative because the quadratic coefficients
are fixed at the base state. Positive definiteness gives h22>0 and det H>0;
the nonzero physical KKS pairing therefore yields the stated positive
inertia, stiffness and frequency without assuming a diagonal example.

The cage is a prescribed part of this microscopic constrained family, and
its q-deformation is tied to the core deformation by that explicit generator.
The theorem does not claim that this constraint is an invariant reduction of
all Euler motions. Within its stated family, the complete same-fluid orbit
energy and symplectic form include its reaction before elimination. No
separate strain reservoir or mechanical mass is introduced.

The small-ratio prescriptions are met through exact analytic inequalities:
the complete Hessian and physical symplectic rank are specified, projection
errors have finite norm bounds, and the positive margin has a derived finite
threshold. No discretized soft eigenvalue, unmeasured numerical floor, fitted
modulus or convergence guess enters the result. Additional mesh or eigenpair
tests would not strengthen this exact oracle.

- Verification: exact analytic proof with symbolic-verified identities;
  no generic EPS integral is asserted to have a numerical value.
- Review: audited; scientific acceptance recommended for this microscopic
  theorem as stated.
- Compatibility: compatible smooth-Euler construction on an actual stationary
  EPS tube, under its declared constrained ensemble.
- Epistemic: the positive direct-EPS angle-action route is established;
  the parent common-rotation, translation and coarse-graining work remains
  active.

Registry incorporation, importable API extraction and affected-consumer
replay belong to the parent's promotion transaction. This review neither
changes the release nor declares those integration steps complete.
Correction check: not needed; no correction was requested.

## Content-addressed boundary

The reviewed SHA-256 identities are:

```
4214e4c0f40f0e6b89698bebe739cbdb89e4ec3eb8168f1df20e7f09447640ec README.md
161d30b976dbb73cfb52ff5860ed59e708eadb359673056a23955258a5d16674 proof.md
04af282f49520dc2f44fd6836f934d1dcf596ca6861e15f599168d28c33cda1c verify_direct_eps.py
8b9ee10b71bb4fca85f47130cfde1dc715f8d4b1aed839da21aa57310fd97f36 stdout.txt
```

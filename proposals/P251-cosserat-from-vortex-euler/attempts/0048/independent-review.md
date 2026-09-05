# Independent review: common rotor on the EPS relative orbit

Reviewer: `/root/smooth_core_review`, distinct from the author of 0048.
Date: 2026-09-05. This is one bounded scientific review under AGENTS.md and
the physics/small-ratio skills. It assesses the uniform common/relative
rotation construction, not the subsequent spatial modulation of its global
generator or the complete continuum theorem.

## Positive decision and precise scope

The relative-orbit common-rotation construction is established as stated.
For an actual EPS field chosen from its finite spherical-harmonic
construction, centered rotation-invariant balls define the relative energy
and angular impulse consistently. The selected global rotation and compact
shape/physical-angle directions have exact nondegenerate KKS blocks, zero
common-rotation energy entries and a positive full compact 3 by 3 Hessian.
Eliminating both compact momenta gives positive joint inertia, the computed
gyroscopic connection and positive relative restoring energy.

The declared time-reversal pair, with coherently tied physical angles and
independent conjugate shapes, cancels the odd connection while retaining
the positive quadratic action. Its field map has a finite positive
denominator by the additional proved threshold, not by positive definiteness
alone. The resulting J_beta is a positive **uniform cage-angle inertia**.
Its interpretation as spatial cage-gradient inertia requires the parent's
separate affine/spatial joining; that interpretation is not independently
proved by this child result.

No blocking finding or scientific correction is requested. The global
direction is explicitly outside the finite-L2 displacement space, and the
relative orbit is a specified constrained ensemble. No unrestricted Euler
invariant finite-dimensional dynamics or spatially varying global rotation
is inferred from this review.

## Frozen evidence and source bound

The reviewed evidence is the README, `direct-rotor.md`, `verify.py`, and the
saved 20/20 receipt. The earlier 18-check output is preserved; the added
checks expose the independent-momentum ensemble ordering. Existing receipts
were inspected and reused. Base provenance is checkpoint 3626fbf and the
campaign's v0.171.0 release; no accepted registry entry changes here.

The proof of [EPS Theorem 8.3](https://arxiv.org/abs/1210.6271), specifically
the finite spherical-Bessel/harmonic sum and its curl-polynomial conversion
on pages 56–57, was inspected directly. Physical angular differentiation
preserves the finite harmonic space and commutes with that conversion.
It therefore retains O(1/r) decay. Moving the rotation center adds a
translation derivative, which has the same decay. This supports the stated
bound on vK; it is stronger than, and is not inferred merely from,
`|grad u|=O(1/r)` multiplied by a growing rotation generator.

## Relative energy, moment and KKS extension

A compact vorticity rearrangement has a compact vorticity difference with
zero integral, obtained by differentiating along its compact isotopy.
Its velocity difference has the O(r^-3) far field. Thus its cross energy
with the O(r^-1) EPS background is absolutely integrable. For pure rotation
about the chosen axis and center, the reference balls are invariant and
the energy integrals cancel at every radius. Combining the two operations
therefore yields the stated finite relative energy, exactly independent of
the common angle. No freely chosen finite counterterm is introduced.

The rotation tangent obeys `curl vK=lambda vK`. The pure common Hessian
vanishes both by relative-energy symmetry and by its pointwise Beltrami
identity. For a compact tangent the mixed energy integrals converge, and
the curl integration boundary term is O(R^-2); moving curl onto vK proves
the mixed Hessian is zero. The KKS pairing of K with a compact generator is
a finite compact integral. Only one global rotation axis is included;
its bracket with a compact generator is compact, so the local closedness
calculation does not require an undefined pairing of two global rotations.
This supplies the local selected symplectic action without treating vK
as an L2 vector.

The relative angular impulse uses the same invariant balls. Its pure
rotation contribution cancels, while its other vorticity contribution is
compact. Since `curl(|r|²e)=-2K`, integration by parts gives the moment-map
sign `dJrel(eta)=Omega(K,eta)`. Also `K cross omega-vK` has zero curl on
simply connected R³ and pairs to zero with compact divergence-free eta.
Consequently `l(eta)=-rho integral vK dot eta`, and the selected
`eta0=curl(chi curl vK)` has the exact strictly nonzero moment
`l0=-rho integral chi |curl vK|²`. Analyticity/nonvanishing permits its
support away from the two physical jets.

## Physical jets and exact compact moment projection

The two compact rotation potentials supply opposite rigid rotation jets
about the same physical axis at the two tracked core points. Common spatial
rotation rotates their material orientation jets together. Hence their
observable angles are B+q and B-q. This statement follows the tracked
material points; it does not confuse fixed-coordinate advection with a
vorticity rotation. Every subsequent moment correction is supported away
from those jets, and projection tails have zero vorticity there.

Using the fixed eta0, rather than a high-frequency body direction, gives

```
r=(1-l(A)/l0)eta0+A,
Q=q0-l(q0)eta0/l0,
S=s0-l(s0)eta0/l0.
```

Disjoint supports remove the cross KKS terms exactly. The resulting blocks
are `Omega(K,r)=l0`, `Omega(Q,S)=c`, with all cross blocks zero. Thus the
common moment is fixed exactly and all four directions are independent.
The source-dependent ratios defining the compact corrections remain
unchanged upon reversing the background circulation sign.

## Full positive Hessian and the finite map threshold

The leading cage matrix is diagonal and grows like
`rho(1+K/|lambda|)diag(A_b,A_i,A_i)`. The projected velocity errors have
the 0045 D/K and bounded-curl estimates. Fixed attachment contributions
are controlled by their velocity and curl norms; moving curl onto the
attachment makes their cross energy bounded independently of K.
The moment corrections have uniformly bounded coefficients because the
oscillatory compact moment integrals admit the stated integrations by parts.
Summing these bounds yields the displayed C_H. The full compact matrix,
including all projected cross entries, is therefore positive above the
finite threshold; energy orthogonality is not assumed from disjoint support.

For the two compact momentum variables, positive P gives
`M=D P^-1 D>0`, and positivity of the complete compact 3 by 3 Hessian gives
`Kq=h-g^T P^-1 g>0`. The exact inverse of P yields
`m01/m00=-(c/b)P12/P22`. The separate carrier threshold bounds this ratio
in magnitude by one-half, since b=l0 is fixed, c has the stated finite
upper bound, P12 is bounded and P22 grows. Hence
`d=m00+m01>m00/2>0`; positive definiteness by itself would not have
established this needed fact.

With beta=B-q, direct substitution into the kinetic quadratic form gives
the mixed coefficient d and q coefficient `e_mass=m00+2m01+m11`.
The map `Psi=beta+(e_mass/d)q` removes that mixed term and yields

```
J_Psi=d²/e_mass,
J_beta=det(M)/e_mass,
K_Psi=Kq d²/e_mass².
```

All are finite and strictly positive, and Psi shifts with unit weight under
common physical rotation. Neither a target mass nor a target stiffness
appears in the selection threshold.

## Gyroscopic term and paired ensemble ordering

Eliminating both momenta from the full first-order action gives
`V^T M V/2-q V^T n-Kq q²/2`, retaining all off-diagonal entries. Only the
q qdot part is a total derivative. The q Bdot term is a real single-sign
gyroscopic connection and is not removed by a change of energy convention.

For the same geometric generators on -u0, lambda and H are unchanged while
Omega, b and c reverse sign. Thus M and Kq are unchanged and n reverses.
The positive/negative realizations must have independent compact shapes:
their extrema are `P^-1(DV-gq)` and `P^-1(-DV-gq)`. Substituting both before
averaging cancels the connection and retains M. Tying those conjugate shapes
instead would cancel the first-order kinetic response, as the verifier's
wrong-ensemble check shows. Equal weights and coherent B/q are explicit
zero-handedness/kinematic premises, not conclusions inferred from Euler.

## Joining boundary, oracle and decision

The established relative orbit contains one uniform global rotation. Its
vK=O(1/r) tangent need not be L2. The compact-base phase derivative estimates
of 0047 do not apply to a nonuniform amplitude B(x) multiplying K. Neither
the two-child combination nor a finite spatial modulus for that global
direction follows automatically. This is the already identified parent
joining problem, not a refutation of the valid uniform relative-orbit
action. In particular the positive J_beta is retained here as the cage-angle
coefficient; a later gradient interpretation needs its own licensed map.

All signs and positivity margins rest on exact symmetry, integrability,
finite-norm bounds and full matrix algebra. No sampled spectrum, arbitrary
precision sign guess or fitted parameter is used. Numerical small-ratio
refinement is not an additional oracle for these exact statements.

- Verification: exact analytic and symbolic-verified relative-orbit action.
- Review: audited; scientific acceptance recommended at the stated uniform
  rotation-plus-compact and independent-momentum paired scope.
- Compatibility: same actual EPS field, explicit invariant-ball relative
  normalization and physical Euler KKS action.
- Epistemic: this common/relative rotation route is established; spatial,
  affine and full-parent joining remain active.

No author files were changed. Correction check: not needed. Canonical
extraction, consumer replay and release promotion remain parent work.

Reviewed SHA-256 identities:

```
af25b78a821985a27b03f0d66093d2ca54b841eb73e795df34f7981697233ce2 direct-rotor.md
1ee23efb0887d062a8bb2abb8fa958612a4450243c8a85fb3413f1bcfc449ace verify.py
```

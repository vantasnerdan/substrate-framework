# Independent review: Bloch angle action and isotropic continuation

Reviewer: `/root/smooth_core_review`, distinct from the author of 0044.
Date: 2026-09-05. This is one bounded scientific review under AGENTS.md and
the physics/small-ratio skills. The previously reviewed 0040 microscopic
sector is an unchanged input; no common-body inertia is added here.

## Decision and strongest positive statement

The full divergence-free Bloch continuation, its exact raw Hessian and KKS
coefficients, and the internal-sector isotropic action are established.
The oriented cell's negative axial optical curvature is real and retained.
After averaging the declared SO(3) ensemble at the action level, the
periodic or compact-support macroscopic gradient action has strictly positive
curl and divergence coefficients. The boundary divergence and physical
local-angle field map are explicit and correct.

No blocking mathematical defect or scientific correction is requested.
This is an internal-sector result before common-angle/body and material
coupling. Its normalized speeds cannot be carried unchanged into a different
joint inertia. The raw H and Omega matrices are the reviewed inputs to that
next parent construction. Full-PDE invariant finite-dimensional dynamics or
all-wavenumber local equality is not added as an acceptance requirement.

The frozen evidence consists of `bloch-action.md`, the README and the three
scripts `bloch_sector.py`, `effective_action.py`, `verify_core_jet.py`, with
their saved 27/27, 15/15 and 8/8 receipts. The initial structural-equality
diagnostic and repaired output preserve the same rational matrices; no
tolerance, geometry or equation was retuned. Existing receipts were inspected
and reused. Base provenance is checkpoint 3626fbf and campaign v0.171.0.

## Physical Bloch representation and raw matrices

The Coulomb potential of each nonzero periodic generator mode is
`i m cross xi(m)/|m|²`. Taking curl at m+kappa gives the correct slow
return terms and exact shifted solenoidality. In the implementation,
convolution with the periodic background is performed before applying
the shifted velocity Leray projection. Thus the projected frequency is the
frequency of the force field, not the original generator frequency. The
finite set of force modes has no zero-frequency singularity near kappa=0,
so the stated Taylor expansions are analytic in that neighborhood.

The conjugate-mode bookkeeping computes the Hermitian energy form and
anti-Hermitian KKS form with the original physical sign. It is exact finite
Fourier algebra, not a numerical Fourier cutoff for an unknown profile.
At zero wave number both matrices return the physical 0040 sector.
Differentiating all entries along the three coordinate axes proves that
the entire first-derivative tensor vanishes: a linear tensor is fixed by
those three directional derivatives. No chiral term is discarded by hand.

For the two-coordinate action, the conjugate elimination gives
`I=B²/s`, not an inertia inferred from H alone. Its derivative coefficient
is `I2=I0(2B2/B0-s2/s0)`. The normalized stiffness coefficient is then
`C=h2-(h0/I0)I2`. The saved raw matrices reproduce the authored entries,
including the positive x/y coefficients and negative
`Cz=-7603/43200`. Since h0 and s0 are strictly positive and B0 nonzero,
negative curvature of this optical branch does not imply instability
arbitrarily close to its positive zero-wave-number gap.

Units agree: H per volume scales as rho b², Omega as rho ell b, and
frequency as b/ell. The gradient stiffness is rho b² ell² C and derivative
inertia rho ell^4 I2. The physical k² speeds therefore have velocity-squared
units without an extra density or microscopic length factor.

## What the local angle observable measures

The precise observable is the infinitesimal polar/skew rotation of the
material core-axis displacement Jacobian: one-half the axial curl of that
displacement. It is not an assertion about the major-axis orientation of
an elliptic contour under arbitrary additional symmetric strain, nor a
rigid rotation of every point in a finite cross-section. This is the
authored local-Jacobian definition, and the review retains it explicitly.

Directly curling the displayed real-space Coulomb potentials recovers the
0040 generators. Applying `curl+i kappa cross` twice and evaluating its
axial component at the two sections gives the exact local rotation
`plus/minus (1+kappa_perp²/2)q`. The shape generator gives zero for this
observable, including its slow-return terms. The return can translate or
strain the local core neighborhood; those effects are not relabeled as an
extra rotational coordinate.

The literal difference between the two section angles uses their own slow
phases. Factoring their sum about the midpoint gives exactly
`2 cos(pi kappa_z/2)(1+kappa_perp²/2)q_mid`, and hence the stated
`-pi² kappa_z²/8` axial correction. If the normalized physical coordinate
is `(1+r kappa²)q`, then I2 and h2 change by `-2r I0` and `-2r h0`.
Their combination C is unchanged. Thus both the finite-section correction
and the local polar-rotation correction are retained consistently rather
than used to hide the negative oriented coefficient.

## Action-level SO(3) average

The declared ensemble rotates the whole cell and its generators together
and sets the scalar amplitude to the axial component `q=n dot Q`.
For each raw second-gradient tensor, averaging rotations about n retains
only its transverse trace and axial entry. Mixed body-frame entries vanish
in this average. Those invariants are determined by the three axis values;
an omitted off-axis value is not needed to fix this SO(3) contraction.

At order zero the averaged inertia and potential are scalar multiples of
the identity, I0/3 and h0/3. Therefore their order-two operator normalization
is linear in the averaged raw derivative terms. This is why averaging
`h2-(h0/I0)I2` gives the correct normalized action in this particular
sector. The derivation does not assume that averaging arbitrary already
normalized dispersion relations commutes with adding another inertia.

Independently performing the spherical contraction gives

```
Ct=(Cx+Cy)/2,
<gradient action> = A ||G||²+B[(tr G)²+tr(G²)],
A=(4Ct+Cz)/15,
B=(Cz-Ct)/15.
```

The stated rational values give
`A=114523/648000>0` and `A+2B=19127/324000>0`.
These are action coefficients. Dividing by the averaged inertia I0/3 only
after this calculation yields the reported squared-speed coefficients
`114523 b²/51840` and `19127 b²/25920`.

## Boundary term and positivity scope

For `G_ij=partial_j Q_i`, the elementary identity
`|curl Q|²=||G||²-tr(G²)` gives the exact difference between the averaged
density and the positive curl/div form:

```
(A+B)[tr(G²)-(tr G)²].
```

Expanding the divergence of the authored flux cancels its second-derivative
terms after exchanging dummy indices, leaving precisely this null
Lagrangian. Its integral vanishes for periodic or compact-support fields.
It is not zero pointwise, and the authored negative pure-trace example
correctly prevents a claim of pointwise positivity on every affine matrix.
A domain with a physical boundary must retain the flux and corresponding
boundary variation.

On a periodic mean-zero field space, the identity
`integral(|curl Q|²+|div Q|²)=integral|grad Q|²` and the two strictly
positive coefficients give the claimed gradient coercivity. The negative
oriented-cell axial curvature therefore does not contradict the established
positive isotropic action under these boundary conditions.

## Oracle, findings and continuation

No numerical soft eigenvalue, energy subtraction, fitted wavenumber or
unmeasured truncation floor is used. The negative coefficient and positive
isotropic coefficients are exact rationals obtained from full matrix
derivatives. The complete action, symplectic convention, physical angle map
and admissible macro boundary conditions are all exposed by the exact
oracles. Additional numerical mesh or eigenpair tests are not requirements
for these statements.

- Verification: exact analytic and symbolic-verified Bloch/action identities.
- Review: audited; scientific acceptance recommended at this internal-sector
  and periodic/compact-support ensemble scope.
- Compatibility: compatible declared smooth-Euler affine ensemble; no
  background superposition, fitted cell ratio or suppressed boundary force.
- Epistemic: the isotropic route is established; the positive-parallel-axial
  route has its stated negative-curvature counterexample; parent joining
  remains active.

Canonical extraction, raw-matrix coupling to common/body coordinates and
the corresponding consumer replay are next parent implementation work. This
review neither performs nor certifies those uncompleted operations.
Correction check: not needed; no scientific correction was requested.

Reviewed SHA-256 identifiers:

```
3af41cea9827a03426f908510b6aefec721ea8e0db48194ff469627112cd428b bloch-action.md
8d59832d97e552d7d272c385ad3ad6b7049eb90edda4b2e12f0fadc313e0132a bloch_sector.py
ec1f2d3839916c2ae80d210216469058bc1c7d5070dc59efc589e160d438228c effective_action.py
7be4bbbe3d4a6674e1895484242f970c236df6436cacca080609f77f8b175aca verify_core_jet.py
```

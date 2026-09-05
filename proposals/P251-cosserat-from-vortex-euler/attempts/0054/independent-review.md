# Independent review: same-EPS finite-core affine shear

Reviewer: `/root/smooth_core_review`, distinct from the author of 0054.
Date: 2026-09-05. One bounded scientific review under AGENTS.md and the
physics/small-ratio skills. The object is the finite-core affine shear
attachment and its conditional full static reduction, not construction of
the parent's joint momentum block or a global assembly theorem.

## Decision and supported result

Established as stated: the actual EPS field supports the five compact
physical symmetric-tracefree deformation jets, the rank-aware fixed
moment projection, and the nine-bond cage attachment. Their complete
Euler Hessian is finite with its exterior velocity contribution retained.
For a fixed positive joint compact momentum block, a finite analytic
carrier threshold makes the full reduced shear matrix strictly positive.
Its simultaneous-rotation ensemble has the stated positive isotropic
shear modulus and no local STF-strain/vector-spin potential cross.

No load-bearing scientific correction is requested. The positive joint
P is an explicit dependency, held fixed while choosing the new shear
carriers. Removing common moment from shear generators does not require
proving a surviving common conjugate; constructing that positive centered
momentum sector is separate. No mass or final continuum assembly is
inferred from this shear attachment.

## Evidence and actual affine core jets

I read the README, complete `affine-shear.md`, verifier and saved 18/18
receipt. Exact curl calculus, finite-norm control, full Schur reduction
and representation algebra are the strongest practical oracles. The
existing receipt was inspected and reused, not redundantly rerun.

For tracefree E, the vector identity gives
`curl[y cross(Ey)]=-3 Ey`. Hence the cutoff potential has exactly the
claimed affine displacement and gradient E wherever its cutoff is one.
The complete return term is present, and the full vector field is a curl,
so its finite flow is volume preserving. A tracked point initially in the
interior constant-cutoff neighborhood remains there for sufficiently small
time; its derivative is then exp(tE), with determinant one. This is a
physical deformation of a smooth core neighborhood, not a substituted
director or singular-filament variable.

The energy uses the same complete second transported-vorticity variation
as the earlier actual-EPS construction. In particular the negative
helicity contribution remains in H; positivity is not asserted from the
velocity norm alone. Outside the compact vorticity-change supports, curl v
vanishes but v need not, leaving the stated positive exterior Gram with
the decaying finite tail. Whole-space Leray and physical exterior motion
are therefore retained throughout this shear matrix.

## Four-function projection and its exact class

The three mean responses use the same ambient projector. The scaled
common response is correctly represented by `-vK/L_ref`, since its
functional is Omega(K,xi)/(rho L_ref). Scaling is a fixed geometric
normalization; it changes neither the response nullspace nor the physical
condition being imposed.

All response curls are analytic in connected D. The positive-on-an-open-
ball cutoff therefore makes the response Gram nullspace equal to the
functional nullspace on curls of compact potentials: zero Gram norm
implies zero analytic curl everywhere, and compact curl integration by
parts gives the desired annihilation. Conversely, testing all compact
potentials detects a nonzero response curl. Thus F(xi) is in range G,
and the pseudoinverse correction annihilates all four moments without
presuming rank four. In particular a possibly dependent common functional
is no obstacle to removing it.

The torus topology is treated correctly: this does not automatically
annihilate arbitrary flux-carrying compact solenoidal directions against
a curl-free response of nonzero period. Every selected field in this
attachment has the required compact potential, so no such extension is
used. The correction is off the physical jets and preserves their affine
map. Zero compact-displacement mean plus zero corrected velocity mean
gives the stated material centroid slice. H and KKS are recomputed there;
the proof assumes no additional symplectic orthogonality.

## Nine bonds and uniform positive full energy

The exact bond sum is
`||E||_F²+(1/2) sum_i E_ii²`. In the displayed orthonormal STF basis,
the bond Gram is diag(3/2,3/2,1,1,1). Thus all five shear components are
resolved with the stated lower and upper bounds. Antisymmetric gradients
have zero extension on every bond. The centered-difference measurement
has precisely this leading affine strain and its stated third-derivative
remainder. The bond-to-cage attachment is an explicit kinematic premise;
its restoring energy is then evaluated from Euler rather than supplied
as an independent constitutive term.

Each cage has positive diagonal energy growing linearly with its carrier
magnitude and a uniformly bounded induced velocity. Disjoint supports
remove the local helicity cross, not the projected kinetic cross. Those
kinetic crosses remain bounded and are included in the remainder. The
fixed-response coefficients are bounded (indeed the compact oscillatory
pairings permit the stated inverse-carrier decay), so projection does not
remove the growing positive diagonal term.

For a unit STF coordinate vector, the bond norm bounds its cage
superposition by M_c. The fixed core/response remainder has velocity and
curl norms bounded by R0,R1. Cauchy--Schwarz and moving curl onto that
fixed remainder give the displayed C_s bound on every cross term. The
frame's unit lower eigenvalue then gives the full five-dimensional
H_shear lower bound, without coordinatewise positivity being mistaken
for positive definiteness.

With fixed conjugate fields, the same integration-by-parts estimate gives
`||N||<=rho F_p(M_c+R0)=N_star`. Static elimination of the complete
momentum block is exactly `C=H_shear-N^T P^-1 N`, including all mixed
entries of P and N. Its worst subtraction is bounded by
`N_star²/lambda_min(P)`. P is fixed and positive here, so all costs are
independent of the newly selected carrier magnitude and the strict
threshold is finite. The argument also permits a larger fixed positive
base block when proving a full joint elastic matrix. It does not claim
that this static elimination removes KKS gradient inertia or gyroscopic
connections from the full frequency-dependent action.

## Isotropic coefficient, separation and verification scope

Rotating the entire field, projection, bonds and physical jets together
preserves the construction and rotates C in the STF representation. Its
Haar average is a scalar identity. The verifier independently constructs
the representation generators and finds the symmetric commutant to be
one-dimensional; the trace fixes that scalar to tr(C)/5. Hence the
energy is `mu ||E||_F²`, with `mu=n_cell tr(C)/10>0`. This convention is
the incompressible shear modulus; no trace strain or bulk modulus was
constructed by these five directions.

The STF and vector representations admit no invariant linear cross
map. The explicit intertwiner calculation has rank 15 on its 15 unknowns,
and equivalently the isotropic epsilon tensor annihilates symmetric E.
Thus there is no local strain/spin potential cross in this ensemble.
Strain/wryness terms are a different gradient-order tensor: their
parity-odd part may be removed by the separately declared reflection
pairing after full reduction, not by deleting even kinetic gradients.

The positivity proof uses explicit finite norm bounds and a fixed positive
input P, not a sampled soft eigenvalue. The exact analytic threshold is
therefore the applicable small-ratio oracle; this review does not turn
the still-needed construction of P into an assumed numerical result.

## Frozen hashes and disposition

- `affine-shear.md`: `bcbb950e072aa76d214f58123c0c925c16039064c57693ae60fc0a95511f879f`
- `verify.py`: `4b42a0dff20feb0f0b5f9c17c8463516ee959ded53aabcf3813ddc86bd7497cb`
- `stdout.txt`: `838ca20e8fc24214c208f90c9c06834835bbb58d8cc1512ecd819ea694d5f885`

Acceptance is recommended at the stated conditional attachment scope.
The parent's centered momentum construction, kinetic normalization,
coherent continuum assembly and completion are separate obligations.

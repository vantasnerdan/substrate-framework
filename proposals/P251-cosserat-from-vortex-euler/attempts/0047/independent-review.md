# Independent review: direct-EPS geometric gradient cage

Reviewer: `/root/smooth_core_review`, distinct from the author of 0047.
Date: 2026-09-05. This is one bounded review under AGENTS.md and the
physics/small-ratio skills. The already reviewed 0045 compact angle sector
is an unchanged input. The review does not combine this two-coordinate
gradient action with the different four-coordinate common-rotation action.

## Positive result and explicit premise

The same-EPS gradient construction is established at its stated conditional
scope. Centered differences of the neighboring macro-angle values drive
explicit additional compact divergence-free fluid cages. Their complete
Euler orbit Hessian supplies a positive long-wave gradient tensor after the
full shape/KKS reduction, at a finite carrier threshold fixed by the stated
norm bounds. The original uniform-angle action and physical core tilt jet
remain unchanged.

The neighboring-angle/cage tie is an additional Cauchy--Born-type kinematic
premise. It prescribes an actual displacement map, not an elastic energy or
a desired coefficient. Unrestricted Euler is not asserted to select that
map. This distinction is load-bearing and must remain in the final ensemble
statement. The conditional result does not prove a universal constitutive
law from Euler without that premise.

No blocking mathematical finding or correction is requested. The failed
minimal-return principal route is correctly kept separate from this positive
replacement and from any verdict on the actual EPS background.

## Frozen evidence and oracle

The reviewed evidence is the README, `principal-symbol.md`,
`gradient-cage-proof.md`, `principal_symbol.py`, and
`verify_gradient_lift.py`. The two saved 16/16 exact receipts were inspected
and reused. Earlier stalled symbolic representations were not numerical
failures or physical counterexamples. The principal calculation is an exact
algebraic diagnostic on a frozen local vorticity jet; that constant jet is
not represented as an exact nonzero-lambda Beltrami field.

The analytic compact-support estimates are the strongest oracle for the
actual-field positive result. The scripts corroborate the vector/matrix
identities and do not replace those estimates with a finite numerical
spectrum. Base provenance is checkpoint 3626fbf and campaign v0.171.0.

## Principal negative route and full action jet

The shifted-wave generator, projected force norm and unprojected helicity
pairing give the authored full principal matrices. Their exact Darboux
transformation retains the diagonal KKS first derivative. It produces a
nonzero drift and negative normalized transverse curvature. Averaging
counter-drifting actions retains the positive `(v q_z)^2` term in the
Lagrangian, which is a negative laboratory-frame stiffness. Averaging
frequency squares would reverse that interpretation. The negative route
verdict therefore concerns the stated principal minimal-return mechanism,
not a sign inferred for the full compact EPS construction.

For the general actual-field two-coordinate jets, independently expanding
the determinant of a symmetric matrix plus an imaginary multiple of J gives

```
det(H-i nu Omega)=D-B²nu²+epsilon nu L
                  +epsilon²(P+nu²Q)+...,
L=2Ba+tr(adj(H0)S),
P=tr(adj(H0)H2)-a²,
Q=det S-2Bb2.
```

The near-identity Darboux map can be chosen to preserve real-field symmetry.
Its determinant-squared factor through order two is
`B²/(B²-epsilon²Q)`. Applying this factor to the characteristic polynomial
removes the varying nu² coefficient and changes the potential coefficient
to `P+nu0²Q`. After the conjugate elimination and kinetic normalization this
gives the authored laboratory coefficient `(P+nu0²Q)/t`. Expanding the roots
also gives `nu1=L/(2B²)` and the identity
`C=I0[coefficient_2(nu²)-2nu1²]`. Thus the chiral terms are included,
and the coefficient is not obtained by blindly averaging squared speeds.

In particular the coefficient of H2_11 is exactly one, while mixed-energy,
momentum-gradient and KKS corrections retain their displayed coefficients.
That unit coefficient is what permits the constructive positive lift.

## Exact support and derivative bounds

The centered difference is `2i sin(d_j kappa_j/2)`: it vanishes at zero,
starts with `i d_j kappa_j`, and has no quadratic amplitude term. Taking
the curl of the full potential supplies every slow-return term. The added
potentials vanish near the physical core. The original rotation potential
has zero value and first derivative there, so its first displacement jet
and resulting physical vorticity tilt retain their normalization.

KKS is local in the displacement generators. Disjoint cage/base supports
therefore remove their cross pairings exactly, including their Bloch
returns. The added self pairing has zero constant term; multiplying it by
the squared first-difference amplitude puts its first contribution at order
three. All KKS jets through order two consequently remain unchanged.

For energy, disjoint support does not remove Leray tails. The proof correctly
retains them. Writing the projected energy with a fixed whole-space P and
compact phase-multiplied force fields makes its kappa derivatives regular:
they differentiate bounded compact phases and potentials, not a singular
multiplier at zero frequency. In a base/cage cross term, self-adjointness of
curl places the derivative on the fixed base force. This yields the stated
M_i0 and M_i1 bounds.

The explicit upper bound on M_i1 follows by differentiating
`(e.x)F_i` and `(e cross A_i) cross omega`: the former gives both F_i and
`(e.x)curl F_i`, and the latter is bounded by the sum of its coordinate
direction curl norms. Compact support makes the coordinate radius finite.
Disjoint added supports give the displayed D_F bound, while their potentials
are O(1/K_j), giving the frequency-independent D_U bound. Thus N_i and
Lstar are genuinely independent of the newly selected carrier magnitudes.

The full matrix correction is

```
delta C=H(Z_e,Z_e)+r11-2g r12/t
        -(2a delta_a+delta_a²)/t.
```

The proof bounds all remaining terms by Rstar with exactly these factors.
No potentially growing momentum-gradient or KKS term is omitted. The
compactness of the fixed base fields is essential to this argument.

## Finite positivity and joining boundary

The positive diagonal cage energies grow linearly with their finite carrier
magnitudes by the 0045 estimate. Off-diagonal helicity terms vanish because
their force supports are disjoint; projected kinetic cross terms remain and
are bounded by the stated force norms. Gershgorin with the bond-length
weights gives the authored uniform lower bound. The original gradient
tensor is allowed to have either sign and is bounded by Cstar. Since every
A_j is strictly positive and all right-hand-side bounds are fixed, finite
carrier choices satisfying the threshold exist. The exact final coefficient
is the full matrix-jet integral, not that lower bound or a fitted stiffness.

The result is positive in every macro direction for this scalar internal
sector. Any parity/orientational removal of its odd time-space drift is an
explicit action-level ensemble step. The proof does not silently remove
that term from a single oriented realization.

The independently identified parent gap remains outside this child theorem:
the global rotation K used in 0048 is not a compact base field, and its
velocity tangent need not be L2. The compact-base derivative bounds here
cannot be applied directly to a spatially varying multiple of K. Likewise
adding the common-body momentum changes the reduction boundary. This review
does not infer that the two child actions already form a complete continuum.

- Verification: exact analytic and symbolic-verified conditional gradient
  construction; no numerical curvature extrapolation.
- Review: audited; scientific acceptance recommended at the declared
  displacement-ensemble scope.
- Compatibility: same actual EPS background and physical Euler orbit,
  with the additional neighboring-angle tie explicitly retained.
- Epistemic: positive gradient-lift route established; parent joining active.

No author files were changed. Correction check: not needed. Canonicalization
and replay of later coupled consumers remain parent implementation work.

Reviewed SHA-256 identities:

```
73fa646df64958a32ddd4b1799c60949e3b909265f6660a4cfea6c801404c367 principal-symbol.md
4cd56483c655903f149de3992e797dfab38d8bca3e043ede113abd9a2a2e8eea gradient-cage-proof.md
0a9b4c2882ae87bd4f19bc9eb6d9749af7b2b64a63a6a96d78b9b6c29261fc48 principal_symbol.py
2fb3d8220b33baee143447994e35fc6185c850a49722dd820e0d5c46bb624130 verify_gradient_lift.py
```

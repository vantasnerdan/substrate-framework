# Independent review: smooth six-core existence and affine angle action

Reviewer: `/root/smooth_core_review`, distinct from the authors of attempts
0036–0038. Date: 2026-09-05. This is one substantive review transaction,
authorized by the parent campaign and the physics skill. The reviewer read
AGENTS.md, the physics and small-ratio skills, the oracle/governance references,
the frozen 0035/0036 contracts, the three mathematical constructions and their
two symbolic scripts. No recursive review or delegation was commissioned.

## Decision and strongest positive statement

The proposed mathematical results are established at their stated scope.
There exists a family of exact planar Euler relative equilibria consisting of
six disjoint compactly supported C-infinity vorticity cores, with fixed equal
positive circulation and sufficiently small but nonzero core size. On the
specified two-triangle isovortical affine trial family, the exact smooth-core
Euler orbit action has a nondegenerate relative canonical pair and a positive
quadratic angle stiffness and inertia. The coefficients are finite-core
interaction integrals, not prescribed point-vortex values.

No load-bearing mathematical error or missing proof step was found in this
transaction. No scope reduction is requested. These results supply a genuine
smooth microscopic base state and its exact restricted action; they are not
merely a conditional statement assuming such a base state exists.

## Frozen boundary and evidence roles

The committed base is 3626fbf, and the campaign's pinned accepted release is
v0.171.0. The reviewed delta is confined to the following new evidence in
attempt 0036. Existing accepted claims are unchanged. Downstream use is the
P251 finite-core molecular action and its declared affine coarse-graining;
the complete parent continuum theorem is outside this child transaction.

- `smooth-polygon-construction.md`: exact analytic proof of the smooth
  polygon branch, using a proved radial kernel calculation, a bordered Green
  equation and the ordinary implicit function theorem.
- `radial-core-gap.md`: exact analytic proof of the radial nontranslation
  gap, including the distinction between weighted and ordinary translation
  gauges; also supplies the angular part of the seed elliptic kernel proof.
- `finite-core-angle-action.md`: exact orbit-form pullback and interaction
  energy estimate proving positivity at sufficiently small finite core size.
- `smooth_core_construction.py` and `radial_core_gap.py`: symbolic
  corroboration of the scaling, Green normalization, solvability and canonical
  factors. Their 13/13 and 18/18 saved receipts are reused. These scripts are
  not substitutes for the functional-analysis proof and are not counted as
  independent proofs of the existence theorem.

The source-only patch/smoothing routes in the attempt history are not needed
as imported existence or spectral theorems by the new direct construction.
The point-vortex Hessian from attempt 0032 is used only as an explicitly
computed limiting quadratic form. No theorem about spectral convergence from
those papers is silently imported.

## Independent analytic check

The following checks address the mechanisms that could produce a false
positive result, rather than increasing the check tally.

1. **Radial existence and the mass border.** For the fourth-power seed,
   positive radial mass gives `U'=-m/r`, which forces a finite transverse zero.
   Scaling gives `Gamma_A=A Gamma_1`; differentiating the scaling solution
   gives the regular radial linearized solution with a nonzero logarithmic
   coefficient. Thus a bounded mass-zero radial kernel is absent. The smooth
   flat nonlinearity converges in C1, sufficient for continuous dependence of
   the profile and its central-value derivative. The transverse zero and the
   nonzero circulation derivative consequently persist. Fixing circulation
   is an actual scalar implicit solve, not a presumed property of smoothing.

2. **Angular kernel and full border.** With `V=F'(U)`, the relation
   `a=V Omega` makes the elliptic Birman--Schwinger kernel exactly
   `sqrt(V) G_m sqrt(V)`. The Green comparison and weighted Schur estimate
   exclude a unit eigenvalue for m greater than one. The strictly positive
   Doob minorization makes the m=1 unit eigenspace simple. Reflection removes
   the second translation. Pairing the bordered equation with `V partial_1 U`
   yields the coefficient `integral y1 V partial_1 U=-Gamma`, so the force
   border vanishes in its kernel. The mass condition removes the global log
   term; the centroid condition then removes the remaining translation. The
   operator has index zero because its Green term is compact and its two
   domain/codomain borders are finite-dimensional of equal size. Injectivity
   therefore supplies surjectivity as well.

3. **Function spaces and smoothing.** On the fixed disk, convolution by the
   logarithmic kernel maps a bounded source into W^{2,p} locally for every
   finite p. Choosing p with `1-2/p>alpha` gives compactness into
   C^{1,alpha} on the smaller fixed observation region; equivalently one may
   extend the source by zero and use a surrounding disk. The core has a
   negative collar, so the nonlinear sources vanish near the disk boundary.
   Uniform convergence of V therefore controls the Green-operator norm on
   the stated space. For fixed positive smoothing parameter, the flat smooth
   nonlinearity has smooth Nemytskii composition on C^{1,alpha}. Neither a
   moving flat edge nor an unproved differentiability of a characteristic
   function enters the implicit function argument.

4. **Nonlinear force solvability and Euler residual.** The polygon Green
   gradient sums to `Gamma (N-1)/(4 pi d)` in the radial direction. The
   first-order external term is purely the bordered linear function, leaving
   zero first-order profile correction. Thus lambda/epsilon extends smoothly
   and has derivative d with respect to the unknown rate. The second implicit
   solve removes lambda exactly. Reassembling the self logarithm and the
   subtracted external constants gives `U=constant-psi+Omega |x|^2/2` on each
   support. Its vorticity is therefore constant on relative streamlines,
   including the smooth zero-vorticity collar. This proves the Euler
   relative-equilibrium statement, not just a small residual.

5. **Physical orbit form.** For a rigidly moved core with centroid c and
   centered material coordinate r, a local orbit potential is
   `-rho/2 integral omega_0 X cross dX`. Integrating
   `X=c+R(theta)r` gives the circulation translation term and the intrinsic
   angular impulse term; the mixed terms vanish by the exact first-moment
   constraint. Summing three cores gives
   `-3 rho (Gamma r_i^2+C_epsilon) dtheta_i/2`. Its relative canonical pair is
   precisely `p=-3 rho Gamma S x/4`, independent of the constant intrinsic
   impulse. Local compactly supported Hamiltonian extensions of the rigid
   motions supply admissible isovortical variations. This is the Euler orbit
   form, not a postulated second-order director kinetic term.

6. **Finite-core positivity.** Differentiating the mutual interaction after
   pulling each density back to its reference core differentiates the smooth
   separated kernel, not an epsilon-dependent sharp profile. The first-moment
   term vanishes and its remainder is uniformly O(epsilon^2/d^2), including
   two collective derivatives. Self energies are constant under the trial
   motions. At fixed total angular impulse the exact rotating equilibrium is
   critical for the restricted Hamiltonian. Dihedral symmetry removes the
   mixed relative-radius/angle Hessian. Independently expanding the point
   expression gives `H/A = constant + x^2/4 + 9 q^2/4 + higher order`, where
   `q=chi-pi/3`; hence the two limiting Hessian entries are A/2 and 9A/2.
   The uniform C2 estimate preserves their strictly positive signs for small
   nonzero epsilon. Eliminating the conjugate momentum then yields exactly
   `I_epsilon=P0^2/H_xx` and `K_epsilon=H_qq`. Both triangles react; neither
   an externally imposed strain nor a fitted elastic coefficient is needed.

## Findings and small-ratio applicability

There are no blocking findings. One non-blocking exposition refinement was
sent to the author: the radial gap document initially says smooth w, whereas
the unsmoothed fourth-power seed is C3 at its outer edge. The displayed kernel,
integration-by-parts and Schur arguments apply to radially regular C2 profiles
and hence to that seed unchanged. Stating this regularity extension makes the proof
dependency explicit without changing a result or adding a premise.

The small-ratio prescriptions are satisfied analytically: full radial
frequency operator and symmetry gauge are specified; the ordinary translation
complement has the proved positive bound `delta Omega_min`; the molecular
Hessian uses a direct mutual-kernel remainder rather than subtraction of
large numerical self energies. There is no discretized soft eigenvalue,
floating-point sign claim, numerical core-radius threshold or quadrature
remainder in the reviewed result. Consequently numerical eigenpair residuals,
mesh/domain extrapolation and roundoff-jitter tests would not strengthen its
present oracle and are not additional requirements.

## Compatibility, promotion and parent frontier

The signs, circulation normalization and dimensions agree with attempt 0032.
The whole-plane nonzero-circulation field is correctly scoped to renormalized
energy changes, not finite total kinetic energy. The result respects smooth
incompressible Euler, fixed circulation and isovortical variations.

The original N3 affine slowly varying ensemble may consume an exact
variational restriction; this review does not invent unrestricted finite-D
Euler invariance or all-wavenumber agreement as a new acceptance criterion.
Conversely, a consumer claiming an isolated full-PDE optical eigenspace would
need its separately stated transport/resolvent construction. Three-dimensional
translation, wryness and stationary EPS compatibility remain parent work,
not defects in these child statements.

- Verification: exact analytic proof, with symbolic-verified corroborating
  identities; no numerical spectrum is claimed.
- Review: audited, scientific acceptance recommended for the statements above.
- Compatibility: compatible extension under the explicitly declared classical
  Euler and standard analytic imports.
- Epistemic: the three route-scoped propositions are established; the parent
  P251 objective remains active.
- Registry/release: unchanged by this review. Canonicalization, consumer replay
  and release incorporation belong to the parent's eventual promotion
  transaction; this review does not authorize an intermediate campaign PR.

## Correction check and content boundary

The author added the C2 applicability clarification and an explicit
near/far-integral estimate proving continuity of the logarithmic Green
operator from L-infinity into C^{1,alpha}. The one correction check verified
these two additions and their direct dependencies. Both are correct and
change no scientific statement. No script inputs or asserted identities
changed, so their existing execution receipts remain applicable.

The final corrected mathematical-document hashes are:

```
5c4683b736b18cd6eddbcccff3857729734b8e36b73e84c03935e81e2cee91de smooth-polygon-construction.md
7be4a625e494c45de455f1789c544b71932c9e06a67bd4674d191b71bb1a4894 radial-core-gap.md
```

At substantive review the SHA-256 identities were:

```
af153bb45ebe09f25696bb8ff356abffd65f397c14e8c83c3ab76f211db07ad4 smooth-polygon-construction.md
1c38b3d58cd3170d9ee6d9054856cc817311ce256c0cc650137fb83a67a91579 finite-core-angle-action.md
da5e13cd7bec4879ffe2ac997967a9d950736986936cd76d5f17c12835deeeaf radial-core-gap.md
71f48d429c80a47e14b16f8aad49c6a3971e0592c30c690563c7158ceb499a57 smooth_core_construction.py
92471a697acff31d663c275f23e0e092db8f44ff27eae3bb5db2a0617be0e5ee radial_core_gap.py
```

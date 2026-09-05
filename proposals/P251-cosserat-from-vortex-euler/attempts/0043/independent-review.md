# Independent review: affine Euler energy of the Beltrami-cell ensemble

Reviewer: `/root/smooth_core_review`, distinct from the author of 0043.
Date: 2026-09-05. This review covers the new affine-energy claim and its
importable API only, under AGENTS.md and the physics/small-ratio skills.
It does not combine this energy with the parent's unfinished joint kinetic
or body-rotation action.

## Positive decision and boundary

The exact volume-preserving affine Biot--Savart energy identity is established,
and its declared isotropic Beltrami-cell ensemble has shear coefficient
`mu_affine=4E0/15>0`. For the actual two-mode tube this is
`2rho(a²+b²)/15`, hence `2rho b²/3` at a=2b. Both the advected vorticity
and the deformed lattice are included. No independent elastic potential or
velocity-only transport approximation is being used.

No blocking defect or scientific correction is requested. The result is an
exact constrained energy calculation. It does not assert that an affine
image remains a stationary Euler field, that independently oriented cells
can be glued or superposed, or that translational inertia and internal
cross-coupling have already been derived. Those exclusions preserve the
actual positive result rather than changing it.

The frozen evidence is the 0043 README and verifier, the new
`src/substrate_framework/euler_affine.py` API and its four tests. The saved
13/13 extraction replay and four-test success were inspected and reused.
The initial 12-check receipt remains valid at its pre-extraction boundary.
Base provenance is checkpoint 3626fbf and the campaign's v0.171.0 release;
no accepted registry statement is changed by this review.

## Independent physical and algebraic check

The affine material map x=Fa sends a covector k to `F^-T k`, while the
vorticity vector is pushed to `F omega_k` when det F=1. Their contraction
is unchanged, so each transformed vorticity mode remains solenoidal.
Inverting physical curl gives
`u_F=i k_F cross omega_F/|k_F|²`, whose curl is exactly omega_F. Parseval
on the transformed lattice then gives the authored energy formula. The
unchanged cell volume is essential and is supplied by det F=1.

For a real Beltrami field each nonzero Fourier pair lies on the curl
eigenvalue shell and is circularly polarized transverse to its direction n.
Its normalized real quadratic covariance is `(I-nn^T)/2`. With C=F^T F,
the numerator is therefore `(tr C-n.C.n)/2`, while the covector denominator
is `n.C^-1.n`. This independently produces the claimed finite-deformation
ratio; it is not inferred from the desired shear coefficient.

For symmetric traceless E and F=exp(tE), direct multiplication gives

```
numerator = 1-ts+t²(T-u),
denominator = 1-2ts+2t²u,
ratio = 1+ts+t²(T-3u+2s²)+O(t³).
```

The complete reciprocal denominator is retained. The exact isotropic
second and fourth moments give `<s>=0`, `<u>=T/3` and `<s²>=2T/15`.
The coefficient is consequently 4T/15 times the unstrained energy, with
the conventional shear-energy normalization `mu tr(E²)`. The verifier
includes all five independent traceless strain components rather than
testing only one deformation. Since the entire cell and its generators
are rotated together, this ensemble operation preserves the internal
core/cage correlations; it does not randomize their relative angle.

Objectivity is exact: a left rigid rotation Q sends both transformed
vorticity and wave covector to their Q-rotates and leaves both norms
unchanged. No spring is assigned to a common rigid rotation. Symmetric
traceless logarithmic strain is an exact admissible SL(3) path, avoiding an
unaccounted second-order volume change or pressure-work term. A bulk modulus
and incompressible longitudinal sound speed are not inferred.

## API and oracle scope

The API implements the full transformed-wave expression and checks matrix
size, determinant one, nonzero wavevectors, solenoidality and decidably
invalid real/positive inputs. Its documented real deformation, physical
wavevector, Fourier normalization and harmonic-velocity hypotheses are the
caller contract. The complex amplitudes are handled with conjugate norms;
real fields use both members of each Fourier pair as documented.

The four tests independently exercise the actual tube, a finite diagonal
deformation, a nondiagonal physical curl inverse, rigid rotation, empty
field and invalid determinant/solenoidality/zero-wave/density inputs. The
verifier's fixed-wavevector mutation exposes the principal physical failure
mode. The empty-amplitude limit removes the derived shear. The unpromoted
sphere-moment utility is used as a mathematical identity provider, not as
an accepted homogenization theorem or an assumed elastic modulus.

This is exact finite algebra with no soft numerical spectrum or energy
subtraction; numerical error-budget or mesh tests are not additional
requirements. The attempt verifier still uses a local tally helper rather
than the repository ledger. That is a non-blocking implementation-cleanup
item, not a defect in the physical API or proof and not a reason to rerun
its unchanged scientific identities.

## Decision and integration frontier

- Verification: exact analytic and symbolic-verified affine energy and
  isotropic shear identity.
- Review: audited; scientific acceptance recommended at this energy scope.
- Compatibility: compatible smooth-Euler constrained ensemble, without an
  imported spring, fitted parameter or change of physical conventions.
- Epistemic: this route is established; parent coupled dynamics remain active.

The new API and named tests have their recorded passing receipts. Full
promotion, cross-sector consumer replay and release updates belong to the
parent transaction. Adding this scalar to an independently reduced internal
sector without evaluating the common action is not licensed by this review.
Correction check: not needed; no scientific correction was requested.

Reviewed SHA-256 identifiers:

```
823de4043b1be8fec7a2749ab2cd5c5ce3d08f9392b754592fdbdab310badeff README.md
eedd6bfb50b4f1bf674ed8d60faf86fda93dea057a52d56e46a167cafcaf6904 verify.py
578613179531f82e640f4b3600e4817e97a3bd759eb36607dbac0be44459a6e8 src/substrate_framework/euler_affine.py
8b05ec1558661d25aa303f8aa7529b1aff455cd0b5707a92794500bbe4fd3264 tests/test_euler_affine.py
```

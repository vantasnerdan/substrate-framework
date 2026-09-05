# Independent review: material pressure and spin balances

Reviewer: `/root/smooth_core_review`, distinct from the author of 0051.
Date: 2026-09-05. One bounded scientific review under AGENTS.md and the
physics skill. The object is an exact weak material-partition balance,
not constitutive closure or identification with a selected orbit impulse.

## Decision

Established as stated: smooth constant-density Euler flow induces the
displayed center-based linear and intrinsic-angular momentum balances
on a material partition. The bond distributions include the complete
noncentral pressure-force moment, with the correct axial-stress sign.
The finite-domain angular-impulse identity retains its generally nonzero
boundary term. No load-bearing scientific correction is requested.

The homogeneous distributional equations apply to the complete partition
with internal shared faces (or with no external resultant). Where there
are external boundary tractions, the explicitly retained external force
and torque enter as source distributions. They are not silently removed
by the internal-pair formulas. This is the natural source convention for
the proof's stated external-traction provision.

## Material moments and pair reactions

I read the README, complete `pressure-balance.md`, verifier and saved
9/9 receipt. Exact Reynolds transport and integration against arbitrary
smooth test functions are the strongest oracle; the polynomial probes
and wrong-sign mutations corroborate those identities. The saved receipt
was reused, not duplicated.

For a material cell, constancy of its mass gives
`M Xdot=integral rho u`. Differentiating once more and applying Euler
gives `M Vdot=-integral p n`. In the derivative of intrinsic spin,
the convective product `(u-V) cross (u-V)` is zero and the center
acceleration term integrates to zero because `integral rho(y-X)=0`.
The surviving acceleration moment is therefore exactly the pressure
torque shown in the proof. Neither cell rigidity nor a fitted inertia
is needed for these kinematic and balance identities.

Opposite pressure tractions on a shared face give opposite forces.
Their intrinsic torques are not opposite: using
`y-X_b=(y-X_a)-R` gives `m_b=-m_a+R cross F` with the displayed sign.
It is this center shift, rather than a microscopic nonsymmetric pressure
stress or an external body torque, that generates the coarse torque term.

## Distributional stress and total angular conservation

For an arbitrary smooth test function phi, the bond force distribution
acts as `-F integral_0^1 d(phi(X_a+sR))/ds ds`. It equals
`F(phi(X_a)-phi(X_b))`, establishing the stated stress divergence.
For the couple bond, integrating the affine moment
`m_a-s(R cross F)` by parts adds the interior term
`-(R cross F) integral_0^1 phi(X_a+sR) ds` to the endpoint torques.
This proves the formula for all smooth tests, not merely degree-five
traces used in the script.

With divergence on the second tensor index,
`epsilon_ijk (F tensor R)_jk=F cross R=-R cross F`.
Consequently the intrinsic torque distribution is `div mu-ax(sigma)`.
Both the verifier's component calculation and its opposite-sign mutation
test this convention. The moving center delta functions supply exactly
the convective fluxes T and C in the displayed time-dependent balances.

The product rule gives
`x cross div sigma=div(x cross sigma)+ax(sigma)`.
The extra term cancels the intrinsic-spin stress torque. The translational
convective tensor T is symmetric, so its own axial contribution vanishes.
Thus orbital plus intrinsic angular momentum is conserved, with any
external force/torque sources retained. Smooth convolution commutes with
these spatial distributional derivatives and preserves the identities;
it does not eliminate the material-cell labels or imply constitutive
closure.

## Angular impulse is not automatically material spin

Integrating
`curl(r² u)=2 r cross u+r² omega`
over a finite domain gives exactly

```
integral rho r cross u
 =-rho/2 integral r² omega
                 +rho/2 integral_boundary r² (n cross u).
```

Subtracting the cell center velocity does not change the left side when
X is its mass center. The orientation of the boundary cross product is
fixed by the integrated curl identity, not a chosen moment convention.

For solid rotation in the unit ball, direct isotropic moments give
spin `8 pi rho Omega/15`, volume impulse `-4 pi rho Omega/5`, and
boundary contribution `4 pi rho Omega/3`. Their sum verifies the identity
and exposes the wrong sign obtained by using volume impulse alone.
No decay condition on a whole-space EPS field removes this boundary term
on a finite material cell. The proof correctly leaves its subsequent
orbit/material-spin matching as a distinct construction.

There is no discretized force, soft Hessian or numerical cancellation in
this claim. Exact weak identities and the explicit counterexample supply
its validation without an artificial numerical precision ladder.

## Frozen hashes and disposition

- `pressure-balance.md`: `679ee33c6aac17fd15511fded76d1c02acbe46ac0a91774392b26e2f7266a7ed`
- `verify.py`: `0b8163ed353d0ef03cd31f0d6098edad3f33fb54df26b6595b299a5a7411d1fb`
- `stdout.txt`: `4ad4ef0d9a120304e83bf7838c1d16223b3f4301aa4a31356af469b5a0108281`

Scientific acceptance is recommended for this exact material-partition
balance theorem. No accepted registry entry is changed by this attachment,
and independent rotational dynamics, constitutive matching and the full
parent continuum are not inferred from these balance identities alone.

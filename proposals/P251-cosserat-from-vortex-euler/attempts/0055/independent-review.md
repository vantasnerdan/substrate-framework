# Independent review: exterior response and centered KKS restoration

Reviewer: `/root/smooth_core_review`, distinct from the author of 0055.
Date: 2026-09-05. One bounded scientific pass under AGENTS.md and the
physics/small-ratio skills. The actual finite-harmonic EPS field and its
relative orbit are the previously reviewed inputs; the new object is the
exterior four-response theorem and resulting centered compact rotor.

## Decision and one precise proof wording

Established: the four response curls are independent on every nonempty
exterior open ball. Four disjoint dual responses consequently give the
exact physical mean projection, mutually isotropic response span and
restored canonical rotor KKS blocks. The complete compact Hessian is
positive above a finite carrier threshold. The selected material cell
has zero centroid response, and the surrounding angular impulse is
retained in the physical current map.

One wording in section 2 needs its precise infinitesimal meaning:
**translation derivatives**, not arbitrary finite translations, preserve
a finite regular spherical-harmonic Helmholtz expansion about a fixed
origin. The derivative statement is exactly what vK needs. A finite
translation can instead create infinitely many angular orders about that
origin. The explicit argument below repairs this sentence without changing
the positive theorem or introducing a new hypothesis. No other load-bearing
correction is requested.

The result centers D, not every parcel of surrounding fluid. Its global
material assembly and finite spatial return remain separate parent work.
No finite-energy whole-space Galilean boost or finished parent continuum
is inferred from this attachment.

## Evidence and far-field argument

I read the README, complete `exterior-centering.md`, verifier and saved
18/18 first receipt. The source-specific finite-harmonic property is the
same EPS Theorem 8.3 construction already inspected in the 0048 review.
The new oracle is analytic far-field independence followed by exact
response and KKS algebra; the finite matrix examples in the verifier are
corroboration, not substitutes for the actual-field Gram proof. The
existing receipt was reused without a redundant run.

In the exterior, the Leray response is minus the gradient of the Newton
potential derivative. Differentiating its leading volume term gives

```
w_a=|D|[3n(n.a)-a]/(4 pi r³)+O(r^-4),
n.w_a=|D|(n.a)/(2 pi r³)+O(r^-4).
```

The coefficient depends on the cell volume, not an approximation of D by
a ball. The corresponding derivative remainder is valid away from the
bounded source. The symbolic dipole sign and radial factor agree with
this direct differentiation.

For the supplied EPS field, regular spherical Bessel asymptotics give
the stated smooth angular sine/cosine amplitudes. If both leading
amplitudes vanish, their spherical-harmonic coefficients vanish separately:
different angular orders cannot cancel each other as functions of n.
Thus a nonzero finite regular Helmholtz sum has a nonzero leading 1/r
coefficient. Curl-polynomial conversion preserves finite angular support.

More explicitly for the rotation tangent, in coordinates about the
original expansion origin,

```
vK=e cross u-(e cross x).grad u+(e cross X).grad u.
```

The first two terms are the infinitesimal vector rotation, preserving
finite angular order. A Cartesian translation derivative of a finite
regular Helmholtz expansion changes angular order by at most one. This
also follows from its sphere integral representation: differentiation
multiplies the finite angular density by i lambda n_j. Hence vK is a
finite regular sum and, since it is nonzero, has a nonzero 1/r amplitude.
This is the precise derivative argument needed in place of the proof's
overbroad finite-translation wording. No translation of the radial origin
is required. Since curl vK=lambda vK, curl f0 cannot decay as O(r^-4).

For the mean response, solenoidality gives
`curl(omega cross w_a)=(w_a.grad)omega-(omega.grad)w_a`.
The second term is O(r^-5). The transverse part of the first is likewise
O(r^-5), whereas its radial part gives exactly the displayed O(r^-4)
oscillatory coefficient. If that curl vanishes identically in the
exterior, taking radial sine/cosine phase subsequences gives
`(n.a)A(n)=(n.a)B(n)=0`. For a nonzero a, its nonzero-direction set is
dense, so continuity forces A=B=0, a contradiction. The three mean
response curls are therefore independent.

Analytic continuation on the connected exterior converts vanishing on
any open ball into exterior vanishing. In a relation involving f0, the
1/r versus 1/r^4 orders first force the f0 coefficient to zero. The mean
argument then forces the remaining three to zero. This establishes the
strict four-response result without a numerical rank or radiation-condition
assumption.

## Disjoint exact duals and complete KKS pairing

For each positive-on-an-open-subset bump, any zero quadratic Gram norm
would give a response relation on that subset and contradict the theorem.
Every local four-by-four Gram is therefore strictly positive definite.
Its inverse defines the stated compact curl dual by integration by parts.
The integral can be small, but it is strictly positive; no size estimate
or floating-point sign is needed to assert a finite exact inverse.

Taking only one dual in each of four disjoint balls is essential. It
gives both F_a(eta_j)=delta_aj and zero mutual KKS entries, the latter
by disjoint supports or antisymmetry. A generic set of duals in one
shared support would not automatically have this isotropy. The script's
nonisotropic mutation correctly exposes that distinction.

All raw fields are disjoint from the response balls. Expanding the
projected r,Q,S pairings therefore cancels every response cross exactly.
The common moment becomes b0 on r and zero on Q,S, all three physical
cell means vanish, and the internal pairing remains precisely the raw
circular cage value c_k. Its existing finite-k bound gives nonzero c_k;
there is no uncontrolled asymptotic symplectic correction. Body/internal
raw supports remove their cross pairings, and all physical core jets
remain unchanged because the corrections are outside D.

The complete selected form is thus the two canonical blocks with
determinant b0² c_k², not merely a full-rank matrix with neglected entries.
All response generators have bounded support inside a finite coherence
region. The use of a far-field theorem to prove their duality does not
put a noncompact tail into their displacement support.

## Full positive Hessian and ensemble convention

The raw cage response coefficients decay under compact oscillatory
integration by parts, while the core response is fixed. Consequently the
projection adds only bounded coefficients of fixed smooth fields.
Their velocity and curl norms are finite. Moving curl to these fixed
fields bounds their energy crosses with high-frequency cages, while
the growing positive principal cage energy is retained. This proves
the stated full three-by-three lower bound, including all projected
kinetic cross entries. A finite carrier exceeds its fixed remainder.
The common energy row remains zero by the relative rotation symmetry.

Every selected compact direction has zero physical D mean, and K has
zero complete moving-domain mean. Hence the finite-cell centroid energy
Gram actually vanishes here; it is not subtracted while leaving old KKS
entries unchanged. This does not remove surrounding parcels' energy
or momentum, which the proof explicitly retains for subsequent assembly.

The positive Hessian and fixed b0 feed the existing complete Routh
reduction. The old finite carrier argument for the nonsingular physical
angle map applies with its newly fixed response norms. For time reversal,
the geometric generators constructed on one field are frozen and reused
on its negative, as in 0048. Their KKS signs then reverse while H is
unchanged. One should not reconstruct sign-dependent normalized duals
independently and silently identify different coordinate conventions.
Equivalently those dual conventions must be transported consistently.
Under the stated frozen-geometry ensemble, independent reaction momenta
give the existing time-even reduction exactly.

## Ambient impulse and physical current

For compact vorticity changes, splitting the full impulse into D and
its exterior gives `J_coh=J_D+J_ext` with both changes finite. The
finite-domain curl identity gives physical spin `delta L_D=J_D+B_surface`.
Combining them therefore requires
`B_effective=B_surface-J_ext`, with precisely the displayed sign.
Outside compact generators can change J_ext while leaving the physical
core jets untouched; their impulse cannot be relabeled as core spin.

The antisymmetric tensor formed from B_effective has axial vector
B_effective and divergence curl(B_effective)/2. The full current
improvement of 0052 therefore applies with this corrected quantity,
retaining the boundary angular flux, force-flux time derivative and
moving-center convective bookkeeping. The positive exterior integral
in t_eff is the derivative of minus J_ext, so its sign is also correct.

For fixed geometric momentum tangents, the surface and exterior impulse
rows reverse under circulation reversal, together with the KKS diagonal.
Independent plus/minus elimination retains exactly
`t_eff P^-1 diag(b0,c_k) (Bdot,qdot)^T`. Its nonzero rate response is an
observable current filter, not an assigned spin modulus or discarded
surrounding-fluid reservoir.

## Frozen hashes and disposition

- `exterior-centering.md`: `e95280606561a10063cd0dc739aef66d40759e86616f8f9b063f11e9d0db917a`
- `verify.py`: `74c7a051fc77e11ba8e3507eb2f70312b4593f31f5b84322e0bbf7e49f5ff21e`
- `stdout.txt`: `aeb7ac3e83fe156cbb8e41c48557866d8e9830b385decc087735f6fd26a492d9`

Acceptance is recommended with the precise translation-derivative
statement recorded above. The strict Gram and Hessian results are
analytic, not numerical small-eigenvalue decisions. Parent material
assembly and spatial-gradient joining remain separate positive tasks.

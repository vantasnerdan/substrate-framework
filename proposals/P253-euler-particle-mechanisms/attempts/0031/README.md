# P253/0031: finite relative observables and action for the exact solitary Euler family

This README preregisters a bounded authoring attempt for the exact solitary
family established in root-owned `0027` and independently reviewed in `0028`.
The positive objective is to determine the strongest exact physical statement
carried by that family about finite relative angular momentum, helicity,
translation/rotation momentum maps, and action. It must also decide whether the
finite-excess phase space contains an actual nontrivial compact symmetry or
internal cyclic coordinate whose finite symplectic action can be derived.

There is no preferred verdict. A finite field integral is valuable even if it
does not generate a compact orbit, and an admissible coadjoint action is
valuable only when its group, phase space, stabilizer, period, surface terms,
and physical normalization all close. Neither outcome is to be adjusted to a
target value. In particular, no value of `hbar`, spin one-half, fitted action
quantum, or electron interpretation enters the attempt.

## Ownership and activation boundary

`particle-balance-review` authors only append-only attempt `0031`. Root owns
the central proposal registration, activation receipt, accepted registry,
memory, source/API changes, validation, and commits. Root-owned `0030`
stability/restoring work and carrier-owned `0029` block-reducibility work are
disjoint concurrent units: their bodies, computations, and conclusions are
not inputs here and cannot inherit from or supply this attempt.

Before activation, this README is the sole `0031` artifact. The author may use
the already completed `0028` review as the inventory of the reviewed `0027`
family, but will not reopen `0027` scientific bodies, start the observable or
action derivation, inspect `0029` or `0030` bodies, create a verifier, or write
any result artifact. Substantive work begins only after root:

1. registers `0031` centrally with this ownership and frozen scope;
2. records the exact repository-schema invocation in `0031`; and
3. writes `attempts/0031/activation-schema.exit` with value exactly zero.

After activation, pin this README, the activation receipt, and every inspected
artifact by SHA-256. The observed preregistration head is
`614cdc63bec394bc61b5d34b928d627133510dcb`; the accepted release is
`v0.183.0`. Shared-worktree files under `0029` and `0030` are expressly outside
the transaction. No central, source, API, test, memory, generated-document, or
commit edit is owned without a later explicit scope extension.

## Frozen supplier and source inventory

The sole scientific supplier at freeze is the `0027` family at the boundary
adjudicated by `0028`:

- `0028/README.md`, SHA-256
  `a90ffcb1d0e05ccdbcf6b97db4941af430eab5c457ba04c96ae1094dd848887a`;
- `0028/review.md`, SHA-256
  `91fa6ffc55c43acec1d87f05b18100e6f0d30f0a18734f44aa6696399601ddc0`;
- `0028/verdicts.yaml`, SHA-256
  `cba564a7e38818d6fa1c156cfb08f0b06ada5d502236373bba8ac66c5076d3c9`.

The underlying supplier hashes recorded by that independent review, and not
reopened during this preregistration, are `0027/README.md`
`6cd66e3acb69c39d4e1347fde56cb8f49c406b421091a6fd11b3534a1fe48d41`,
`0027/solitary-wave-construction.md`
`c18d394a757f9fd967975ba23a486e62087c285c8fda4f08d8aeaa0159c6f850`,
`0027/exterior-construction.md`
`0c2cda195dd5bea9163a7a6c582bfcff746a608a2e48d0a9cd7a007c1b260fc3`,
`0027/result.yaml`
`50c022db1a6a9f542b10e732f3aa043b1443b2832a18e5ae5ff241ad0df43c29`,
and `0027/validation.md`
`fbfd4c387d1e1a5b39d82b2edf0eb64b2fc6fd97b4559cabb3453992f3ee18cc`.
Only artifacts at that boundary inherit `0028`'s review; any later hash drift
must be reconciled before it is used.

That review establishes a smooth nonzero axisymmetric traveling excitation

    u(t,r,z)=U(r)+v_c(r,z-c t),
    U(r)=(L(r)/r)e_theta,

of constant-density incompressible Euler on `R^3`. The background has compact
transverse axial vorticity but a nonzero `1/r` circulation tail and infinite
total kinetic energy. The perturbation is smooth, belongs to every finite
Sobolev order, has an exact irrotational same-fluid exterior, and has absolutely
convergent literal kinetic excess. The branch exists for `c>c_0` sufficiently
close to `c_0`, with small parameter `mu` and large axial scale `L_mu` satisfying

    mu L_mu^2=f_0(R)^2 log L_mu.

These facts are supplier inputs, not conclusions about moment, helicity,
coadjoint geometry, or stability. The `0028` review explicitly withholds
physical spin. After activation, the exact `0027` construction, exterior,
result, and validation artifacts may be reopened only to derive the frozen
observables from their actual fields and estimates. The executable
`euler_column_wave` identities may corroborate local algebra but do not prove
integrability, cutoff independence, a global group action, or a symplectic
period.

No inaccessible article is imported as an action or momentum-map theorem.
Standard vector calculus, cylindrical-coordinate integration, Noether's
theorem, Lie--Poisson reduction, and KKS coadjoint-orbit calculus are permitted
mathematical tools only after their hypotheses, boundary terms, and sign
conventions are checked on this affine infinite-background phase space.
Earlier P253 Euler-action attempts are not frozen dependencies. All
load-bearing identities are to be rederived for the actual `0027` field; any
later source expansion must be recorded append-only before reliance.

## Frozen phase space and observable questions

The starting physical phase space is an affine finite-excess class

    P_U={U+v : v has the regularity, exterior behavior, and integrability
                 required by the reviewed solitary family},

with the background `U` fixed at infinity. The derivation must sharpen this
placeholder to a precise function space and tangent class. In particular,
finite `||v||_2` alone does not make the cross term `integral U dot v`, a
moment, helicity, or KKS pairing finite. Every observable must be defined first
on finite cylinders or another explicit exhaustion and proved absolutely
convergent or cutoff-independent before it is called literal or physical.

The direct candidate for literal relative angular momentum is

    J^rel[u;U]=rho_m lim integral_C x cross (u-U) d^3x.            (F0)

Its axial component is

    J_z^rel[u;U]
      =rho_m lim integral_C (x cross (u-U))_z d^3x
      =2 pi rho_m lim integral r^2 (u_theta-U_theta) dr dz,       (F1)

where the limits and cylinder `C` will be specified after activation. If the
stream/swirl representation is used, the factor `u_theta=F(psi)/r` and the
cylindrical Jacobian must be propagated rather than folded into an unnamed
normalization. The attempt must distinguish (F0)--(F1) from vorticity impulse or
moment formulas. Any conversion between them must display the axis and
infinity surface terms and prove that they vanish or retain them explicitly.
All three components and their absolute convergence must be addressed;
azimuthal cancellation on a symmetric cutoff is not a substitute for a
well-defined vector integral.

The origin dependence must also be explicit. Translating the origin mixes
angular and linear momentum, so an axial value can be called intrinsic or
spin-like only after an actual center/translation reduction or the relevant
zero-momentum condition is derived. A same-density velocity perturbation does
not by itself define a localized mass density or physical centroid.

The kinematic relative helicity candidate is

    H^rel[u;U]
      =lim integral_C [(U+v) dot curl(U+v)-U dot curl U] d^3x
      =lim integral_C [U dot curl v+v dot curl U+v dot curl v] d^3x.  (F2)

The calculation must separately state the convention without density, the
density-weighted physical quantity `rho_m H^rel`, and any optional factor
`1/2`; these conventions cannot be interchanged. Integration by parts between
the two cross-helicity terms must retain the signed boundary flux involving
`U` and `v`. The background's `1/r` tail makes omission of that term
inadmissible without a proved limit. Pointwise zero background helicity, if it
holds, is not by itself convergence of (F2).

If translation momentum is used to construct the traveling-wave functional,
the attempt must distinguish literal relative linear momentum, hydrodynamic
impulse, and the momentum map that generates axial translations. Their equality
is a proposition requiring its surface terms, not a naming convention. The
same applies to rotation: a finite value of (F1) is not proof that it is the
Hamiltonian generator of an admissible compact orbit.

## Route A: exact fixed-background relative moment and KKS construction

Route A constructs the strongest possible exact observable/action theorem in
the affine space `P_U`:

1. Define the volume-preserving diffeomorphism group or local regular cover
   that preserves the allowed asymptotic background and circulation class.
   Characterize its Lie algebra and show that every generator used is tangent
   to `P_U` with the required decay and regularity.
2. Derive the momentum one-form density, Lie--Poisson bracket, KKS two-form,
   Hamiltonian sign, and any cocycle or affine correction with the physical
   density `rho_m`. Track pressure-exact, harmonic/closed-nonexact, circulation,
   and stabilizer rows rather than recovering them only after taking an
   exterior derivative.
3. Prove finiteness and cutoff independence of the KKS pairing and of the
   symplectic-potential action on the actual orbit. A compactly supported
   tangent or path lift supplies only that tangent/path unless its brackets and
   global composition really integrate to the claimed group action.
4. Derive the finite relative Hamiltonian from the already established literal
   kinetic excess, and derive rather than assume the translation or rotation
   momentum term in a relative traveling-wave functional such as
   `E^rel-c P_z^rel` with the convention-dependent sign fixed by Euler.
5. Evaluate (F1), (F2), every required surface flux, and the action on the exact
   solitary family. State which quantities are exact integrals, which have only
   a leading small-`mu` asymptotic, and which are merely candidate expressions.

The admissibility audit must treat rotations geometrically. An axial rotation
may preserve the axisymmetric background but can lie in the stabilizer and act
trivially on an axisymmetric wave. A rotation that tilts or otherwise changes
the infinite vortex background does not automatically act within `P_U` and
cannot be assigned a finite relative action without an explicit enlarged
phase space and renormalization. Axial translations preserve `U` but form a
noncompact `R` action; traveling motion is not thereby a cyclic `S^1` degree of
freedom.

Route A earns `established` if it supplies exact finite relative observables
and the precise admissible action/orbit statement, including a rigorous
negative conclusion when a proposed compact orbit is absent from the fixed
phase space. It earns `refuted` only when its specified construction is
contradictory, or `blocked` with the exact missing convergence, tangent-group,
or global-action step. A lack of compact spin orbit does not erase valid finite
moment or helicity results.

## Route B: actual internal mode or cyclic coordinate

Route B searches for a genuine compact degree of freedom of the same exact
Euler family, rather than declaring one from a scalar integral. A candidate
must provide:

1. a family `u_alpha` of actual Euler states in the same finite-excess phase
   space, with a globally defined group parameter and derived least period;
2. a nontrivial tangent after quotienting the stabilizer, with the same
   background, circulation class, exterior matching, and regularity;
3. invariance of the full physical Hamiltonian and a finite KKS form or
   symplectic potential on the orbit;
4. the exact Noether momentum paired with that coordinate, including all
   density, `2 pi`, sign, and surface factors; and
5. a derived action over one orbit, with chart transitions or cocycles included
   and no assigned normalization.

The wave amplitude, speed, axial position, time phase, streamline label, or an
arbitrary rotation of coordinates is not a cyclic internal coordinate unless
these five conditions are proved. A parameter interval is not a circle, an
axisymmetric stabilizer motion has zero orbit dimension, and a noncompact
translation has no finite period merely because one can choose an interval.
Route B receives its own `established`, `refuted`, or `blocked` verdict. If no
such mode is present in the actual `0027` family, name the missing construction
and preserve Route A's strongest result rather than inferring a global Euler
no-go.

## Exact evaluation and observable-transfer checks

For each claimed quantity, the attempt will freeze the physical convention and
then establish the chain

    exact Euler field -> convergent cutoff formula -> exact relative integral
      -> controlled small-parameter expansion -> physical interpretation.

The following distinctions are mandatory:

- `Q(mu)=C mu^p+o(mu^p)` with `C!=0` proves eventual nonzero only after the
  remainder and parameter side are controlled; a formal coefficient alone is
  not an exact nonzero integral.
- A leading-order formula is not the literal full observable. If no closed form
  is available, the exact convergent integral remains the theorem and the
  asymptotic is labeled separately.
- Cancellation between core, exterior, background cross, and surface terms is
  checked before assigning a sign. Positivity of kinetic excess does not imply
  nonzero angular momentum or helicity.
- All cylindrical reductions retain `2 pi`; physical energy and momentum/action
  formulas retain `rho_m`; kinematic helicity conventions are stated
  separately. Dimensional analysis must agree across KKS, Hamiltonian,
  momentum, and period factors.
- Translation of the localized profile through the fixed column may make a
  convergent integral time independent, but this is not a general conservation
  theorem until the appropriate Euler flux at infinity is controlled.

The strongest oracle is direct analytic evaluation from the exact `0027`
streamfunction/swirl construction, together with an independently derived
cutoff/surface calculation and group-action/KKS calculus. Symbolic algebra may
check cylindrical identities, expansions, signs, and factors, but cannot prove
absolute convergence, function-space tangency, global group integration, or
existence of an internal mode. Representative exposing mutations are removal
of `rho_m` or `2 pi`, reversal of the helicity boundary-flux sign, replacement
of an exact integral by its leading term, use of a tilt that changes `U`,
identification of a stabilizer with an orbit, or assignment of a period to an
`R` translation.

No production numerical observable, small force, soft eigenvalue, stability
edge, or energy splitting is registered. The small-ratio-numerics skill does
not bind at freeze. If exact analysis leaves an irreducible small numerical
remainder, its quantity, scale, error floor, and verifier design require a new
append-only registration before computation.

## Selection, verdict, and completion contract

The two routes compete by exact Euler closure, finite and cutoff-independent
observables, natural preservation of the fixed background, globally defined
symmetry and stabilizer, finite physical KKS/action, parameter economy,
dimensional consistency, and observable reach. They are not ranked by proximity
to a desired spin or action value. Route generation remains open to a
failure-derived representation change within this bounded objective, but each
attempted route receives exactly one scoped verdict.

After activation, `0031` will contain a derivation/source audit, a structured
result with separate Route A and Route B verdicts, and a validation receipt.
Reusable executable algebra is added only if it materially exposes the exact
claim and only under separately authorized source/test ownership. Completion
requires the strongest true relative-observable statement, an exact decision
about the admissible finite action/orbit claim, all surface and normalization
terms, and one named next dependency. Report artifact hashes, verdicts, and
that dependency to `herdr agent prompt w3:p1`.

This attempt does not license stability, restoring dynamics, interaction
carriers, an all-time open neighborhood, finite total background energy,
particle identity, a mechanical rigid-body spin, quantum spin or statistics,
`hbar`, spin one-half, electron/neutrino completion, parent-campaign
completion, or a global no-go for Euler internal modes.

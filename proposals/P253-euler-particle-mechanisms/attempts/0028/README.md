# P253/0028: independent exact-existence review of the unbounded-background Euler solitary wave

This README preregisters one fixed-theorem, non-author review of root-owned
attempt `0027`. The positive object is the strongest actual smooth localized
traveling-wave solution of the full constant-density incompressible Euler
equations on the unbounded background and spatial domain that `0027` proves.
The review will audit the complete operator, implicit-function or equivalent
existence argument, symmetry axis, exterior matching, and finite relative or
excess-energy statement as one coherent existence unit.

This is not a stability review and does not ask the existence proof to supply
a restoring spectrum, nonlinear orbital control, a persistent particle, spin,
or electron/neutrino completion. A correct exact existence theorem is retained
at its strongest proved hypotheses even when those later mechanisms remain
open.

## Independence, ownership, and activation boundary

`particle-balance-review` owns only append-only attempt `0028` review
artifacts. Root owns `0027`, central proposal and registry changes, source/API
materialization, tests, memory, validation, and commits. The reviewer did not
author or implement `0027` and will not edit its source, equations, API, tests,
or receipts.

Before this freeze, the only `0027` scientific artifact opened was
`attempts/0027/README.md`, SHA-256
`6cd66e3acb69c39d4e1347fde56cb8f49c406b421091a6fd11b3534a1fe48d41`.
Shared-worktree status exposed names of later `0027` files but none of their
contents. The central proposal currently registers `0027` as root-owned with
scope “P2/P3/P4 exact solitary Euler waves on vortex backgrounds with physical
exterior matching.” No result, equation, source extract, implementation, test,
or output from a `0027` body has been inspected.

No further `0027` artifact may be opened until root:

1. centrally registers `0028` with this independent ownership and frozen
   target;
2. places the exact repository-schema invocation and its exit in `0028`; and
3. `attempts/0028/activation-schema.exit` is exactly zero.

After activation, pin every reviewed artifact and source receipt by SHA-256,
perform one substantive pass, and allow at most one bounded correction check
limited to an evidence-driven repair. Root may materialize the `0027` API in
parallel; only the frozen hashes actually reviewed can support the verdict.
The accepted base is `v0.183.0`; the observed preregistration head is
`bd61f489484701c849bc309f4c1df361b85b7f13`.

## Frozen target and maximum verdict

The target kind is `fixed_theorem`. No empirical comparator or mechanism
selection is involved, so no artificial competing construction is required.
The maximum review verdict is an exact existence theorem of the following
form, with every placeholder fixed by the actual `0027` proof:

    u(t,r,theta,z)=u_bg(r,theta,z)+v(r,theta,z-c t),
    p(t,r,theta,z)=p_bg(r,theta,z)+pi(r,theta,z-c t),

where `(u,p)` solves incompressible Euler pointwise or in an explicitly stated
class on its entire unbounded physical domain, `v` is nonzero and localized in
the claimed directions, the field is regular at the symmetry axis and across
every support/interface boundary, and its declared relative/excess energy is
finite. Axisymmetry and swirl/no-swirl content, background asymptotics, wave
speed, parameter range, decay, regularity, and uniqueness or branch
quantifiers must be those actually proved.

A rigid-cylinder theorem, formal Bragg--Hawthorne ansatz, local core solution,
wall eigenmode, asymptotic KdV/NLS profile, numerical continuation, or a
solution of a truncated radial problem does not by itself establish this
unbounded-domain target. Conversely, if `0027` honestly proves an exact
wall-bounded or conditional theorem rather than the full target, preserve that
strongest useful statement and name the precise exterior construction needed
to restore the stronger form.

## Frozen source and dependency inventory

The only preactivation scientific source inventory is the one stated in
`0027/README.md`:

- S. M. Sun, *On the existence of solitary waves in rotating fluids*, Proc.
  R. Soc. Edinburgh A 125 (1995), 1105--1129,
  DOI `10.1017/S0308210500022678`. Only publisher metadata or an extract was
  accessible at the `0027` freeze and it names a rigid cylinder and equilibrium
  hypotheses. The inaccessible article is **not an imported theorem**. It can
  provide provenance or a route name only unless an exact accessible statement
  and hypotheses are later pinned and audited.
- T. B. Benjamin, *Theory of the vortex breakdown phenomenon*, J. Fluid Mech.
  14 (1962), 593--629, and *Some developments in the theory of vortex
  breakdown*, J. Fluid Mech. 28 (1967), 65--84. These may support only the
  exact equations or spectral facts verified from accessible pinned bodies;
  historical motivation is not an existence theorem.
- S. Leibovich and A. Kribus, *Large-amplitude wavetrains and solitary waves in
  vortices* (1990; exact journal/version/DOI to be verified). Numerical
  continuation may corroborate a branch but cannot prove the exact continuum
  existence claim.

After activation, the review may open `0027`'s actual construction,
source-access record, result/verdict/validation artifacts, exact verifier and
captured output, and each importable module/test explicitly inventoried by
those artifacts. Any additional primary source must have title, version,
stable locator, hash, exact theorem/page, hypotheses, role, and scope limit
recorded before it is credited.

Accepted Euler-column or compact-vorticity claims may supply only the exact
background fields, modes, tails, or conventions in their accepted statements.
In particular, an accepted stationary column or linear mode does not supply a
nonlinear solitary branch, an unbounded exterior match, or relative-energy
convergence. Active `0022` and other P253 attempts are proposal evidence, not
accepted theorem imports. Standard Sobolev/weighted elliptic estimates,
Sturm--Liouville theory, Fredholm alternatives, and Banach-space implicit
function or Lyapunov--Schmidt theorems are permitted mathematical imports only
after their hypotheses are checked explicitly on the constructed operator and
domain.

## Exact Euler and traveling-frame checks

1. **Physical equation and frame.** Substitute the final velocity and pressure
   into
   `partial_t u+(u dot grad)u+grad(p/rho_m)=0`, `div u=0`, with one consistent
   material density and traveling coordinate. Check the sign of `c`, every
   background-wave cross term, and whether a Galilean frame change is being
   used. A steady traveling-frame equation must reconstruct an exact global-in-
   time translating Euler solution, not a prescribed transport history.
2. **Axisymmetric representation.** Derive the streamfunction, azimuthal
   velocity/circulation, Bernoulli function, vorticity, and
   Grad--Shafranov/Bragg--Hawthorne equation from Euler using the artifact's
   stated conventions. Retain all `1/r` and `1/r^2` terms, background terms,
   pressure integration constants, and density factors. Verify the converse:
   a solution of the reduced scalar equation plus its constitutive functions
   really reconstructs Euler wherever claimed.
3. **Background state.** Verify that `u_bg,p_bg` themselves solve Euler on the
   declared unbounded domain and have the stated axis/exterior behavior. State
   whether the background has finite total energy, finite energy per unit
   length, or infinite energy, and whether its circulation, axial flux, swirl,
   or ambient translation is fixed. No inaccessible-source statement may
   replace this direct check.

## Operator, kernel, and exact existence checks

4. **Linear operator and domain.** Starting from the nonlinear traveling-wave
   map, derive its full Fréchet derivative at the background. Identify the real
   weighted Sobolev/Hölder spaces, axis conditions, decay or radiation
   conditions, dense/operator domain, codomain, and parameter dependence.
   Check formal and actual adjoints in the correct measure. If variables or
   weights conjugate the operator, track every boundary term and singular
   coefficient.
5. **Spectrum and Fredholm property.** Prove the exact simple trapped mode or
   kernel used for bifurcation, its normalization and nodal class, the cokernel,
   closed range and index, and separation from essential/continuous spectrum.
   A radial eigenvalue obtained with a rigid wall cannot be carried to the
   unbounded exterior without its Dirichlet-to-Neumann or matching condition.
   Translation, scaling, flux, circulation, and pressure-gauge zero modes must
   be identified and either retained as parameters or removed by explicit
   phase/constraint conditions.
6. **Nonlinear map and IFT.** Verify that the complete nonlinear residual maps
   the declared Banach neighborhood into the stated codomain with the
   differentiability required by the invoked theorem. Check invertibility on
   the complement, transversality/nondegeneracy, nonlinear remainder estimates,
   constraint preservation, and the exact order of quantifiers among
   background, wave speed/eigenvalue, amplitude, localization scale, and
   regularity. Establish nontriviality of the resulting branch and distinguish
   existence from formal series or approximate residual cancellation.
7. **Homoclinic/spatial-dynamics alternative.** If the proof uses spatial
   dynamics rather than a direct IFT, audit the phase space, reversible or
   Hamiltonian structure, center/stable/unstable splitting, exponential
   dichotomy, homoclinic intersection, and reconstruction with the same rigor.
   A reduced amplitude homoclinic is insufficient until its remainder is solved
   in the full Euler operator.

## Axis, exterior, localization, and energy checks

8. **Symmetry axis.** In cylindrical variables verify Cartesian smoothness at
   `r=0`, not merely bounded scalar coefficients. The streamfunction and swirl
   must vanish to the orders that make `u_r`, `u_theta`, `u_z`, vorticity, and
   pressure smooth and single valued. Integration by parts must have no hidden
   axis boundary term.
9. **Unbounded exterior.** Verify velocity and pressure matching across each
   vorticity/support interface and solve the actual exterior elliptic problem.
   Audit the sign and mapping domain of every exterior Dirichlet-to-Neumann
   operator, far-field multipoles, and decay in both radial and axial
   directions. A cutoff, wall, imposed zero trace, or compact scalar source
   generally produces a noncompact velocity/pressure tail and cannot be called
   vacuum matching without the exact cancellation.
10. **Free interfaces.** If compact vorticity rather than a smooth global
    profile is used, check the material-interface kinematic condition, pressure
    continuity, velocity regularity, jump laws, and the function space in which
    the free boundary is solved. If all cutoffs are flat, verify flatness of the
    reconstructed Euler fields rather than only of the scalar ansatz.
11. **Localization and relative energy.** State exactly what tends to zero and
    in which norm as `|z|` and `r` tend to infinity. Independently derive the
    finite quantity called energy. For constant density this may be the full
    excess

        (rho_m/2) integral (|u_bg+v|^2-|u_bg|^2) d^3x,

    the positive perturbation norm `(rho_m/2) integral |v|^2`, or a renormalized
    Hamiltonian with flux/impulse terms; these are not interchangeable. Prove
    convergence of every cross term and include the cylindrical `2 pi r`
    measure. Record dependence on Galilean frame and background subtraction.
    Finite relative energy does not imply finite total energy of the background.

## API, oracle, and scope checks

12. **Implementation surface.** Inspect the actual `0027` source inventory
    after activation. A module name, metadata row, or test count is not evidence
    until the implementation is opened and matched to the theorem. Reusable
    functions must encode the same operator, boundary maps, branch parameters,
    and domains as the proof; tests may corroborate exact identities but cannot
    manufacture existence by checking literals or prescribed profiles.
13. **Strongest oracle.** The load-bearing oracle is the complete analytic
    existence proof with direct Euler substitution, applicable operator
    theorem, and global axis/exterior estimates. Exact symbolic checks are
    regression or exposing subclaim evidence for signs, conjugations, series,
    matching identities, and residuals. Floating-point spectra, shooting, or
    sampled profiles are not exact existence evidence. No soft numerical
    quantity is registered, so the small-ratio-numerics prescriptions do not
    bind this frozen review.
14. **Sensitivity.** Representative exposing changes are: reverse the
    traveling-frame sign; omit a background-wave cross term; change the
    cylindrical `1/r` coefficient; replace the exterior DtN map by a wall
    condition; violate the axis parity/order; include a nonintegrable
    background cross term in the energy; retain an unremoved kernel in the IFT
    inverse; or claim exact existence from a truncated series. Each chosen
    mutation must be caught by direct calculus, an operator-domain argument, or
    a sensitive executable assertion.
15. **Claim boundary.** Record whether the result is on `R^3`, a half-space,
    an infinite cylinder, a periodic axial domain, or a rigid vessel; whether
    it is axisymmetric and has swirl; and whether localization is absolute or
    relative to an unbounded background. Exact existence gives a translating
    solution for all time but does not imply spectral, orbital, asymptotic, or
    nonlinear stability of nearby data.

## Verdict and completion contract

The coherent `0027` route receives exactly one review verdict:

- `established` for the strongest exact existence statement whose full
  operator, axis/exterior reconstruction, quantifiers, and energy definition
  close;
- `refuted` only if a stated solution fails Euler or an explicit contradiction
  holds under its hypotheses; or
- `blocked` with the exact missing operator, exterior, or existence step when
  the candidate is not yet proved.

A missing accessible copy of Sun's paper is a source-availability limit, not a
refutation. If the self-contained `0027` proof closes the theorem, Sun is not
needed. If only a wall-bounded theorem, formal branch, or exact reduced
subclaim survives, preserve it and identify the minimum extension that would
recover the unbounded solitary wave.

After activation, write `0028/review.md` and `0028/verdicts.yaml` with frozen
hashes, the strongest supported statement, evidence roles, findings and minimum
corrections, four-axis status, exact remaining dependency, and exclusions.
Send artifact/verdict/dependency to `herdr agent prompt w3:p1`. Do not edit
`0027`, central files, source modules, tests, memory, generated documents, or
commits.

Nothing in this review can license solitary-wave stability, an open all-time
neighborhood of initial data, all-time localization of perturbations, a
mechanical restoring force, particle identity, physical or quantum spin,
statistics, electron/neutrino completion, parent-campaign completion, or a
global Euler no-go.

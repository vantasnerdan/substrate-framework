# P253/0018: independent review of the field/action map and global twisted carrier

This is a non-author, fixed-theorem review with **two separately adjudicated
units**. Unit A reviews root-owned attempt 0013's Lamb/Madelung/spinor
field-and-action map and its flat-Heaviside distributional counterexample.
Unit B reviews root-owned attempt 0017's global Hodge-corrected twisted initial
carrier, its same-field invariants and tail, its local-time Euler lift, and its
radial stationary/translation obstruction. Success or failure of one unit does
not propagate to the other.

The accepted base is `v0.183.0` at
`b6fc902a0942d07996f12a81028fbd3f7c909a43`. The review is exact/source-level
and has no empirical comparator or mechanism-selection contest. Root owns the
central proposal, activation receipts, source/API/test changes, memory, and
commits. `particle-balance-review` owns only append-only attempt 0018 review
artifacts and did not author or implement attempts 0013 or 0017, their modules,
tests, or verifier receipts.

## Preregistration and activation boundary

No 0013/0017 scientific body, module, test, output, or result is an admissible
review input until root registers 0018 centrally and
`attempts/0018/activation-schema.exit` is exactly zero. After activation, pin
the frozen README and every reviewed artifact by SHA-256, then perform one
substantive pass and at most one bounded correction check.

During preregistration, a repository citation search intended to exclude both
target directories used a glob that did not match paths prefixed by `./`. Its
printed output exposed isolated line hits from the 0013/0017 directories,
including their README summaries and brief snippets naming the already assigned
constructions. Neither body was opened sequentially; no target equation,
verifier, output, or result was inspected. Because the user had already fixed
both theorem statements and there is no comparator choice, this does not alter
the criteria below. It is nevertheless recorded: the reviewer is independent
and non-author, but not perfectly body-blind.

## Frozen source and dependency inventory

After activation, the review may open the following target evidence exactly at
its claimed role:

- attempt 0013's frozen README, boundary hashes, field/action derivations,
  flat-interface verifier and receipt, `src/substrate_framework/euler_field_maps.py`,
  and `tests/test_euler_field_maps.py`;
- attempt 0017's frozen README, global-carrier derivation and receipts,
  `src/substrate_framework/euler_twisted_carrier.py`, and
  `tests/test_euler_twisted_carrier.py`;
- attempt 0005 and the completed 0011 correction review, only for the exact
  Euler--Poincare/coadjoint-orbit action, Hodge reconstruction with harmonic
  row, and the already-scoped local Clebsch pullback;
- attempts 0007/0008, only for independently checked compact-vorticity
  Biot--Savart tail conventions and finite angular-momentum/helicity
  definitions. These proposal attempts are evidence, not accepted authority.

The frozen external primary inventory is:

- H. Marmanis, `hep-th/9602081`, and the associated *Physics of Fluids* paper,
  DOI `10.1063/1.869762`, for the precise Euler/Lamb-vector Maxwell analogy,
  including any erratum affecting the current definition;
- E. Madelung, *Quantentheorie in hydrodynamischer Form*, DOI
  `10.1007/BF01400372`, for the historical scalar transform, with a verified
  translation if the original text is inaccessible;
- D. Fusca, *The Madelung transform as a momentum map*, DOI
  `10.3934/jgm.2017006`, for the modern momentum-map domain and symplectic
  statement;
- the primary spinor-hydrodynamic source actually invoked by 0013, including
  Meng--Yang `arXiv:2302.09741` if used, only for its stated map and not as an
  automatic incompressible-Euler or quantum interpretation;
- N. D. Mermin and T.-L. Ho, *Circulation and angular momentum in the A phase
  of superfluid helium-3*, DOI `10.1103/PhysRevLett.36.594`, only if a
  Mermin--Ho/Berry-curvature identity is imported rather than directly derived;
- T. Kato, *Nonstationary flows of viscous and ideal fluids in R3*, *Journal
  of Functional Analysis* 9 (1972), for any local classical incompressible
  Euler existence/uniqueness transfer. The exact Sobolev index, decay,
  divergence-free hypothesis, and persistence conclusion must be checked
  against the theorem actually used.

Hodge/Leray projection, Newton-potential asymptotics, distributional
differentiation of a flat Heaviside interface, and the displayed action
pullbacks should be derived directly. A target citation not listed above may be
inventoried after activation only at its narrow declared role; it cannot
silently strengthen either unit or import a quantum/stability conclusion.

## Unit A: frozen field/action-map checks

Unit A's strongest possible result is an exact coordinate/action statement,
not a universal equivalence between Euler flow and quantum mechanics.

1. **Lamb/Maxwell conventions.** Starting from incompressible Euler with fixed
   material density, derive the Lamb vector, Bernoulli scalar, vorticity
   evolution, and every proposed electric/magnetic charge/current analogue.
   Check all curl/divergence signs, density and `4*pi`/unit conventions,
   boundary terms, full-space integrals, and whether Ampere-like closure or a
   propagation speed is an identity, an added constitutive law, or absent.
2. **Scalar Madelung domain and action.** On the nonvanishing scalar-wavefunction
   domain, independently expand
   `psi=sqrt(n) exp(iS/hbar)` in the Schrödinger action. Check the continuity
   term, phase kinetic term, density-gradient energy and quantum-pressure
   variation with exact coefficients. Keep probability density `n` distinct
   from Euler's constant material mass density. Determine precisely where
   nodes, multivalued phase, circulation, or nonzero Euler vorticity obstruct a
   global scalar chart.
3. **Normalized spinor map.** For the claimed unit spinor and phase variables,
   rederive the Berry one-form/velocity, curvature/vorticity, texture-energy
   identity, normalization constraint, gauge redundancy, and local chart
   transition. Audit whether the proposed action is the pullback of the actual
   Euler kinetic Hamiltonian and coadjoint symplectic structure, with pressure
   and Hodge reconstruction retained, rather than the Pauli/Schrödinger
   Hamiltonian substituted because the coordinates are spinorial.
4. **Actual dynamics.** Verify any advected-spinor or label equation by
   substituting the reconstructed velocity into full nonlinear Euler and by
   varying the complete action. A representation identity at one time does not
   prove that an arbitrary Schrödinger/Pauli evolution stays on the Euler
   solution manifold.
5. **Flat-Heaviside counterexample.** Compute distributional first and second
   derivatives at the flat interface, including the `delta` and `delta-prime`
   terms and any undefined nonlinear distribution products. Check the exact
   residual, sign, coefficient, matching conditions, and test-function action.
   The counterexample may refute the claimed sharp-interface transfer; it does
   not by itself refute smooth interfaces, every multicomponent chart, or all
   possible quantum completions.
6. **Interpretive boundary.** Separate invertible classical coordinates,
   canonical action, and conserved observables from Hilbert-space
   superposition, a linear quantum Hamiltonian, Born probabilities,
   measurement, statistics, or a selected value of `hbar`.

Unit A receives its own route verdict—`established`, `refuted with mechanism
named`, or `blocked with the missing construction named`—and its own four-axis
verification/review/compatibility/epistemic statuses.

## Unit B: frozen global-carrier checks

Unit B's strongest possible result is a smooth finite-energy initial carrier
with an actual short-time Euler solution and exact initial invariants. It is not
a stationary, orbitally stable, or particle theorem.

1. **Global spinor and regularity.** Check the radial `SU(2)`/unit-spinor profile
   at the origin, every transition radius, and infinity. Verify normalization,
   single-valuedness, differentiability to the Sobolev order used for Euler,
   compact support of the derived vorticity, and absence of hidden axis,
   chart, or distributional singularities.
2. **Actual Hodge correction.** Derive the raw Berry one-form, its divergence,
   and the decaying whole-space Leray/Hodge phase solve. Confirm the corrected
   velocity is divergence-free, has exactly the claimed curl, lies in `L2`,
   has finite kinetic energy/action, and uses the unique decaying harmonic
   convention. A compact raw profile does not license compact velocity.
3. **Same-field invariants.** Compute helicity and axial intrinsic angular
   momentum from the same corrected physical velocity and vorticity, with the
   same origin, mass density, orientation, and sign conventions. Establish
   convergence and nonzero values analytically; do not combine the helicity of
   one representative with the moment of another gauge/Hodge representative.
4. **Far-field tail.** Expand the corrected Biot--Savart/Newton field from its
   actual compact vorticity moments. Check the leading power, tensor/angular
   coefficient, impulse normalization, remainder hypotheses, and whether a
   vanishing leading moment changes the asserted tail. No Coulomb, force, or
   compact-velocity inference follows.
5. **Local-time Euler lift.** Map the initial field into the exact natural
   local classical Euler class: divergence-free `H^s(R3)` for an explicitly
   sufficient `s`, finite energy, and the stated decay. Check the imported
   theorem's interval and uniqueness. If spinor/label variables are claimed to
   persist, verify their advective reconstruction and gauge/Hodge solve for the
   same Euler solution; ordinary local well-posedness alone evolves `u`, not an
   unproved global spinor chart.
6. **Action and leaf accounting.** Compare the global construction with the
   corrected 0005 Euler orbit action and 0013's scoped spinor pullback. Track
   helicity, circulation, impulse, and any topological integer under permitted
   variations. A dilation that changes conserved leaf data is a family of
   initial states, not automatically a stability direction or a scale
   selection mechanism.
7. **Stationary/translation obstruction.** Substitute the radial candidate
   into stationary Euler and into a translating-frame relative-equilibrium
   equation, retaining pressure and all Hodge terms. State the exact ansatz and
   variation class touched by any contradiction. An obstruction for that
   radial profile is not a universal nonexistence or instability theorem for
   twisted/swirl carriers.

Unit B receives a separate route verdict and separate four-axis statuses. A
local-time solution preserves existence of the candidate over that interval;
it does not establish nonlinear stability, persistent localization, a compact
rotation orbit, exchange dynamics, or all-time regularity.

## Frozen oracle and sensitivity plan

Exact differential-form/vector-calculus derivation, distribution theory,
direct invariant integration, full Euler residuals, and applicable primary
theorems are the strongest oracles. Existing tests may corroborate these
statements but cannot prove global chart coverage, action descent, source
applicability, or stability merely by passing.

After activation, audit the recorded processes and check-call sensitivity once.
Representative exposing mutations are:

- reverse one Lamb-vector, Berry-curvature, or Hodge-projection sign;
- remove or alter the scalar/spinor gradient-energy coefficient;
- suppress the flat-interface `delta-prime` contribution;
- omit the Leray correction when the raw carrier has nonzero divergence;
- change one helicity/angular-momentum density or tail normalization; and
- weaken the local Euler Sobolev/domain hypotheses or promote its finite
  interval to a global one.

Each mutation must be caught by the relevant exact derivation, test, or source
applicability check. No production numerical remainder, soft eigenvalue,
small force, or energy splitting is registered, so the small-ratio numerical
prescriptions do not bind at freeze. If a target artifact makes such a claim,
that evidence is out of the frozen review scope until separately preregistered.

## Completion and non-inheritance boundary

Complete 0018 with one review artifact and one verdict artifact containing two
independent unit decisions, the strongest supported statement for each,
minimum corrections supported by concrete evidence, hashes, oracle roles, and
the exact remaining dependency. Do not reopen corrected 0011 or unrelated P253
attempts.

Neither unit may license a universal Euler--quantum equivalence, quantum
Hilbert space, Born rule, measurement theory, fermionic statistics, physical
spin one-half, a selected action/length/mass scale, charge, an invariant
propagation speed, nonlinear or orbital stability, all-time localization,
generic three-dimensional regularity, an electron or neutrino, parent campaign
completion, or a global no-go. Write only attempt 0018 review artifacts and
send artifact/verdict/dependency to `herdr agent prompt w3:p1`.

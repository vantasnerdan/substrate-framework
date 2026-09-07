# P253/0085 — joint charged-Cao Hessian and full-3D persistence

Owner: `root`

Status: **centrally registered and schema-activated; one bounded
postactivation spectral-scope repair is recorded and replayed before claim
bodies.**  The exact helper/test skeleton for the field identities was created
in the delivery race before that repair arrived; no derivation, result,
validation, verifier run, or numerical execution existed.

## Frozen objective and authority

The full Issue #203 electron and neutrino campaign remains active.  This
attempt advances the same-carrier P2/P3/P5 rung after the existence work in
P253/0077 and P253/0080.  It asks whether a small-charge Cao relative
equilibrium has a useful constrained restoring structure and whether any such
structure survives in the full three-dimensional, nonaxisymmetric joint
Euler--Maxwell--tag dynamics.

Reviewed P253/0068/0071 supplies the exact constrained `U(1)` extension,
joint action and balances, local constraint-preserving flow, and retained
history.  Reviewed P253/0077/0081 supplies the traveling Lorenz-gauge Maxwell
elimination and exact axisymmetric no-swirl steady equations, but leaves the
charged branch conditional.  P253/0080 may be used only at final author scope
pending independent review: its claimed output is a finite-subluminal-window
charged branch from two bordered derivatives, not a stability theorem.
Reviewed P253/0074/0078 supplies two positive-Krein uncharged Cao modes and a
fixed-fiber graph/Riesz transfer; it does not supply global constrained
coercivity.

This is a `candidate_comparison` and exact Hessian/spectral classification
attempt.  It neither adopts the `U(1)` extension into bare Euler nor selects
charge, action, Planck's constant, statistics, a Lorentz cone for the coupled
pressure system, an electron, or a neutrino.

## Exact joint phase space and relative functional

Fix one charged branch member with `|c_g|<c_EM`, fixed coupling constants
`g, epsilon_EM, mu_EM`, fixed circulation and tag leaf, and total charge
`Q_g=g integral chi_g`.  The tangent space must separately declare:

1. dynamically accessible fluid variations on the Cao vorticity orbit,
   including the stabilizer quotient, circulation, impulse/center, and
   translation rows;
2. material-tag variations induced by the same displacement,
   `delta chi=-xi dot grad chi_g`, with any additional semidirect-product
   relabeling directions stated explicitly;
3. Maxwell variations satisfying the instantaneous phase-space constraints
   `div delta B=0` and
   `epsilon_EM div delta E=g delta chi`, with gauge and Coulomb-affine rows
   quotiented or fixed; and
4. the exact axis, interface, exterior-decay, and weighted norms in which the
   quadratic form and linearized generator are closed.

Use the total relative functional, not fluid kinetic energy alone,

    H_rel = H_fluid + H_EM - c_g P_total,z
            - (circulation/tag/center multipliers).                 (1)

Ampere/current equations determine the equilibrium and evolution.  They are
not added to the instantaneous Hessian constraint set.

The sign convention for total field momentum is fixed so that the exact
field quadratic completes as

    E_EM-c_g P_EM,z
      = (epsilon_EM/2)||E+c_g e_z cross B||_2^2
        + (1/(2 mu_EM))(1-c_g^2/c_EM^2)||B_perp||_2^2
        + (1/(2 mu_EM))||B_z||_2^2,                                (2)

where `c_EM^2=1/(epsilon_EM mu_EM)`.  Equation (2) must be rederived from the
joint action and the adopted momentum sign.  With a fixed positive
subluminal margin it is a positive field block; it is not by itself a fluid
or tag coercivity theorem.

Any field elimination must be the explicit minimization/Schur complement of
the closed quadratic form over the Gauss and `div B` constraint fiber on the
affine weighted spaces.  The resulting charge block must be the anisotropic
comoving `H^{-1}` form, with its sign and two-sided norm bounds derived.
Matter-current terms then enter the full reduced Hessian through the declared
joint tangent and equilibrium map.  The
induced fluid/tag quadratic correction is `O(g^2)`, but its sign, kernel, and
operator norm are outputs.  Static Coulomb intuition may check the answer and
may not replace this derivation.

## Route A — axisymmetric semidirect-product Hessian

On the exact charged steady branch, derive the second variation of (1) in the
axisymmetric no-swirl constrained tangent.  Carry the dependence
`chi_g=F(I_g)` and the band-local stabilizer condition where both `zeta` and
`zeta'` are nonzero.  Prove or refute:

- closedness and finite-codimension slicing of the joint quadratic form;
- positivity, finite Morse index, or an explicit negative direction after
  quotienting true translations and gauge rows;
- the exact Maxwell Schur correction relative to the uncharged Cao
  constrained form; and
- control of `delta chi` in the topology actually supplied by the energy.

At `g=0` the material tag is passive and has infinitely many neutral
directions.  For `g!=0`, the natural electric control is an
`H^{-1}`-type charge norm with coefficient `O(g^2)`.  A tag Casimir defines a
leaf but contributes no physical restoring energy.  Route A must state
whether its coercivity constant degenerates as `g->0` and must not upgrade
this control to an unweighted `L2` or Sobolev tag norm without a separate
estimate.

If the form has a finite negative or neutral sector, compute that sector and
continue by a finite-dimensional KKS/Feshbach reduction rather than treating
failure of coercivity as failure of the carrier.

## Route B — full three-dimensional spectral persistence

Axisymmetric Cao variational control does not answer P2.  Construct the full
three-dimensional linearized joint generator on one common constraint
domain.  Separate toroidal Fourier sectors, the passive-tag block, Maxwell
wave branches, whole-space pressure/Hodge coupling, and the finite symmetry
rows.  Establish the exact Hamiltonian/Krein identity on this domain.

At `g=0`, the full decoupled joint generator already contains whole-space
Maxwell radiation with spectrum `i R` and the passive-tag transport block,
which generally also has continuous or essential spectrum.  Therefore a Cao
internal eigenvalue `i omega` is embedded in the joint continuum even when it
is isolated in one fluid fiber.  Ordinary common-domain Kato/Riesz persistence
may be used directly only for an unstable eigenvalue with nonzero real part.
Neutral and positive/negative-Krein imaginary modes require one of:

- an actual invariant symmetry sector that decouples every open channel;
- a proved full-joint gap;
- a weighted limiting-absorption/Feshbach resonance construction with its
  radiation condition; or
- a radiation-damping/Fermi-golden-rule calculation.

Test off-axis unstable modes by the direct perturbative route and embedded
imaginary modes by those open-channel alternatives.  If an axisymmetric
positive block is found, identify what additional estimates are needed before
it controls nonaxisymmetric sectors.  A full-3D positive result requires
either a constrained coercive estimate on an open real same-leaf neighborhood
or a spectral/resonance decomposition with all nonpositive and radiating
sectors controlled.  The two positive modes from 0074/0078 are inputs to this
classification, not isolated eigenvalues of the decoupled joint operator.

The strongest Route-B success is a full joint spectral/restoring or controlled
resonance theorem for the charged carrier.  Persistence of one off-axis
unstable mode refutes small charge as its cure.  A nonzero on-shell radiation
coupling for an embedded imaginary mode instead gives a resonance/damping
verdict and activates Route C without being mislabeled as eigenvalue
persistence.

## Route C — failure-derived persistence mechanisms

If passive-tag degeneracy or a full-3D negative sector defeats ordinary
coercivity, execute at least these materially different continuations:

1. keep the same charged relative equilibrium and derive a center-stable,
   orbital-modulation, or finite-time metastability theorem with radiation
   retained;
2. test whether the reviewed positive-Krein internal modes admit a finite
   invariant/approximately invariant symplectic block whose complement has a
   signed energy estimate; and
3. preserve the bare-Euler weighted-tail and Cao-mode routes independently,
   because success or failure of the declared `U(1)` extension does not close
   the original Euler objective.

Any nonlinear persistence statement must track charge localization and tag
filamentation in the norm it claims.  A small-`g` existence IFT, a tag
stabilizer identity, or an axisymmetric Hessian does not by itself establish
all-time particle persistence.

## Evidence and verdict boundary

The analytic proof carries the functional domains, constraint reduction,
Schur complement, and full-3D operator result.  Exact symbolic checks may
rederive (2), block-Hessian parity in `g`, finite-dimensional Schur signs, and
Krein identities, but cannot establish closedness, coercivity, or spectral
completeness.  No production numerics are preregistered.  If a soft spectral
quantity later becomes load-bearing, the small-ratio-numerics contract must
be appended before its verifier is designed.

Each route receives exactly one scoped verdict.  The parent campaign remains
active after this attempt unless the complete electron and neutrino success
contract or a separately reviewed scientific-exhaustion certificate is
earned.

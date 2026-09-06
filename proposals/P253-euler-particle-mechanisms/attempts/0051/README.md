# P253/0051: physical Schwinger--Hopf analyzer on an Euler mode doublet

Root preregisters a materially different P4 continuation after `0043/0047`
and `0046` exposed the distinction between one classical KKS oscillator and
a physical quantum two-state sector.  The algebraic candidate is the
Schwinger construction: two positive canonical modes form `z in C^2`; fixing
their total action and quotienting the common phase gives `CP1 = S2`, while
the Stokes bilinears close `su(2)`.  The scientific target is stricter than
that identity.  This attempt asks whether the two modes, their mixing, their
action constraint, and their readout are all realized by one autonomous Euler
carrier, with a controlled complement.

The accepted base is `v0.183.0`.  `0030/0035` may supply only their reviewed
positive axisymmetric column metric and mode-domain statements.
`0043/0047` supplies the reviewed quantization/nonselection boundary.
`0046/0050` may be used after review only at its exact compact-pair KKS and
Heisenberg scope.  The active `0048` ring construction is neither authority
nor evidence for a doublet until independently reviewed.  No electron,
neutrino, measured constant, `hbar`, Born rule, or target history is an input.

## Frozen positive object

On one persistent localized Euler carrier, construct all of the following:

1. two linearly independent, positive-energy, dynamically accessible normal
   modes with the same physical frequency and exact KKS normalization;
2. an invariant or quantitatively controlled four-real-dimensional mode
   bundle for the autonomous full Euler evolution, including pressure and a
   bound on leakage into the complement;
3. a physical family of Euler/analyzer operations inducing noncommuting
   `U(2)` mixing on that bundle, rather than arbitrary coordinate changes;
4. conservation or physical selection of the total mode action used in the
   common-`U(1)` reduction;
5. an actual same-substrate observation map for the three Stokes components,
   with the retained unresolved state and complete pressure dependence; and
6. a mechanism selecting the action unit, event probabilities, composition,
   exchange character, and a finite propagation band, or a precise ledger
   identifying which of these require a new law.

Items 1--5 would establish a classical physical `CP1` analyzer sector.  LP4
requires item 6 as well.  The quotient identity by itself does not select a
quantum state, spin one-half, fermionic statistics, or Born probabilities.

## Route A: two internal modes of one carrier

Search the actual positive spectrum of the same persistent carrier for a
symmetry-protected two-mode degeneracy.  Construct physical modes `e_1,e_2`
with

    Omega(e_a,J e_b) = B delta_ab,
    H_2 = omega B (|z_1|^2+|z_2|^2),

on a declared real energy/graph domain.  Derive the full pressure-coupled
spectral projection and prove that geometric symmetries or an exact internal
operator, rather than tuned initial data, preserve the doublet.  Calculate
all complement and symmetry-zero couplings.  If degeneracy splits, retain the
actual two-frequency dynamics and determine whether a controlled analyzer
operation can still mix the modes.

## Route B: two identical localized copies

Use two separated copies only if their product is an actual same-field Euler
sector with controlled interaction.  Derive the two KKS blocks, total action,
frequency splitting, exchange topology, and leakage.  A direct product at
infinite separation is not a single-particle two-state space, and a label
permutation is not a physical exchange.  This route succeeds only if a
localized operation coherently mixes the copy amplitudes while preserving a
common carrier and controlled complement.

## Route C: prepared C-CST polarization channels

Test whether the accepted two-polarization finite-window response can supply
an exact or controlled mode bundle.  Preserve its preparation, nonzero-wave,
finite-window, and source-cost hypotheses.  An invertible source-to-current
matrix earns a classical analyzer map; it becomes autonomous state dynamics
only if the prepared controls are replaced by one invariant Euler evolution
and the source/complement leakage is bounded.

## Exact Schwinger--Hopf calculation

For any physically earned canonical doublet, first use the positive Hessian
and its compatible complex structure to put each mode in physical
energy-action coordinates.  KKS normalization alone is insufficient because
`q->c q`, `p->p/c` preserves the symplectic form.  For a diagonal normal form
use `Q_a=sqrt(nu_a)q_a`, `P_a=p_a/sqrt(nu_a)` with every density and unit
factor retained.  Then independently derive

    J = |z_1|^2+|z_2|^2,
    S_i = z^dagger sigma_i z/2,
    {S_i,S_j} = epsilon_ijk S_k,
    |S|^2 = (J/2)^2.

Fix `J=J_0>0` and divide by the common phase only after showing that both are
physical invariants/gauges on the actual Euler bundle.  Record every density,
`2*pi`, KKS, and sign factor before rescaling to the displayed convention.
The strongest exposing checks are: unequal frequencies, a nonzero
mode-complement block, failure of exact total-action conservation, and the
finite-dimensional CCR trace obstruction.

For a compressed real Hamiltonian perturbation `V` on the doublet, derive

    V_C=(V-J V J)/2,     V_A=(V+J V J)/2.

The first part commutes with the compatible complex structure and is
number-preserving; the second anticommutes and produces squeezing/action
drift.  A physical analyzer needs `V_A=0` or a bound over its actual gate
time, plus two number-preserving traceless compressions with noncollinear
Pauli vectors.  Static isolation or one avoided-crossing axis is insufficient.

## Measurement, scale, and propagation continuation

If the physical `CP1` sector is earned, compose it with the separately
reviewed `0049/0050` detector theorem only after deriving a same-field
analyzer-to-port flux.  A deterministic first-exit construction must state
its invariant measure or sampling ensemble and derive channel cross-sections;
an exponential clock remains an added hypothesis unless Euler supplies its
hazard law.  Reversible measure-preserving Euler cannot furnish an attracting
autonomous reset on a finite invariant phase volume; scattering to radiation
or external repreparation must be distinguished.

Action normalization and propagation remain independent tests.  Continuous
Euler scaling must be confronted directly.  A common rate constant that
cancels from event probabilities does not select an action.  A finite-speed
band must come from the same carrier/background spectrum with a controlled
regime; incompressible elliptic pressure and Galilean boosts are retained in
the comparison.

## Verification and verdict boundary

The analytic oracle is the exact KKS/Hessian/spectral projection and the
physical mixing/observation map.  A small matrix is evidence only after its
full Euler invariant-bundle and leakage estimate is established.  The
Schwinger algebra will be implemented in an importable API with mutation and
domain tests, but code cannot establish the carrier realization.

Every route receives one verdict.  A failure activates a new physical
two-state candidate or an explicit minimal-foundation comparison in-run.  No
formal `CP1`, conditional quantization, stochastic clock, or finite-order
carrier expansion earns LP4, particle identity, or campaign completion.

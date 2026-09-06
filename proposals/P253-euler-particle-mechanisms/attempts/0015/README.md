# P253/0015: periodic-pair variational or Floquet restoring construction

## Activation state

This README preregisters the 0015 route only. No García–Hassainia–Hmidi
source body, new derivation, verifier, or numerical result may be opened or
produced in 0015 until root centrally registers this attempt and records a
repository-schema exit of zero here. Root owns central proposal/API work and
commits; 0015 owns only its append-only attempt artifacts after activation.

Target kind: fixed periodic-orbit stability/restoring construction with two
competing proof representations. Route A tests an autonomous constrained
energy–impulse maximizer. Route B is the preregistered representation change
if the actual non-equilibrium periodic motion defeats Route A.

## Frozen source inventory before bodies

No source body was opened while writing this contract. Inventory is frozen
from the user-supplied lead and the already durable 0002 source receipt:

| ID | Source | Frozen version/provenance | Planned post-activation audit |
|---|---|---|---|
| GarciaEtAl2603.21644 | García–Hassainia–Hmidi, *Time-periodic leapfrogging vortex rings in the 3D Euler equations* | arXiv:2603.21644v1; cached-PDF SHA-256 7ba22ef57ba9453dd2631c33b1ad094a35ce2be171426cecce80fba8a5cdfac6; 150 pages | Exact Theorem 1.1 quantifiers, patch amplitudes/areas, label transport, translating speed, period, regularity, symmetry quotient, parameter range, and whether the solution is nonconstant in the translating frame |
| CaoEtAl2206.10165 | Cao–Lai–Qin–Zhan–Zou, *Uniqueness and stability of steady vortex rings for 3D incompressible Euler equation* | arXiv:2206.10165v2; cached-PDF SHA-256 6d90be6b7aab2ea7d92dba843b06e72e319cad50b0a2e39cfe5589cf6aa528ca | Exact rearrangement/maximizer hypotheses that can transfer to a labeled pair, without transferring steady-orbit conclusions to a periodic orbit |
| P253/0002 | Published-source transfer ledger | attempts/0002/source-provenance.yaml | Verify the García page map: Hamiltonian formulation pp. 3–4, translating-frame period pp. 5–6, Theorem 1.1/Remarks pp. 6–7, function spaces from p. 48 |
| P253/0012 | Nearby same-leaf deficit construction | attempts/0012/construction.md | Reuse only exact label-rearrangement, total-energy/impulse, and cross-energy identities; do not assume its one-ring modulus is a pair-orbit modulus |

Any version change, corrected theorem, or source-family mismatch discovered
after activation is recorded before transfer. Full papers remain in the
existing /tmp/p253-primary-source-cache, not in campaign artifacts.

## Frozen actual claim

Let the published solution, after the exact source normalization is audited,
be a labeled pair of axisymmetric vortex patches

\[
 q_*(t)=(\zeta_{1,*}(t),\zeta_{2,*}(t)),
\qquad q_*(t+T)=q_*(t)
\]

in a frame translating axially at the source speed \(c\). The laboratory
solution may therefore be relative-periodic, while \(q_*\) is periodic in the
translating frame. The proposed stability object is the full orbit set

\[
 \Gamma=\left\{
 (\zeta_{1,*}(t,\cdot+a e_z),\zeta_{2,*}(t,\cdot+a e_z)):
 t\in\mathbb R/T\mathbb Z,\ a\in\mathbb R
 \right\},
\]

not one frozen phase and not a prescribed center history.

The strongest requested 0015 claim is:

> On an explicitly defined labeled pair-rearrangement class containing the
> actual published patches and invariant under the common Euler flow, either
> (A) the entire phase/translation orbit \(\Gamma\) is a compactness-controlled
> maximizing set of an actual conserved translating-frame functional, or
> (B) the actual periodic orbit is a critical loop of the Euler
> coadjoint/KKS action and has a constraint-preserving modulated
> action/Floquet estimate that yields a stated orbital restoring conclusion
> modulo time phase and axial translation.

“Restoring” must be stated at the strength actually earned: nonlinear orbital
stability of \(\Gamma\), a coercive modulated action, or a spectrally stable
Floquet return map are distinct statements. Spectral Floquet stability alone
does not become nonlinear orbital stability.

## Common exact objects to construct after activation

The first post-activation work freezes source conventions and then defines
the labeled pair leaf, provisionally

\[
 \mathcal A_{\rm pair}
 =\mathcal R(\zeta_{1,*}(0))
  \times\mathcal R(\zeta_{2,*}(0)),
\]

with the exact weighted measure, patch amplitudes, areas/circulations,
regularity, disjointness or topology conditions, and weak closure dictated by
the source. It must be proved—not assumed—that the common Euler transport
preserves each factor and the additional patch conditions used by the
argument.

On that class, construct from actual Euler invariants

\[
 \mathscr H_c(\zeta_1,\zeta_2)
 =E(\zeta_1+\zeta_2)-cP(\zeta_1+\zeta_2),
\]

including the full self and cross energies. Verify conservation of total
energy and total axial impulse under the common velocity while retaining
individual impulse exchange. No separate conservation of \(P_1\), \(P_2\),
self energies, or cross energy may be inserted.

The first-variation test uses independent source-admissible, labelwise
measure-preserving perturbations and removes only symmetries actually present
in the source. Any weak closure used for compactness must be distinguished
from the exact transported patch class.

## Route A — autonomous pair maximizer

Route A asks whether

\[
 \Gamma\stackrel{?}{=}
 \operatorname*{arg\,max}_{q\in\mathcal A_{\rm pair}^{\,w}}
 \mathscr H_c(q)
\]

or is at least one connected component of that maximizing set, with every
maximizing sequence compact modulo the joint phase and axial-translation
actions.

The audit order is:

1. Prove the exact pair leaf and conservation law.
2. Compute the complete first variation of \(\mathscr H_c\) under independent
   patch-boundary/isovortical variations at an arbitrary phase of \(q_*\).
3. Determine whether every phase is critical on the full invariant leaf.
4. Only if criticality holds, determine maximum versus saddle and seek
   compactness/coercivity modulo the actual neutral directions.
5. Derive the strongest orbit-set stability statement from that structure.

Relative periodicity is not relative equilibrium. A time-periodic Euler
solution merely keeps every conserved functional constant along its
trajectory; it does not follow that the functional derivative vanishes in
all leaf directions. In the translating frame, if \(\mathscr H_c\) is the
Hamiltonian on a nondegenerate coadjoint leaf, an actual critical point
generates equilibrium. Therefore a nonconstant \(q_*(t)\) can be a family of
autonomous maximizers only if the exact constrained geometry supplies a
degeneracy or larger maximizing manifold consistent with its nonzero
Hamiltonian motion. This implication must be checked from the source
variables, not assumed either way.

Route-A success requires all of:

- source-admissible criticality at every phase;
- a maximum characterization rather than constancy along the orbit;
- compactness of maximizing sequences modulo exactly the time-phase and
  axial-translation directions;
- a distance and perturbation class that preserve both actual patch labels;
- an orbital conclusion for \(\Gamma\), not for a frozen surrogate.

If first variation is nonzero in one admissible direction, Route A receives
the route-scoped verdict “refuted with the noncritical direction named,” and
Route B starts in the same activated attempt. Failure of a Hessian estimate
without a first-variation obstruction is instead “blocked with the missing
coercivity/compactness construction named.” No uniqueness-to-Hessian shortcut
is allowed.

## Route B — periodic-orbit action and Floquet/modulation

Route B is the preregistered representation change, not a weakened relabeling
of Route A. It preserves the actual published patch orbit and works on
periodic histories in the same coadjoint leaf.

Construct the translating-frame loop action

\[
 \mathscr A_T[q]
 =\int_0^T\left(\Theta_q(\dot q)
 -\mathscr H_c(q)\right)\,dt,
\]

where \(\Theta\), or an equivalent boundary-contour/KKS primitive, must be
derived with its gauge and patch-label domain. Show that the published
periodic pair is a critical loop with periodic endpoint conditions. If no
global primitive exists, formulate the action by symplectic area with the
extension dependence proved harmless, or use the exact contour-dynamics
linearization directly.

Then construct the linearized period map

\[
 \mathcal M=D\Phi_c^T(q_*(0))
\]

on the source-admissible tangent space. Account explicitly for neutral
directions from time phase and axial translation, any circulation/area
constraints, and any additional source symmetry. A gauge or symplectic slice
must remove neutral directions without deleting physical patch-shape modes.

The analytic ladder is:

1. exact linearized contour/Euler equation and invariant tangent space;
2. symplectic or conserved quadratic structure of \(\mathcal M\);
3. identification and multiplicity of neutral Floquet multipliers;
4. second variation of \(\mathscr A_T\) or a modulated energy controlling the
   complementary tangent space;
5. only then, a nonlinear iteration/normal-form estimate if nonlinear orbital
   stability is claimed.

A unit-circle Floquet spectrum with unresolved Jordan blocks earns only
spectral evidence. A coercive loop-action second variation must still handle
phase, translation, patch-boundary regularity, and the continuum spectrum.
Any soft multiplier, Hessian edge, force, or small action splitting proposed
for numerical evaluation first triggers the small-ratio-numerics skill and
a separate analytic-closure receipt; no production numerics are licensed by
this README.

## Frozen route comparison

Selection criteria, fixed before the source body is opened:

1. exact match to the published patch family and quantifiers;
2. invariance of the perturbation class under the full common Euler flow;
3. use only actual total Euler invariants and source symmetries;
4. honest treatment of periodic phase versus equilibrium;
5. retention of cross interaction, pressure, and all patch-shape modes;
6. strongest restoring conclusion with the fewest added assumptions;
7. an analytic path to compactness or nonlinear control, not a finite-mode
   surrogate.

Route A is preferred only if the criticality/maximizer test passes. Route B
is preferred after a genuine noncriticality or maximizer obstruction because
periodic-loop criticality is the natural variational object for a
non-equilibrium periodic orbit.

## Verdict and evidence contracts

Each tested route receives exactly one verdict:

- established: the route's complete stated claim is proved;
- refuted: an exact admissible variation or structural implication excludes
  that route, with the mechanism named;
- blocked: the exact missing construction is named without promoting an
  absence of proof to a no-go.

Required evidence after activation:

- source receipt with URL, version, hash, and exact theorem/page mapping;
- exact field/measure/energy/impulse and patch-label convention table;
- pair-leaf preservation proof;
- first-variation calculation for Route A;
- if Route A fails, loop-action or contour-linearization derivation for Route
  B and its strongest justified restoring statement;
- exact symbolic checks for tractable identities;
- one explicit next executable construction for any remaining all-time or
  nonlinear-stability dependency.

Maximum 0015 verdict is an axisymmetric patch-class pair-orbit restoring or
stability theorem at the exact strength proved. This attempt cannot license
full three-dimensional/swirl stability, a smooth Cao-family theorem, a
restoring frequency without a physical mode and kinetic metric, helicity,
intrinsic rotation, quantum spin or spin-\(\tfrac12\), electron/neutrino
identification, parent completion, or a global no-go.

## Activation handoff

Root should append 0015 to the central candidate universe/obligation graph,
validate the repository, and place the schema receipts in this directory.
Only after activation-schema.exit records 0 does substantive 0015 work begin.

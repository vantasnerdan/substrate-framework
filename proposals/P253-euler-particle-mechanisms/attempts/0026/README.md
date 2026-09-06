# P253/0026: independent operator-transfer review of the Gavrilov return and GHH monodromy

This is one non-author fixed-theorem operator-review transaction with **two
separately adjudicated units**. Unit A reviews carrier-foundations attempt
`0019`'s actual Gavrilov whole-space axisymmetric Euler field, Leray/Kelvin
return construction, and finite-time packet transfer. Unit B reviews
carrier-owned attempt `0021`'s actual García--Hassainia--Hmidi two-label contour
initial-value/Jacobi propagator and monodromy construction. Neither unit may
inherit evidence, a verdict, an operator domain, or parent-object completion
from the other.

The accepted base is `v0.183.0` at
`b6fc902a0942d07996f12a81028fbd3f7c909a43`. The review has no empirical
comparator or mechanism-selection contest. Root owns the central proposal,
activation receipts, source/API/test edits, memory, and commits.
`particle-balance-review` owns only append-only attempt `0026` review artifacts
and did not author or implement either target unit, its module, tests, verifier,
or receipts.

## Preregistration and activation boundary

No `0019` or `0021` README, scientific body, equation, verifier, receipt,
module, test, output, or result is an admissible review input until root
registers `0026` centrally and `attempts/0026/activation-schema.exit` is exactly
zero. After activation, pin this README and each reviewed artifact by SHA-256,
perform one substantive pass over both units, and permit at most one bounded
correction check per unit when concrete evidence requires it. Until activation,
this README is the only file owned by `0026`.

The target directories were queried only for filenames before this freeze.
During the preceding `0023` task, repository searches exposed isolated
`0021/README.md` snippets, and a root-pane status message exposed a high-level
claim that `0021` had constructed a bounded one-period propagator and a reduced
symplectic map. No `0021` body, equation, verifier, or result artifact was
opened. The user's present fixed scope already names the actual two-label IVP
and monodromy, so this exposure selects no comparator or proof route; it is
nevertheless recorded. No `0019` body or result has been exposed.

Attempts `0024` and `0025`, including any equations, evidence, or conclusions,
are excluded. This review does not expand into root's concurrent `0022` work.

## Frozen Unit A source and target inventory

After activation, Unit A may open only `0019`'s frozen README,
`access-inventory.md`, derivation, result/validation records, exact return
verifier and captured first-run receipt, together with the actual public module
and tests it names. Directory metadata identifies the possible reusable surface
as `src/substrate_framework/euler_core_packet.py` and
`tests/test_euler_core_packet.py`; ownership and exact claim role must be
verified rather than inferred from the filename.

The frozen primary inventory for Unit A is:

- A. V. Gavrilov, *A steady Euler flow with compact support*,
  [`arXiv:1810.08020v1`](https://arxiv.org/abs/1810.08020v1), seven pages,
  submitted 2018-10-18. The arXiv metadata states existence of a nontrivial
  smooth steady incompressible three-dimensional Euler flow with compact
  support. Review the theorem and construction actually used: whole-space
  domain, axisymmetry, swirl, smoothness, velocity and pressure/Bernoulli
  support, Grad--Shafranov/localization conditions, nontrivial region, and
  whether any streamline period or nondegeneracy is proved.
- P. Constantin, J. La, and V. Vicol, *Remarks on a paper by Gavrilov:
  Grad-Shafranov equations, steady solutions of the three dimensional
  incompressible Euler equations with compactly supported velocities, and
  applications*, author PDF
  [`CLV1.pdf`](https://cims.nyu.edu/~vicol/CLV1.pdf), locally retrieved PDF
  SHA-256
  `1c945acecb8242fd686d546f4d0876522d39122f1666e2d265473304d3af275b`.
  Use it only for the precise reproduced/extended Gavrilov construction and
  its hypotheses; its steady existence theorem does not supply perturbation
  stability or packet transport.
- T. Kato, *Nonstationary flows of viscous and ideal fluids in R3*, *Journal
  of Functional Analysis* 9 (1972), for any classical whole-space local Euler
  existence, uniqueness, differentiability, or perturbation estimate. The
  exact Sobolev index, divergence-free/decay hypotheses, derivative loss,
  amplitude dependence, and interval quantifiers must be checked against the
  proposition actually invoked.

Whole-space Leray projection, the pressure Poisson equation, Kelvin's theorem,
linearized Euler, Lagrangian-flow differentiation, and any finite-time
nonlinear remainder estimate should be derived directly. An additional source
named by `0019` may be inventoried after activation only at its declared role;
it cannot silently add an all-time return or stability theorem.

## Frozen Unit A checks

Unit A's strongest possible result is an actual whole-space Euler
operator/transport theorem on a stated finite interval. It is not an orbital-
stability or persistent-particle theorem.

1. **Actual Gavrilov state.** Verify the displayed velocity and pressure solve
   steady constant-density incompressible Euler pointwise on all of `R^3`,
   including the symmetry axis, localization interfaces, and exterior. Check
   axisymmetry/swirl conventions, smoothness, divergence, vorticity,
   Bernoulli/pressure relation, nontriviality, finite energy, and exactly which
   fields are compactly supported. A formal Grad--Shafranov ansatz or local
   core chart is not yet the global source solution.
2. **Whole-space Leray operator.** Audit the precise Fourier/Newton definition,
   sign and normalization of `P=I-grad Delta^{-1} div`, treatment of the zero
   mode, Sobolev mapping property, symmetry preservation, commutation used by
   the evolution, and the induced nonlocal tail. A compact seed generally does
   not remain compact after Leray projection; any exception needs an exact
   cancellation.
3. **Linearized and nonlinear equations.** Derive the full linearized Euler
   operator about the steady Gavrilov field after Leray projection, retaining
   both `u_G dot grad v` and `v dot grad u_G`, pressure, and domain. For the
   nonlinear packet, substitute the actual divergence-free initial datum into
   full Euler and distinguish exact evolution from the linearized tangent and
   from a prescribed transported profile.
4. **Kelvin/Lagrangian return.** Identify an actual closed material streamline
   or return map of the Gavrilov field, its period and parameter dependence,
   and every nonvanishing/regularity hypothesis. Derive the circulation or
   packet return using the material loop and deformation gradient actually
   generated by the flow. A closed base streamline does not make an arbitrary
   Eulerian vector packet periodic, and a scalar phase return does not establish
   vector/polarization return.
5. **Packet admissibility and locality.** Check that the seed and its Leray
   image are nonzero, finite-energy, in the stated `H^s` tangent space,
   axisymmetric or fully three-dimensional exactly as claimed, and not merely
   translation/rotation/gauge. Track support versus decay and whether any
   weight, streamline tube, or cutoff is fixed independently of future
   evolution.
6. **Finite-time transfer.** Audit the quantifier order among base field,
   packet, amplitude, Sobolev index, and time endpoint. Re-derive the norm and
   order of the difference between the nonlinear perturbed Euler solution and
   the base plus linearized packet, including constants and derivative loss.
   The result must use the same whole-space solution and remain explicitly
   finite-time unless a primary theorem supplies more.
7. **Scope.** Separate exact steady global existence of the base from
   finite-time perturbation transport. Neither a Kelvin return nor local
   continuous dependence establishes all-time recurrence, localization,
   monodromy power bounds, nonlinear orbital stability, or a restoring force.

Unit A receives its own route verdict—`established`, `refuted with the
mechanism named`, or `blocked with the missing construction named`—and its own
verification/review/compatibility/epistemic statuses.

## Frozen Unit B source and dependency inventory

After activation, Unit B may open only `0021`'s frozen README, construction,
source receipt, verdict/validation records, exact verifier, and captured
output, plus every source/API/test artifact explicitly named there at its
declared role.

The controlling primary source is C. García, Z. Hassainia, and T. Hmidi,
*Time-periodic leapfrogging vortex rings in the 3D Euler equations*,
[`arXiv:2603.21644v1`](https://arxiv.org/abs/2603.21644v1), cached PDF SHA-256
`7ba22ef57ba9453dd2631c33b1ad094a35ce2be171426cecce80fba8a5cdfac6`.
Audit the exact finite-core two-patch theorem, contour equation, self/cross
kernels, periodic Sobolev chart, label-delay/reversibility restrictions,
first-mode modulation, parameter/Cantor quantifiers, and mapping estimates
actually transferred. The source's periodic Nash--Moser inverse is existence
machinery, not automatically a full initial-value propagator.

Corrected attempt `0015` at construction hash
`90eec03cd419ca92e6796e6e4d3839b6c6d2152ffa2b04ef49048eacd3a5bff9`
and independent review `0023` may be used only for the actual relative-periodic
source orbit, physical `2*pi*rho_m` action normalization, compact leaf
realization, full formal two-label Jacobi equation, phase/translation modes,
and the already-scoped distinction between source periodic inversion and a
full IVP. They do not supply the operator constructed by `0021` or prejudge
Unit B. Accepted `C-FLO-001` supplies only its conditional finite-dimensional
Jordan-sensitive matrix theorem; it does not license an unbounded contour
operator, its domain, or its powers.

Any additional evolution-family or contour theorem cited by `0021` must be
inventoried and checked after activation at its exact role. It cannot replace
the source-specific two-label construction with an abstract bounded-generator
assertion.

## Frozen Unit B checks

Unit B's strongest possible result is a bounded invertible one-period operator
for the actual full two-label linearized contour IVP, with exact invariant and
symmetry structure on a declared source Sobolev space. It is not a spectral,
power-bounded, nonlinear, or orbital-stability theorem.

1. **Actual two-label contour equation.** Starting from the source patches,
   derive two independent boundary-normal graph equations before imposing the
   existence proof's half-period label relation or reversibility. Retain the
   exact self and cross Green kernels, common translating speed, elliptic base
   charts, area constraints, tangential reparametrization gauge, and all first
   Fourier modes. Check signs and the slow/physical-time conversion.
2. **Linearized IVP operator and domain.** Identify the precise real Sobolev
   spaces, zero-area tangent conditions, regularity index, dense domain, and
   time continuity of the nonautonomous generator. Check its principal
   transport/Hilbert terms, lower-order self/cross blocks, derivative count,
   and parameter uniformity. A periodic boundary-value inverse or a formal
   Jacobi expression does not by itself generate an IVP evolution family.
3. **Existence, boundedness, and inverse.** Verify the theorem or energy
   estimate constructing `U(t,s)` for every allowed source parameter, its
   composition law, uniqueness, boundedness, derivative loss or zero-loss
   claim, and inverse `U(s,t)`. Constants must be uniform only over the exact
   declared epsilon/lambda/time range. Check that the operator advances actual
   arbitrary two-label tangents rather than the reversible delayed subspace.
4. **One-period monodromy.** Define `M=U(mathcal_T,0)` with the exact physical
   period and source phase. Prove boundedness/invertibility on the same space
   and distinguish a map between moving tangent fibers from an endomorphism
   after the periodic identification. Check dependence on contour gauge and
   label ordering.
5. **KKS, invariants, and modes.** With the corrected `2*pi*rho_m`
   normalization, verify symplectic conservation and the exact phase and common
   axial-translation unit modes. Derive linearized energy/impulse invariant
   covectors and every invariant flag or quotient used. Do not delete their
   conjugate/generalized companions merely by calling the eigenvectors neutral;
   area, gauge, phase, translation, energy, and impulse reductions are
   distinct.
6. **Reduced operator.** If a quotient or transverse monodromy is asserted,
   prove that the relevant kernel/intersection is invariant, identify the
   induced norm and domain, and establish boundedness/invertibility of the
   induced map. A quotient by symmetry eigenvectors does not automatically
   remove Jordan shear or make the complement symplectic.
7. **Scope.** One-period boundedness is not power boundedness. Symplecticity,
   invertibility, zero derivative loss, or unit symmetry multipliers do not
   locate the remaining spectrum, establish semisimplicity, control powers,
   close nonlinear remainders, or prove orbital stability.

Unit B receives a route verdict independent of Unit A and its own four-axis
statuses. A failure in one operator realization narrows only that route; it is
not a full carrier or restoring no-go.

## Frozen oracle and sensitivity plan

Primary-source applicability, exact PDE/contour derivation, energy estimates,
operator-domain analysis, and direct invariant identities are the strongest
oracles. Existing exact scripts/tests and captured receipts may corroborate the
defined predicates. Inspect their source and sensitivity once; reuse an
unchanged passing receipt unless a scientific predicate is invalidated.

Representative Unit A exposing mutations are: reverse the whole-space Leray
sign; omit one linearized advection term; assert compact support after a generic
projection; replace a material Kelvin loop by a fixed Eulerian curve; or remove
the fixed finite-time/amplitude/Sobolev quantifier from the nonlinear estimate.

Representative Unit B mutations are: drop a self/cross contour block; impose
the source half-period/reversibility restriction on the claimed full domain;
project away a first mode; reverse the translating or corrected KKS sign;
claim invertibility without the reverse evolution; or promote boundedness of
one `M` to bounded powers. Each selected mutation must be exposed by a direct
derivation, source applicability check, or genuinely sensitive assertion.

No production small force, energy splitting, soft Hessian eigenvalue, or
stability-window edge is registered. The small-ratio numerical prescriptions
do not bind at freeze. Any numerical spectrum or near-unit multiplier found in
a target artifact is outside this review until separately preregistered.

## Completion and non-inheritance boundary

After central activation, complete `0026` with one review artifact and one
verdict artifact containing two separate decisions, the strongest supported
statement for each unit, any minimum evidence-driven correction, boundary
hashes, oracle roles, and each exact remaining dependency. Send
artifact/verdict/dependency to `herdr agent prompt w3:p1`.

Neither unit licenses the other's object or operator. Together they do not
license all-time packet recurrence, a power-bounded or spectrally stable GHH
monodromy, nonlinear/orbital stability, a generic three-dimensional regularity
theorem, global contour propagation beyond the exact Unit B result, a
mechanical force, particle identity, physical spin, quantum statistics,
electron/neutrino completion, parent-campaign completion, or a global no-go.
No source, central, API, test, memory, or commit edits are owned by `0026`.

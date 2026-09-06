# Fresh predecessor claim/oracle audit: P251 original route and C-CST-008--016

Reviewer: `/root/predecessor_audit`, fresh non-author/non-implementer

Boundary: `substrate-framework@b6fc902a0942d07996f12a81028fbd3f7c909a43`, release `v0.183.0`

Audit obligation: P252/A3

Date: 2026-09-06

## Verdict

The accepted predecessor claims C-CST-008--016 have meaningful positive
content at their written scopes. I found no equation-level counterexample to
those scoped statements. Their strongest common meaning is a collection of
exact static ensemble identities together with analytic **quadratic,
linear-response and deliberately prepared** Euler constructions. They do not derive an autonomous finite-amplitude
nonlinear Cosserat material law from unrestricted Euler, and the source files
usually say so precisely.

The original C-CST-001--007 route remains unresolved and must not be read as a
promoted predecessor theorem. Its frame kinematics and target-system algebra
survive, but its advertised straight-tube self-energy-to-locking chain does
not. The durable campaign record itself calls the original route unsupported
and records `alpha=L_v*T/6` as unsupported
(`attempts/0271/closure-map.md`, lines 8--18). The active historical graph keeps
N2--N7 open; later C-CST-008--016 use different compact-response, prepared
history, or column-mode mechanisms rather than proving those original nodes.

The main missing bridge to Federico Comparsi's broader nonlinear objective is
concrete, not a request for more generic caution: C-CST-009's finite phase
Cauchy--Born action has the reconstruction residual

```text
R = Udot - V + E zdot - A E z,
```

which remains explicitly nonzero/uncontrolled for unrestricted Euler
trajectories (`attempts/0102/canonical-coefficient-join.md`, lines 320--337).
The later prepared-history claims replace that missing invariance by
input-dependent controls and ordered limits. In particular C-CST-014 chooses a
finite preparation curvature `B_*`, C-CST-015 chooses a preparation parameter
to make energy equal response, and C-CST-016 changes the microscopic control
with the accuracy and then chooses the macro wave number after its norm is
known. Those are valid existence/preparation results because the statements
declare them. They are not an Euler-derived intrinsic nonlinear constitutive
coefficient or an autonomous nonlinear invariant manifold.

## Scope and coverage honesty

I read the registry statements and cited claim-level source/proof bodies for
each C-CST-008--016, the original frozen N1--N7 graph and its append-only repair
ledger, the archived Gaussian-Beltrami source applicability audit, and the
critical exact/symbolic verifier surfaces. I used prior independent reviews
only after forming source-level conclusions and only as navigation/correction
history, never as proof. I did not read all 271 attempt directories line by
line. This is complete claim-boundary coverage with deeper inspection of the
proofs that supply coefficients, kinetic forms, physical observations and
history claims, especially C-CST-009/015/016 consumed by C-CST-017 and then
C-CST-018.

## Original C-CST-001--007 route

| Historical claim/node | Strongest source-supported result | Oracle and evidence role | Present route verdict / exact missing step |
|---|---|---|---|
| C-CST-001 / N1 | Frenet/Bishop/material-frame transport identities, including the gauge split of axial frame rate, are exact kinematics. | `verify_cst001.py` performs symbolic substitutions and sign mutations; attempt 0003 is an exact algebra regression. | **Established for filament-frame kinematics only.** It does not identify the resulting angle with a globally defined Euler microrotation, symmetry, charge or nonlinear material director. |
| C-CST-002 / N2 | The corrected Rankine-column branch equations contain the m=1 logarithmic branch with constant `1/4-EulerGamma` and the m=2 coefficient `-1/6`; the historical ledger records the corrected Bessel derivative and signed-root corroboration. | Exact asymptotic algebra is primary; 40/60-digit Bessel roots are corroboration. Earlier numerical eigen-solves are explicitly void/representation-scoped. | **Blocked as the original node.** The promised straight-tube bending and core-twist moduli from regularized Biot--Savart self-energy were never closed. Branch dispersion is not the promised static modulus. |
| C-CST-003 / N3 | Sphere moments and coefficient matching for a supplied micropolar energy are exact. | `verify_cst003.py` checks tensor identities and target-energy matching. | **Refuted for the original free-director line-stretch locking mechanism.** Exact relative Green--Lagrange line strain is independent of free frame rotation, so it cannot yield the claimed `alpha`; the ledger states this at `historical-route-state.yaml`, lines 496--498. The later compact finite-core response is a different route. |
| C-CST-004 / N4 | Varying the stated quadratic micropolar action gives the displayed linear and angular momentum balances; the repaired longitudinal spin coefficient is `2(c_s+c_tr)`. | SymPy Euler--Lagrange/operator algebra with sign and dropped-term mutations. | **Established only as target-system algebra conditional on supplied positive kinetic and potential coefficients.** It does not derive those coefficients from the original Euler tube ensemble. |
| C-CST-005 / N5 | The two-branch dispersion and reduced linear-system energy conservation follow from the C-CST-004 matrices. | Exact determinant algebra is primary; DOP853 is only reduced-equation regression. | **Established only for the conditional target operator.** It is not an Euler spectrum or evidence for an invariant nonlinear Euler sector. |
| C-CST-006 / N6 | Zero coherent first moment follows for independent orientation averaging, but quadratic fluctuations can remain. | Original Monte Carlo is finite-sample corroboration. C-CST-010's exact product-Haar change of variables is the repaired theorem. | **Original inference repaired, not established as originally broad.** Uniform one-point marginals alone are insufficient; correlated uniform phases are a counterexample. A no-angular-current statement needs its separate dynamical closure. |
| C-CST-007 / N7 | EPS sources establish appropriate steady Beltrami fields/invariant knotted tubes; the elementary periodic Beltrami field correctly checks `curl u=lambda u =>` stationary Euler. | Primary-source import plus exact equation regression and file provenance. | **Established as an existence/applicability import only.** It supplies neither the original Biot--Savart moduli nor the coarse-grained action/history bridge. |

The historical defects are real and repaired in the later route rather than
silently erased. The most important examples are: the original interior
Poincare wavenumber carried `1/omega` instead of `1/wtilde^2`; the m=1 constant
was later corrected; a rate-quadratic coefficient was found to have inertia
units and could not serve as angle stiffness; the original mutual-energy
kernel used the wrong pairing; and zero signed response was initially treated
too strongly. These are recorded around
`attempts/0271/historical-route-state.yaml`, lines 330--514. None is a present
counterexample to C-CST-008--016, whose statements changed mechanisms and keep
their scopes explicit.

## Accepted predecessor inventory

### C-CST-008 -- compact stationary Euler angular cells

**Claim object.** A positive-density stationary isotropic Gaussian Beltrami
marked law with actual invariant tubes; compact isovortical fields `Q,S`; a
positive complete two-by-two orbit Hessian; nonzero KKS and independently
normalized generator/spin moments; and arbitrarily stiff additional compact
pairs.

**Source/proof chain.** The external existence input is the
Enciso--Peralta-Salas--Romaniega Gaussian Beltrami support and knot-density
result as audited in `attempts/0060/source-review.md`. The isotropic, ergodic
Gaussian law is constructed from twelve independent rotated copies in
`attempts/0071/README.md`; this correctly averages covariances before invoking
Gaussian determination and does not confuse a Haar mixture with a Gaussian
law. `attempts/0085/compact-spin-theorem.md` constructs the exact differential
syzygy and characteristic-zero rank certificate. `0085/positive-compact-pair.md`
uses the full Hessian and KKS forms (lines 17--18), retains the subprincipal
return, and proves a finite positive carrier choice. `0103/six-moment-construction.md`
separates generator moment `G` from mechanical-spin moment `L`.
`0098/stationary-compact-assembly.md` supplies the measurable disjoint marked
selection and Palm averages.

**Oracle/evidence roles.** The proof-bearing oracle is analytic construction
plus exact rational/mod-101 rank, with a second exact integer-rank replay.
Symbolic scripts expose row conflation and rank loss. The cited Gaussian and
EPS papers are applicability/import evidence. There is no load-bearing small
eigenvalue or floating-point sign.

**Verdict.** **Supported as written.** The statement is about tangent fields,
the complete quadratic orbit Hessian and a local registered core-angle
observation. It correctly denies a global compact action or unrestricted
dynamics. That limitation is essential for downstream use.

### C-CST-009 -- conditional one-action Euler-to-Cosserat continuum

**Claim object.** Pullback of one quadratic Euler material phase action to an
explicit phase Cauchy--Born family; common macro momentum and separate
microscopic reactions; exact Schur coefficients `j`, `kappa`, `alpha=kappa/4`;
positive selected shear/spin-gradient coefficients; physical observation maps;
and linear micropolar balance/dispersion identities.

**Source/proof chain.** `0097/canonical-phase-pullback.md` derives the canonical
one-form, Hamiltonian, signed KKS pullback and all macro/internal crosses.
`0097/material-cotangent-bridge.md` explicitly distinguishes phase velocity
from material reconstruction. `0102/canonical-coefficient-join.md` performs
the six-moment triangular completion, common-mean reduction, time-reversal
reaction split and kinetic/potential normalization. The construction enforces
`C=J` only after `J=D P^-1 D*` has been computed (lines 160--169); it does not
insert an external rotor inertia. Positive `mu,C_T,C_L` are obtained by finite
compact attachments whose complete action integrals and inertia costs remain
in the formulas (lines 233--279).

**Euler-derived versus assigned.** Within the declared restricted action, the
coefficients are Euler integrals/Schur complements of explicit fields. They are
not copied target numbers. The geometry, carrier and attachment amplitudes are
freely selected construction inputs, so the result is an existence theorem,
not a unique prediction from `rho,lambda` alone. The older
`alpha=L_v*T/6` is not used.

**Oracle/evidence roles.** The analytic phase/action derivation is primary.
`0102/verify.py` checks finite-block algebra largely by comparing constructed
expressions with the displayed expected normal forms; it is exposing regression,
not independent proof of the Euler reduction or PDE realization. The DOP853
run proves only the already-derived reduced ODE behavior.

**Verdict.** **Supported as the explicitly conditional linear second-jet
Cauchy--Born action.** It supplies C-CST-017/018 only with action/current and
branch-substitution algebra. It does not supply actual Euler invariance or
nonlinear histories. The exact missing bridge is the retained reconstruction
residual quoted above.

### C-CST-010 -- product-Haar conservative contrast

**Claim object.** For a complete unresolved local product-Haar orientation law,
independent coherent right shifts leave the averaged interaction energy
constant, so conservative coherent torque and stiffness vanish; fluctuations
need not vanish, and angular current is separate.

**Source/proof chain and oracle.** `0058/haar-contrast.md`, lines 18--31, is a
complete exact change-of-variables proof and provides the correlated uniform-
marginal counterexample. The 200000-sample Monte Carlo is narrower first-moment
corroboration only.

**Verdict.** **Supported as written.** The product law and complete-state
qualification are load-bearing. It cannot be weakened to one-point isotropy,
and it cannot be strengthened to a dynamical no-current theorem without the
separate closure stated at line 44.

### C-CST-011 -- finite optical packet on an actual stationary tube

**Claim object.** A finite-carrier, finite-action, Kelvin-prepared optical packet
with actual material tag angle and literal tagged angular momentum, controlled
carrier curvature, positive two-phase action and periodic stationary geometry,
all on fixed finite optical windows.

**Source/proof chain.** `0147/finite-packet.md` derives the Laguerre packet,
full Plancherel KKS, physical sheet tag, finite moment/IFT construction and
moving action. It explicitly evolves the preparation by the full Euler/Lin
system (line 52) and chooses nonlinear disturbance amplitude only after the
fixed geometry and time window (line 265). `0147/toroidal-transfer.md` builds
the common-circle finite CK background, proves a local twist margin, retains
the full pressure through weighted angular-harmonic estimates, transfers the
packet to the torus and then to a periodic Beltrami approximation.

**Oracle/evidence roles.** Analytic asymptotics, finite-dimensional IFT and the
cited KAM/Beltrami approximation theorems are primary. `0147/verify.py` checks
Laguerre identities, finite row ranks and several scaling-loss calculations;
it does not independently prove the weighted pressure estimate, KAM transfer
or nonlinear PDE existence.

**Verdict.** **Supported at the fixed-window, prepared linear-response/small-
amplitude scope stated.** It does not establish an autonomous optical mode, a
globally conserved tag spin or a nonlinear coupled continuum. Those exclusions
prevent its finite-action construction from being misused as such.

### C-CST-012 -- triangular-array acoustic window

**Claim object.** An exact triangular periodic stationary Euler array, a
uniformly bounded zero-axial horizontal Euler group on the declared vector
isotypic sector, and a slow-time acoustic estimate for common-velocity and
Kelvin-displacement phases with the full material action.

**Source/proof chain.** `0161/planar-acoustic-array.md` derives the field,
separatrix tiling, weak transport primitives, first-shell kernel and exact
Fourier complement gap. Its adjoint rows control the translation kernel and
yield the uniform group bound (line 105). The normal-form proof then gives the
`0 <= t <= T/|k|` estimate (line 144), while the initial action calculation
derives the positive phase energy (line 215).

**Oracle/evidence roles.** The functional-analytic proof is primary. The
symbolic verifier checks field/pressure identities, the lattice gap, boundary
term, covariance and Kelvin initial action, but not the long-time semigroup
estimate itself.

**Verdict.** **Supported on the stated symmetry sector and acoustic window.**
The nonlinear amplitude can depend on `k`; the source explicitly denies a
uniform nonlinear amplitude and infinite-time stability. It gives no nonaxial,
lifted-background or optical joining theorem.

### C-CST-013 -- prepared mean/action transfer

**Claim object.** Exact first- and second-directional mean jets for a general
stationary periodic Euler field, plus transfer of the fixed-time observed
coefficients from a one-wave reference to a selected finite constant-curl
insertion carrying the C-CST-011 optical object.

**Source/proof chain.** `0170/prepared-mean-jet.md` differentiates the complete
Jacobi equation, retains both pressure rows, eliminates the unknown second cell
by divergence and integration by parts, and distinguishes material mean from
physical velocity mean. It proves period-uniform first-cell estimates using
normalized `H2` defect and a fixed global `C3` bound, then derives positive
`K_N/(rho k^2)` only after selecting a finite cell and then small `k`
(line 223).

**Oracle/evidence roles.** Analytic PDE estimates are primary. Exact Fourier
tests expose the projector, pressure-sign and current-factor errors. They are
corroboration rather than a numerical convergence proof.

**Verdict.** **Supported as a fixed-time prepared mean/action transfer.** The
source itself denies full Bloch-operator `C2` convergence, a period-uniform
analytic radius and acoustic-time homogenization (lines 238 onward).

### C-CST-014 -- correlated smooth-column optical response

**Claim object.** Two genuine smooth compact-vorticity column-mode families,
positive physical tilt actions and current rows, combined by a positive mixture
and signed state preparation to cancel frequency variance and realize positive
transverse/longitudinal second-jet curvatures.

**Source/proof chain.** `0181/nodal-domain-mode.md` gives the Rankine determinant,
monotonic simple nodal branch, positive KKS/action, a positive fixed material
tag, carrier-jet current identities and smooth-vorticity continuation.
`0187/nodal-family-realization.md` obtains common frequency by Euler amplitude
scaling and unequal slopes by geometric scaling. `0187/correlated-family-action.md`
derives the variance and inherited-energy conditions; energy equality holds
only when its separate scalar condition `D=0` is imposed (line 72). The two
preparation slopes solve that condition together with `B=B_*` (line 92).
`0191/transverse-family-correction.md` repairs the earlier wrong Haar tensor and
constructs a pointwise-positive curvature representative.

**Euler-derived versus assigned.** Mode frequencies, slopes, phase masses and
currents are actual Euler quantities. The positive `B_*` is an explicit freely
chosen preparation coefficient, and the source calls the resulting curvature
preparation-dependent (`correlated-family-action.md`, line 101). Therefore this
is not an intrinsic modulus prediction.

**Oracle/evidence roles.** Exact Sturm/Bessel calculus, KKS sign and IFT are
primary. Symbolic correction checks test the finite mixture algebra and Haar
factor, not Euler mode existence by themselves.

**Verdict.** **Supported as a prepared positive-mixture column response per
axial period.** It gives neither a single same-field stationary ensemble nor a
finite-bandwidth three-dimensional mode or nonlinear material law.

### C-CST-015 -- positive displacement preparation on a fixed Euler cell

**Claim object.** One explicit periodic Beltrami field and a correlated finite
initial displacement preparation with positive exact second-order inherited
energy and matching physical restoring response.

**Source/proof chain.** `0188/positive-correlated-preparation.md` solves the two
first-shell moment equations by exact rational Fourier algebra, constructs the
stationary forced corrector and full first-shell Euler-kernel return, and
computes energy `h2` and mean response `R_D` as different contractions. It then
chooses the explicit preparation parameter `t` so that `h2=a=-R_D>0`
(lines 95--103). `0188/initial-phase-energy.md` keeps the actual initial
circulation/cotangent distinction.

**Euler-derived versus assigned.** The positive value `a` is algebraically
derived after the finite cell/corrector is computed, and is not an empirical
fit. The choice of `t` is nevertheless target-directed state preparation. It
therefore establishes existence of a prepared restoring column, not spontaneous
or universal Euler elasticity.

**Oracle/evidence roles.** The finite Fourier calculation in
`verify_chebyshev.py` is a strong exact oracle: it constructs the field,
projectors, generator, kernel return, energy and current separately before
checking their equality. Some review scripts merely restate final identities,
but they are not needed for the load-bearing result.

**Verdict.** **Supported as written.** Its own final paragraph correctly leaves
the common-velocity history, optical/current join and autonomous acoustic
Hamiltonian open.

### C-CST-016 -- common acoustic action by passive controls and a diagonal limit

**Claim object.** On the C-CST-015 field, finite smooth velocity/configuration
controls in an exact passive Euler sector approximate the common-velocity
physical output, cancel the added second-order energy and phase moment, and
produce a positive fixed-time second-order acoustic action along an ordered
`accuracy -> k` diagonal.

**Source/proof chain.** `0196/actual-velocity-output.md` derives the complete
three-component source and physical current. `0196/passive-physical-control.md`
uses actual low-frequency wrapped streamlines and a signed Vandermonde
combination of smooth positive band densities. `0196/configuration-phase-energy.md`
sets `h=T^-1 g/2` on nonzero harmonics, cancels the full inherited energy, and
adds one negative-frequency moment to restore the phase. `0196/fixed-cell-diagonal-limit.md`
records the control norm `N_n`, obtains a finite cubic bound
`C_n=O(N_n^2)`, constructs the controls at `delta_n=2^-n`, and only then chooses
`k_n` (lines 109--130).

**Oracle/evidence roles.** The transport-sector identity, coarea formulas,
Vandermonde sign, functional approximation and fixed-cell Euler/Lin estimates
are analytic. `0207/verify_review.py` checks a small representative moment
system and exact transport identities; it does not prove the infinite sequence
or a uniform control norm. The proof does not need such a uniform norm because
the diagonal explicitly makes `k_n` depend on `N_n,C_n`.

**Verdict.** **Supported as a diagonal prepared-history theorem.** The source
expressly notes that the full radial family need not be `C2` at preparation
switches and does not yield an acoustic-time window. Consequently C-CST-016
cannot by itself license a fixed microscopic preparation with a differentiable
Bloch second derivative, an autonomous constitutive law, or finite-amplitude
nonlinear closure.

## Kinetic metric, phase, current and history audit

The source record now distinguishes the objects that earlier attempts
conflated:

* C-CST-008 uses the complete isovortical Hessian and KKS, not a positive
  kinetic norm alone.
* C-CST-009 pulls both the one-form and Hamiltonian from one material phase
  action, averages the common macro momentum before eliminating independent
  microscopic reactions, and keeps generator moment `G` distinct from induced-
  velocity spin `L`.
* C-CST-011/014 normalize phase mass from the actual KKS and measured angle;
  their physical spin/current rows are additional observations rather than
  definitions of the canonical momentum.
* C-CST-015 changes circulation data explicitly and C-CST-016 repairs both the
  full initial energy and the remaining phase moment. The raw packet current is
  retained even when the reduced physical-chart connection becomes small.

This is substantial repair of the historical assigned-inertia pattern. It also
sets the exact ceiling: these are quadratic phase restrictions or prepared
linear histories. None proves that the same finite coordinate family is closed
by the nonlinear Euler vector field.

## Symbolic and numerical oracle audit

No accepted C-CST-008--016 conclusion depends on a soft Hessian eigenvalue,
small force obtained from a difference of large energies, sampled stability
window, or a floating-point energy splitting. The small-ratio prescriptions
therefore do not reveal a present precision defect in these claims:

* C-CST-008's positivity uses strict analytic high-carrier inequalities and
  exact rank certificates.
* C-CST-009's DOP853 run is labeled reduced-equation regression.
* C-CST-010's Monte Carlo is labeled narrower corroboration.
* C-CST-011--016 use exact/analytic asymptotics, finite-dimensional existence,
  functional approximation and ordered bounds; their scripts exercise finite
  identities rather than making a numerical sign call.
* The original C-CST-002 Bessel-root calculations at 40/60 digits corroborate
  branch asymptotics only. The historical numerical eigen-solve was rejected
  as representation-polluted and does not support an accepted claim.

Several symbolic scripts compare a hand-constructed expression to the same
displayed expected normal form. Those checks are useful mutations/regressions,
but they would be false-green if cited as independent proof of Euler existence,
semigroup bounds, KAM transfer or action reduction. The registry assumptions
for C-CST-009 and C-CST-011--016 correctly name the analytic construction/PDE
proof as primary and the symbolic checks as exposing corroboration. I found no
accepted claim whose only load-bearing evidence is one of those copied-normal-
form assertions.

## Dependency effect on C-CST-017/018

C-CST-012--014 do not appear in the accepted dependency lists of C-CST-017 or
C-CST-018. They are parallel predecessor routes and cannot be silently combined
into terminal closure.

C-CST-017 consumes C-CST-009, C-CST-015 and C-CST-016. That consumption is
sound only at these exact interfaces:

1. C-CST-009: canonical action/current and branch-substitution algebra for the
   declared conditional Cauchy--Born restriction;
2. C-CST-015: the explicit fixed-cell positive displacement preparation and
   its second-order phase/energy/current coefficient;
3. C-CST-016: common-input fixed-time prepared Euler/Lin histories along the
   norm-aware diagonal.

C-CST-018 in turn lists C-CST-009 and C-CST-017 and explicitly says those
dependencies supply action/current and branch-substitution algebra only. No
predecessor reviewed here licenses unrestricted Euler invariance,
acoustic-time-uniform nonlinear dynamics, an intrinsic preparation-independent
modulus, or a finite-amplitude nonlinear rotational continuum. A terminal
claim that stays within prepared quadratic second variation can legitimately
reuse these results; a statement of Federico's broader nonlinear objective
cannot.

## Minimum continuation for the broader nonlinear objective

The shortest scientifically meaningful next route is to replace preparation-
by-preparation output synthesis with one fixed nonlinear reduction. Freeze a
specific Euler stationary ensemble and a globally defined finite-rotation
orientation observable, derive the full nonlinear pulled-back action and
kinetic metric, and then prove one of:

1. an exact invariant manifold under Euler; or
2. a normally controlled approximate slow manifold with the reconstruction
   residual and all eliminated complement modes bounded on the target
   macroscopic time scale.

The proof must derive the finite-amplitude stress/couple-stress and its
coefficients from that fixed ensemble before selecting histories to match a
desired law. It must preserve the physical/canonical current distinction and
show the Legendre map remains positive/nondegenerate beyond the tangent point.
If only a well-prepared homogenization theorem is intended, the minimum repair
is to state that objective explicitly and provide uniform estimates for one
fixed preparation family rather than the C-CST-016 switching diagonal.

## Route verdict

`route_verdict: established` for the bounded predecessor audit: the support
boundary and the concrete missing nonlinear bridge are identified claim by
claim.

`evidence_scope: SOURCE_AND_ORACLE_AUDIT` -- no new claim promotion, no
counterexample to the accepted scoped statements, and no assertion that all
271 attempts were read line by line.

## Reconciliation with `target-and-routes.md`

This is a bounded source-level check of `target-and-routes.md` and
`root_nonlinear_separation.py`; it is not a new proof or promotion gate.

The finite-angle nonuniqueness example is exact.  For a relative rotation
`Q=R_1^T R_2`, the objective scalar

    s(Q)=(3-tr Q)/2=1-cos(theta)

realizes the two displayed potentials as `k s` and `k s+beta s^2`.  Taking
`k>0` and `beta!=0` makes the intended counterexample explicit: the Hessians
at the identity agree and the first difference is cubic in the potential
gradient.  One wording correction is needed at
`target-and-routes.md:70-71` and `root_nonlinear_separation.py:9-12`:
the checked positive expression is
`partial_theta(V2-V1)=2 beta(1-cos theta)sin theta`.  If "torque" means the
usual mechanical generalized force `-partial_theta V`, its difference has
the opposite sign.  Calling the checked object the potential-gradient or
declaring the torque convention removes the ambiguity; the nonuniqueness
conclusion is unchanged.

The affine Euler covariance example is also exact: divergence and steady
Euler residual vanish, the flow-map determinant is one, and the covariance
eigenvalues vary.  Its scope needs to remain attached to the result.  This
velocity and pressure are affine/quadratic on all of `R^3`, hence unbounded,
nondecaying, nonperiodic and not finite-energy.  It refutes rigid-orientation
closure as a consequence of the local Euler equations alone; it is not an
admissible counterexample inside a compact, periodic or finite-energy class
without a separate localization construction.  Lines 88-91 already exclude
P251's prepared compact background, so this is a discussion qualification,
not a defect in any accepted C-CST claim.

The target distinction is supported by the cited primary text.  Böhmer,
Downes and Vassiliev explicitly describe a special Cosserat theory with
geometrically nonlinear rotation kinematics and physically linear quadratic
energy.  Erofeev, Pavlov and Vikulin section 5 explicitly describes Vikulin's
alternative as rigid blocks with intrinsic moments, a surrounding elastic
stress field, symmetric stress, and the one-dimensional sine-Gordon equation
(3).  These facts distinguish nonlinear rotational kinematics from nonlinear
constitutive response, and distinguish the Vikulin block/host model from a
generic nonlinear Cosserat target.  Neither source supplies an Euler
derivation; `target-and-routes.md` says so correctly.

The three proposed routes do not silently assume the missing closure:

* Route A retains the hierarchy/memory and asks for decay, a scale limit or an
  invariant/slaved family before localization.  Its commuting filter is an
  explicit hypothesis; a moving-cell realization must retain the stated
  boundary and pressure fluxes.
* Route B expressly says restricted action is insufficient and requires a
  normal-residual estimate.  The claimed variation formula is exact with the
  advecting velocity held fixed, or in a material-coordinate variation; a
  simultaneous Eulerian variation of the advecting field must include its
  convective contribution.  The proposed family must also enforce
  `det F=1`, compatibility and its boundary conditions rather than assigning
  them after restriction.
* Route C requires the finite-angle interaction, inertia and full angular
  current to be derived from Euler and rejects matching only the linearized
  stiffness.  A symmetric effective stress is therefore an output to derive
  with its spin/couple-current split, not a licensed import from Vikulin.

No further source-level correction to these two P252 artifacts is indicated
by this bounded check.

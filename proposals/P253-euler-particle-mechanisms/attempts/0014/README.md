# P253/0014 preregistration: dynamically accessible swirl stability and restoring operator

Status: **README-only preregistration awaiting central registration and schema activation.** No
new external source body, calculation, numerical design, or comparator value may be opened in
0014 before the central activation receipt is zero. This attempt owns only its eventual `0014`
artifacts; source modules, tests, governance, the proposal manifest, and other attempts remain
outside its write surface.

## Parent obligation and positive deliverable

This is the failure-derived P2/P4 continuation from 0005/0010. Attempt 0010 supplies two actual
localized intrinsic-swirl candidates and the exact joint Cao--Zhan Hessian, but it does not yet
supply LP2 persistence. The positive deliverable is an exact, carrier-specific stability or
restoring statement on the **dynamically accessible Euler leaf**, not on an arbitrary pair of
Eulerian functions `(eta,chi)`:

1. derive the physical coadjoint-orbit tangent at one Cao--Zhan/Turkington-generator ring and at
   the compact Gavrilov carrier from `delta m=-ad*_v m`, including its boundary/decay class and
   all circulation constraints;
2. pull the complete augmented Euler Hamiltonian to that tangent, remove translation and genuine
   stabilizer zero modes, and determine the exact constrained second-variation signature;
3. prove a coercive modulated deficit and nonlinear remainder estimate in the declared
   perturbation class, or identify the exact neutral/negative mode and then test the corresponding
   linearized Hamiltonian mechanism; and
4. state the physical restoring observable licensed by the result. A Hessian eigenvalue is not
   silently renamed a temporal oscillation frequency: a frequency requires the restricted KKS
   operator and the linearized Hamiltonian flow.

Success licenses an LP2 axisymmetric restoring/persistence statement for the same finite-`j`
carrier and its controlled KKS sector. It does not license full three-dimensional stability,
scale selection, exchange statistics, quantization, relativity, or electron/neutrino identity.

## Frozen object, ensemble, and exact primary target

Work first in decaying axisymmetric `R^3`, where 0005 representation A is complete and no
unproved arbitrary-`H^1` Clebsch completion is needed. Let

    T_DA(m_*) = closure{ -ad*_v m_* : div v=0,
                         v axisymmetric, decaying, regular },

with closure taken in the perturbation energy space. Only after deriving its image may a tangent
be written as `(eta,chi)=(delta zeta,delta xi)`. Preserving the linearized Casimirs
`int C(xi)` and `int zeta F(xi)` is necessary but is not assumed sufficient for membership in
`T_DA`. In particular, a freely chosen pure `(0,chi)` direction cannot be used to infer a saddle
or instability unless it is produced by an admissible displacement.

For the Cao--Zhan generator `H(psi)=psi_+`, `-B'(psi)=alpha` at fixed thickness `epsilon`, the
already-derived conserved augmented functional has Hessian

    Q_e[eta,chi]
      = int [eta K eta-epsilon^2 r^2 eta^2] dnu
        +int (chi/r-epsilon r eta)^2 dnu.                 (1)

Freeze the perturbation norm

    ||(eta,chi)||_X^2=int [eta K eta+chi^2/r^2]dnu,       (2)

and let `Z` contain the exact axial-translation tangent and any further zero modes proved from
the genuine stabilizer. The one primary exact second-variation/spectral target is

    gamma_e = max_{s in {+1,-1}}
                inf { s Q_e[q] : q in T_DA(m_e) intersect Z^perp,
                                  ||q||_X=1 }.            (3)

Equation (3) is fixed before any spectrum is inspected. `gamma_e>0` means that one globally
consistent sign of the Hessian has a nonzero constrained gap; the sign `s` is reported, not
selected profile-by-profile. `gamma_e=0` identifies an exact/essential neutral obstruction.
`gamma_e<0` means the accessible Hessian has mixed signature. Mixed signature alone is not a
claim of dynamical instability: that route must next analyze the actual linearized Hamiltonian
operator `J_orbit Hess(A_e)` and exhibit a growing mode or a rigorously neutral mechanism.

The expected maximum branch suggested by the reduced Cao--Zhan variational principle is tested,
not assumed. A positive (3) becomes a physical restoring theorem only after conservation of the
augmented functional, modulation fixing `Z`, a controlled higher-order remainder, and exclusion
of support escape are proved in the same norm. The thin limit `epsilon->0` is a separate uniformity
question; fixed-`epsilon` coercivity does not imply an epsilon-uniform particle scale.

## Competing registered routes

### Route A — Cao--Zhan joint constrained maximum

Use the exact energy--impulse--angular-momentum--helicity functional from 0010. Derive the orbit
tangent directly from representation A, reconcile it with the source rearrangement class, and
evaluate (3) analytically. Attempt compactness/nondegeneracy and a modulation estimate around the
axial translate family. This route has the clearest reduced functional but must prove that its
allowed joint swirl directions are the physical leaf directions.

### Route B — Gavrilov compact pressure-shell orbit

Use the smooth compact carrier and its real meridional-plus-swirl field. Derive the second
variation directly from the Euler--Poincare/coadjoint-orbit action so the flat cutoff boundary is
handled without dividing by a vanishing profile derivative. Determine the analogue of (3), with
translations and rotations modulated. This route has exact compact support and finite nonzero
`j`, but no published variational extremum and potentially has pressure-shell modulation zero
modes across different conserved-label leaves.

### Route C — Hamiltonian spectral/representation continuation

If neither energy--Casimir route yields a sign-definite gap, retain the same carrier and derive
the exact linearized orbit generator rather than declaring instability from an unconstrained
Hessian. Seek an analytic growing eigenmode, continuous-spectrum obstruction, or bounded neutral
dynamics. A materially different compact helical carrier may be registered append-only if the
failure mechanism is carrier-specific; absence of a gap for A or B is not a global Euler no-go.

Routes are compared by: exact membership in the physical Euler orbit; finite energy and finite
nonzero `j`; correct zero-mode quotient; regularity of the variational functional at support
boundaries; analytic control before numerics; perturbation breadth; and whether the result yields
a conserved modulated distance or an actual Hamiltonian time rate without imported restoring
potentials.

## Frozen source-access inventory

No new body is opened by this preregistration. After activation, access status is recorded before
any additional download; full papers remain in `/tmp/primary-source-cache/P253-0014`.

| Source/input | Existing access and intended role | Import boundary |
| --- | --- | --- |
| Cao--Zhan, arXiv:2009.13210v2 | Full primary body already audited with hash/theorem locations in `0010/source-audit.md` | Exact swirling-ring family and reduced maximization. The paper explicitly does not supply stability/nonstability or local uniqueness. |
| Turkington, DOI 10.1137/0520005 | SIAM metadata/abstract only in 0010; full body was inaccessible | No unseen theorem is imported. Its generator is used only through Cao--Zhan's explicit accessible equations. |
| Gavrilov, arXiv:1810.08020v1 | Full primary body already audited in 0010 | Smooth compact steady carrier and pressure cutoff; no stability theorem. |
| Constantin--La--Vicol, arXiv:1903.11699v1 | Full primary body already audited in 0010 | Independent Grad--Shafranov reconstruction and regularity boundary; no stability theorem. |
| P253/0005 plus independent P253/0011 review | Local repository derivation/review, not accepted canon | Exact full Euler orbit action and corrected general-`H^1` row. Local Clebsch is corroboration only in its proved coverage. |
| Choi, arXiv:2011.06808v2; Cao--Lai--Qin--Zhan--Zou, arXiv:2206.10165v2 | Primary bodies and exact no-swirl stability scopes already audited in `0002/source-transfer.md` | Method comparators for modulation, compactness, and rearrangement. Their no-swirl perturbation classes do not transfer automatically to intrinsic swirl. |

Any additional dynamically-accessible-variation or Hamiltonian-stability source is inventoried
after activation with URL, version, lawful access result, hash, and exact theorem/equation
location before use.

## Analytic ladder and oracle boundary

The attempt earns its strongest result by calculus first:

1. derive the `(eta,chi)` image of `-ad*_v m_*` and prove constraint closure;
2. rederive the first and second variations before imposing the equilibrium relations;
3. identify translation/rotation/stabilizer kernels and choose modulation conditions;
4. determine compactness, essential spectrum, boundary degeneracy, scaling, and fixed-thickness
   limits of the constrained operator;
5. decide (3) analytically where possible and derive the nonlinear remainder needed for a
   modulated conserved deficit; and
6. only if a named spectral remainder survives, freeze a numerical representation in a later
   activated artifact.

No production numerical verifier is preregistered here. If a soft eigenvalue, spectral edge, or
small restoring rate remains, the small-ratio-numerics skill must be read before freezing mesh,
domain, zero-mode gauge, precision, tolerances, or acceptance thresholds. Required evidence would
then include the continuum operator, residuals, `lambda_min/lambda_2`, crossed mesh/domain or
basis refinement, observed-order/error budget, and jitter/sign sensitivity. A boxed eigenvalue
cannot establish open-space coercivity.

Maximum verdict before the nonlinear remainder closes is an exact constrained-Hessian signature
for a named carrier and perturbation class. Full LP2 stability requires the conserved modulated
estimate; unrestricted 3-D stability requires a separate non-axisymmetric continuation. The
parent P253 electron-and-neutrino objective remains active in every outcome.

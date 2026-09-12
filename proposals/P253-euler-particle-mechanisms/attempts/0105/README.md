# P253/0105: normalized Cao lock response and admissible source lift

## Frozen objective

This attempt advances the still-open P6 same-carrier supplier after independent
P253/0103. It has two coupled but separately verdictable obligations:

1. evaluate the physically normalized straight-column response
   `Q_P[chi_0' Phi_hat-chi_0 H_hat]` for the actual Lane--Emden/Cao profile,
   orientation and tagged band; and
2. construct the exact bounded derivative from admissible fixed-leaf
   `(delta P,delta zeta,delta chi)` variations to the full zero-mean response
   function over the tagged band together with the annular matching moments
   `(A,B)`, including the carrier/free-boundary and finite constraint rows.

A source-space determinant, a copied sign, or an unnormalized contour value
cannot answer the second obligation. A zero response activates curvature and
profile-control continuations; it does not close constant locking. A nonzero
response decides only the frozen fixed-profile route and activates the
admissible-control problem.

This is a fixed analytic construction with no empirical comparator. The
Fourier, cylindrical, Green, Lorenz, Ampere, density, charge and orientation
conventions remain those independently checked in P253/0077--0084 and
P253/0102--0103.

## Authority and exact inherited boundary

The accepted authority base is release `v0.183.0`. P253 remains an active
proposal and supplies no accepted particle claim.

This attempt consumes these independently reviewed active-proposal atoms only
at their exact scopes:

- P253/0077 with review 0081: traveling Lorenz elimination, toroidal Ampere
  primitive, modified Grad--Shafranov/Bernoulli equations and strict
  regular-band tag geometry;
- P253/0080 with review 0084: projected QAQ HSE, two bordered derivatives and
  the finite uniformly subluminal charged Cao branch;
- P253/0083 with review 0089: the connected fixed-`(kappa,R)` thin Cao path and
  its uniform common-domain carrier topology;
- P253/0088 with review 0096: physical per-member KKS/energy mode normalization
  and the source-specific off-diagonal response bridge, without a uniform Y4
  control-seed theorem;
- P253/0100 with review 0101: the exact forced Ertel/current rows and the
  distinction between zero density and a possible magnetization flux; and
- final P253/0102 with independent review 0103: the exact contour reduction,
  explicit nonlocalized Liouville cancellation, and unconstrained annular
  source determinant `-1/4`.

The final P253/0103 review and verdict SHA-256 values are
`c506400fac1377c545cb9e94666a36964b1a802f8b4f35753779674ae0c2485f`
and
`03fa9b87bef653f8bc11f04e00c5a812f20fb5dbc1061e7042b29ad02e167293`.
The frozen P253/0102 post-Liouville manifest and completion receipt are
`ffb741bc9850ded93e1a0cb984fc68226dd3c51c2b23b457e81e2f84fdde19b1`
and
`3945247c073a267c0fac5b92e82a74c69ef46bfcfecb5fe1b2396773473083ae`.

P253/0095 is excluded. Its body, evidence and conclusions will not be used.

## Domains and pairing

Work first on one compact regular positive-core band on which the uncharged
Cao action coordinate `I`, vorticity derivative, streamline period and tag are
smooth and the required denominators stay uniformly nonzero. Every division by
`chi_0` is band-local unless a new full-core tag branch is constructed.

The response is paired in the actual physical energy/KKS normalization. The
Maxwell solve uses the affine Coulomb spaces and common Green normalization of
reviewed 0077/0080. The hydrodynamic derivative uses the fixed-core source
space and bordered complement topology of reviewed 0080. Boundary traces,
zero extensions, free-boundary motion and the axis/exterior Green rows are part
of the map; they are not encoded as independent formal sources.

Let `X_leaf` be the closed tangent space satisfying the linearized Cao profile,
`delta chi=delta(zeta/lambda)` on the nonzero-tag set, circulation, mean radius,
charge, center, stabilizer/Casimir and Hodge rows. Let `Y_band` be the declared
zero-mean contour-function space over the whole tagged band in the common
physical response norm, and define

    L_resp : X_leaf -> Y_band

as the derivative of the complete interior Maxwell lock response. The
baseline target is the actual fixed-profile function
`S_Cao in Y_band`, not one scalar or two boundary coefficients. Also let

    L_AB : X_leaf -> R^2,    h |-> (delta A,delta B)

be the boundary/matching trace obtained only after solving the induced
streamfunction and Maxwell equations. The unconstrained source determinant is
only an exposing finite block of this trace; it cannot replace `L_resp` for a
generic Lane--Emden/Cao response. Let `m_0=(A_0,B_0)` denote the boundary
mismatch of the physically normalized fixed-profile response. Define

    L_full=(L_resp,L_AB): X_leaf -> Y_band x R^2,
    T_0=(S_Cao,m_0).

The actual affine cancellation condition is

    T_0 in Ran L_full,

while `<ell,T_0>=0` for every `ell in ker L_full^*` is necessary. It is
sufficient only after proving closed range/Fredholm control or constructing a
bounded right inverse on the relevant complement. Without closed range, the
adjoint condition characterizes membership in the closure of the range, so
approximate cancellation and exact cancellation must be reported separately.
Surjectivity of `L_AB` is neither necessary nor sufficient for the full-band
equation. Even within its finite trace, rank two is sufficient while rank one
can cancel `m_0` when the target lies in its range. A nonzero adjoint pairing
refutes exact and approximate cancellation; zero pairings alone establish
neither.

Before invoking any profile inverse, derive the fixed-material-leaf tangent
itself. In the straight circular limit,

    delta zeta=-xi dot grad zeta,    div xi=0,
    partial_s(s <xi_s>_alpha)=0.

Regularity at the center and compact support make
`s <xi_s>_alpha=0`. Consequently the axisymmetric radial mean of
`delta zeta` vanishes; when `chi=chi(zeta)` on the frozen leaf, the same
restriction applies to the corresponding `delta chi`. This is an exposing
rank obstruction for the formal annular source controls. On a finite torus
the proof must use the exact `r dr dz` measure and retain boundary
displacement, circulation, impulse, tag-distribution and center rows.

## Candidate routes

### Route A — actual normalized fixed-profile response

Derive the Lane--Emden/Cao straight-column limit from the same physical
orientation as the finite rings. Solve the signed Lorenz and Ampere radial
problems with one global matching constant, then evaluate the contour
zero-mean projection in the physical pairing. The proof must control the tag
cutoff and its boundary terms and must state whether the result is nonzero,
zero, or sign-changing as a function of the declared band/profile parameters.

A route verdict is earned only by an exact identity, a sign/zero theorem, or a
rigorous enclosure whose analytic remainder has already been isolated. A
formula regression or arbitrary normalization is supporting evidence only.

### Route B — bordered admissible-leaf lift

Differentiate the fixed-domain Green/free-boundary formulation with respect to
admissible vorticity-profile and tag variables. Use the reviewed QAQ complement
inverse and add only the actual finite rows. Prove the resulting derivative is
bounded on `X_leaf`, compute the full interior operator `L_resp` and its
matching trace `L_AB`, and test the target class
`[T_0] in (Y_band x R^2)/Ran L_full`. If its cokernel is nonzero, first decide
whether the exact mean-flux identity above forces range loss on the entire
declared coadjoint leaf and compute `<ell,T_0>` for its adjoint generators. A
nonzero pairing earns this cancellation route a `refuted` verdict with that
mechanism; finite-block rank loss alone does not. It then activates an
independent admissible label degree in a redesigned full-core tag/defect, a
second-order nonaxisymmetric mean response, and a materially different steady
carrier leaf. Identify the exact adjoint compatibility functional in every
case rather than closing the parent obligation.

Prove a closed-range/Fredholm estimate for `L_full` in the declared topologies
or build a bounded right inverse before using adjoint orthogonality as a
sufficiency statement. If the range is nonclosed, give separate verdicts for
`T_0 in closure(Ran L_full)` and `T_0 in Ran L_full`; a sequence of controls
with unbounded norm is not an exact admissible cancellation.

The unconstrained matrix `[[0,-1/2],[-1/2,0]]` is an exposing comparison. It
cannot be substituted for either the constrained trace `L_AB` or the full
operator `L_full`; positivity, support, charge, tag-leaf and free-boundary rows
may reduce or collapse its range.

The full-band construction can close in either of two honest ways: solve the
infinite-dimensional response equation together with its `A/B` and
free-boundary traces, or first construct an admissible exponential-law band
for which the interior response is proved identically zero and then solve the
remaining matching rows. The latter must still be a compact Cao/material-leaf
construction; the noncompact Liouville model is only its analytic seed.

### Route C — tag-domain representation change

If the existing strict-band tag prevents the desired ratio globally, compare
without conflation:

- a band-local lock with its explicit boundary/extension current;
- a redesigned full-core tag whose vanishing order matches `zeta`; and
- a distinct compact moving phase defect.

Any redesigned tag must re-enter the charged-carrier existence map and prove
its stabilizer, regularity, charge and affine Maxwell rows. A compact defect
must construct its phase map and finite action; `Theta=theta` on the
noncompact removed axis is not such a construction.

The first exposing comparator for a redesigned full-core tag is analytic. For
the Lane--Emden law `zeta=C P^p`, `P_s<0`, `P(a)=0`, and
`chi=zeta/lambda_0`, the common-normalization Volterra formula predicts on a
small punctured center annulus

    S_1(s)=-C^2 P(0)^(2p) s
           /[2(p+1)lambda_0^2 epsilon_EM c_EM^2 R]+O(s^3).

Execution must derive this from the total source moment, local electric field
and Ampere response rather than hard-code it. If confirmed, it separates the
full-core power-law tag from the Liouville cancellation. It does not evaluate
the existing strict-band cutoff, whose global Volterra moment includes its
transition regions; those are distinct routes.

### Route D — failure-derived curvature continuation

Whenever the residual class
`[T_0] in (Y_band x R^2)/Ran L_full` vanishes—because the fixed response itself
is zero or because an admissible leaf tangent cancels it—compute
the first finite-curvature coefficient on the exact Cao path with a common
domain and matched left/right normalization. If that coefficient also
vanishes, continue to the next harmonic or a materially different profile
control. A column cancellation is a branch point rather than an obligation
no-go.

## Verification contract

Reusable exact formulas belong in an importable `substrate_framework` module
with focused tests. The strongest practical oracle must independently derive:

1. the signed radial Lorenz/Ampere response and contour projection;
2. the Lane--Emden/profile identity actually used in Route A;
3. the domain, band topology and finite traces of `L_resp` and `L_AB`;
4. the range or explicit adjoint cokernel of `L_full` together with its pairing
   against the actual full target `T_0`; and
5. the closed-range/Fredholm estimate or bounded right inverse used to turn
   adjoint orthogonality into exact solvability; and
6. every claimed curvature coefficient if Route D activates.

Symbolic substitutions and formula regressions will be labeled accordingly.
No production numerical claim is preregistered. If a small singular value or
rank decision later requires numerics, a new append-only attempt will load the
small-ratio prescription before freezing that verifier.

The normalized fixed-profile response in Route A is evaluated independently
of Route B. An inadmissible formal profile/source variation cannot be used to
cancel a nonzero physical response.

## Route verdicts and exclusions

Each of A, B, C and D receives exactly one route-scoped verdict: established,
refuted with mechanism, or blocked at its named missing construction. A route
verdict does not decide P6 or the parent campaign.

This attempt cannot establish P2 persistence, a normalized universal P4
action, Lorentz chirality, weak coupling, flavor oscillation, Born/reset,
fermionic exchange, an electron, a neutrino, LP2, LP4, LP5, LP6 or parent
completion. Those remain separate live obligations.

## Activation and ownership

Root owns P253/0105 and will not edit P253/0095. Before substantive 0105 body
work, the coordinator will centrally register this attempt and capture the
repository-schema command, stdout, stderr and exact zero exit under this
directory. At preregistration this directory contains only `README.md`.

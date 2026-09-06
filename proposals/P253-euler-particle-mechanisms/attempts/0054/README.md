# P253/0054: common-domain moving-boundary Euler graph-jet theorem

Status: **README-only preregistration awaiting central registration and schema
activation.** Until `0054/activation-schema.exit` exists and contains exactly
`0`, this README is the sole 0054 artifact. Root owns `proposal.yaml`, the
central attempt registry, activation receipts, commits, review, and promotion.
No derivation body, verifier, source download, numerical run, API, or claim is
authorized before activation.

## Parent boundary and positive theorem target

Attempt 0052 established the corrected crossing law

    n_*(delta)=A/(c delta L)-A d/(c^2 delta L^2)
       +(B/A-beta_0)+o(1/(delta L^2))+O(1),
    L=log(1/delta),                                      (1)

the exact weighted action--angle/Piola representation, the physical
KKS/energy/frequency normalization, and the full-sector limiting Grushin
reduction. It proved the exact-Cao graph-Riesz implication conditional on the
named common-domain, modewise global-Euler graph jet `HJ2`. The leading scales

    g_(n_*)=O(delta^3 L),
    Delta sigma_(n_*)=Theta(delta^3 L^2),
    g_(n_*)/Delta sigma_(n_*)=O(L^-1)                    (2)

remain fixed.

The positive target of 0054 is the missing **HJ2 theorem** for one actual
polynomial-profile Cao vortex-ring subfamily. After the physical dilation is
factored and the 0052 chart is applied, twice differentiate the complete
moving-boundary, whole-space Euler generator and prove either:

1. a strong common-graph polyhomogeneous expansion with uniform order-minus-one
   and reduced order-minus-two modewise bounds; or
2. a common-form, weighted norm-resolvent estimate that implies the same
   contour and Riesz conclusions without asserting strong graph
   differentiability.

Passing either route discharges the explicit 0052 `HJ2` dependency at exactly
its stated scope. It does not itself construct a nonlinear rotating branch,
prove nonlinear stability, or license a quantum, electron, or neutrino
interpretation.

## Fixed physical object and quantifiers

Fix an integer Cao profile exponent `p>=6`, circulation `Gamma>0`, density
`rho_0>0`, one toroidal Fourier number `ell>=2`, and a sufficiently small
thinness interval `0<delta<delta_0`. Let

    delta=a_delta/R_delta,       k_delta=ell delta,
    Omega_delta=s_delta e_theta,
    A_delta eta=curl(u_delta cross eta
                     +B eta cross Omega_delta),           (3)

where `B=curl(-Delta_R3)^-1` is the decaying whole-space Biot--Savart/Hodge
operator and (3) is interpreted on compact-core dynamically accessible (DA)
vorticity. The physical velocity, pressure, collar traces, and exterior
harmonic field are always those obtained from the global operator; they are
not replaced by a core-only elliptic surrogate.

The regularity threshold `p>=6` is frozen to supply two shape derivatives,
one generator derivative, and the required interface multipliers. The proof
must derive the actual trace threshold from the positive-part profile. A
better threshold may be recorded as a corollary, but failure at a smaller `p`
cannot refute the frozen `p>=6` theorem.

## Common DA quotient and fixed chart

Let `K_*` be the fixed solid torus of 0052 with polar measure

    dmu_*=s ds d beta d theta.                             (4)

The cross-`delta` Hanzawa/action--angle map

    Phi_delta:K_* -> K_delta                              (5)

is orientation preserving, not volume preserving between unequal physical
carrier volumes. Its core Jacobian relative to (4) is the explicit constant
`Jbar_delta=2 I_a(delta)` after normalized Moser correction. Its collar
extension has a variable positive Jacobian `J_delta(y)` with uniform relative
bounds after physical scaling. Use the full contravariant Piola transform

    P_delta eta=J_delta^-1 D Phi_delta eta o Phi_delta^-1, (6)

including every determinant and cofactor factor. Only a displacement within
one fixed carrier is a determinant-one physical coadjoint motion.

Do not identify the DA ranges for different `delta` by declaration. Start
from the fixed displacement space

    Y_*={xi: div_* xi=0, xi has the declared collar trace}/Stab_*,  (7)

and the exact tangent maps

    C_delta xi=P_delta^-1 curl[(Phi_delta)_*xi
                               cross Omega_delta].         (8)

Construct finite-dimensional impulse, centering, and stabilizer duals and
prove that the gauge-fixed `C_delta` have a common closed quotient domain and
uniformly equivalent ranges. This earns a single fixed DA energy space `X_*`;
no unproved global Clebsch completion is allowed.

Let `S_delta=L^-1/2 P_c+(1-P_c)` be the 0052 logarithmic centerline rescaling
and set

    U_delta=P_delta S_delta,                              (9)
    Ahat_delta=U_delta^-1 A_delta U_delta.                (10)

The exact kinetic metric is pulled back from
`rho_0 integral_R3 |B(P_delta eta)|^2 dx`. Define

    D_*={eta in X_*: T_0 eta+K_0 eta in X_*},             (11)

with the nonzero streamline sectors carrying one `beta` derivative and the
`m=0` sector carrying the compact Kelvin graph norm. Prove uniform energy and
graph equivalence for (9)--(11), including both directions and the finite KKS
quotient.

## Obligation A: two shape derivatives with interface closure

Differentiate the exact Cao free-boundary/profile equations, circulation,
impulse, center gauge, and action inversion to construct

    Phi_delta=Phi_0+delta H_1+delta^2 H_2+o_(C^2)(delta^2),
    J_delta=J_0+delta J_1+delta^2 J_2+o_(C^1)(delta^2),    (12)

after the physical dilation has been separated. If the chart itself has a
source-forced logarithmic cell, isolate it rather than placing it in the
remainder. Differentiate (6) explicitly through second order and verify the
Piola divergence and DA curl covariance at each order.

On the core interface derive the vanishing order of `Omega_delta` and
`C_delta xi` from `(Psi_delta)_+^p`. Prove that two parameter derivatives
produce no unaccounted vortex sheet, identify the exact normal/tangential
trace spaces, and show that the global Hodge velocity has the correct
interior/exterior trace and pressure transmission. The collar diffeomorphism
is fixed to the identity outside a declared compact set; its variable
Jacobian is retained throughout.

## Obligation B: twice differentiated whole-space Hodge/Leray operator

Write the conjugated Biot--Savart kernel before differentiating:

    B_delta(y,y')=P_delta^-1 curl_x(-Delta_x)^-1
                   P_delta |_(x=Phi_delta(y),x'=Phi_delta(y')).   (13)

Split (13) into near-diagonal core, interface/collar, and smooth exterior
pieces. For the near kernel, subtract the homogeneous singular kernel and
derive the first two commutator/shape derivatives, including all derivatives
of `J_delta`, `D Phi_delta`, and the Euclidean distance. For the velocity
formulation also derive

    P_Leray,delta=I-grad_(g_delta)
       Delta_(g_delta)^-1 div_(g_delta)                   (14)

through second order on the same domain. Equations (13)--(14) are two
representations of the same global Hodge response and must agree on the
vorticity-velocity bridge.

The direct route must prove, for a fixed admissible Sobolev interval
`s in [s_0,s_1]`, uniform parameter-symbol estimates

    ||Lambda^(s+1) partial_delta^j B_delta Lambda^(-s)||<=C_s,
       j=0,1,2,                                           (15)

and the corresponding trace/collar/exterior estimates. Here `Lambda` is the
fixed action--angle/Kelvin elliptic weight, not a boxed Fourier surrogate.
The local `delta^2 L` kernel coefficient is extracted explicitly as
`C_(2,log)` before (15) is applied to the finite remainder.

## Obligation C: common-domain Euler expansion

Combine the differentiated profile, transport/stretching, Piola, Hodge, and
Leray blocks without imposing the background equation before all derivatives
are taken. Establish on the common graph domain

    Ahat_delta=A_col(k_delta)+delta C_1
      +delta^2 L C_(2,log)+delta^2 C_2+R_delta.            (16)

The longitudinal `i ell delta` terms already present in
`A_col(k_delta)` are not counted again as curvature. The exact streamline
chart must show which first-order blocks preserve `m` and which map
`m -> m+/-1`; a coordinate harmonic alone is not the selection rule.

Freeze the following strong-route target. For every `N` needed by the Schur
summation, the order-minus-one coefficient blocks obey

    |<e_(m,n),B_j e_(m',n')>_E|
      <=C_N <m-m'>^-N <n-n'>^-N
        /<|m|+|m'|+n+n'>,                                (17)

and the nonzero-`m` Schur correction to the `m=0` block obeys the same bound
with the last denominator squared. After the diagonal order-minus-one symbol
is absorbed into the exact renormalized Kelvin frequencies, require

    |<e_(m,n),R_delta e_(m',n')>_E|
      <=epsilon_delta delta^2 <m-m'>^-N <n-n'>^-N
        /<|m|+|m'|+n+n'>,
    epsilon_delta ->0,                                   (18)

with the analogous translation/Kelvin estimate

    |<t,R_delta e_(0,n)>_E|
       <=epsilon_delta delta^2/(1+n).                     (19)

Here `<q>=1+|q|`; zero denominators are interpreted through the finite
translation/stabilizer projections. The proof must state the precise
orthonormal or KKS-biorthogonal basis and cannot infer (17)--(19) from
unstructured form convergence.

At (1), (18)--(19) give

    epsilon_delta delta^2/n_*
       =O(epsilon_delta delta^3 L)
       =o(delta^3 L^2),                                  (20)

while an order-minus-two second Schur block gives

    delta^2/n_*^2=O(delta^4 L^2)
       =o(delta^3 L^2).                                  (21)

Equations (20)--(21), not merely `o(delta^2)` in an unidentified norm, are the
bridge to the 0052 contour radius.

## Route A: direct kernel and graph theorem

Route A earns the theorem by:

1. constructing (12) with uniform inverse/Jacobian bounds;
2. proving the differentiated singular-kernel cancellations in (13), the
   metric Leray identities in (14), and the order-minus-one estimates (15);
3. deriving the interface, collar, and exterior transmission estimates;
4. proving common-domain equivalence and the expansion (16);
5. deriving the modewise bounds (17)--(19) by integration by parts in the
   actual action/Kelvin variables; and
6. inserting (20)--(21) into the full-sector Grushin inverse and Riesz contour
   of 0052.

The strongest practical exact oracle will independently expose a missing
Jacobian/cofactor term, a wrong shape-commutator sign, double counting of the
longitudinal mode, loss of one pseudodifferential order, and comparison of an
unweighted remainder directly with the resonant gap.

## Route B: common-form and norm-resolvent fallback

If the positive-part interface has the required form regularity but not two
strong graph derivatives, keep the same carrier and construct a common closed
Hamiltonian/KKS form on a fixed form domain `V_*`. Derive the exact form
coefficients through `delta^2 L` and `delta^2`, including the exterior Hodge
energy. Let `A_app,delta` denote the operator represented by those forms.

Route B must prove directly on the 0052 contour `Gamma_delta`, whose distance
from the complement is

    r_delta=gamma delta^3 L^2,                            (22)

the weighted relative-resolvent estimates

    sup_(z in Gamma_delta)
      ||(Ahat_delta-A_app,delta)
          (A_app,delta-z)^-1||_(X_* -> X_*)
       <=epsilon_delta<1/2,                               (23)

    sup_(z in Gamma_delta)
      ||(Ahat_delta-z)^-1-(A_app,delta-z)^-1||_(X_* -> D_*)
       <=C epsilon_delta/r_delta.                         (24)

The form proof must derive (23) from order-minus-one/order-minus-two weighted
bounds equivalent to (17)--(21); abstract Mosco or strong-resolvent
convergence alone is insufficient because `n_*(delta)` diverges. Contour
integration of (24) must give an `o(1)` graph-bounded Riesz-projection
difference and preserve the real/KKS rank. Success of Route B discharges HJ2
without claiming the unavailable strong Taylor map.

## Failure-derived continuation inside 0054

If Route A loses derivatives only at the moving interface, first repair it by
good-unknown/paradifferential conjugation and repeat the kernel estimate. If
that still yields only form control, execute Route B. If a global collar chart
prevents uniform estimates, switch to a core boundary-integral plus exterior
Dirichlet-to-Neumann representation and prove its equivalence to (13). A route
verdict applies only to that representation; none of these failures implies
nonexistence of the Cao ring or rotating branch.

An unresolved scalar coefficient is not licensed for numerical production by
this README. Exact source recurrence and kernel calculus are exhausted first.
If a genuinely irreducible numerical remainder is later found, a separate
append-only expansion must freeze its object, scales, error enclosure, and
maximum verdict before execution.

## Source inventory and applicability boundary

Full papers remain in `/tmp/primary-source-cache`, not in campaign artifacts.

- Cao--Lai--Qin--Zhan--Zou, arXiv:2206.10165v2: exact weighted Green kernel
  (2.1)--(2.4), steady equation (3.1), rescaling (3.11)--(3.12), parameter
  system (3.36), free-boundary graph Lemma A.2, and odd inverse Theorem C.2.
  These define the carrier, shape cells, and kernel. Their printed elliptic
  remainder is not treated as the HJ2 graph theorem.
- Gallay--Smets, arXiv:1805.05064v3: exact column Fourier equations,
  Biot--Savart relations, and limiting energy estimates. Its ambient
  full-enstrophy essential spectrum is not transferred to the closed-core DA
  space.
- Attempts 0044, 0048, and 0052: append-only derivations at their reviewed or
  explicitly conditional scopes. The universal filament matrix is not used as
  the finite-core operator, and the 0048 formal residual is not imported as an
  exact solution.

Standard Hanzawa, Moser, Piola, Calderon--Zygmund, shape-commutator, and closed-
form facts may organize the proof only after their hypotheses are verified for
the source-defined Cao interface and the exact domains above. A theorem name
does not supply (17)--(24).

## Obligation and verdict ledger

| Node | Positive intent | Pass licenses | Does not license | Maximum verdict | Failure scope | Unlocks |
| --- | --- | --- | --- | --- | --- | --- |
| HJ2-A | Twice differentiated fixed-domain profile/Piola/interface chart | Common DA quotient and coefficient jets | Riesz control without global Hodge estimates | Exact analytic | This chart/regularity route | HJ2-B, HJ2-C |
| HJ2-B | Whole-space Hodge/Leray order-minus-one derivatives with collar/exterior closure | Global velocity and pressure blocks in (16) | Modewise decay by itself | Exact analytic | Direct kernel representation | HJ2-C |
| HJ2-C | Strong graph bounds (17)--(21) or fallback (23)--(24) | Exact 0052 full-sector Cao Riesz family | Nonlinear relative-periodic branch | Exact analytic or conditional with named remaining construction | Selected graph/form representation | P253 Cao rotating-wave continuation |

Every executed route receives exactly one route-scoped verdict: established,
refuted with its mechanism, or blocked with its missing construction. A support
gap keeps HJ2 active and triggers the continuation ladder; it is not a physical
no-go or campaign completion.

## Frozen success contract

0054 succeeds only if one of the following is established for the actual Cao
carrier and global Euler operator:

1. equations (12)--(21), including common graph domain, interface/collar/
   exterior closure, and full modewise bounds; or
2. the common-form construction and equations (22)--(24), with the same
   fixed-domain physical reconstruction and contour consequence.

The completion record must preserve exact source locations, chart and domain
definitions, all determinant/cofactor and density factors, first execution
commands and outputs for exposing symbolic checks, route verdicts, and the
strongest supported theorem. No production numerics are part of this attempt.
The parent P253 objective remains active, and no nonlinear branch claim is made
by either success alternative.

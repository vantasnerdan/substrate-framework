# P253/0025 preregistration: repeated moving-fibre Kelvin--Leray return

Status: **README-only preregistration awaiting central registration and schema
activation.**  Until `0025/activation-schema.exit` exists and contains `0`,
this attempt opens no new source body, performs no new carrier calculation,
executes no symbolic or numerical verifier, and writes no file other than this
README.  Root owns the central append-only registration and activation.  This
worker owns only eventual `0025` artifacts; source modules, tests, governance,
the proposal manifest, memory, and all other attempts remain outside its write
surface.

## Parent result retained and exact open dependency

Attempt 0019 established a deliberately finite-time statement for the same
Euler family.  For every sufficiently small fixed smooth compact finite-`j`
Gavrilov pressure-shell member, the exact whole-space Leray linearisation has
finite-energy dynamically accessible axisymmetric data with one-circuit
kinetic-energy gain greater than `1+pi^2/2`.  Its solvable small-shell limit has

    M_0=[[1,0],[-sqrt(2)pi,1]],       |M_0(1,0)|^2=1+2pi^2,          (1)

but 0019 did not iterate this matrix on the actual carrier.  Finite-shell twist
moves the wave covector between polarization fibres, and the fixed-time WKB
remainder used there gives no estimate at a circuit count tending to infinity.
That exact transient-gain scope is preserved: 0025 neither weakens it nor
silently promotes it to instability.

The positive deliverable here is the actual repeated moving-fibre
Kelvin--Leray skew product, transferred to the returned **packet profile** and
then to the exact global Euler semigroup in the physical energy metric.  It
must decide between two scientifically different outcomes:

1. an explicit polynomial or modulated bound for repeated returns, with the
   global pressure, annulus collars, support boundary, exterior velocity tail,
   and symmetry quotient controlled; or
2. genuine unbounded amplification by exact dynamically accessible solutions,
   first as an unbounded semigroup-norm sequence and, if the required
   almost-orthogonality estimates close, as one fixed finite-energy solution
   with an unbounded time sequence.

These outcomes are competitors, not success/failure labels for a sampled
matrix.  A bounded modulated cocycle may license linear LP2 persistence in the
declared sector; a polynomial upper bound licenses only that growth ceiling;
an unbounded packet sequence refutes uniform linear boundedness in the declared
norm.  None by itself supplies nonlinear particle stability, quantum mechanics,
relativity, or electron/neutrino identity.

| Field | 0025 frozen value |
| --- | --- |
| `positive_intent` | Repeated same-carrier global Kelvin--Leray propagation, with an exact polynomial/modulated bound or a genuine unbounded accessible-packet construction |
| `requires` | 0019 global generator/domain, exact coadjoint tangent, regular Gavrilov annulus, one-circuit transient gain, action sign, and moving-fibre warning; 0005/0010 orbit action and finite-`j` carrier |
| `pass_licenses` | A time-unbounded axisymmetric linear propagation verdict for the named fixed carrier and physical quotient norm |
| `does_not_license` | Nonlinear LP2; unrestricted 3-D stability; a restoring particle frequency; scale selection; exchange/Born rules; relativity; electron or neutrino identification |
| `maximum_verdict` | A theorem-backed global semigroup bound with explicit growth law, or an exact dynamically accessible unbounded-growth construction on one fixed carrier |
| `failure_scope` | Only the fixed Gavrilov member, perturbation sector, profile class, time horizon, and norm for which the full remainder is closed |
| `unlocks` | Nonlinear modulation if uniformly controlled; a materially different stable carrier if uniform boundedness is refuted; nonaxisymmetric completion after the axisymmetric verdict |

The parent campaign and LP2 obligation remain active regardless of either
single-route verdict.  No exhaustion claim is registered.

## Fixed carrier, physical domain, and three structures that may not be merged

After activation select **one fixed** sufficiently small `sigma>0` satisfying
0019's annulus estimates and retain it throughout the theorem.  The carrier is

    u_sigma=W(a/sigma)u,
    A_sigma={5sigma/4<=a<=7sigma/4},                                 (2)

with `W=1` on the displayed annulus and flat support collars.  Sending `sigma`
to zero may establish constants before this choice; it may not replace the
final fixed-carrier statement.  No annulus wall is introduced.

Let `P` be the whole-space Leray projector and

    G_sigma v=-P[(u_sigma dot grad)v+(v dot grad)u_sigma],
    S_sigma(t)=exp(t G_sigma).                                      (3)

Use 0019's skew-adjoint projected-transport realization plus bounded shear,
not an unproved assertion that the maximal domain equals `H^1` in the flat
exterior.  The invariant dynamically accessible space and graph domain are

    X_DA=closure_X{-ad_w^*m_sigma:w smooth, div w=0},
    D_DA={q in D(G_sigma) intersect X_DA:G_sigma q in X_DA},         (4)

with the Euclidean symmetry directions removed by quotient distance, not by
assuming an arbitrary orthogonal complement is invariant.

Three distinct geometric structures are frozen and separately audited:

1. **Orbit symplectic form.**  On the smooth stabilizer quotient use the KKS
   form `Omega_sigma`, convention `i_X Omega=dH`, and quadratic action

       S2=int[-Omega_sigma(q,qdot)/2-Q_sigma(q,q)/2]dt.              (5)

   Thus `qdot=J_sigma L_sigma q`.  Prove the descent/domain on which
   `S_sigma(t)^*Omega_sigma S_sigma(t)=Omega_sigma`; do not infer it from a
   two-by-two determinant.  The indefinite conserved Hamiltonian form
   `Q_sigma` is not the positive kinetic metric.
2. **Ray symplectic form.**  The characteristic map is a cotangent lift and
   preserves `omega_can=sum_i dk_i wedge dx_i`.  This phase-space fact does
   not identify `omega_can` with the KKS form or quantize the orbit.
3. **Physical energy metric.**  The exact observable is

       ||q||_X^2=rho0 int_R3 |delta u[q]|^2 dx
         =2pi rho0 int_Pi[eta K eta+chi^2/r^2]dnu.                  (6)

   On a compact regular annulus its principal fibre metric is

       h_E,x,k=2pi rho0[r^2|eta|^2/|k|^2+|chi|^2/r^2],              (7)

   equivalently Euclidean in `Y=(r eta/|k|,chi/r)`.  The full theorem uses
   (6), while (7) fixes the energy-normalized symbol and its connection terms.
   Growth means growth in (6) after symmetry modulation, not merely growth of
   vorticity derivatives or a coordinate amplitude.

## Exact continuous characteristic and Poincare return

Use the exact travel coordinate on the regular annulus,

    Phi_t(p,vartheta)=(p,vartheta+Omega(p)t),
    dnu=[T(p)/(2pi)]dp dvartheta,       Omega(p)=2pi/T(p).           (8)

For an advected covector `k=k_p dp+k_vartheta dvartheta`, the uniform-time
cotangent lift is

    F_t:(p,vartheta,k_p,k_vartheta)
      ->(p,vartheta+Omega(p)t,
          k_p-t Omega'(p)k_vartheta,k_vartheta).                    (9)

Equation (9), including its sign and the density in (8), is rederived from the
actual flow.  It is not replaced by powers of (1).  In three dimensions retain
the actual azimuthal advance

    DeltaTheta(p)=int_0^T(p) u_theta(x_p(t))/r_p(t) dt.              (10)

On a meridional section the base Poincare map and its cotangent lift are

    Pi(p,theta)=(p,theta+DeltaTheta(p)),
    (k_p,k_theta)->(k_p-DeltaTheta'(p)k_theta,k_theta).              (11)

The axisymmetric sector has `k_theta=0`, but it still has the meridional twist
in (9).  Nonaxisymmetric sectors remain visible in (10)--(11) and may not be
declared controlled by an axisymmetric estimate.

Along the actual three-dimensional trajectory retain the complete Kelvin--
Leray system

    xdot=u_sigma(x),
    kdot=-(grad u_sigma(x))^T k,
    Adot=-(grad u_sigma)A
       +2k[k dot (grad u_sigma)A]/|k|^2,
    k dot A=0.                                                       (12)

The second term in the amplitude equation is the order-zero pressure
projection.  Let `C(z)` be the induced operator on `k^perp`, in a declared
parallel frame, for phase point `z=(x,k)`, and define the actual one-step map

    M(z)=Pexp[int_0^T C(F_t z)dt].                                  (13)

Its repeated object is the skew product

    z_{j+1}=F z_j,
    M^(n)(z_0)=M(z_{n-1})...M(z_1)M(z_0).                           (14)

The fibrewise physical gain is

    Gamma_n(z,A)=h_E,z_n(M^(n)A,M^(n)A)/h_E,z(A,A).                 (15)

All frame identifications and metric factors in (15) are explicit.  A
determinant, eigenvalue, or Jordan form of one `M(z)` cannot decide (14)--(15).
If the polarization cocycle has a symplectic form, derive its transported form
and sign; do not equate it with either symplectic structure above by notation.

## The actual returned profile, not only a central ray

For a smooth compact annular profile `a` and integral rapid phase, let

    I_N(a)=BiotSavart[J_sigma(0,N^-1 Re(a exp(iN phi_0)))]           (16)

denote the exact finite-energy dynamically accessible injection, including
the divergence-free reconstruction and symmetry modulation.  Let
`phi_t=phi_0 composed Phi_{-t}`.  The profile-return target is the physical-
norm identity

    S_sigma(t)I_N(a)
      =I_{N,t}[C_t(a composed Phi_{-t})]+E_N(t)a,                    (17)

where `C_t` is the full ordered cocycle from (12), `I_{N,t}` uses the advected
phase and energy normalization, and `E_N(t)` includes subprincipal transport,
Leray/Biot--Savart pseudodifferential remainders, chart transitions, collars,
and the exterior tail.  Equation (17) is proved uniformly over the whole
profile support.  Evaluating (13) on one selected streamline is only the
leading-ray input to (17).

At reference streamline times `t_n=nT(p_*)`, the returned profile generally
does not reproduce itself: neighbouring streamlines acquire the exact phase

    N[Omega(p)t_n-2pi n].                                          (18)

Two regimes are kept distinct.  A coherent-return construction must prove its
chosen width `delta_N` makes the phase variation in (18) small.  A phase-mixing
construction retains (18), integrates its cancellations in the physical
profile norm or in the nonlocal kernel, and may not call dephasing an error.

## Finite-time versus growing-time WKB license

The fixed-time estimate used in 0019 has the form `||E_N(T)||=O_T(N^-1)`.
It licenses no conclusion at `t=t_N->infinity`.  Before choosing a growing
time, derive from (3), (9), (12), and the exact profile derivatives an explicit
bound of the form

    ||E_N(t)a||_X
      <= N^-1 delta_N^(-s) P_m(t)
         exp[int_0^t Lambda_sigma(tau)d tau] ||a||_Y,                (19)

with derived integers `s,m`, polynomial `P_m`, coefficient
`Lambda_sigma`, and graph norm `Y`.  Every derivative loss caused by twist,
the moving frame, energy normalization, and the global projector is counted.

A proposed circuit sequence `n_N` and `t_N=n_N T(p_*)` earns the WKB bridge
only if

    ||E_N(t_N)a||_X=o(||I_{N,t_N}C_t a||_X)                         (20)

and its tube/profile remains inside the declared regular charts.  If (19) is
exponential, a logarithmic horizon such as `t_N=c log N` is used only after the
admissible `c` is derived.  If analytic cancellations reduce (19) to
polynomial growth, a power-law horizon is used only within its derived
exponent.  Neither `N` nor `delta_N` may be chosen after observing a sampled
gain.  For coherent packets also require the exact phase-width condition from
(18); for mixed packets replace that condition by a proved nonstationary-phase
or kernel estimate.

## Global annulus, collar, and exterior ledger

Choose a smooth identity partition

    1=chi_A+chi_C+chi_E,                                            (21)

where `chi_A` lies in the regular plateau annulus, `chi_C` contains both flat
cutoff collars and the remaining compact carrier, and `chi_E` is exterior.
This is an operator decomposition, not a boundary condition.  In vorticity
variables retain every block

    K_ij=chi_i K chi_j,       i,j in {A,C,E},                       (22)

and the commutators generated by moving the cutoffs through `K`; in velocity
variables retain the corresponding `[P,chi_i]` pressure blocks.  The exact
interaction-frame identity is

    S_sigma(t)-U_A(t)
      =int_0^t S_sigma(t-s)R_global U_A(s)ds,                        (23)

and repeated-return estimates use the full accumulated integral or discrete
Duhamel sum, not `n` times an unnamed one-circuit constant.

The analytic ledger separates:

- adjacent annulus--collar commutators, controlled in the declared graph norm;
- separated annulus--exterior Green kernels, with their actual smoothing and
  high-frequency integration-by-parts power;
- the global pressure/velocity tail, including which vorticity moments make
  its leading multipoles vanish and which remain;
- the flat support boundary, where base coefficients vanish to all orders but
  the perturbation velocity remains nonlocal;
- transport of exact dynamically accessible vorticity, which has no invented
  exterior wall; and
- symmetry-mode leakage and modulation in the quotient metric.

Every block receives a time-dependent bound compatible with (19)--(20).
Smallness for one circuit is not summability over infinitely many circuits.

## Competing routes and predeclared verdicts

### Route A -- twist-normal-form and polynomial/modulated control

Conjugate the full operator by the exact characteristic/profile transport and
derive the repeated cocycle in the energy variables (7).  Compute the circuit
average and resonant part before estimating it.  If an explicit modulator
`N_res(t)` removes the nonintegrable nilpotent or phase term, prove

    sup_t ||N_res(t)^(-1)U_A(-t)S_sigma(t)||_(Y->X)<infinity         (24)

and state the physical growth of `U_A(t)N_res(t)`.  Direct boundedness gives
axisymmetric linear persistence in the named quotient norm.  A bound
`C(1+|t|)^alpha` with `alpha>0` is still a useful global propagation theorem,
but it does not supply LP2 stability.  Bounded physical velocity with growing
vorticity gradients is recorded as phase mixing, not particle destruction.

### Route B -- genuine unbounded dynamically accessible packets

Use the actual product (14), profile identity (17), and growing-time license
(19)--(20) to seek normalized exact data `q_j in D_DA` and times `t_j->infinity`
such that

    ||q_j||_X=1,       ||S_sigma(t_j)q_j||_X->infinity.              (25)

Equation (25) proves the semigroup norm is unbounded on the declared accessible
space.  It is stronger than repeated leading-ray gain.  The uniform
boundedness principle then gives some `q in X_DA` with
`sup_t||S_sigma(t)q||_X=infinity`, but does not by itself give a smooth or
graph-domain packet with controlled localization.  The stronger constructive
route chooses a summable almost-orthogonal sequence of packets and proves that

    q=sum_j c_j q_j in X_DA,
    sup_j ||S_sigma(t_j)q||_X=infinity,                              (26)

with cross-packet evolution controlled, rather than assuming orthogonality at
later times.  If only (25) closes, report the operator-norm result and the
Banach--Steinhaus `X_DA` consequence separately from the stronger `D_DA`
profile construction.
If a single exponential multiplier is found, it still needs (17)--(20) and
the global ledger before it contributes to (25).

### Route C -- exact resolvent or defect-measure continuation

If derivative growth makes every available WKB horizon shorter than the
circuit count required by Route B, that is a representation limit, not a
stability result.  Continue on the same carrier by either an exact
characteristic resolvent for `lambda-G_sigma` or a semiclassical defect measure
that is propagated for fixed times and concatenated with a proved uniform
compactness estimate.  A resolvent claim must connect its half-plane bound to
the time-domain norm through a stated semigroup theorem.  A defect measure
that sees ray amplification without norm recovery licenses only a necessary
mechanism.

Routes are compared in this order: exact fixed-carrier and open-space domain;
dynamical accessibility; physical kinetic metric and symmetry quotient;
retention of moving covectors and pressure; uniform profile rather than
central-ray control; explicit time horizon; collar/exterior closure; strength
of the time-domain verdict; assumption and parameter economy; then
implementability.  No empirical or sampled growth value selects a route.

## Source-access boundary

This preregistration opens no new external body.  After activation, reused
0019 sources retain their URL/version/hash/theorem-location receipts.  Any new
Kelvin-wave, geometric-optics, inviscid-damping, semigroup, or resolvent source
is inventoried before its body is opened; full papers are cached only under
`/tmp/primary-source-cache/P253-0025`, while the attempt retains URL, version,
hash, and exact theorem/equation location.  A paywall or absent source is not a
physical verdict because (3), (9), (12), and (17) remain direct derivation
targets.

Accepted abstract APIs may check an already-derived finite product or
semigroup implication.  They may not be treated as the Euler realization,
physical energy map, accessible profile, or global-pressure closure.

## Analytic ladder and production-oracle boundary

After central activation, the work order is frozen:

1. pin one fixed Gavrilov member and derive its exact `T`, `Omega`,
   `DeltaTheta`, twist signs, frame, and full profile support;
2. derive (9)--(15), including all metric connections and the relationship--or
   nonrelationship--among the three symplectic/metric structures;
3. construct the exact accessible injection and profile comparison (16)--(18);
4. derive the growing-time remainder (19), freeze `delta_N,n_N,t_N`, and prove
   the bridge condition (20);
5. close every block of the global ledger (21)--(23) with time dependence;
6. execute Route A and Route B against the same carrier and metric, continuing
   through Route C if the WKB representation rather than the physics is the
   obstruction; and
7. assign one route-scoped verdict to each attempted route and preserve the
   strongest positive propagation statement.

No production numerical verifier is registered or permitted by this README.
Exact algebra, calculus, profile transport, asymptotic remainder, and operator
bounds must first leave one named irreducible proposition and state what its
resolution would license.  If that proposition is later numerical--especially
a near-unit multiplier, small growth exponent, pseudospectral boundary, or
long-time cancellation--a separately activated design must load the
small-ratio-numerics skill and freeze precision, conditioning, domain/tail and
time refinement, residuals, nonnormal sensitivity, analytic limiting case,
and maximum representation-scoped verdict before its first production run.

The strongest oracle presently frozen is an exact repeated-characteristic and
profile-return derivation plus a continuum physical-norm estimate.  A boxed
eigenvalue, powers of (1), fixed-time `O(N^-1)` remainder used at `t_N`, or a
sampled unit multiplier cannot earn the 0025 verdict.

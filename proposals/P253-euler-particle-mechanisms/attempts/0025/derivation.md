# Exact moving-fibre reduction and fixed-time diagonal bridge

## 1. Fixed carrier and source-derived action chart

Fix one sufficiently small smooth compact Gavrilov member whose velocity
cutoff is identically one on a closed pressure interval.  Let
`B=[I_-,I_+]` be a compact action interval inside that plateau.  This carrier
is retained throughout; no limiting carrier is substituted into a theorem.
It is smooth on `R^3`, flat at its support boundary, and finite energy.

Baldi's Theorem 1.1 and equations (3.52)--(3.68) give an analytic
diffeomorphism

    x=Phi(sigma,beta,I),       (sigma,beta) in T^2, I in (0,I*),   (1)

in which the actual particle flow is

    sigma_dot=Omega_1(I),  beta_dot=Omega_2(I),  I_dot=0.          (2)

On the selected cutoff plateau the source expansions are

    K(I)=I+(1065/1024)I^3+O(I^4),
    Omega_1=K'=1+(3195/1024)I^2+O(I^3),
    Omega_2=sqrt(I) R(I) Omega_1,
    R(I)=1+(7/4)I+O(I^2).                                         (3)

Writing `a=Omega_1'` and `b=Omega_2'` therefore gives

    a(I)=(3195/512)I+O(I^2),
    b(I)=1/(2sqrt(I))+(21/8)sqrt(I)+O(I^(3/2)),
    D(I):=b(I)/a(I)
        =(256/3195)I^(-3/2)(1+O(I)),
    D'(I)=-(128/1065)I^(-5/2)(1+O(I)).                             (4)

For a sufficiently small plateau, `a,b>0` and `D` is nonconstant and strictly
decreasing.  Its image contains rationals.  Thus there are an exact torus
`I=I_0` and coprime positive integers `(n,m)` with

    D(I_0)=n/m,       n a(I_0)-m b(I_0)=0.                         (5)

This is selection from a source-derived analytic range, not fitting a growth
rate.  Condition (5) concerns derivatives of the frequencies.  It neither
assumes nor implies rationality of `Omega_2/Omega_1`.

The chart flow and its derivative are

    f_t(sigma,beta,I)
      =(sigma+Omega_1(I)t,beta+Omega_2(I)t,I),
    F_t=Df_t=[[1,0,t a],[0,1,t b],[0,0,1]],  det F_t=1.             (6)

The physical deformation is

    Dvarphi_t=D Phi(f_t z_0) F_t D Phi(z_0)^(-1).                  (7)

The chart and its inverse have uniform singular-value bounds on
`T^2 x B`.  Those bounds are used whenever a chart component is compared with
physical length; the chart is not assigned an artificial Euclidean metric.

## 2. Exact covector resonance

The globally single-valued integral phase

    theta_0=n sigma-m beta                                          (8)

has transported eikonal

    theta_t=n sigma-m beta-[n Omega_1(I)-m Omega_2(I)]t.            (9)

Its chart covector is the exact cotangent lift

    k_z(t)=F_t^(-T)(n,-m,0)
          =(n,-m,-t[n a(I)-m b(I)]).                              (10)

On `I=I_0`, (5) makes `k_z(t)=(n,-m,0)` for every time.  The physical
covector `k_x=D Phi(f_tz_0)^(-T)k_z` stays bounded above and away from zero.
Because the carrier is axisymmetric, shifting `beta` is a rigid rotation.
After that physical rotation identifies the endpoint frame, all coefficients
on the central torus return after the meridional period

    T_1(I_0)=2pi/Omega_1(I_0).                                    (11)

Thus (5) constructs an actual repeated moving-fibre problem even when the
particle orbit itself is quasi-periodic.

The six-dimensional map `(F_t,F_t^-T)` preserves
`sum_i dk_i wedge dx_i`.  This is the canonical ray symplectic form.  It is
neither the KKS form nor a physical-energy metric.

## 3. The full-pressure amplitude and the exposing correction

Along the physical trajectory let `L=grad u_*`, `varpi_*=curl u_*`, and retain
the complete Kelvin--Leray system

    k_dot=-L^T k,
    A_dot=-L A+2k[k dot L A]/|k|^2,
    k dot A=0.                                                     (12)

The coefficient two is the whole-space pressure projection.  It is required
for `d(k dot A)/dt=0`; deleting it leaves `-2 k dot L A` in general.

An initially tempting shortcut was to put `c=k cross A` and treat `c` as a
passive Cauchy vector.  The exposing calculation in the verifier instead gives

    c_dot=L c+mu A,
    mu=varpi_* dot k=constant,
    A=c cross k/|k|^2.                                            (13)

The constancy of `mu` follows from
`varpi_dot=L varpi` and `k_dot=-L^T k`.  Equation (13) is exact.  The omitted
term is not small on the dynamically accessible sector: for an oscillatory
divergence-free displacement with leading amplitude `d`,

    delta varpi_pr
       =i k cross(d cross varpi_*)
       =i[(k dot varpi_*)d-(k dot d)varpi_*].                      (14)

With `k dot d=0`, realizing an arbitrary transverse vorticity amplitude
requires `mu!=0`, exactly the case in which the second term of (13) survives.
Therefore (5) freezes the covector but does **not** by itself establish
`c(t)=Dvarphi_t c(0)` or velocity growth.  The three early verifier runs that
encoded that false shortcut are retained as failed receipts; the repaired
oracle detects the missing term directly.

This is a route-scoped refutation with a useful continuation, not a stability
claim.  Substituting (10) into (12), using the source map (1), gives a fully
specified periodic two-dimensional polarization equation on the resonant
torus.

## 4. Exact returned monodromy and its one scalar decision

Let `E(t):R^2->k(t)^perp` be a smooth physical orthonormal frame transported
through one meridional circuit and returned by the rigid axial rotation.  In
that frame, (12) becomes

    y_dot=C_res(t;I_0,n,m)y,       C_res(t+T_1)=C_res(t),           (15)

where every entry is the explicit analytic expression obtained from
`D Phi`, `Omega_1`, `Omega_2`, and their first derivatives.  In particular,
(15) retains the physical metric connection and the full pressure row; it is
not `F_t` acting on a velocity vector.

For two solutions of (12), direct differentiation gives the invariant

    [(A_1 cross A_2) dot k](t)
      =[(A_1 cross A_2) dot k](0).                                 (16)

Since the physical covector and frame return after (11), the normalized
one-circuit map

    M_res(I_0;n,m)
      =Pexp integral_0^T1 C_res(t;I_0,n,m)dt                       (17)

lies in `SL(2,R)`.  This transported polarization area is not identified with
the globally normalized KKS form.  The exact repeated ray object is
`M_res^j`, because the base shift in `beta` is removed only by the genuine
rigid rotation symmetry.

Consequently the complete central-fibre alternative is decided by

    Delta_res=[tr M_res]^2-4.                                     (18)

- `Delta_res>0` gives a real reciprocal pair and exponential ray growth.
- `Delta_res=0` with `M_res!=+/-I` gives parabolic polynomial growth.
- `Delta_res<0` gives bounded powers on this central fibre, although it says
  nothing by itself about nearby profiles or other covectors.
- `M_res=+/-I` is the exceptional bounded return.

The source small-shell limit does constrain but does not decide (18).  Since
`m/n=1/D(I_0)=O(I_0^(3/2))`, the resonant covector approaches the
axisymmetric meridional covector as `I_0->0`.  Coefficient continuity gives

    M_res=M_0+o(1),
    M_0=[[1,0],[-sqrt(2)pi,1]],                                   (19)

in the returned frame of 0019.  Thus the off-diagonal transient survives and
`tr M_res->2`, but the sign of the first nonzero correction to (18) is not
fixed by the published frequency expansions.  Treating closeness to a Jordan
matrix as proof that its powers grow would repeat the error this attempt was
designed to remove.

Equation (18), with (15)--(17), is the precise remaining same-carrier
calculation.  It requires either a higher source expansion of `Phi` and the
ordered amplitude integral, or an independently licensed interval enclosure
of that exact integral.  It is a route gap, not evidence of physical
impossibility.

## 5. An exactly bounded invariant-phase subroute

There is one nontrivial sector in which the amplitude can be solved without
evaluating (18).  Take the advected phase `theta=P`, so

    k=grad P.                                                      (20)

Steady Euler and Gavrilov localizability give

    L u_*=-grad P=-k,       u_* dot k=0.                           (21)

Substitution in (12) shows that `A_1=u_*` is an exact amplitude solution.
Moreover, steady Euler and `u_* dot grad P=0` imply
`u_* dot grad |u_*|^2=0`.  Baldi proves that irrational trajectories are
dense on their pressure tori for a full-measure set of actions.  Continuity
then gives

    |u_*|^2=S(P)                                                   (22)

throughout the regular action region.

Define the normalized transverse complement

    A_2=(k cross u_*)/(|k|^2 |u_*|^2).                             (23)

The invariant (16) makes the defect
`A_2dot-C_A A_2` proportional to `u_*`.  Its coefficient is

    gamma=(1/(2|u_*|^2)) A_2 dot grad |u_*|^2=0                   (24)

by (22), so `A_2` also solves (12).  Both fields are bounded above and away
from degeneracy on a compact regular annulus.  The normal invariant-phase
return is therefore the identity in the moving basis `(A_1,A_2)` and is
uniformly bounded at principal order.

This exact positive result is narrow but physical: pressure-normal high-
frequency perturbations do not inherit the 0019 lift-up.  It does not bound
angle phases, subprincipal profile coupling, or the global semigroup.

## 6. Dynamically accessible packet injection and normalization

Suppose a chosen solution of (15) has a strict gain at a finite target time.
Let `chi_delta(I)` be supported in a band around its torus and let
`theta_0` be the corresponding integral phase.  On a subband where
`mu=k dot varpi_*` is bounded away from zero, any desired transverse leading
vorticity `c_0=k cross A_0` is realized by

    d_0=c_0/(i mu),       k dot d_0=0.                             (25)

Choose a compact vector-potential amplitude `B_delta` with
`i k cross B_delta=d_0` on the middle band, and set

    w_{N,delta}=N^-2 Re curl[chi_delta B_delta exp(iN theta_0)],
    q_{N,delta}=-ad^*_{w_{N,delta}}m_*.
                                                                        (26)

Then `w` is exactly divergence free and compactly supported, `q` is an exact
smooth coadjoint-orbit tangent, and (14) supplies the prescribed principal
polarization.  Its raw Biot--Savart velocity norm is comparable to `N^-1`.
Carrying the accepted 0019 correction, the only unit data used below are

    qhat_{N,delta}=q_{N,delta}/||delta u[q_{N,delta}]||_X.          (27)

The normalization is scalar, hence preserves accessibility and graph-domain
membership.  All WKB and Duhamel errors below are relative to
`||qhat(0)||_X`, never absolute errors compared with an order-`N^-1` raw
velocity.

If `mu=0`, (14) has smaller rank and (25) is unavailable for a general
polarization.  That is recorded as an accessibility restriction rather than
silently filled with an arbitrary `chi` perturbation.

## 7. Exact fixed-time diagonal transfer theorem

The growing-time estimate frozen in the README is sufficient but not
necessary for an operator-norm conclusion.  The following quantifier is
strictly weaker and closes the continuum bridge whenever the exact ray
cocycle itself has unbounded gains.

Let `t_j` be any finite time sequence and suppose the exact solution of (15)
on one fixed carrier has velocity gains `G_j->infinity`.  For each fixed `j`:

1. continuity over the compact angle torus permits a width `delta_j>0` for
   which the complete leading profile gain is at least `G_j/2`;
2. the global two-term pseudodifferential construction on `[0,t_j]` gives a
   finite constant `C(j,delta_j)` and the **relative** bound

       sup_{0<=t<=t_j}
       ||S_*(t)qhat_{N,delta_j}-v^WKB_{j,N}(t)||_X
          <= C(j,delta_j)N^-1 ||qhat_{N,delta_j}(0)||_X;            (28)

3. after `t_j` and `delta_j` are fixed, choose `N_j` so that (28), the initial
   symbol error, and symmetry modulation consume less than one quarter of the
   leading gain.

The resulting exact dynamically accessible data obey

    ||q_j||_X=1,        ||S_*(t_j)q_j||_X >= G_j/4.                 (29)

The quantifiers are

    for every finite j, there exist delta_j and then N_j.          (30)

No estimate uniform for `t_N->infinity` is used or claimed; the constants in
(28) may grow arbitrarily with `j`.  This is enough because an operator norm
at each time may be tested on a different unit datum.  Thus either the
hyperbolic or nontrivial parabolic alternative in (18) implies

    sup_t ||S_*(t)||_(X_DA->X_DA)=infinity.                         (31)

Banach--Steinhaus then gives some `q in X_DA` with an unbounded time sequence,
but does not make that datum smooth, localized, or graph-domain.  An explicit
single smooth packet still requires an almost-orthogonal superposition and
cross-time estimates.

## 8. Whole-space pressure, collar, exterior, and quotient ledger

The estimate (28) is a whole-space statement, not an annulus-wall model.

1. Each action band is invariant under the base flow and, for fixed `j`, is a
   positive distance from both flat collars and the support boundary.
2. The Leray projector is the global order-zero multiplier.  Its commutator
   with a smooth cutoff is order minus one, contributing to
   `C(j,delta_j)/N`; no local pressure closure is substituted.
3. Biot--Savart is order minus one.  Its local symbol reconstructs physical
   velocity from (12)--(14), and explains both the raw `N^-1` scale and the
   mandatory normalization (27).
4. Between separated annulus and exterior cutoffs the Green kernel is smooth.
   Since the integral phase has nonzero covector, repeated integration by
   parts gives `C_{M,j}N^-M` for arbitrary fixed `M`.  The same cancellation
   controls the oscillatory curl source's exterior multipoles.  The
   perturbation velocity itself is not compactly supported.
5. Flatness removes a base-coefficient jump, not the nonlocal perturbation
   tail.  That tail and every `[P,chi]` block are retained in (28).
6. The entire transported profile is used.  Choosing `delta_j` first controls
   its exact covector variation; choosing `N_j` last controls the PDE
   remainder.  There is no infinite-circuit summability assumption.

Euclidean symmetry tangents form a finite-dimensional smooth subspace.
Fixed-time oscillatory pairings are `O(N_j^-M)` for every fixed `M`; taking
quotient distance, rather than declaring an orthogonal complement invariant,
preserves (29) after increasing `N_j`.

## 9. Orbit action and separation of structures

On the smooth stabilizer quotient the KKS form descends because stabilizer
generators are null representatives of the orbit tangent.  With
`i_X Omega=dH`, the quadratic action is

    S_2=int[-Omega(q,q_dot)/2-Q(q,q)/2]dt.                         (32)

Euler equivariance makes the smooth linearized flow symplectic and preserves
the indefinite quadratic Hamiltonian `Q`.  No strong extension of KKS to the
kinetic completion is asserted.  Conditional growth in (31) concerns the
positive physical kinetic norm and is compatible with indefinite-`Q`
conservation.  The KKS form (32), ray cotangent form from (6), polarization
area (16), and physical kinetic metric are four distinct objects.

## 10. Route verdicts and the concrete missing mechanism

- **Normal invariant-phase route:** established bounded at principal order by
  (20)--(24).
- **Naive Cauchy/deformation route:** refuted by the exact extra term in (13).
  Derivative resonance alone does not produce velocity growth.
- **Resonant returned-fibre route:** reduced exactly to the physical periodic
  `SL(2,R)` monodromy (15)--(18), including metric connection and pressure.
  The sign of `Delta_res` is not fixed by the available source expansion.
- **PDE transfer route:** established conditionally in the sharp quantifier
  form (28)--(31).  A growing ray sequence needs no single WKB expansion
  uniform on growing times.
- **Axisymmetric angle-phase route:** the moving covector prevents powers of
  the 0019 Jordan matrix; a global upper bound remains open.

The strongest exact verdict is therefore not an unbounded-semigroup theorem
for this carrier.  It is the exact construction of the returned physical
monodromy, the exposure of the term that invalidates the tempting positive
shortcut, an exactly bounded invariant-phase subfamily, and a complete
fixed-time diagonal bridge that will promote either a hyperbolic or parabolic
monodromy finding directly to exact global accessible packets.

For LP2, the next physical mechanism is concrete: compute or control
`Delta_res` for the source-defined carrier, then, if it is elliptic, control
the nearby profile/subprincipal nonlocal cocycle rather than declaring
stability from one fibre.  A stable particle carrier must eliminate these
resonances or supply a restoring sector that changes the physical monodromy.
This remains classical linear Euler work.  Nothing here supplies `hbar`, Born
probabilities, exchange statistics, relativistic dispersion, an electron
scale, or a neutrino mechanism.

# P253/0032: source-defined resonant Kelvin cocycle discriminant

This README preregisters a new fixed-carrier analytic attempt.  The parent
objective remains a robust localized Euler carrier capable of supplying the
electron mechanism (and subsequently the neutrino mechanism); 0032 is not an
instability target and cannot close that parent by finding growth.  Its
positive purpose is to decide the exact physical propagation metric and
returned cocycle on the same source-defined Gavrilov carrier, so that a robust
carrier is assessed with the full pressure and accessible perturbation space.

## Ownership and activation boundary

Root owns central registration, activation, source/API/test changes, memory,
commits, and promotion.  This worker owns only append-only files in
`attempts/0032`.  Before root registers 0032 and writes
`activation-schema.exit` containing exactly `0`, this README is the sole
artifact: no source body, derivation, verifier, numerical design, or result is
opened or produced.  After activation, every inspected input and generated
receipt is pinned by SHA-256 in this directory.  Attempts 0019 and 0025 are
suppliers only; no conclusion from them is silently promoted.

Accepted base: release `v0.183.0`.  Parent campaign, LP2, electron-first and
neutrino obligations remain active.  A route verdict here is scoped to the
fixed carrier, cocycle, perturbation sector, metric, and domain stated below.

## Frozen object and exact positive target

Fix one smooth compact Gavrilov carrier `u_*` with a closed regular pressure
plateau and no artificial annulus wall.  On the plateau use Baldi's exact
action-angle chart

    x=Phi(sigma,beta,I),
    sigma_dot=Omega_1(I),  beta_dot=Omega_2(I),  I_dot=0.          (1)

Choose a compact interval `B=[I_-,I_+]` in that plateau.  The source expansions
are frozen as

    K(I)=I+(1065/1024)I^3+O(I^4),
    Omega_1=K',
    Omega_2=sqrt(I)(1+7I/4+O(I^2))Omega_1.                        (2)

Write `a=Omega_1'`, `b=Omega_2'`.  Since

    b/a=(256/3195)I^(-3/2)(1+O(I)),

select one exact `I_0 in B` and coprime positive integers `(n,m)` satisfying

    n a(I_0)-m b(I_0)=0.                                           (3)

The source-defined integral phase and exact chart cotangent lift are

    theta_0=n sigma-m beta,
    F_t=[[1,0,t a],[0,1,t b],[0,0,1]],
    k_z(t)=F_t^(-T)(n,-m,0)
          =(n,-m,-t[n a-m b]).                                   (4)

Thus the covector returns on `I=I_0`; the beta displacement is identified only
by the genuine axial rigid rotation.  The physical covector is
`k_x=D Phi^{-T}k_z` and is never replaced by a chart Euclidean covector.

The target is the exact physical polarization cocycle on this returned
torus, including metric connection, pressure, and accessible-sector descent:

    k_dot=-L^T k,
    A_dot=-L A+2 k(k dot L A)/|k|^2,
    k dot A=0,                                                     (5)

where `L=grad u_*` on the actual trajectory.  Let `E(t)` be a physical
orthonormal frame of `k(t)^perp`, returned after one meridional period

    T_1(I_0)=2pi/Omega_1(I_0)

by the axial symmetry identification.  Define the exact ordered-integral
representation

    C_res(t)=E(t)^*[ -L+2 k \otimes(k^T L)/|k|^2 ]E(t)-E(t)^*dot E(t),
    M_res=Pexp integral_0^T1 C_res(t) dt.                         (6)

Every symbol in (6) is source-defined; the frame connection is part of the
operator.  The invariant physical polarization area is

    [(A_1 cross A_2) dot k](t)=[(A_1 cross A_2) dot k](0),           (7)

so the returned map in the returned physical frame is in `SL(2,R)`.

The primary success contract is an analytic determination of

    Delta_res=(tr M_res)^2-4,                                     (8)

with the following physically distinct alternatives retained:

| Exact outcome | Meaning licensed here |
|---|---|
| `Delta_res>0` | hyperbolic reciprocal multipliers and a genuine ray-growth route |
| `Delta_res=0`, `M_res != +I,-I` | nontrivial Jordan/parabolic return and polynomial ray growth |
| `M_res=+I` or `-I` | scalar identity return; no Jordan conclusion is allowed |
| `Delta_res<0` | elliptic central return; seek a positive periodic metric and profile control |

The symbols `+/-I` are tested separately from the discriminant.  A repeated
unit eigenvalue is not called a Jordan matrix without a rank test.

## Analytic trace route (first and preferred)

Before any numerical design, derive the trace from pressure-localizability,
conserved quantities, and reducible-coordinate algebra.

1. Steady Euler and localizability give

       L u_*=-grad P,       u_* dot grad P=0,
       varpi_dot=L varpi,   k_dot=-L^T k,
       mu:=varpi dot k=constant.                                  (9)

   For `c=k cross A`, the exact full-pressure identity is

       c_dot=L c+mu A,
       A=c cross k/|k|^2.                                         (10)

   The `mu A` term is retained: dynamic accessibility of a general transverse
   polarization uses

       delta varpi_pr=i[(k dot varpi)d-(k dot d)varpi],             (11)

   and cannot erase the factor `mu` by relabelling.

2. The normal invariant phase `theta=P` is an exact check on reducibility.
   Steady Bernoulli plus localizability gives `u dot grad|u|^2=0`; Baldi's
   dense irrational tori and continuity imply `|u_*|^2=S(P)`.  Therefore

       A_1=u_*,
       A_2=(grad P cross u_*)/(|grad P|^2 |u_*|^2)                 (12)

   are two exact bounded amplitudes.  This sector supplies a positive bounded
   periodic metric candidate, but does not decide angle phases.

3. For the resonant angle phase, pull (10) and (5) into the action chart,
   retain the physical metric `g=D Phi^T D Phi`, and reduce the returned
   equation to a scalar second-order equation or an explicit ordered integral.
   Use (7) to prove `det M_res=1`, then calculate the trace and the rank of
   `M_res-I` and `M_res+I`.  The small-shell limit may be compared with

       M_0=[[1,0],[-sqrt(2)pi,1]],                                  (13)

   but closeness to (13) cannot decide the sign of (8).

4. In parallel, attempt a positive periodic energy metric `H(t)>0` solving

       H_dot+C_res^T H+H C_res=0,
       H(t+T_1)=H(t).                                              (14)

   Existence of (14) proves bounded central propagation and distinguishes an
   elliptic return from a merely sampled unit multiplier.  It is not inferred
   from `det M_res=1`, and it is not identified with KKS.

## Competing representations and route ledger

The frozen candidates are:

- **A: reducible-coordinate trace.**  Use (9)--(13), pressure-normal neutral
  amplitudes, conserved `mu` and polarization area, and exact quadratures to
  determine (8) symbolically.
- **B: invariant-metric construction.**  Solve (14) on the physical energy
  metric, then extend the metric to a profile band and control connection,
  Leray, collar, and exterior terms.
- **C: direct ordered-integral enclosure.**  If A and B leave one scalar
  integral or discriminant sign, freeze that remainder for a separately
  activated rigorous/convergent computation.  C cannot select a carrier or
  reinterpret an elliptic return as instability.

Selection is by exact pressure retention, physical metric and accessibility,
global-domain closure, symmetry descent, and predictive reach.  A candidate
that only computes a boxed eigenvalue or a fixed Jordan power is not a route
to the physical object.

## Profile, global pressure, and accessible-domain bridge

The exact perturbation domain is the kinetic-norm closure of smooth coadjoint
orbit tangents, modulo Euclidean symmetry distance.  For a finite target time
`t_j=jT_1`, a smooth action-band profile is first narrowed to `delta_j`; only
then is its oscillation frequency `N_j` chosen.  The required continuum bridge
has the quantified form

    ||E_{N_j}(t)a||_X <= C(j,delta_j) N_j^(-1)||a||_X,              (15)

relative to the normalized physical velocity, where the raw
Biot--Savart-normalized generator has velocity scale `N_j^(-1)`.  The whole
space Leray projector, pressure tail, `[P,chi]` collar blocks, separated
annulus/exterior smoothing, flat support boundary, and symmetry leakage all
remain in `C(j,delta_j)`.  This is `for every j exists delta_j exists N_j`,
not a single expansion uniform for `t_N -> infinity`.

If the exact cocycle gains `G_j -> infinity`, (15) transfers that gain to
exact normalized dynamically accessible packets and hence to an unbounded
semigroup norm.  If (14) and the profile ledger close instead, the result is a
positive propagation metric for the declared sector.  Neither result alone
creates nonlinear particle stability, a restoring frequency, quantum
statistics, or an electron/neutrino identity.

## Numerical boundary (no production run at preregistration)

No production numerical verifier is licensed by this README.  If analytic
routes A and B leave only the scalar discriminant (8), route C must first load
`small-ratio-numerics/SKILL.md` and freeze:

- the dimensionless action interval and the small parameter controlling the
  near-Jordan limit;
- the exact ODE (15)/(6), initial frame, period, and quadrature/IVP method;
- interval or observed-order error bounds for `M_res`, `tr M_res`, and both
  singular values, with a residual check for the ordered-integral solution;
- parameter-box refinement in `I_0`, `(n,m)`, frame gauge, and period;
- separate tests for `M_res=+/-I`, nontrivial Jordan rank, and elliptic versus
  hyperbolic sign, including perturbation/conditioning sensitivity; and
- the exact bridge from an enclosed sign to either (14) or the fixed-time
  accessible-profile theorem (15).

Until that design is activated, no sampled trace, multiplier, or approximate
metric is evidence.  A numerically unresolved sign earns
`evidence_scope: NUMERICALLY_UNRESOLVED`, not a physical no-go.

## Source and dependency inventory

- Baldi, *The dynamics of Gavrilov's compactly supported Euler flows*,
  [`arXiv:2302.02982v1`](https://arxiv.org/abs/2302.02982v1), cached outside
  the repository at `/tmp/primary-source-cache/P253-0025/2302.02982v1.pdf`,
  SHA-256 `4fe82be0db15a7f117017e505e6682b3438ed14ebcf56225bd7e74bc759136f3`.
  Theorem 1.1, pp. 4--6; exact action-angle construction
  (3.52)--(3.68), pp. 16--18; coefficients (4.48)--(4.61), pp. 26--28.
- P253/0019, fixed-time global Gavrilov carrier and one-circuit pressure
  return, with accepted normalization repair.  It supplies only finite-time
  transient scope.
- P253/0025, exact moving-fibre reduction and normalization convention.  Its
  corrected exposing identity (10) is imported as a route input, not as a
  proof of the discriminant sign.
- Whole-space Euler linearization, Leray projection, KKS action convention
  `i_X Omega=dH`, coarea, smooth-flow dependence, and standard ordered-ODE
  facts are declared mathematical imports whose hypotheses must be checked on
  the fixed carrier.

Accepted abstract APIs may validate a finite matrix or ordered-integral
implication after the physical object is derived; they are not Euler
realizations and do not supply the pressure, accessible profile, or metric.

## Verdict and continuation contract

Each route receives exactly one scoped verdict: established, refuted with the
mechanism named, or blocked with the missing trace/metric construction named.
If A fails algebraically, continue through B and C; if C is numerically
unresolved, preserve the exact ordered-integral object and continue with a
different carrier or profile representation.  No finite discriminant failure
closes the parent robust-carrier objective, and no positive metric on one
phase sector is promoted to full particle stability without the global
profile/domain replay.

# P253/0039: global accessible WKB lift of the resonant Gavrilov cocycle

Status: **README-only preregistration awaiting central registration and schema
activation.** Until `0039/activation-schema.exit` exists and contains exactly
`0`, this README is the sole 0039 artifact. No source body, derivation,
verifier, numerical run, API, test, result, or receipt is produced. Root owns
central files, commits, activation, review, and promotion. This worker owns
only eventual append-only artifacts in `attempts/0039`.

## Parent objective and inherited result

The parent P253 objective remains an actual robust localized Euler mechanism
for the electron first and the neutrino subsequently. Attempt 0039 addresses
one fixed LP2 dependency only: whether the source-derived hyperbolic Kelvin
cocycle from 0032 lifts to the exact whole-space linearized Euler semigroup on
finite-energy dynamically accessible perturbations. Neither success nor
failure of this fixed-carrier route decides the parent particle objective.

Fix the same smooth compact Gavrilov carrier `u_*` and one sufficiently small
derivative-resonant regular pressure shell `I=I_*` wholly inside a flat cutoff
plateau. Attempt 0032 supplies, with `epsilon=sqrt(2I_*)`,

    tr M_*=2+22*pi^2*I_*+O(I_*^(3/2)),
    det M_*=1,
    Delta_*=88*pi^2*I_*+O(I_*^(3/2))>0.                    (1)

Thus its returned physical polarization map has reciprocal real multipliers
`lambda_+>1>lambda_->0`. This is an input at active-proposal authority, not an
accepted global instability theorem. The exact meridional circuit time is

    T_*=2*pi/Omega_1(I_*),   gamma_*=(log lambda_+)/T_*>0. (2)

The resonant integral phase is retained with coprime `(n,m)` satisfying
`n Omega_1'(I_*)-m Omega_2'(I_*)=0`; the physical covector is always the
cotangent lift `D Phi^{-T}(n,-m,0)`.

| Frozen field | 0039 value |
| --- | --- |
| `positive_intent` | Lift the corrected hyperbolic same-carrier cocycle to normalized finite-energy dynamically accessible packets and an exact Euler-semigroup lower bound for every circuit count |
| `requires` | 0032 corrected source recurrence and hyperbolic return; 0019 whole-space domain/normalization; 0025 moving-fibre, pressure, collar, and finite-time quantifier ledger |
| `pass_licenses` | A linear semigroup/essential-norm lower bound on the named fixed carrier, energy space, accessible closure, quotient, and circuit times |
| `does_not_license` | Nonlinear Euler instability; orbital breakup; LP2 stable-particle embedding; quantum statistics; Born rule; relativity; electron or neutrino identity |
| `maximum_verdict` | An exact all-circuit semigroup lower bound, and optionally a generator approximate-spectrum/essential-spectrum corollary if its domain bridge closes |
| `failure_scope` | The selected carrier, resonant shell, packet class, norm, and one of the two frozen microlocal constructions only |
| `unlocks` | Nonlinear modulation/instability analysis if growth is transferred; otherwise a different stable-carrier or invariant-sector construction for LP2 |

## Exact Euler space, domain, pressure, and accessibility

Let `P_L` be the whole-space Leray projector and define the closed linearized
Euler generator on the kinetic-energy space `X=L^2_sigma(R^3)` by

    Gv=-P_L[(u_* dot grad)v+(v dot grad)u_*],
    S(t)=exp(tG).                                             (3)

Use the skew-adjoint transport plus bounded-shear realization from 0019; do
not replace the maximal transport domain by an unproved global `H^1` domain.
The dynamically accessible subspace and graph domain are

    X_DA=closure_X{-ad_xi^* u_* : xi smooth, compact, div xi=0},
    D_DA={v in D(G) intersect X_DA : Gv in X_DA}.             (4)

Euclidean symmetry tangents are removed by quotient distance. Invariance of
`X_DA` and the quotient under `S(t)` is proved from the coadjoint evolution,
not assumed from the ray system.

For a compact packet amplitude `a` and an advected integral angle `phi`, take
an integer `N>=1`, put `h=N^-1`, and construct the injection at vorticity
level,

    delta omega_N=curl(xi_N cross omega_*),
    v_N=BiotSavart(delta omega_N)=P_L V_N,
    xi_N=N^-1 Re[a exp(i N phi)].                          (5)

The integer frequency makes the oscillation globally single-valued on the
action-angle torus.  Include any needed lower-order divergence correction and
normalize the actual velocity, not the formal amplitude:

    q_N=v_N/||v_N||_X,   ||q_N||_X=1,   q_N in D_DA.          (6)

Record separately the graph size
`||q_N||_D=||q_N||_X+||Gq_N||_X=O(N)`; no bounded graph-norm claim is
inferred from energy normalization. The WKB remainder is always relative to
`||v_N(0)||_X` before (6).

The proof retains the complete pressure reconstruction: the order-zero
principal Leray term in the Kelvin amplitude, subprincipal Leray symbols,
the nonlocal Hodge/Biot--Savart tail, commutators with packet cutoffs, the
flat carrier-support collar, and the exterior harmonic velocity. No compact
vorticity support is identified with compact velocity support. Packet support
begins in a tubular action-angle chart at positive distance `d_0` from every
chart and cutoff boundary.

## Frozen all-circuit statement and quantifiers

Let `e_+` be a returned unit polarization eigenvector for `lambda_+`. For
every integer `j>=1`, define `t_j=jT_*`. The primary theorem target is

    ||S(t_j)||_{X_DA/sym -> X_DA/sym}
       >= (1-r_j) lambda_+^j,   0<=r_j<=1/2,                 (7)

and consequently

    limsup_{j->infinity} t_j^-1 log ||S(t_j)|| >= gamma_*.  (8)

The quantifiers are

    for every j, choose delta_j>0, then choose an integer
    N_j>N_0(j,delta_j), and construct q_{N_j,j}
    with ||q_{N_j,j}||_X=1.                                 (9)

This proves one fixed semigroup has arbitrarily large operator norms; it does
not assert one frequency expansion uniform as `j->infinity`, nor one fixed
initial datum with an unbounded orbit. A stronger single-packet or uniform-time
result is recorded only if a derived remainder closes that stronger order of
quantifiers.

The support width is chosen so the exact profile cocycle and base return vary
by at most `r_j lambda_+^j/4` on the tube. Only afterward is `N_j` chosen so
the finite-time microlocal remainder is at most the same amount. Every
constant may depend explicitly on `j` and `delta_j`; it may not depend on
`N_j`. This is sufficient for (7)--(8) because the semigroup norm permits a
different normalized input at each time.

## Route A: direct Gaussian beam and Egorov propagation

Construct a divergence-free/coadjoint Gaussian beam centered on the exact
closed resonant bicharacteristic `z_*(t)=(x_*(t),k_*(t))`, with transverse
width `delta_j` and integer semiclassical scale `N^-1`. Prove a finite-time Euler Egorov
formula on `[0,t_j]`,

    S(t) I_N(a)=I_{N,t}(C(t)a)+R_N(t)a,                      (10)

where `C(t)` is the full pressure-retaining Kelvin cocycle and `I_{h,t}` uses
the advected phase, physical energy metric, transported density, and Hodge
normalization. The required estimate is

    sup_{0<=t<=t_j} ||R_N(t)a||_X
      <= N^-1 C_j delta_j^-s_j exp(B_j t_j)||a||_{Y_j}.      (11)

Derive `s_j`, `B_j`, and `Y_j` from symbol derivatives; do not insert them as
fit parameters. Propagate support through every chart, show the central tube
stays inside the flat plateau, estimate collar/exterior leakage, and prove
that the final physical frame and density return. Apply (10)--(11) to `e_+`
and then normalize as in (6) to obtain (7).

The route may alternatively invoke a published Euler essential-norm/Egorov
theorem only after mapping its manifold, Sobolev index, pressure projector,
flow regularity, cocycle convention, compactness remainder, boundary, and
accessible-subspace hypotheses to (3)--(6). An abstract principal-symbol API
does not supply that map.

## Route B: packet pseudomode and essential-spectrum construction

Independently construct a Weyl packet sequence localized on the same returned
orbit. With Floquet exponent `gamma_*`, freeze candidate generator parameters

    z_l=gamma_*+i(2*pi*l+theta_F)/T_*,   l in Z,              (12)

where `theta_F` includes the returned phase/Maslov and polarization sign and
is derived rather than chosen. Seek integer-frequency `q_N in D_DA`,
`||q_N||_X=1`, weakly
convergent to zero modulo symmetries, such that

    ||(G-z_l)q_N||_X -> 0.                                   (13)

The longitudinal cutoff is closed around the genuine periodic orbit; the
transverse localization, Hodge tail, Leray commutators, and graph-domain
membership are estimated explicitly. If (13) closes, determine whether it
proves approximate point spectrum, essential spectrum under the repository's
declared Fredholm convention, or only a pseudospectral lower bound. A small
residual for a nonnormal generator is not by itself a spectral theorem.

The time-domain alternative is an orthogonal packet family whose propagated
images retain `lambda_+^j` gain modulo compact remainders. Prove weak-nullness
and asymptotic orthogonality in the physical kinetic metric, then obtain

    ||S(t_j)||_ess >= lambda_+^j                              (14)

or the explicitly derived constant multiple. Route B competes with Route A
by removing profile-uniformity burdens; it does not replace the exact
dynamical-accessibility and pressure bridges.

## Route comparison and continuation ladder

Route A is preferred if the finite-time Euler Fourier-integral propagation
can be proved directly with transparent support control. Route B is preferred
if compact-remainder/weak-packet arguments yield the sharper essential norm
or generator statement. Compare them by:

1. exact preservation of coadjoint accessibility and quotient descent;
2. retention of Leray/Hodge pressure and physical energy normalization;
3. explicit packet support and chart/collar/exterior control;
4. honest time/frequency quantifiers and graph-domain closure;
5. strength of the earned semigroup statement with fewer imported hypotheses.

An implementation fault is repaired within its route. A chart or phase
failure triggers the alternate packet representation. Failure of an imported
Egorov theorem's hypotheses activates a direct parametrix. Failure of a Weyl
sequence to control the nonlocal tail leaves Route A active. Neither route's
failure is a fixed-carrier no-go until the other route and these representation
repairs have verdicts.

## Nonlinear and particle separation

Equations (7)--(8) are linear operator-norm growth. They do not imply nonlinear
Lyapunov instability without a differentiable Euler flow map in the stated
topology and a controlled nonlinear remainder on a growth-compatible time
scale. They also do not imply carrier destruction: nonnormal linear growth,
an unstable Floquet packet, orbital instability, and failure of an LP2
stable-particle embedding are distinct propositions. Any nonlinear route must
be separately preregistered with its topology, solution class, modulation,
and time scale. No quantum, exchange, action-normalization, relativity,
electron, or neutrino claim is licensed here.

## Analytic and numerical boundary

0039 is an exact/microlocal construction. No production numerics are planned
or authorized. Symbolic scripts after activation may verify source recurrences,
signs, normalization, cocycle/frame covariance, and exact model remainders;
they may not substitute sampled packets or matrices for (10)--(14). If an
irreducible numerical remainder unexpectedly appears, it requires an
append-only contract expansion and the applicable numerical skill before any
design is frozen.

## Source and dependency inventory

- P253/0032 corrected derivation and `verify_c2_source_recurrence.py`: exact
  same-carrier hyperbolic cocycle (1), source recurrence, physical frame, and
  flat-plateau period scaling. It supplies no global packet lift.
- P253/0019: whole-space linearized Euler realization, dynamically accessible
  injection, finite-energy normalization repair, regular Gavrilov annulus, and
  fixed-time WKB scope.
- P253/0025: repeated moving-fibre object, physical energy metric, profile
  return, global Leray/Hodge/collar ledger, and the distinction between fixed
  circuit count and a uniform growing-time expansion.
- Gavrilov, *A steady Euler flow with compact support*, arXiv:1810.08020v1,
  cached at `/tmp/primary-source-cache/P253-0010/1810.08020.pdf`, SHA-256
  `fcaca85faa77e3876b11d16718037169fc026112dfd3e4248e03e963a0ebc3c9`;
  theorem p. 1, defining ODE/PDE pp. 2--4, field and cutoff pp. 5--6.
- Baldi, *The dynamics of Gavrilov's compactly supported Euler flows*,
  arXiv:2302.02982v1, cached at
  `/tmp/primary-source-cache/P253-0025/2302.02982v1.pdf`, SHA-256
  `4fe82be0db15a7f117017e505e6682b3438ed14ebcf56225bd7e74bc759136f3`;
  action-angle map (3.52)--(3.68), pp. 16--18, and exact source recurrences
  (4.29)--(4.61), pp. 24--28.

Potential published Egorov, essential-spectrum, or Euler short-wave theorems
are inventory candidates only until their exact statements and hypotheses are
opened after activation and pinned with URL/version/hash/theorem location.

## Route verdict and completion contract

Each route terminates with exactly one route-scoped verdict: established,
refuted with its mechanism, or blocked with its missing construction. The
attempt earns its positive result only when (7) or a stronger correctly bridged
version of (14) is proved for the exact semigroup (3), accessible domain (4),
normalization (6), and every circuit count. A central-ray matrix, formal beam,
principal symbol, fixed finite time, or generator pseudomode without the
pressure/domain bridge does not meet that contract. The next dependency after
success is a separately licensed nonlinear modulation/instability analysis;
after a route failure the alternate construction remains active.

# P253/0019 preregistration: global streamline return and Euler semigroup control

Status: **README-only preregistration awaiting central registration and schema
activation.**  Until `0019/activation-schema.exit` exists and contains `0`, this
attempt opens no new source body, derives no new carrier coefficient, executes no
symbolic or numerical verifier, and writes no file other than this README.  It owns
only eventual `0019` artifacts; source modules, tests, governance, the central
proposal, and every other attempt remain outside its write surface.

## Parent obligation and exact positive deliverable

This is the failure-derived LP2 continuation from 0014, not a repetition of its
Hessian calculation.  Attempt 0014 supplies the exact axisymmetric coadjoint tangent,
mixed constrained Hessian signature, and the regular-core `J Hess(A)` expression for
the same finite-`j` Cao--Zhan and Gavrilov carriers.  Independent 0020 review fixes
two bounded convention/domain points that 0019 consumes: the KKS quadratic-action
sign is set by a canonical `dq wedge dp` test, and the global free-boundary
linearization retains `(delta J)dA` wherever `dA=0` is unavailable.  Its raw frozen
symbol excludes an `O(|k|)` exponential characteristic rate only.  In the physical
energy weighting, the retained `K` coupling can still produce an `O(1)` short-wave
rate.  Neither the raw principal symbol nor the mixed Hessian decides the global
initial-value problem.

The positive deliverable in 0019 is an actual global linear propagation theorem for
one of those same carriers.  It will:

1. construct a closed, densely defined linearized Euler generator on the physical
   dynamically accessible kinetic-energy space, with exact symmetry modulation and
   no artificial wall;
2. solve the characteristic/return problem on a proved regular invariant streamline
   annulus and compute its circuit monodromy, including the possible cancellation or
   accumulation of 0014's local Jordan shear;
3. retain the inverse-elliptic Biot--Savart operator `K` and control its cross-streamline,
   annulus-collar, and exterior-tail contributions by an exact Duhamel, resolvent, or
   modified-scattering estimate; and
4. prove a global semigroup/phase-mixing bound with a stated physical observable, or
   construct a genuine growing mode or wave packet of the full operator and identify
   its mechanism.

A positive no-growth/scattering result licenses an **axisymmetric linear LP2
persistence statement in the declared norm** for the named carrier.  A full growing
mode licenses linear instability of that carrier and activates a different LP2
candidate.  A local multiplier or high-frequency cocycle alone licenses only its
route-scoped necessary mechanism.  Neither outcome by itself supplies nonlinear or
three-dimensional stability, a restoring frequency, quantum mechanics, relativity,
or electron/neutrino identity.

The obligation record is frozen as follows:

| Field | 0019 value |
| --- | --- |
| `positive_intent` | Global same-carrier linear Euler propagation with the physical orbit, nonlocal pressure/Biot--Savart reconstruction, and a time-domain persistence or instability verdict |
| `requires` | 0010 finite-`j` carriers and conserved functional; 0005 Euler orbit/KKS action; corrected 0014 exact tangent, regular-core operator, and raw order-`|k|` characteristic; 0020 bounded review corrections |
| `pass_licenses` | Axisymmetric linear LP2 control for the named carrier and a justified next nonlinear modulation problem |
| `does_not_license` | Nonlinear Lyapunov stability; unrestricted 3-D stability; discrete particle spectrum; quantum probabilities/statistics; relativity |
| `maximum_verdict` | Exact or theorem-backed semigroup/scattering theorem for the continuum operator, or a full-operator growing solution in the declared accessible space |
| `failure_scope` | Only the named carrier, perturbation sector, annulus hypotheses, and operator representation actually closed in this attempt |
| `unlocks` | Nonlinear modulation if controlled; materially different stable carrier if unstable; nonaxisymmetric completion after either axisymmetric verdict |

## Carrier, domain, orbit space, and physical norm

The preferred analytic carrier is the compact Gavrilov pressure-shell field from
0010.  Its meridional streamfunction is

    psi_G(p)=-int^p omega_c(s) ds,

and its pressure contours are invariant under the actual meridional flow.  After
activation, the first construction must exhibit numbers `I_-<I_+` for the **same
fixed cutoff/carrier** such that

    A_G={x in Pi : I_- <= psi_G(x) <= I_+}                            (1)

is a compact regular annulus strictly inside the positive cutoff shell,
`grad psi_G!=0` on `A_G`, and every level `Gamma_I={psi_G=I}` is a smooth closed
meridional streamline.  The annulus is an analytic chart, not a boxed physical
domain: no boundary condition is imposed at `I=I_+/-`.  Its two collars, the rest of
the carrier, and all Green-kernel coupling across them remain in the global operator
and error ledger.

The Cao--Zhan ring is the second same-family route.  Its annulus

    A_e={I_- <= psi_e <= I_+}                                        (2)

must lie strictly inside the active topological-disc core, away from its free
boundary and the axis.  Theorem 1.3(iii)'s `partial_z psi_e<0` may be used only after
the already-audited primary receipt is rechecked.  The result may not upgrade an
interior annulus theorem to the whole free-boundary ring without controlling the
omitted collars and exterior.

For either carrier, let `u_*` be the actual steady or translating-frame velocity,
`m_*=u_*^flat`, and define

    T_DA(m_*)=closure{-ad^*_v m_* :
                       div v=0, v smooth, regular and decaying}.      (3)

The closure uses the physical perturbation kinetic norm

    ||q||_X^2=rho0 int_R3 |delta u[q]|^2 dx.                          (4)

In the axisymmetric variables of 0014 this is exactly

    ||(eta,chi)||_X^2
      =2 pi rho0 int_Pi [eta K eta+chi^2/r^2] dnu.                   (5)

Let `Z_ax=span{q_Z}` be the proved axial-translation tangent, with any additional
nonzero symmetry vector admitted only after it is derived.  Work on

    X_ax=closure_X(T_DA(m_*) intersect Z_ax^perp).                    (6)

The generator domain is not inferred from a formal differential expression:

    D(G_*)=closure_graph{q=-ad^*_v m_* smooth and decaying,
                         q perpendicular Z_ax,
                         G_*q in X_ax},                              (7)

where the graph norm is dimensionally normalized by a proved reference streamline
period `T_ref`,

    ||q||_Y^2=||q||_X^2+T_ref^2 ||G_*q||_X^2.                        (8)

The construction must prove density, closability, invariance of the orbit tangent,
and generation of a strongly continuous group or semigroup before any spectral or
time-growth verdict is stated.

## Frozen global operator and nonlocal split

The KKS convention is frozen before using an action: with `i_X Omega=dH` and
`dot q=J grad H`, the quadratic kinetic term is
`-Omega(q,dot q)/2`.  The opposite sign produces `dot q=-J grad H` and is an
exposing canonical failure, not an alternative convention used later.

For Cao--Zhan, corrected 0014 gives the regular-active-core Hamiltonian
expression

    G_e,core(eta,chi)
      =({zeta,K eta-e chi}+{xi,chi/r^2-e eta},
        {xi,K eta-e chi}).                                           (9)

It equals `J_e Hess(A_e)` only where `dA_e=0` is valid.  The global
linearization frozen for 0019 is

    G_e,global q=J_e Hess(A_e)q+(delta J_e[q])dA_e,                   (10)

interpreted through direct linearization of Euler across the active-set/free-boundary
transition.  The second term, including any interface or weak contribution, is
derived and retained rather than set to zero by the interior constitutive relation.

The exact regular-core split used in 0019 is

    G_e=G_e,loc+R_K,
    R_K(eta,chi)=({zeta,K eta},{xi,K eta}),                           (11)

with every remaining term, including derivatives of `r^-2`, retained in
`G_e,loc`.  For a regular Gavrilov energy--Casimir patch, the same split is made from

    L_G(eta,chi)=(K eta+A'(xi)chi, A'(xi)eta+D chi),
    G_G=J_G L_G,                                                      (12)

and `R_K` is again exactly the pair in (11).  At the flat cutoff boundary, where a
global regular Casimir has not been constructed, `G_G,global` is derived directly
from the physical linearized Euler equation rather than extrapolated from (12).
Equation (11) is not discarded as “lower order”: it is the global pressure/velocity
reconstruction whose control is the main achievement of this attempt.

Use a smooth partition subordinate to the regular annulus and its physical collars
only as an identity decomposition.  Commutators `[K,varphi_A]`, Green-kernel transfer
between annulus and complement, free/support-boundary behavior, and the dipolar Cao
exterior are explicit terms in `R_global`; no cutoff-supported subproblem may be
renamed the open-space carrier.

## Exact characteristic and return-map object

On a regular annulus, define travel time from a fixed transverse section by

    T(I)=oint_{Gamma_I} r ds/|grad psi_*|,
    Omega(I)=2 pi/T(I),
    theta=Omega(I) tau mod 2 pi.                                     (13)

Then the meridional characteristic map is exactly

    Phi_t(I,theta)=(I,theta+Omega(I)t),                               (14)

up to the sign fixed by substitution into (9), and coarea gives the physical measure

    dnu=[T(I)/(2 pi)] dI dtheta.                                     (15)

Equation (15), the section orientation, and the sign of (14) must be rederived rather
than copied from canonical action--angle notation.

Let `U_0(t)` be the exact propagator of `G_*,loc` on the global orbit variables after
transport through (14), with annulus localization used only inside the partition of
identity.  The primary characteristic object is the circuit operator

    M_0(I)=U_0(T(I);I)                                                (16)

on the actual accessible generator/amplitude fibre.  It must include the integrated
off-diagonal shear, not merely the eigenvalues of the frozen matrix.  Record its full
Jordan structure and the streamline twist

    tau_tw(I)=Omega'(I).                                             (17)

For Gavrilov, first derive `T(I)` and the first nonzero twist coefficient from the
published local analytic series and the selected cutoff.  An identically isochronous
annulus is a route result, not permission to assume mixing.  For Cao--Zhan, prove any
twist or circuit cancellation from its PDE/monotonicity data; do not fit a period
profile.

The full nonlocal return over the reference period is

    P_*=S_*(T_ref),                                                   (18)

where `S_*` is generated by the closed full `G_*`.  Its exact comparison with (16) is
the Duhamel identity

    P_*-U_0(T_ref)
      =int_0^Tref S_*(T_ref-s) R_global U_0(s) ds.                    (19)

This operator identity, not a matrix obtained by truncating `K`, is the frozen
return-map target.

## Primary semigroup/phase-mixing target

Conjugate the retained remainder into the exactly solved characteristic frame,

    V(t)=U_0(-t) R_global U_0(t).                                    (20)

The first target is to compute the exact circuit average of `V(t)` and decide which
of the following predeclared normal forms applies:

- **direct scattering:** the twist and circuit average give
  `int_0^infinity ||V(t)||_(Y->X) dt<infinity`;
- **modified scattering:** a nonzero explicit resonant `t^-1 B_res` term is removed
  by its derived normal form `N_res(t)`, and the conjugated remainder is integrable;
  or
- **critical/isochronous obstruction:** the nonintegrable term is retained in the
  return operator and tested by the resolvent route below rather than hidden in an
  error constant.

The concrete positive propagation target is

    sup_{t in R} ||N_res(t)^(-1) U_0(-t) S_*(t)||_(Y->X) < infinity,  (21)

with `N_res=I` in the direct-scattering case, together with the explicit growth or
decay of `U_0(t)N_res(t)` in the physical kinetic norm.  Equation (21) must control
the carrier complement and Green tail as well as the regular annulus.  It earns no
more than controlled linear persistence unless the resulting physical velocity norm
is uniformly bounded after translation modulation.

If time-domain integrability fails, solve the full continuum resolvent equation

    (lambda-G_*)q=f,
    q=R_0(lambda)f+R_0(lambda)R_global q,                             (22)

by characteristic integration on (13)--(16).  Freeze the exposing spectral target

    ker(lambda-G_*)={0} for Re(lambda)>0                              (23)

plus a derived resolvent bound strong enough, through a stated semigroup theorem, to
imply the claimed time growth.  A formal Fredholm or Neumann series is not a result
until compactness/smallness, critical layers, limiting absorption, and the operator
domain are proved.  Conversely, a nonzero kernel vector in (23) is a global growing
mode only when it belongs to `D(G_*)`, is dynamically accessible, satisfies the
open-space conditions, and reconstructs finite perturbation kinetic energy.

The physical observables reported with (21) are the modulated perturbation kinetic
energy, the Biot--Savart velocity on and outside the core, and the displacement of a
declared vorticity/impulse centre.  Vorticity-gradient growth with bounded or decaying
velocity is phase mixing, not automatically particle destruction; unbounded kinetic
energy or loss of core localization is a distinct verdict.

## Registered routes and failure-derived continuation

### Route A -- Gavrilov regular pressure annulus

Use the exact compact field, pressure foliation, and flat exterior.  Prove (13)--(17)
first on a strictly interior annulus, solve the local return exactly, then estimate
the collars using cutoff flatness and retain all cross-annulus `K` terms in
`R_global`.  This route is preferred because its invariant pressure shells and zero
exterior are geometrically explicit.  Its main risks are an isochronous leading
centre, degenerating period at the flat cutoff, and absence of a global regular
Casimir; none prevents deriving the physical linearized Euler operator directly.

### Route B -- Cao--Zhan active-core annulus and limiting absorption

Use the exact translating-frame ring and its monotone `z` structure.  Establish a
regular interior foliation, compute the return object, and treat `K` with the full
open-space Green operator.  If twist estimates are unavailable, use (22) as an
operator-valued characteristic equation and separate discrete growing modes from
continuous convective spectrum.  The free boundary and exterior tail remain part of
the theorem rather than a box truncation.

### Route C -- energy-weighted axisymmetric Kelvin/trajectory return

Before any raw-symbol conclusion is used, derive the meridional Kelvin ray

    dot x=u_mer,       dot k=-(grad u_mer)^T k,       kappa=|k|.      (24)

The kinetic-energy symbol `K#=r^2/kappa^2` makes the natural compact-annulus
variables equivalent to

    y=(eta/kappa,chi).                                                (25)

After removing the common `O(kappa)` convective phase, the Cao active-core
system contains, among all coefficient-transport and wave-vector-normalisation
terms, the paired entries

    y1_dot <- i [2 xi k_z/(r^4 kappa)] y2,
    y2_dot <- i [r^2 b/kappa] y1,
    b=(xi_r k_z-xi_z k_r)/r.                                        (26)

Both entries are `O(1)` along a homogeneous Kelvin ray and their product can
generate a finite real rate.  The actual trajectory system must therefore
evolve `k(t)` by (24), include the `dot kappa/kappa` connection created by
(25), retain the complete subprincipal coefficient matrix and the
global/interface term in (10), and compute its circuit monodromy.  Equation
(26) is an exposing pair, not a two-entry truncation licensed to decide
stability.  This route can find an `O(1)` short-wave instability that the raw
first-differential-order Jordan matrix cannot see.

### Route D -- nonaxisymmetric helical Kelvin/Floquet diagnostic

Axisymmetric control cannot by itself establish physical particle persistence.
Decompose a full dynamically accessible perturbation into azimuthal sectors
`exp(i m theta)`, `m in Z`, and retain the base swirl in the three-dimensional
Lagrangian trajectory.  For a meridional circuit `Gamma_I`, define its actual
azimuthal advance

    DeltaTheta(I)=int_0^T(I) xi_*(x_I(t))/r_I(t)^2 dt.                (27)

Starting from the full linearized incompressible Euler equation and Leray pressure
projection, derive rather than import the Kelvin bicharacteristic system for

    dot x=u_*(x),
    dot k=-(grad u_*)^T k,
    k dot b=0,                                                        (28)

and its divergence-free amplitude `b`.  The frozen nonaxisymmetric return object is

    M_m(I)=exp(-i m DeltaTheta(I))
           Texp int_0^T(I) A_Kelvin(x_I(t),k_I(t))dt.                (29)

The sign in (29), pressure projection, density, accessible polarization, and full
Jordan structure are derived explicitly.  Accepted claim C-FLO-001 may check the
finite-matrix power criterion only **after** (29) has been physically derived; it
does not supply the Euler cocycle, operator domain, or continuum stability.

A multiplier with modulus greater than one, or a non-semisimple unit multiplier,
identifies a short-wave/Floquet instability mechanism.  Before promoting that to a
full linear Euler instability, construct finite-energy dynamically accessible wave
packets whose residual is controlled over repeated circuits.  Unit semisimple
multipliers only refute this local Kelvin mechanism; they do not prove the global
three-dimensional semigroup bounded.  Finite `m` sectors still require their exact
nonlocal Biot--Savart operator.

### Route E -- representation repair if action--angle coordinates degenerate

If a critical streamline or flat boundary makes `(I,theta)` singular, retain the
same carrier and pass to travel-time/coarea charts on overlapping annuli, with
transition and Green-coupling terms explicit.  If no regular annulus exists for one
carrier, that carrier/representation route is refuted by the named mechanism and the
other same-carrier route continues.  Coordinate failure is not a global Euler
stability verdict.

Routes are selected in this order: exact physical orbit and open-space domain;
analytic solvability of the characteristic return; retention and controllability of
`K`; finite energy and finite nonzero `j`; assumption/parameter economy; correct
symmetry, dimensions and limiting behavior; strength of the time-domain physical
observable; and only then implementability.  No empirical comparator selects a
route.

## Frozen source-access inventory

No new external body is opened by this preregistration.  After activation, access is
recorded before each new body and full papers are cached only under
`/tmp/primary-source-cache/P253-0019`.

| Source/input | Pre-activation access | Permitted role after activation |
| --- | --- | --- |
| Corrected P253/0014 derivation/result and independent `0020/review.md` | Local active-attempt artifacts; formal review is complete and source/result hashes are pinned at its review boundary | Exact axisymmetric tangent and mixed sign; corrected action sign, regular-core/global scope, and raw order-`|k|` inference ceiling |
| P253/0010 and 0005 | Local proposal artifacts with primary receipts | Same finite-`j` Cao--Zhan/Gavrilov carriers, conserved functional, full orbit/KKS action |
| Gavrilov arXiv:1810.08020v1 and Constantin--La--Vicol arXiv:1903.11699v1 | Full bodies already hashed and theorem-located in `0010/source-audit.md` | Recheck local analytic pressure series, cutoff identities, and regular annulus; no stability theorem imported |
| Cao--Zhan arXiv:2009.13210v2 | Full body already hashed and theorem-located in `0010/source-audit.md` | Recheck regularity, monotonicity, core/free-boundary facts and exact carrier equations |
| Lifschitz--Suters--Beale DOI 10.1051/proc:1996017 | Metadata only; prior full-text attempts failed | Stability-route inventory only unless a lawful full body is obtained; no result inferred from title/abstract |
| C-FLO-001 and `rotating_stability.py` | Accepted finite-dimensional algebraic API | Jordan/power criterion for an already-derived finite cocycle only; never a field realization or semigroup theorem |
| Standard coarea, elliptic regularity, semigroup and resolvent theorems | Declared mathematical imports only after exact hypotheses are sourced | Measure, mapping/generation, and implication steps; hypotheses must be checked for the actual weighted open-space operator |

Any geometric-optics/Kelvin-wave, inviscid-damping, limiting-absorption, or Euler
semigroup paper found after activation receives URL/version, lawful access result,
hash, and exact theorem/equation location before use.  Literature absence or a
paywall is not a physical no-go; the full linearized Euler and characteristic
equations remain directly derivable.

## Analytic ladder and oracle boundary

After central activation, the work order is frozen:

1. fix the canonical action sign, then derive the full operator, `(delta J)dA`
   and interface rows, orbit range, symmetry quotient, and graph domain;
2. prove a concrete regular annulus for the fixed Gavrilov carrier and derive
   (13)--(15), including orientation and measure;
3. solve `G_loc` along characteristics, compute `M_0(I)`, its Jordan structure,
   circuit shear, period and twist, then derive the energy-weighted system
   (24)--(26) with evolving covector and every order-zero amplitude term;
4. derive exact mapping, commutator and tail estimates for `K`, including annulus
   leakage and boundary collars;
5. compute the interaction-frame leading term and execute direct or modified
   scattering in-run; if it is nonintegrable, continue through the resolvent equation
   rather than label mixed signature an instability;
6. repeat the decisive analytic object on Cao--Zhan if it materially changes carrier
   coverage; and
7. derive the nonaxisymmetric Kelvin/Floquet cocycle and, when it exposes growth,
   construct a full finite-energy accessible wave packet before any instability claim.

The strongest practical oracle is an exact characteristic/return derivation plus a
continuum operator estimate in the physical norm.  No production numerical verifier
is registered in this README.  If an analytically irreducible Floquet multiplier,
resolvent boundary, or soft growth rate later requires numerics, a new activated
design must first load the small-ratio-numerics skill and freeze precision, residual,
domain/tail, refinement, pseudospectral sensitivity, and error-budget requirements.
A boxed eigenvalue, a finite matrix truncation of `K`, or unit-modulus sampled
multipliers can earn only representation-scoped exploratory evidence.

Every route receives one route-scoped verdict and triggers method repair,
representation change, and materially different candidate continuation when needed.
No result in 0019 closes the electron-and-neutrino campaign unless the full parent
closure map separately proves that implication.

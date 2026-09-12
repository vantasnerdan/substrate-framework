# P253/0052: fixed-domain, full-sector Cao graph-Riesz continuation

Status: **README-only preregistration awaiting central registration and schema
activation.** Until `0052/activation-schema.exit` exists and contains exactly
`0`, this README is the sole 0052 artifact. Root owns `proposal.yaml`, the
central attempt registry, activation receipts, commits, review, and promotion.
No derivation body, verifier, source download, numerical run, API, or claim is
authorized before activation.

## Failure-derived objective and preserved result

Attempt 0048 found the correct leading finite-window scale on the compact-core
dynamically accessible (DA) Cao domain:

    c_n^E=O(k/n),             g_n=delta c_n^E=O(delta^2/n),
    |g_n|/Delta sigma_n=O(1/log(1/delta))                 (1)

at the bending/Kelvin crossing. It follows that the interaction window is
finite for each fixed small `delta`, its high-index width is smaller than one
local Kelvin spacing, and the dyadic bad-parameter fraction is
`O(1/log(1/delta))`. Those leading scaling statements are retained.

They do **not** yet establish the full graph-domain Riesz theorem. Attempt
0052 targets the six constructions identified by independent supervision:

1. the corrected crossing-index expansion, including its two distinct error
   scales;
2. an explicit fixed-domain transport with uniform energy/graph equivalence
   and a common-domain conjugated Euler expansion;
3. a complement resolvent over every cross-sectional Fourier sector, not
   only the `m=0` tail;
4. exact and mutually consistent KKS, energy, generator, Fourier-measure, and
   physical-frequency normalizations;
5. a proved full moving-boundary/global-Green Euler graph jet, or an explicit
   conditional theorem whose hypothesis is not attributed to Cao; and
6. an honest nonlinear boundary: the 0048 formulas corresponding to
   `(52)--(53)` remain approximate finite-order residual expansions, not an
   exact branch.

The positive target is a physical fixed-`delta` Riesz family for the exact Cao
ring, with the two `m=1` translation polarizations and every `m=0` Kelvin mode
inside the declared coupling window, reconstructed to compact vorticity and
whole-space finite-energy velocity. The universal thin-filament matrix remains
only a leading compression and is never substituted for the finite-core
operator.

## Fixed carrier, parameters, and exact crossing target

Fix one polynomial-profile no-swirl Cao ring with integer `p>=3`, circulation
`Gamma>0`, major radius `R_delta`, core radius `a_delta`, and one ring Fourier
number `l>=2`:

    delta=a_delta/R_delta,          k_delta=l delta,
    L_delta=log(1/delta).                                  (2)

For the limiting column put `Phi=2 Omega W` and

    L_Phi=integral_0^a sqrt(Phi(r)) dr,
    A_l=|l| L_Phi/pi.                                     (3)

Let the renormalized Kelvin and bending frequencies in one common core-time
unit satisfy

    sigma_n(delta)=A_l delta/(n+beta_0)
                    +O(delta/n^2+delta^3/n),
    tau_l(delta)=c_l delta^2 L_delta+d_l delta^2
                    +o(delta^2),                           (4)

where `beta_0`, `c_l>0`, and `d_l` are derived from the source-defined
Sturm--Liouville, KKS, and graph-jet problems. Inverting (4), rather than
dropping its finite terms, must give

    n_*(delta)=A_l/(c_l delta L_delta)
       +O(1/(delta L_delta^2)+1).                          (5)

The first error in (5) is the finite bending-frequency correction; the second
is the Sturm--Liouville phase. A verifier must expose the incorrect smaller
`O(L_delta^-1)` remainder used in 0048. The finite interaction window and
measure estimate are recomputed using (5), while preserving the leading
conclusions in (1).

## Obligation A: explicit fixed-domain transport

Let `K_delta` be the exact compact toroidal vorticity support in physical
space, and choose a fixed solid torus

    K_*={(s,alpha,theta):0<=s<=1}.                         (6)

Construct, rather than merely name, a volume-preserving Hanzawa/Piola map

    Phi_delta:K_* -> K_delta                              (7)

from the Cao free-boundary graph and its center/scale gauges. It must extend
to a collar and the identity at infinity without changing circulation or
impulse. Define the vorticity Piola transport and its whole-space velocity
completion explicitly:

    U_delta eta_*=det(DPhi_delta)^-1
                    DPhi_delta eta_* o Phi_delta^-1,
    v_delta=BiotSavart(U_delta eta_*).                     (8)

If exact volume preservation makes the determinant one, prove it after the
collar correction. Show that (8) maps the quotient of fixed-domain DA
tangents onto the physical DA closure and preserves distributional support.

On the fixed space define

    ||eta||_(E,*)^2=||BiotSavart_* eta||_2^2+||eta||_H^-1^2,
    D_*={eta: A_* eta in X_(E,*)}.                         (9)

Earn constants independent of sufficiently small `delta` such that

    C^-1||eta||_(E,*)<=||U_delta eta||_(E,delta)
                         <=C||eta||_(E,*),
    C^-1||eta||_(G,*)<=||U_delta eta||_(G,delta)
                         <=C||eta||_(G,*).                 (10)

The conjugated generator

    Ahat_delta=U_delta^-1 A_delta U_delta                 (11)

must have the **common domain** `D_*`. Derive its coefficients from the exact
metric, base velocity/vorticity, moving boundary, whole-space Hodge operator,
and field-dependent Leray projector. The required expansion is

    Ahat_delta=A_col(k_delta)+delta C_1
      +delta^2 L_delta C_(2,log)+delta^2 C_2+R_delta,
    ||R_delta||_(D_* -> X_(E,*))=o(delta^2).               (12)

No comparison of operators on moving, unidentified domains counts as (12).

## Obligation B: exact physical normalizations

Carry the density `rho_0`, the physical toroidal Jacobian, the `2*pi` Fourier
factors, and the core-to-laboratory time conversion through every formula.
On the right-reduced orbit use

    Omega_KKS(X_xi,X_zeta)
      =-<m,[xi,zeta]>
      =rho_0 integral Omega_delta dot (xi cross zeta) dx,
    i_(X_H) Omega_KKS=dH.                                 (13)

For a complex eigenvector `q` satisfying `Aq=i nu q`, the convention in (13)
requires the exact identity

    L(q,conjugate(q))=i nu
       Omega_KKS(q,conjugate(q)).                          (14)

The canonical test `Omega=dq wedge dp`,
`H=(p^2+nu^2q^2)/2`, `q_+=(1,i nu)` gives

    Omega(q_+,conjugate(q_+))=-2 i nu,
    L(q_+,conjugate(q_+))=2 nu^2.                          (15)

Equations (14)--(15) are the sign/factor oracle. Derive the exact energy-unit
and KKS-unit scale of every Kelvin mode and both real translation modes, the
biorthogonal duals, and the matrix conversion. Close the relation among
`lambda_n`, column `sigma_n`, the eigenvalue of `Ahat_delta`, core time, and
laboratory frequency. The earlier proportionality `sqrt(frequency)` is not
used until all factors in (13)--(15) have been evaluated.

## Obligation C: all-sector complement resolvent

Let `P_delta` be the KKS-biorthogonal projection onto the bending translation
pair plus all real/Hamiltonian partners of `m=0` modes within four times their
actual coupling. Let `Q_delta=1-P_delta`. Decompose the fixed-domain operator
by cross-sectional Fourier number before curvature coupling:

    X_(E,*)=direct_sum_(m in Z) X_m.                       (16)

Construct estimates for every part of `Q_delta`:

- `m=0`: use the exact weighted Sturm--Liouville eigenbasis, corrected index
  (5), high-`n` coefficient decay, the retained finite window, and the tame
  tail;
- `m=1`: remove both translations with their exact KKS duals and use the
  finite-core ground-state gap on the common domain;
- `|m|>=2`: combine the positive radial form increment with large-`|m|`
  coercivity and the whole-space Biot--Savart estimate;
- negative `m`, real conjugates, axial/centering modes, collar fields, and
  exterior harmonic velocities: include them explicitly rather than invoking
  symmetry without a domain map; and
- curvature: control every `m -> m+/-1` block, the longitudinal terms that
  preserve `m`, and the second-order infinite matrix in the weighted graph
  norm.

Choose a contour `Gamma_delta` around the complete enlarged cluster and prove

    sup_(z in Gamma_delta)
      ||[Q_delta(Ahat_delta-z)Q_delta]^-1||_(X_* -> D_*)
        <=C/r_delta,                                      (17)

with an explicit `r_delta` compatible with the local Kelvin spacing and all
remainders in (12). Large-`m` summation and finite low-sector estimates must
close (17); a scalar `m=0` Schur sum does not. Then evaluate the exact Feshbach
map

    K_delta(z)=P_delta(Ahat_delta-z)P_delta
      -P_delta Ahat_delta Q_delta
       [Q_delta(Ahat_delta-z)Q_delta]^-1
       Q_delta Ahat_delta P_delta                         (18)

and reconstruct

    Pi_delta=(2*pi*i)^-1 integral_(Gamma_delta)
                 (z-Ahat_delta)^-1 dz.                    (19)

Prove that (19) is bounded on both `X_*` and `D_*`, is KKS nondegenerate, has
the declared real rank, and transports through (8) to compact DA vorticity
with its exact exterior velocity.

## Obligation D: the full Euler graph jet

The exact rescaled Cao elliptic equation determines the local profile cells,
but the published `O(epsilon^2 |log epsilon|)` remainder does not by itself
prove (12). Execute both of the following analytic rungs:

### Route D1: direct moving-boundary/Green proof

Fix the free boundary with (7), differentiate the exact Cao equation and all
circulation/impulse/center constraints, and construct the first and second
profile/shape cells. Expand the exact Green kernel (not only its printed big-
`O` estimate) into singular logarithmic, finite core, collar, and exterior
pieces. Differentiate whole-space Biot--Savart and Leray/Hodge twice, prove
the boundary trace and commutator bounds in (9)--(10), and establish (12).
Determine rather than assume whether the profile has a
`delta^2 log(delta)` cell; keep it separate from the universal Green/local-
induction term.

### Route D2: form/resolvent replacement

If a strong graph Taylor jet fails at finite `p`, construct a common closed
quadratic/Hamiltonian form with norm-resolvent control strong enough for
(17)--(19). State the exact form-domain convergence rate, compactness input,
and contour consequence. This is a representation change, not permission to
assert the strong jet.

If neither route closes, state a named hypothesis `HJ2` containing precisely
(7), (10), (12), and the needed KKS differentiability. Then prove the Riesz
theorem **conditional on `HJ2`** and keep `HJ2` unproved; do not attribute it
to Cao or report the unconditional graph-Riesz claim.

## Route E: nonresonance versus infinite-block control

Recompute the finite-window theorem from (5), the exact normalizations, and
the full all-sector resolvent. Route E1 proves a dyadic parameter estimate
with explicit transversality constants and retains every exact resonance in
`P_delta`. Route E2 keeps the complete `m=0` tail and constructs a KKS-
symplectic tame block transform on a declared weighted graph sequence space.
At an exact resonance, evaluate the enlarged finite KKS and energy matrices;
their Krein signs and discriminant follow from (14), not from a proportional
normalization or from the filament comparator.

The target alternatives are:

1. an unconditional full-sector Riesz family on an asymptotically full-
   measure set of exact Cao carriers;
2. an unconditional tame infinite-block Riesz family including resonant
   exchanges; or
3. a conditional Riesz theorem with the sole unproved input `HJ2`, plus the
   strongest independently established fixed-domain/resolvent subclaims.

Each alternative is a route verdict, not a parent no-go.

## Nonlinear boundary and next bridge

The 0048 polynomial formulas denoted there by `(52)--(53)` are retained only
as an **approximate fixed-order expansion**:

    V^(N)(a)=V_0+sum_(r=1)^N a^r V_r,
    ||F(V^(N),c^(N),Omega^(N))||<=C_(N,delta)|a|^(N+1).    (20)

Equation (20) neither produces an exact solution nor licenses a physical
period, KKS action, stability, or particle interpretation. Record its Sobolev
index loss and the bound `N<=p-2` (or the sharper earned finite-regularity
limit). Do not call it a branch.

Only after (17)--(19) are established may this attempt formulate the exact
nonlinear dependency: a range/transparency construction for the generated
critical divisors

    m Omega(r)-j l Omega_pattern,                           (21)

followed by convergence in the same fixed solid-torus chart. The linear Riesz
theorem is not equated with that nonlinear construction.

## Sources, oracle, and verdict boundary

Primary scope begins with Cao--Lai--Qin--Zhan--Zou,
arXiv:2206.10165v2, for the exact carrier, Green kernel, rescaled elliptic
equation, finite regularity, and published remainder; and Gallay--Smets,
arXiv:1805.05064v3, for the source-scoped column equations and Biot--Savart
estimates. Existing 0040/0044/0048 receipts are reused only after checking
their exact equations. Full papers remain in `/tmp/primary-source-cache`, with
URL, version, SHA-256, and theorem/equation locations recorded. Abstract APIs
may supply conventions but not the physical domain map or resolvent.

The strongest oracle consists of:

1. an exact symbolic/asymptotic check of (5) and (14)--(15);
2. an independent pullback/pushforward check of (7)--(12), including a test
   field crossing the collar and a whole-space pressure field;
3. an all-sector resolvent ledger whose finite low sectors, high-`m` tail,
   `m=0` window, and exterior reconstruction sum to (17); and
4. a Riesz rank/KKS-sign check performed on the reconstructed physical field.

No production numerics are authorized. If an irreducible small spectral
separation survives the analytic ladder, the small-ratio-numerics skill is
read before a later numerical design is frozen.

Every obligation and route receives exactly one of `established`,
`refuted_with_mechanism`, or `blocked_by_named_construction`. Success requires
the actual full-sector physical graph-Riesz family, or an explicitly
conditional theorem whose one unproved graph-jet hypothesis is visible. It
does not require or imply the nonlinear solid-torus branch. P253, LP2, and P4
remain active, and no global impossibility is inferred from a failed graph
representation.

## Post-execution correction boundary (P253/0057)

The corrected resonance law, exact local Cao cells, unequal-volume
Hanzawa/Piola representation, toroidal distance, and KKS normalization remain
established. The former limiting/polyhomogeneous all-sector complement claim
is reduced to its scalar and sectorwise calculations: a nonnormal coupled
`X_* -> D_*` graph-resolvent/Grushin estimate `GR` is independently required.
The exact physical graph-Riesz statement is therefore the conditional
implication `HJ2 + GR => Riesz transfer`; `HJ2` alone is insufficient.

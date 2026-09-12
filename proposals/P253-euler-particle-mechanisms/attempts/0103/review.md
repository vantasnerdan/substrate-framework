# P253/0103 independent final review of P253/0102

Reviewer: `particle-balance-review`

Target owner: `root`

Activated review README SHA-256: `d28e70bc7d738c85d91da53cce88541a76040fba8b99e0e7ff2ab4afd363c1c2`

Frozen target manifest SHA-256: `ffb741bc9850ded93e1a0cb984fc68226dd3c51c2b23b457e81e2f84fdde19b1`

## Boundary and provenance

This was a content-blind, independent, non-author review of final
post-Liouville P253/0102. The activated README and all frozen target selectors
matched, and `activation-schema.exit` contained exactly `0`. P253/0095 was not
opened or used. No P4, P2, chiral/flavor, particle or parent conclusion was
adjudicated.

The manifest uses repository-relative paths. All seventy immutable target and
receipt entries verify. Its one mutable central selector,
`proposals/P253-euler-particle-mechanisms/proposal.yaml`, changed from the
historical manifest hash
`9622c808412f71ecf08b309bdd81e476b3f4abd104798840bfca7bfbcddd122f` to
`1f7388df8633b251e144afec34a941d3f16c156d38de47b60576a440bc3cedd4`
when P253/0103 was centrally registered; this expected post-target registry
change does not alter a 0102 artifact or scientific predicate. The final
exact-v2 receipt has thirteen derived PASS predicates, empty stderr and exit
zero. The focused-v2 receipt separately has fourteen passing tests, empty
stderr and exit zero. Static-v3 and repository-v3 exit zero. These receipts
were reused without rerunning unchanged evidence.

The post-Liouville v2 completion receipt supersedes the older completion
receipt scientifically while preserving it as chronology. The expansion
receipt transparently records removal of the failed universal-sign argument
and addition of the Liouville and annular source routes before independent body
review.

## Source and authority audit

The reviewed inputs supply only the transported electromagnetic tag, charged
axisymmetric no-swirl Cao equilibrium, modified Grad--Shafranov and Maxwell
conventions, Ertel current, and conditional two-state algebra at their stated
scopes. P253/0102 directly derives the new quotient, contour, Green, Liouville
and source-moment identities. None of its sources supplies the missing
normalized Cao response, finite-carrier transfer, admissible profile lift,
moving compact defect, shared P4 action or P2 persistence.

## A — variable two-label forced equilibrium lock and support

Let `chi_c` be transported and let the distinct circle phase `Theta` have an
advected closed one-form. With

    q_A=omega dot grad Theta=lambda chi_c,

the forced Ertel equation gives the undivided compatibility law

    chi_c D_t lambda=(curl f) dot grad Theta.

Because the labels are distinct, contraction along a vorticity line does not
differentiate `chi_c`; the old same-label exponential monodromy inference does
not apply.

For

    f=-(g/rho_m)chi_c(grad Phi+H grad P),

direct cylindrical curl calculation gives

    curl f=(g/rho_m)grad(chi_c' Phi-chi_c H) cross grad P,
    (curl f) dot grad theta
      =(g/rho_m)W dot grad(chi_c' Phi-chi_c H).

The sign and factor of `r` are correct for
`W=r^(-1)grad P cross e_theta`. At the no-swirl equilibrium,
`q_A=zeta`; the modified Grad--Shafranov row then makes
`lambda=zeta/chi_c` satisfy the forced compatibility equation wherever
`chi_c` is nonzero. This is an exact variable equilibrium identity, not force
selection or dynamical persistence.

The support obstruction is equally exact. The existing `chi_c` vanishes
outside its strict regular band, while `zeta` is nonzero elsewhere in the Cao
core, so the undivided residual `zeta-lambda_0 chi_c` cannot vanish globally.
The quotient is therefore only a local/band object unless a new full-core tag
and stabilizer are constructed. Any band boundary current remains part of that
extension problem.

**A verdict: established at the nonzero-tagged equilibrium band, with the
global constant-lock route refuted for the existing support geometry.**

## B0 — exact contour reduction and Maxwell normalization

At `g=0`, constant locking imposes
`chi_0=zeta_0/lambda_0`. With `tau=g^2`, the order-`tau` profile, tag and
constant-ratio corrections are functions only of the streamline label. A
period/action contour projection `Q_P` therefore leaves the necessary row

    Q_P[chi_0' Phi_hat-chi_0 H_hat]=0.

This is a necessary compatibility equation, not a nonvanishing theorem. The
finite API correctly models only already-integrated contour sectors; the
continuum action-period measure and mapping properties remain analytic inputs.

The signed-charge-factored Lorenz source is

    L_c Phi_hat=(chi_0/epsilon_EM)
      (1-c^2/c_EM^2-c W_z/c_EM^2),

with no residual factor of `g`. The Ampere primitive is

    epsilon_EM a c_EM^2 r H_hat
      +epsilon_EM c partial_r Phi_hat=K(P)/r,
    a=1-c^2/c_EM^2,

so elimination gives

    chi_0' Phi_hat
      -chi_0 K(P)/(epsilon_EM a c_EM^2 r^2)
      +chi_0 c partial_r Phi_hat/(a c_EM^2 r).

The signs, `epsilon_EM`, `c_EM`, and radius factors are consistent. The curved
`r^(-2)` term cannot be projected away before a straight-column replacement.
All global potentials must retain the same Green normalization and their
affine/logarithmic tails; the target says so explicitly.

**B0 verdict: established as an exact contour compatibility reduction with
the full Lorenz/Ampere normalization.**

## B1 — normalized straight-column Cao response

For a fixed straight profile and compact radial tag, the first-speed scalar
dipole satisfies

    L_1 v=f,
    v(s)=-1/2[s^(-1) integral_0^s t^2 f(t)dt
                   +s integral_s^infinity f(t)dt].

The Ampere derivative is
`H_1=-Phi_0'/(c_EM^2 R)`, so the first harmonic of the lock response is

    S_1=chi_0' v-chi_0 H_1.

The Green formula and relative normalization are exact. Positive Cao
vorticity does not fix the sign because the physical orientation of `P_s`
makes the electric and magnetic terms capable of opposing each other. The
target neither evaluates `S_1` for the normalized Lane--Emden/Cao profile and
tag nor proves a normalized response lower bound.

**B1 verdict: blocked at evaluation of the actual normalized Cao profile,
physical orientation, major-radius scale and tag cutoff.**

## B2 — finite-Cao transfer

Even a nonzero straight-column coefficient would require uniform coupled
Maxwell Green convergence on an actual subluminal Cao path, retention of the
curved `r^(-2)` Ampere row, a common contour domain and control of the tagged
band boundary. None is constructed, and none is inferred from the exact
column calculation.

**B2 verdict: blocked at the uniform finite-curvature Green and band-boundary
transfer.**

## C0 — exact Liouville profile, cancellation and energy boundary

For `a_0,b>0`, the proposed radial pair is

    P_b=P_0-(2/a_0)log(1+b s^2),
    zeta_b=8b/[a_0 R^2(1+b s^2)^2].

Direct differentiation gives

    -Delta_2 P_b=R^2 zeta_b,
    d zeta_b/dP_b=a_0 zeta_b.

It is regular at `s=0` and has finite cross-sectional vorticity:

    2 pi integral_0^infinity zeta_b s ds=8 pi/(a_0 R^2).

With `chi_0=zeta_0/lambda_0`, the globally matched Poisson derivatives obey
`Phi_0'=P'/(lambda_0 epsilon_EM R^2)`. Differentiating the fluid equation in
the first harmonic gives the unique regular-at-zero, decaying-dipole solution

    v=-P'/(lambda_0 epsilon_EM c_EM^2 R^3 a_0).

There is no additional homogeneous `s` term when both fields use the same
global Green/matching normalization. With
`H_1=-P'/(lambda_0 epsilon_EM c_EM^2 R^3)`, the two response terms cancel:

    chi_0' v-chi_0 H_1=0.

The scope boundary is load-bearing. As `s` tends to infinity,
`P_b'` is proportional to `1/s`, so the straight-column velocity has a
`1/s` tail and its cross-sectional kinetic energy per unit axial length
contains `integral ds/s`, which diverges logarithmically. The vorticity decays
as `s^(-4)` but is noncompact. Thus this is a genuine smooth global Liouville
column and exact representation-level kernel, not a localized finite-energy
Cao/free-boundary carrier. It also disproves only the abandoned universal-sign
argument; it does not determine the sign for arbitrary profiles.

**C0 verdict: established exactly at global straight-column first-speed
scope, with finite cross-sectional vorticity and logarithmically divergent
energy per unit length.**

## C1 — annular unconstrained Maxwell-source range

On an annulus where the exponential law holds, the homogeneous mismatch is
`A s+B/s`. The dipole Green formula shows that an inner source changes only

    delta B=-(1/2) integral t^2 delta f(t)dt,

whereas an outer source changes only

    delta A=-(1/2) integral delta f(t)dt.

Using these two moments as source coordinates, the map to `(A,B)` is

    [[0,-1/2],[-1/2,0]],

whose determinant is exactly `-1/4`. Smooth bumps supported strictly inside
and outside can independently realize the two unconstrained moments.

This is only the range of radial Maxwell source data. The constant-lock rows
tie the tag and vorticity through
`delta chi=delta(zeta/lambda)`; charge, positivity, free-boundary, Casimir and
material-leaf conditions may reduce or collapse the apparent two-dimensional
range. The target consistently makes that distinction and claims no universal
physical response sign.

**C1 verdict: established at unconstrained annular Maxwell-source range
scope.**

## C2 — physical Cao profile and material-leaf lift

No bounded map from admissible fixed-leaf
`(delta P,delta zeta,delta chi)` variations to the two moments `(A,B)` is
constructed. In particular, no proof simultaneously enforces
`delta chi=delta(zeta/lambda)`, the profile equation, free boundary,
circulation, mean radius, charge/tag distribution, stabilizer, center and Hodge
rows. Recomputing the unconstrained determinant cannot provide that lift.

**C2 verdict: blocked at the joint profile/free-boundary/material-leaf map and
its range.**

## D — noncompact axis phase versus compact defect

For `Theta=theta`, the puncture is the entire symmetry axis. Its global closed
one-form has norm `|dtheta|=1/r`; the unweighted phase-gradient energy diverges
at the axis and also needs explicit outer treatment. This noncompact geometric
phase is not a compact defect carried by the toroidal vortex. Replacing it by a
ring defect requires a new circle map, and neither `Theta=theta` nor
`q_A=zeta` transfers automatically. No moving-domain, branch-independent,
finite-action construction is supplied.

**D verdict: blocked at a compact advected phase/defect with controlled inner
and outer energy/flux.**

## E — supplier boundary and evidence attribution

The exact verifier derives the cylindrical sign, undivided quotient and
support residual, finite weighted projection identities, Maxwell factors,
radial Green identity, Liouville equations, cancellation and source determinant.
The finite contour projector and algebraic helper functions are not proofs of
continuum mapping properties or physical normalization. The focused tests are
API regression evidence. Static and repository validation are agreement and
workflow evidence.

The final claim-bearing artifacts consistently deny a universal sign,
normalized Cao response, finite-Cao transfer, admissible material-leaf source
lift, compact defect, P2 persistence, shared P4 action, Lorentz chiral current,
flavor interaction, particle, electron, neutrino or parent conclusion.

**E verdict: the exclusions are established and the flavor/P4/P2 bridge
remains blocked at its named constructions.**

## Combined verdict

P253/0102 is **established at its exact equilibrium, contour, Liouville and
source-range supplier boundary**. It proves the variable two-label equilibrium
lock on the nonzero tag band, the support obstruction to the existing global
constant lock, the normalized Maxwell contour reduction, an exact nonlocalized
Liouville first-speed kernel, and independent unconstrained annular source
moments. The actual normalized Cao response, finite-carrier Green transfer,
admissible fixed-leaf profile lift and compact moving phase remain open.

No bounded correction is required. The strongest next achievement is the
actual fixed-leaf map from admissible `(delta P,delta zeta,delta chi)` to the
two annular moments, followed by normalized Cao response evaluation and a
uniform finite-curvature transfer; the compact moving phase and P4/P2 bridge
remain independent downstream constructions.

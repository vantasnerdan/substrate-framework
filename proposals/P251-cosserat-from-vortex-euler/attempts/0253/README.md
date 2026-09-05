# 0253 — outer-only localization and non-Lipschitz compact-core candidates

Owner: `herdr geometry-review pane w3:p4`, continuing independently from the
issue200 coordinator and the fixed-ring response pane.  Write scope is new
attempt 0253 only.  Base checkpoint is `d9d8bbb`; the parent remains the full
issue200/P251 objective.  No source module, prior attempt, central proposal
entry, registry, release, generated file, or commit is changed here.

This contract is frozen before substantive predecessor/source-body opening.
The coordinator must add the central manifest entry and a successful schema
receipt before either candidate is analyzed beyond the equations and
failure-derived facts supplied by the user.  The construction is analytic;
no numerical solver or long script is currently licensed.

## Four frozen levels

The parent objective is unchanged: one stationary smooth constant-density
Euler ensemble supports actual coupled histories, inherited physical action
and current, and compact Euclidean nonzero finite-core tubes at usable
positive measured density.  The completed periodic join and the direct fixed
exact Beltrami-background plus 0250 response route remain active independent
inputs.  No reduction of the parent to compact geometry alone is permitted.

The active obligation is an exact compact-**velocity** stationary Euler
carrier which preserves a nonzero nondegenerate compact Euclidean core while
avoiding the reviewed 0252 incompatibility.  Localizability may hold in an
outer cutoff collar, but not in a neighborhood of the nonzero core.  A field
with merely compact vorticity, an exterior uniform velocity, a distributional
free boundary, or nonsmooth pressure does not meet this obligation.

The two candidate routes are:

1. match the nonlocalizable core to an outer shell on which pressure is an
   actual first integral, then apply the exact pressure cutoff only there;
2. construct the compact velocity directly from a genuinely non-Lipschitz
   Grad--Shafranov label law, without pressure localization.

This concrete attempt will source-audit novelty and applicability after
activation, derive the edge/core/pressure equations as far as exact analysis
permits, and give one route-scoped verdict to each route.  No finite route
outcome establishes scientific exhaustion.

## Pre-source filename and memory inventory

The archive was searched by filenames and durable memory before freezing the
candidate details.  The nearest records are:

| Located record | Pre-source role | 0253 novelty boundary |
| --- | --- | --- |
| `0143/global-localizable.md` | prior global/localizable construction route | Does not by filename supply an outer-only collar preserving a nonlocalizable core |
| `0169/localization-and-interaction.md` | exact pressure localization and disjoint compact-support assembly | Reuse requires proving localizability precisely where the cutoff changes; it cannot be extended across the 0252 core |
| `0198/localizable-thin-ring.md` | actual localizable annular thin-ring supplier | Its annular/quiet-core geometry is not the required nonzero core |
| `0193/compact-vorticity-alternative.md`, `0186/global-ring-and-tag.md`, `0195/swirling-ring.md` | compact-vorticity/global-ring alternatives | Compact vorticity is not compact velocity or an exact disjoint assembly license |
| `0211/same-ring-beltrami-torus.md` and `0252/center-compatibility.md` | reviewed literal nonzero inner core and reviewed local incompatibility | The 0211 core is preserved; only the core-localizable matching subclass is excluded |
| `0251/source-transfer.md` | new source-license boundaries | The 2026 axisymmetry theorem is guidance, and Baldi is a particle action-angle supplier; neither is the matching or compacton theorem |
| `0250/fixed-ring-response.md` | fixed-background response route | Remains an independent active route and an eventual dynamical consumer, not evidence that either stationary candidate exists |
| durable memory `codex/efforts/pr199-completion.md` | records Gaussian positive-density, compact-vorticity, localizable-shell and ring continuations | No memory entry located the split localizable collar or the non-Lipschitz compact-velocity profile proposed here |

This table is a novelty inventory, not authority.  Each proposition actually
used after activation must be checked in its source body and applicable review.

## Shared axisymmetric Euler equations and full pressure

Work in physical cylindrical coordinates `(r,theta,z)`, away from the axis
unless and until axis regularity is separately proved.  Use

\[
 u=-\frac{\psi_z}{r}e_r+\frac{F(\psi)}{r}e_\theta
       +\frac{\psi_r}{r}e_z,
 \qquad
 \Delta^*\psi=\psi_{rr}-\frac1r\psi_r+\psi_{zz}.          \tag{1}
\]

The stationary Grad--Shafranov equation is frozen with the sign convention
specified by the user:

\[
 -\Delta^*\psi=F(\psi)F'(\psi)-r^2B'(\psi).              \tag{2}
\]

Direct cylindrical curl gives

\[
 \omega_r=F'(\psi)u_r,\qquad
 \omega_z=F'(\psi)u_z,\qquad
 \omega_\theta=-\frac{\Delta^*\psi}{r}.                 \tag{3}
\]

Equations (2)--(3) imply

\[
 u\times\omega=\nabla B(\psi),\qquad
 p=B(\psi)-\frac12|u|^2,                                \tag{4}
\]

so (1)--(2), not a later pressure patch, must produce the full stationary
Euler pair.  Smoothness of `B` as an abstract label function at `psi=0` is
not required a priori in Candidate B, but the physical composition in (4)
must extend smoothly across the support edge.

The weighted variational functional associated with (2) on a meridional
domain is

\[
 \mathcal J[\psi]
 =\int\!\left[
    \frac{|\nabla\psi|^2}{2r}
   -\frac{F(\psi)^2}{2r}+rB(\psi)
   \right]dr\,dz.                                       \tag{5}
\]

Its first variation is

\[
 \delta\mathcal J
 =\int\frac{-\Delta^*\psi-F F'+r^2B'}{r}\,
        \delta\psi\,dr\,dz                             \tag{6}
\]

after the actual boundary term is shown to vanish.  Candidate B must specify
an admissible weighted space in which (5) is finite, the non-Lipschitz
potential is differentiable in the needed weak sense, and a weak critical
point becomes the required smooth flat-edge solution.  A formal Euler--
Lagrange equation alone is not an existence proof.

## Candidate A — nonlocalizable core with an outer localizable collar

Start from a smooth axisymmetric carrier on nested stream surfaces.  On an
inner neighborhood `N_core` it must agree with, or preserve by a proved open
perturbation, a nonzero nondegenerate constant-curl core.  In a smaller inner
region require

\[
 F(\psi)=\lambda\psi,\qquad B'(\psi)=0,
 \qquad -\Delta^*\psi=\lambda^2\psi,                    \tag{7}
\]

with a positive interior critical value `psi_c`, definite meridional Hessian,
and `F(psi_c)/r0 != 0`.  Pressure need not and, by 0252, cannot be a first
integral near this core when its poloidal derivative is invertible.

Choose an interface stream surface `Sigma_1` outside `N_core`.  Across it,
match the full Cauchy and label data needed for a globally smooth solution of
(2): `psi`, its normal derivatives to all required orders, `F`, `B`, and the
physical velocity, vorticity and pressure jets derived from them.  In an
outer collar `S_loc`, impose the additional speed first-integral equation

\[
 \frac{|\nabla\psi|^2+F(\psi)^2}{r^2}=A(\psi).           \tag{8}
\]

Because `u dot grad B(psi)=0`, equation (8) is exactly
`u dot grad p=0` on that collar.  It is not imposed in `N_core`.

Select a smooth scalar cutoff `chi(p)` with `chi=1` on all pressure values
attained by the core and inner matching zone, with `supp chi'` contained
strictly in pressure values realized inside `S_loc`, and with `chi=0` flat
before the outer boundary.  Define

\[
 U=\chi(p)u,\qquad P(p)=P_0+\int_{p_0}^{p}\chi(s)^2ds.    \tag{9}
\]

The exact localization residuals are part of the candidate specification:

\[
 \operatorname{div}U=\chi'(p)u\mathbin\cdot\nabla p,
 \qquad
 (U\mathbin\cdot\nabla)U+\nabla P
 =\chi\chi'(u\mathbin\cdot\nabla p)u.                  \tag{10}
\]

Thus (8) only has to hold on `supp chi'`; the nonlocalizable core is unchanged
because `chi=1` there.  Success for A requires an actual global solution of
the mixed system (2), (7), (8) with a pressure-range separation that makes
(9) possible.  Elliptic Cauchy matching at `Sigma_1` is potentially
overdetermined and is a named load-bearing construction, not an assumption.
The outer zero extension must be `C-infinity`, including (4), and the support
must stay a positive distance from `r=0` unless full axis regularity is proved.

## Candidate B — a genuinely nonlocalizable Grad--Shafranov compacton

Seek a bounded smooth meridional domain

\[
 D\Subset\{(r,z):r>r_->0\}
\]

and a nonnegative physical streamfunction `psi` such that `psi>0` in `D`,
`psi` and all of its spatial derivatives vanish on `partial D`, and its zero
extension solves (2) distributionally and smoothly.  Take

\[
 F(0)=0,\qquad B(0)=0,
\]

and retain (7) on an inner stream-value interval.  Near `psi=0`, the proposed
absorption law is

\[
 B'(s)\sim b\,s\,[\log(s_*/s)]^p,
 \qquad b>0,\quad p>2,                                  \tag{11}
\]

while `F F'` is lower order; the simple choice `F(s)=lambda s` has
`F F'=lambda^2 s`.  The quotient `B'(s)/s` diverges, so the usual reduction
of a locally Lipschitz semilinear equation to one with bounded potential is
unavailable.  This only evades that particular unique-continuation mechanism;
it does not prove a compactly supported solution.

The one-dimensional edge model fixes the relevant exponent rather than merely
suggesting `exp(-1/d)`.  Put

\[
 \alpha=\frac{2}{p-2},\qquad
 \psi(d,y)\sim s_*\exp[-c(y)d^{-\alpha}]                 \tag{12}
\]

for inward normal distance `d` and boundary coordinate `y`.  With
`L=log(s_*/psi)=c d^{-alpha}`,

\[
 \partial_{dd}\psi
 \sim \alpha^2c(y)^{2-p}\psi L^p.                     \tag{13}
\]

The leading balance of (2) and (11) therefore requires

\[
 \alpha^2c(y)^{2-p}=b\,r(y)^2,                         \tag{14}
\]

with the signs in (2), because both leading terms are negative in that
equation.  For `p=4`, `alpha=1` and (12) reduces to the motivating
`exp[-c(y)/d]` profile.  The
radial metric term, tangential derivatives, boundary curvature, `F F'`, and
all subleading logarithmic orders remain to be solved, not discarded.

Candidate B earns a stationary compact field only if it constructs all of
the following rows:

1. a global positive solution of (2) on `D` with the free boundary and flat
   jet (12)--(14), rather than a one-dimensional formal asymptotic;
2. a positive nondegenerate interior critical point `(r0,z0)`, with definite
   meridional Hessian and `F(psi_c)/r0 != 0`;
3. the literal inner constant-curl identities (7), connected through smooth
   positive-label profiles to the non-Lipschitz edge law;
4. `C-infinity` zero extension of `u`, `omega`, and the physical pressure
   (4).  In particular, although abstract `B''(s)` may diverge at zero,
   `B(psi(x))` and every required spatial derivative must be flat;
5. no distributional vortex sheet or pressure traction at `partial D`;
6. support separated from the axis, or the exact cylindrical parity and
   vanishing conditions needed for smooth axis extension;
7. an admissible variational/free-boundary construction for (5), including
   compactness, nontriviality, regularity and exclusion of collapse to the
   zero critical point;
8. persistence of a positive-radius invariant tube and its nonzero flux
   twist on the constructed physical field.

Ordinary smooth or locally Lipschitz label laws with `F(0)=B'(0)=0` remain a
negative control: where their right side divided by `psi` is bounded, a
solution that vanishes on the exterior open set falls into the standard
bounded-potential unique-continuation class.  The non-Lipschitz law (11) is
a deliberate representation change, not permission to identify compact
vorticity with compact velocity.

## Compact assembly, density and physical costs

If either candidate supplies one `C-infinity` compact velocity and pressure
template supported in a bounded toroidal region, separated rigid copies have
zero pointwise cross velocity and pressure gradients.  Their sum is therefore
an exact stationary Euler field.  For `N` copies in a cell of physical volume
`V_cell`, with exact tube volume `V_T`, total measured ring coefficient
`J_T>0`, positive tag fraction `f_tag`, and whole-law detector Gram
`g_law>0`, the required densities are

\[
 \phi_T=\frac{NV_T}{V_{\rm cell}}>0,\qquad
 j_{\rm vol}=g_{\rm law}f_{\rm tag}
              \frac{NJ_T}{V_{\rm cell}}>0,              \tag{15}
\]

at one fixed finite packing.  The full constant material mass density remains
`rho` in the quiescent exterior.  The exact kinetic-energy density is

\[
 e_{\rm kin}=\frac{N}{V_{\rm cell}}
       \frac{\rho}{2}\int_{\operatorname{supp}U}|U|^2dx, \tag{16}
\]

with no Biot--Savart or pressure cross term omitted.  A positive tube count
without positive `J_T`, or a phase mass divided by a separately chosen volume,
does not meet the normalization requirement.

The actual 0250/full-parent history, action and current transfer remains a
separate subsequent row.  It must be performed on this exact compact
background with the full Euler/Lin pressure projection and the same density
normalization (15)--(16).  Baldi's particle action-angle chart can guide a
fixed compact carrier but does not supply the field-changing response.

## Frozen selection criteria

After central/schema activation, compare A and B in this order:

1. exact stationary Euler and pressure closure, with the 0252 core hypothesis
   avoided rather than contradicted;
2. a genuinely nonzero nondegenerate compact Euclidean core and positive
   radius/twist margins;
3. `C-infinity` compact **velocity**, vorticity and pressure, including edge
   and axis regularity and no distributional sheet;
4. a complete construction theorem for the mixed collar or non-Lipschitz free
   boundary, not a formal edge profile or bare existence citation;
5. exact disjoint assembly and fixed positive physical densities
   (15)--(16), including all ambient mass;
6. compatibility with 0250's fixed-background response, action and current
   targets with finite source/observation costs;
7. fewer new label functions, matching parameters and unproved inverses, with
   natural compatibility with accepted Euler invariants.

Candidate A is provisionally tried first because the exact localization
identity would make compact assembly immediate after one successful mixed
matching theorem.  Candidate B is the materially different representation
when the outer speed constraint makes that matching overdetermined.  This
preference is not a verdict.  The direct fixed exact Beltrami compact-core
background plus 0250 response remains active independently and is not made to
compete for permission to continue.

## Typed licenses and analytic oracle

`requires`: the reviewed 0211/0215 nonzero-core geometry at its exact scope;
the reviewed 0252 exclusion only for a pressure-localizable nondegenerate
core; the exact stationary/localization propositions verified after schema
activation; and the unchanged 0250 response objective as a later consumer.

`pass_licenses`: if A succeeds, an exact compact-velocity stationary Euler
template whose nonlocalizable core is untouched and whose cutoff operates
only in a matched localizable collar; if B succeeds, an exact nonlocalizable
Grad--Shafranov compacton with smooth pressure, nonzero inner constant-curl
core, exact disjoint packing and positive density.  Either license remains a
stationary geometry/density supplier until the separate dynamical transfer is
earned.

`does_not_license`: treating (11)--(14) as global existence; calling compact
vorticity compact velocity; weakening 0252 or 0211; importing passive
particle frequencies as Euler/Lin action; arbitrary knot types; nonlinear
finite-amplitude stability; parent completion; claim promotion; or scientific
exhaustion.

`maximum_verdict`: exact analytic construction of one of the two compact
stationary carriers, with all PDE, pressure, edge, core, regularity,
variational, assembly and density rows above, plus a precise interface for
0250.

`failure_scope`: the outer-only localizable matching route or the stated
non-Lipschitz Grad--Shafranov compacton class.  A failure of either route does
not refute the other, the fixed Beltrami-background route, or the parent
obligation.

`unlocks`: independent review of the new stationary compact supplier, then
0250/full-parent response, action/current and joint normalization on that same
field.

The strongest oracle is exact PDE and variational analysis: cylindrical
substitution of (1)--(4); a proved mixed-system or free-boundary existence and
regularity theorem; explicit boundary asymptotics with all lower orders;
smooth zero-extension and distributional residual checks; nondegenerate-core
and twist bounds; and direct density/action normalization.  A symbolic script
is justified only for a concrete sign, edge exponent, variational derivative,
or residual uncertainty.  No production numerical remainder is currently
specified, so the small-ratio numerical prescriptions do not bind this frozen
design.

`route_verdict: pending central/schema activation and substantive source/PDE
audit`.

`evidence_scope: preregistered split-localization and non-Lipschitz compacton
candidate universe, exact equations, variational target, edge scaling,
regularity obligations and density contract only`.

Next coordination achievement: the coordinator in `w3:p1` adds the central
0253 manifest entry and schema receipt.  Only then may this pane open the
located source bodies or advance either candidate beyond this frozen contract.

## Append-only activated result

The central entry and successful schema receipt (`WORKFLOW VALID: 270 claims`,
exit zero) activated the substantive work.  See `source-transfer.md` for the
primary-source scope and `construction.md` for the derivation.

Candidate B now has an exact physically `C-infinity` straight compacton with a
nonzero nondegenerate constant-curl core.  Its logarithmic label derivative is
nonsmooth at label zero while velocity, vorticity and full pressure are flat in
space.  The exact finite-radius equation has been isolated, together with a
nonzero curvature residual that rules out literal circular bending and a
regular defining-function formulation of the remaining degenerate
free-boundary inverse.  This does not yet give compact Euclidean support.

For Candidate A, analytic continuation from a localizable collar to the
retained core is incompatible with 0252; only a genuinely nonanalytic
`C-infinity` mixed matching remains.  This is an analytic-subroute refutation,
not a refutation of the smooth route.

`route_verdict: Candidate A blocked at nonanalytic mixed collar matching;
Candidate B blocked at the finite-radius degenerate free-boundary inverse after
an exact straight compacton and all edge/pressure identities were established`.

`evidence_scope: exact analytic precursor, exact finite-curvature residual,
and precise minimum repair; no compact Euclidean carrier, density supplier,
response/action/current join, parent completion, or exhaustion conclusion`.

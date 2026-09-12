# P253/0102 — two-label forced axial lock and charged flavor bridge

## Frozen objective and authority

This root-owned supplier attempt continues P6 after the independent final
P253/0101 review. It asks whether the same charged Cao carrier can support a
material axial density whose ratio to the electromagnetic charge tag is a
single constant, and whether that structure can supply an actual input to the
conditional two-state flavor kinematics reviewed in P253/0099.

The authority boundary is accepted release `v0.183.0` plus these individually
reviewed inputs at their exact scopes:

- P253/0068 with P253/0071: the transported real charge tag, conserved
  Euler--Maxwell current, action, constraints and local smooth flow;
- P253/0077 with P253/0081 and P253/0080 with P253/0084: the subluminal charged
  axisymmetric no-swirl Cao branch, Lorenz-gauge Maxwell elimination, modified
  Grad--Shafranov first integral, regular-band material tag and exact local
  branch. No dynamical persistence is imported;
- P253/0097 with P253/0099: the exact Euler helicity axial Galilean current and
  conditional two-state propagation/coherence/shared-action kinematics. No
  physical Cao flavor mechanism, Lorentz chiral current or P6 result is
  imported;
- P253/0100 with P253/0101: the unforced and forced Ertel continuity currents,
  the exact forced-lock compatibility equation, zero regular Cao density,
  Cao magnetization superpotential, and punctured-domain equilibrium phase.
  In particular, the independently reviewed moving material defect and force
  selection remain open and are not treated as premises; and
- P253/0085 with P253/0086: the strong-saddle and Maxwell-continuum boundary.
  It prevents a steady identity or prepared mode from being relabeled as P2.

No empirical comparator is used. The target is an exact theorem/construction
competition, so no numerical selection criterion is manufactured. Exact
algebra or symbolic checks may expose signs and projection identities, but
they cannot prove the moving-domain, free-boundary, Fredholm, P2, P4 or P6
claims.

## Geometry, fields and conventions

Work in the comoving frame of a subluminal charged Cao relative equilibrium.
On its regular axisymmetric no-swirl core,

    W=r^(-1) grad P cross e_theta,
    omega=r zeta e_theta=zeta partial_theta,
    chi_c=F(P),
    f=-(g/rho_m) chi_c(grad Phi+H grad P),                 (1)

and the reviewed modified Grad--Shafranov first integral is

    zeta=h(P)+(g/rho_m)[chi_c'(P)Phi-chi_c(P)H].          (2)

Introduce a distinct circle-valued material phase `Theta` on a time-dependent
punctured domain `M_t=R^3\Gamma_t`. Its closed one-form `a=dTheta` must be a
globally defined advected one-form even when no real global lift of `Theta`
exists. Locally,

    D_t Theta=0,              q_A=omega dot grad Theta.   (3)

At the base axisymmetric no-swirl equilibrium, `Theta=theta` and `q_A=zeta`.
The geometric azimuth is not silently used for nonaxisymmetric perturbations:
`D_t theta=u_theta/r`, so beyond the base the independently advected phase and
moving defect are part of the state and domain.

The electromagnetic tag `chi_c` and the axial phase `Theta` are different
objects. Therefore the P253/0101 closed-line obstruction for
`omega dot grad chi=lambda chi` does **not** apply to
`omega dot grad Theta=lambda chi_c`; no differential equation
`d chi_c/dtau=lambda chi_c` follows.

Both labels are transported. The forced Ertel identity gives

    D_t q_A=(curl f) dot grad Theta,
    chi_c D_t lambda=(curl f) dot grad Theta,
    lambda=q_A/chi_c                                      (4)

only where `chi_c` is nonzero. On the equilibrium, direct cylindrical
calculus must rederive

    (curl f) dot grad theta
      =(g/rho_m) W dot grad[chi_c' Phi-chi_c H].          (5)

Equations (2) and (5) imply that the variable ratio

    lambda(P,alpha)=zeta/chi_c                            (6)

satisfies (4) automatically wherever the quotient is defined. This is a
steady forced-lock identity, not a selection theorem.

## Load-bearing support boundary

The reviewed charge tag is supported strictly inside a regular band where the
vorticity stabilizer locks it. The axial density `q_A=zeta` occupies the full
Cao core. Hence `zeta/chi_c` is undefined off the tagged band, and the existing
construction cannot satisfy a global pointwise constant lock.

Every result must choose one of two honest domains:

1. a **band-local lock**, with the tag bounded away from zero on the declared
   sub-band and with the induced boundary/extension current explicitly
   included; or
2. a **redesigned global tag**, accompanied by a new stabilizer, center, edge,
   zero-set and free-boundary proof. Extending the tag into isochronous or
   vanishing-vorticity regions is not inherited from the reviewed branch.

The axial density is never identified with electromagnetic charge
`g chi_c`. Its integral remains circulation-derived and continuous unless a
separate action/character mechanism is proved. Since `q_A` is a pseudoscalar
for a true scalar phase, a nonzero constant `lambda_0` is itself a
pseudoscalar datum: its sign needs the oriented defect/carrier handedness and
cannot be presented as an `O(3)`-scalar coupling.

## Candidate routes and one-verdict contract

### Route A — variable two-label equilibrium lock

Rederive (4)--(6), including cylindrical factors, parity, the zero set of the
tag and band-boundary flux. Determine the axial continuity current

    (q_A, q_A u-f cross grad Theta)                       (7)

and its relation to the separately reviewed magnetization current. The route
earns `established` only for the exact local/band equilibrium identity it
proves. It does not earn a constant vector/axial ratio, Lorentz chirality,
flavor, quantization or persistence.

### Route B — fixed-profile constant-ratio obstruction

A constant ratio `lambda_0` requires the pointwise equation

    h(P)+(g/rho_m)[chi_c' Phi-chi_c H]=lambda_0 chi_c(P)  (8)

on the tagged domain. Merely showing the bracketed correction is a streamline
function gives a variable `lambda(P)`, not a constant.

Use `tau=g^2` and the reviewed branch parity. At `g=0`, choose on the declared
band

    chi_0(P)=zeta_0(P)/lambda_0.                          (9)

Write `Phi=g Phi_hat+O(g^3)`, `H=g H_hat+O(g^3)`. For each regular contour of
`P`, define the action-period average `Pi_P` and zero-mean projection
`Q_P=I-Pi_P`. The first fixed-profile obstruction is

    Q_P[chi_0'(P) Phi_hat-chi_0(P) H_hat]=0              (10)

on every tagged contour. Streamline-average rows may be absorbed by the
declared `h_2(P)`, `chi_2(P)` and `lambda_2`; they cannot cancel (10).

Derive `Phi_hat,H_hat` from the actual comoving Maxwell Green system with its
affine Coulomb tail, and make (10) a reusable exact or rigorously bounded API.
If one contour has a nonzero projected response, the fixed-profile
constant-lock route is `refuted` with that Maxwell mechanism. A missing
evaluation is `blocked` at that named response, not a refutation.

### Route C0 — exponential straight-column cancellation candidate

If Route B exposes a nonzero response, first test whether a materially
different vorticity law can cancel it without treating a profile as a fitting
knob. In the global straight-column reduction, impose
`zeta_0'(P)=a zeta_0(P)` with a nonzero constant `a` and
`chi_0=zeta_0/lambda_0`. Derive the fluid-Poisson, electric-dipole and Ampere
identities in the same physical normalization and test whether their first
speed responses cancel. Audit circulation, decay, energy and localization:
an everywhere positive exponential law cannot be silently relabeled as the
compact Cao free-boundary profile. This route receives its own verdict and,
if it cancels, supplies a kernel candidate for the compact-profile matching
problem rather than a charged-particle result.

### Route C1 — annular source-matching control

On a tagged annulus where the exponential law holds, derive the homogeneous
matching defect `A s+B/s`. Test whether source variations supported strictly
inside and outside the annulus independently control its two coefficients
through the exact radial Green moments. This is a source-level range theorem;
it does not by itself realize those variations as admissible Cao profiles or
place them on one material leaf.

### Route C2 — compact profile/carrier Fredholm control

If Route B exposes a nonzero response, vary the profile as an explicit
competing control. The zeroth-order constant lock already fixes
`chi_0=zeta_0/lambda_0`, so an order-`tau` tag correction cannot change the
zero-mean obstruction (10). A valid control must therefore vary the
uncharged vorticity/profile and tag together at order one, preserving
`chi_0=zeta_0/lambda_0`, or introduce a separately derived phase geometry that
changes `q_A`. Derive the corresponding linearized map to the contourwise
zero-mean part of `chi_0' Phi_hat-chi_0 H`, with circulation, mean radius,
tag integral/leaf, stabilizer, charge and center rows. Test its range, kernel,
compactness/Fredholm character and compatibility with the uncharged Cao and
Maxwell/Grad--Shafranov/free-boundary maps. A profile function is not a hidden
fitting knob: each allowed variation needs a conserved distribution or
Casimir row and must remain on the same declared material leaf.

This route receives one verdict independent of Routes B, C0 and C1. Success constructs a
constant-ratio tag/carrier pair on the declared band or redesigned global
domain; failure must identify the cokernel/range obstruction and activate a
materially different label or defect coupling.

### Route D — global advected circle phase and moving defect

Construct `Theta`, `a=dTheta` and the moving puncture `Gamma_t` as an actual
Euler semidirect-product state. Prove its advection, branch cut independence,
inner/outer flux law, finite action/energy topology, regularization of the
defect core, and compatibility with nonaxisymmetric perturbations. Decide
whether a band-local axial current extends globally or terminates in a
physical boundary current.

For the base identification `Theta=theta`, `Gamma_0` is the entire noncompact
symmetry axis. It is not a compact defect localized on the toroidal carrier,
and an unweighted phase-gradient energy for `dtheta` diverges without a
declared core and outer treatment. Replacing the axis by a compact ring defect
requires a new phase map; neither `Theta=theta` nor `q_A=zeta` then follows
automatically.

This route cannot inherit P2 from stationarity. It must state its dependence
on the active P253/0095 persistence problem and cannot call an equilibrium
defect an enduring particle carrier.

### Route E — same-carrier flavor bridge

Only after a nontrivial axial lock and moving material phase exist, couple the
vector electromagnetic current and axial current to the conditional
P253/0099 two-state sector. Require:

- a shared, physically normalized P4 action space;
- a derived noncommuting propagation and interaction basis;
- momentum-dependent branches and the reviewed coherence/relativistic phase
  ledger;
- an electron-linked endpoint current on the same carrier;
- parity and axis behavior appropriate to a Lorentz chiral representation,
  or an explicit declaration and pricing of the required extension; and
- compatibility with P2 rather than a prepared finite-time beat.

The axial lock alone is not V--A, flavor, a weak current, a spinor projector,
Born statistics, reset, fermionic exchange, an electron or a neutrino.

## Exact checks and continuation ladder

The first exact verifier must derive, rather than print:

1. the forced two-label quotient identity (4);
2. the cylindrical sign and factor in (5);
3. the automatic variable equilibrium lock from (2);
4. the distinction between (6) and the constant equation (8);
5. the contour projector identities `Pi_P Q_P=0`, `Q_P^2=Q_P` and the
   perturbative necessity of (10); and
6. the tag-support counterexample showing why division through `chi_c=0` is
   invalid.

After any route verdict, continue in-run. A nonzero fixed-profile response
activates Routes C0--C2. A compact profile-map cokernel activates a redesigned full-core tag,
defect-localized interaction, or materially different carrier. A successful
constant lock activates Route D and then Route E; it does not complete P6 by
itself. A moving-defect obstruction activates an equivalent globally defined
compact one-form/higher-form representation only if its action coefficient and
new structure are explicitly priced.

## Maximum verdict and exclusions

The maximum attempt verdict is a same-carrier, domain-correct two-label axial
lock with an evaluated fixed-profile Maxwell response and, if needed, a proved
profile-control construction, plus an exact statement of what remains for the
moving material defect and flavor bridge.

This attempt cannot establish LP2, LP4, LP5 or LP6; a Lorentz chiral weak
current; universal or quantized action/charge; a nonzero physical weak
coupling; neutrino masses or empirical mixing; Born outcomes; reset;
fermionic exchange; an electron; a neutrino; or parent completion. P253/0095
remains an active independent P2 route and is not consumed as solved.

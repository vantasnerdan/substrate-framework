# Independent fixed-theorem review of P253/0007 particle-facing calculus

Reviewer: `particle-balance-review`, a non-author reviewer coordinated through
Herdr pane `w3:p1`.  I did not author or implement the reviewed theorem,
module, or tests.  Prior review work in `0004` did not contribute to `0007`.

## Frozen transaction

This is one equation-level substantive pass over:

- `0007/particle-facing-calculus.md`:
  `fdb47e5129c797ae02a0c6da11509161e7863181845bbe17c77668fc0a5aedc8`;
- `src/substrate_framework/euler_impulse.py`:
  `e194baa950dba46140bce4f94c52bef621f6162bb61e0ea61b1b27d075caf85e`;
- `tests/test_euler_impulse.py`:
  `39df617a91fb4b804424550c1d14dc08a1cf2323e4fd53820c276fb7e7f46ffd`;
- the recorded focused-test output and exit token:
  `ea9106f92819092b5d3e8ea219b43794355bb75051142f4e0c91ab4c03bf17f9`
  and
  `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

The accepted base is `v0.183.0` at
`b6fc902a0942d07996f12a81028fbd3f7c909a43`, and the active campaign
checkpoint is `4620bf0e50d09722a0202143938c6f29f03c6a54`.  No claim ID, registry,
release, migration, or existing API consumer changes in this transaction.
The three-test receipt is regression evidence, not the source of the review
verdict.

## Verdict and strongest supported statement

`route_verdict: established as stated at compact-vorticity asymptotic,
finite-tag no-swirl, orientation-orbit, and bare-background symbol scope,
with two non-load-bearing wording clarifications`

For a smooth compactly supported divergence-free vorticity on `R^3`, the
decaying Hodge velocity has the impulse dipole (3), with the stated uniform
differentiated Newton-kernel remainder.  For two fixed translated cores at
sufficient separation, the positive Euler kinetic cross energy has the exact
leading sign and anisotropic `d^-3` coefficient (5), with the moment remainder
(6).  This is an energy asymptotic, not a mechanical force law.

For a centered axisymmetric no-swirl Euler field, helicity density vanishes
pointwise and centered internal angular momentum vanishes on every finite
axisymmetric tag.  Whole-space angular momentum is zero only when its integral
is actually defined; finite kinetic energy alone does not supply the needed
weighted absolute convergence.

For a directed axisymmetric carrier with nonzero impulse, rigid rotations have
stabilizer `SO(2)` and the orientation-only orbit is `SO(3)/SO(2)=S^2`.
Consequently its orientation-only `2pi` rotation loop is contractible.  This
does not compute the topology of the full Euler configuration orbit or produce
quantum spin/statistics.

The exact transverse linear Euler symbol about rest is zero.  About a uniform
background it is only the Galilean convective shift `omega=U_0 dot k`.  For
`U_0 != 0`, this statement belongs either to a periodic domain or to perturbation
theory relative to a homogeneous infinite-total-energy `R^3` background; it is
not an additional finite-total-energy Euclidean carrier.

These results falsify only the named direct identifications: the bare fixed-core
kinetic cross term is not a Coulomb `d^-1` potential, the no-swirl orientation
orbit alone does not supply a nontrivial spinorial loop, and the rest/uniform
bare Euler background does not supply a finite-speed transverse wave cone.
They do not establish or refute an actual particle mechanism globally.

## Objective bridge

The checked predicates are exact moment cancellations, a Newton-kernel Taylor
coefficient with an explicit remainder, an exact finite-tag symmetry integral,
an exact homogeneous-space calculation, and a Fourier-symbol elimination.
They imply the four scoped mathematical statements above for the declared
objects.  Their upstream licenses are smooth compact divergence-free vorticity
and decaying Hodge recovery for the impulse calculation; a centered
axisymmetric no-swirl field and finite axisymmetric tag for the angular
calculation; a nonzero directed impulse for the rotation orbit; and a decaying,
periodic, or relative-background Fourier perturbation for the linear symbol.

No numerical remainder, fitted parameter, empirical comparator, force
observable, collective action, full-orbit topology, or quantum Hilbert space is
present.  The maximum verdict is therefore exact symbolic verification of
these predicates and direct-identification tests.  The parent particle
obligations receive constraints and next routes, not completion licenses.

## Impulse moment and epsilon audit

Let

    I_k=(1/2) epsilon_klm integral y_l omega_m dy.

Compact support and `div omega=0` give, by integrating
`div(y_i omega)` and `div(y_i y_j omega)`,

    integral omega_i dy=0,
    integral(y_i omega_j+y_j omega_i)dy=0.

Thus the first-moment matrix is antisymmetric.  Contracting its antisymmetric
part with the definition of `I` gives the sign convention

    M_ji=integral y_j omega_i dy=epsilon_jik I_k.

For `G(x)=1/(4pi|x|)`, the first nonzero term of
`A_i(x)=integral G(x-y)omega_i(y)dy` is

    -partial_j G(x) M_ji=(I cross x)_i/(4pi|x|^3).

Taking the curl independently yields

    u(x)=[3n(n dot I)-I]/(4pi|x|^3).

This fixes both the epsilon sign and the `4pi` normalization.  A circular
vortex loop with `I_z=pi Gamma R^2` then gives the familiar axial leading term
`Gamma R^2/(2d^3)`, agreeing with the test fixture.  The filament is used only
for normalization; it is not substituted for the smooth carrier hypothesis.

## Differentiated Newton remainder and boundary convergence

For `|x|>2a`, every point on the Taylor segment from `x` to `x-y` stays at
least `|x|/2` from the origin.  The Newton derivatives obey

    |D^m G| <= C_m |x|^(-1-m)

uniformly on that segment.  Taylor expansion through first order and
`integral omega=0` therefore gives

    |A-A_dip| <= C M2/|x|^3,
    |curl(A-A_dip)| <= C M2/|x|^4,

where `M2=integral |y|^2|omega(y)|dy`.  Differentiating the integral is
legitimate because the source is smooth and compactly supported.  This is the
claimed uniform remainder in (3), not an unevaluated asymptotic symbol.

The same estimates give `A_1=O(r^-2)` and `u_2=O(r^-3)`.  Hence the boundary
term in

    integral curl(A_1) dot u_2
      =integral A_1 dot curl(u_2)+boundary

is `O(r^2 r^-2 r^-3)=O(r^-3)` and vanishes at infinity.  The volume cross
energy is absolutely convergent at infinity, while smoothness removes inner
boundary terms.  Equation (4) consequently has the stated sign.

## Mixed Newton expansion and cross-energy sign

Place the two core coordinates at `y` and `z` and let `d` join their centers.
The exact double-integral kernel is `G(d+z-y)`.  Zeroth and linear Taylor terms
vanish because each vorticity has zero zeroth moment.  Pure quadratic terms
also multiply a zero zeroth moment.  The surviving mixed term is

    -partial_jk G(d) M1_ji M2_ki.

Using `M_ji=epsilon_jia I_a` and `Delta G(d)=0` gives

    -partial_jk G M1_ji M2_ki
      =I2_j partial_jk G I1_k
      =[3(I1 dot n)(I2 dot n)-I1 dot I2]/(4pi d^3).

Multiplication by positive physical density proves (5) as the cross term in

    (rho/2) integral |u_1+u_2|^2
      =E_1+E_2+rho integral u_1 dot u_2.

This sign is therefore the positive kinetic-energy convention, not the sign
of a maintained-current magnetic potential.  Parallel impulses on-axis give a
positive coefficient, whereas parallel impulses in an equatorial placement
give a negative coefficient, matching the exposing test.

If the supports lie in balls of radii `a_1,a_2` and
`d>2(a_1+a_2)`, the third derivatives of `G` are uniformly `O(d^-4)` on the
full Taylor segment.  The integral remainder is bounded by

    C rho d^-4 integral integral |omega_1(y)||omega_2(z)|
                              (|y|+|z|)^3 dy dz,

which is (6).  No derivative of this remainder with respect to a collective
separation coordinate is used.  Promoting `-grad_d E_cross` to force would
first require a carrier manifold, its symplectic/kinetic reduction, a statement
of what is held fixed, and control of shape and ambient contributions.  None is
silently supplied here.

## No-swirl helicity and angular momentum

For (7), `u` lies in the meridional `e_r,e_z` plane while `omega` is purely
azimuthal, so `u dot omega=0` pointwise.  Direct cylindrical cross products
give

    x cross u=(z u_r-r u_z)e_theta.

On a finite axisymmetric tag, the coefficient is independent of the azimuth
and the full integral of `e_theta` is zero.  Axisymmetry also puts both the tag
centroid and centroid velocity on the symmetry axis.  Replacing `x,u` by their
centered versions therefore leaves an azimuthal integrand and proves zero
centered internal angular momentum exactly.  Density only multiplies this
zero and no density factor is missing.

For a compact-vorticity field with nonzero impulse, the generic velocity tail
is `O(r^-3)`.  Kinetic energy converges because its radial tail behaves as
`integral r^2 r^-6 dr`, but the absolute angular-momentum integral can behave
as `integral r^2 r r^-3 dr=integral dr`.  This directly confirms the source's
warning: finite energy does not license a whole-space angular-momentum value.
The finite-tag theorem and pointwise helicity need no such convergence.

## Rigid orientation orbit

Impulse transforms as a spatial vector under rigid rotations.  If a rotation
stabilizes a field with nonzero impulse `I`, it must fix the directed unit axis
`n=I/|I|`, so its stabilizer is contained in rotations about `n`.  Axisymmetry
supplies every such rotation.  The rigid stabilizer is therefore exactly
`SO(2)` and the orbit is `SO(3)/SO(2)=S^2`.

Since `pi_1(S^2)=0`, every orientation-only loop produced by a `2pi` spatial
rotation is contractible (rotation about the carrier axis is already a
constant orbit loop).  This is an exact statement about the rigid orientation
quotient only.  Internal phase, swirl, relabeling, vorticity-knot configuration
space, Berry data, prequantum bundles, and the topology of the full Euler orbit
are absent.  The calculation therefore blocks an orientation-only shortcut to
spin half but neither proves integer quantum spin nor globally excludes a
spinorial collective sector.

## Rest and uniform Euler transverse symbol

Linearize physical-pressure Euler about a constant background `U_0`:

    (partial_t+U_0 dot grad)v=-grad(pi)/rho,
    div v=0.

For a Fourier mode with `k != 0`, taking the divergence gives
`|k|^2 pi_hat=0` under the standard periodic, decaying-perturbation, or fixed
mean pressure normalization.  Hence `pi_hat=0` and the two-dimensional
transverse symbol is

    partial_t v_hat=-i(U_0 dot k)v_hat.

With the convention `exp(i(k dot x-omega t))`, this is
`omega=U_0 dot k`.  At rest the generator is exactly zero; for uniform flow it
is removed by a Galilean frame change.  It has no nonzero rest-frame wave cone
and does not reproduce `omega^2=c^2|k|^2`.

The domain qualifier matters.  `U_0=0` is compatible with the campaign's
decaying finite-energy Euclidean class.  A nonzero uniform `U_0` has infinite
total energy on `R^3`; its exact symbol is instead a periodic-domain statement
or a statement for finite-energy perturbations measured relative to that
homogeneous background.  The source's mathematics is correct, but this typing
should remain explicit in any promoted wording.

## API and focused-test audit

The module is a pure conditional symbolic API.  It implements the exact dipole
velocity and the positive kinetic cross coefficient `rho I_b dot u_a(d)`.
It rejects known zero separation and known nonpositive/nonfinite density, while
leaving undecidable symbolic reality, nonzero separation, support separation,
and remainder bounds as caller hypotheses.  Its documentation correctly says
that algebraic evaluation constructs no carrier, stability theorem, potential,
or collective mechanics.  The cross-energy expression is symmetric under
exchanging the impulses and reversing the separation, as required by the
Hessian kernel.

The three focused tests have bounded, exposing roles:

- differentiating `I cross x/(4pi|x|^3)` independently checks the velocity and
  circular-loop normalization;
- differentiating the Newton kernel twice, constructing both epsilon moment
  matrices, and retaining the mixed Taylor minus sign checks the cross-energy
  coefficient independently of the public helper;
- axial/equatorial sign probes expose a magnetic-potential sign substitution;
- zero impulse and representative invalid domains check the leading-term and
  API boundaries.

The no-swirl, homogeneous-space, and transverse-symbol results are analytic
theorem steps rather than module entrypoints; the calculations above audit them
directly.  The recorded output is `3 passed in 2.26s` with exit zero.  I did not
rerun an unchanged exact regression merely to duplicate that count.  No
small-ratio numerical skill prescription binds: the cross energy is an exact
asymptotic coefficient with an analytic remainder, no mechanical force is
computed, and no numerical floor or energy splitting is used.

## Findings and minimum correction

No load-bearing scientific defect was found.  Two non-blocking wording repairs
would make the already-supported scope maximally precise:

1. Type the `U_0 != 0` calculation explicitly as periodic or relative-background
   perturbation theory, because a nonzero uniform flow is not a finite-total-
   energy field on `R^3`.  The exact convective symbol and direct-identification
   verdict are unchanged.
2. Replace “two propagating photon helicities with `omega=+/-c|k|`” by “the two
   transverse Maxwell polarizations with wave cone
   `omega^2=c^2|k|^2` (each having positive/negative temporal branches).”
   Helicity labels and frequency signs are distinct.  This is terminology only;
   the zero/convective Euler symbol still excludes the direct bare-background
   identification.

Neither repair introduces a new hypothesis into the impulse, no-swirl, or
orientation calculations.  No additional theorem route or production test is
needed to sustain the fixed result.

## Four-axis decision

- Verification: `symbolic_verified` for the exact identities and controlled
  asymptotic; the three tests are regression evidence.
- Review: `audited` in this non-author pass.
- Compatibility: `native` for compact-vorticity, finite-tag, and quiescent
  Euler; the nonzero uniform background is `compatible_extension` only in its
  explicitly periodic or relative-energy domain.
- Epistemic: `proposed` campaign evidence; no claim is accepted or promoted.
- Relationship: bounded P3 interaction asymptotic and P4/P5 direct-route tests,
  without upward inheritance to any particle obligation.

## Result and remaining dependency

The fixed theorem route is established with the two bounded wording repairs
above.  The next P3 achievement is to construct a persistent same-carrier
interaction reduction with its actual symplectic/kinetic structure, controlled
shape and ambient terms, and a force observable if one exists.  The next
P4/P5 achievement is a physically nontrivial internal/swirl or collective
quantum sector on the same carrier together with a structured-background
autonomous propagation/current construction.  Full-orbit topology, quantum
statistics, mechanical force, all-time dynamics, electron/neutrino mechanisms,
and P2/P4--P7 completion remain open.

## Bounded correction check

The requested two wording repairs are closed at the amended
`0007/particle-facing-calculus.md` hash
`8f35441acdd64b0738a1dac0241db32ce82e28529de6970add9e2df710dd00ef`.
This correction check inspected only section 3 and its directly affected route
summary.

The amended statement now types constant `U_0 != 0` as a periodic-cell result
or perturbation theory relative to a uniform background, explicitly excluding
finite-total-energy `U_0` on `R^3`.  It also describes the absent target as a
transverse wave equation with `omega=+/-c|k|` and states separately that the
two transverse polarizations are distinct from the positive/negative temporal
branches.  These are exactly the two minimum corrections requested; neither
changes the zero/convective Euler symbol or broadens its verdict.

Equations (9)--(10) are unchanged.  The module remains at
`e194baa950dba46140bce4f94c52bef621f6162bb61e0ea61b1b27d075caf85e`
and the tests remain at
`39df617a91fb4b804424550c1d14dc08a1cf2323e4fd53820c276fb7e7f46ffd`.
The existing three-test receipt therefore remains applicable and was not
rerun.  No other finding, dependency edge, or scientific boundary was
reopened.

`post_correction_verdict: established as corrected at the exact scopes stated
above; both requested wording refinements are closed`

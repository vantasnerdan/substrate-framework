# Independent fixed-theorem review of P253/0001 material balances

Reviewer: `particle-balance-review`, Herdr worker pane associated with the
coordinator at `w3:p1`.  I did not author or implement the reviewed source,
module, or tests.

## Frozen transaction

This one substantive pass reviews the exact theorem and its additive API at
the following working-tree boundary:

- `0001/material-balances.md`, including the later-frozen section 5:
  `9171a613dc7626fa59360d73cda709e17f445bba31422c3ceee12125898dc5dc`;
- `src/substrate_framework/euler_material_balance.py`:
  `5ded2b458ccf46d7f32fe87618b96f414477c9069e4115aad342ffec7a84876b`;
- `tests/test_euler_material_balance.py`:
  `cb2ba5ba2c2a78ace4b3c75bd1891f03e807626b1f1c9500fd8a9cde118b0a51`;
- the scoped receipt and captured nine-test output:
  `8db82442d246837ff2fd9c805ab736aa23b9f5854f81d4b8279540d0a6eace23`
  and
  `2e34549c8e79426cd89c4c1dcae864debe7149dacd2995fe25ec6ec984648d34`,
  with recorded exit-zero token
  `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

The base release is `v0.183.0` at
`b6fc902a0942d07996f12a81028fbd3f7c909a43`; the campaign checkpoint is
`4620bf0e50d09722a0202143938c6f29f03c6a54`.  There is no proposed claim ID,
registry/release delta, or existing API consumer beyond the focused tests.
The test receipt is regression provenance, not the basis of this review.

## Verdict and strongest supported statement

`route_verdict: established as stated at exact smooth-local-Euler and
initial-response scope`

For constant physical density `rho>0`, every smooth incompressible Euler
solution on `R^3` and every bounded smooth material tag obey the centered
mass, centroid, first-momentum, second-mass, spin, kinetic-covariance, and
kinetic-energy balances (1)--(5) throughout the common smooth existence
interval, under the stated decay or finite-tag boundary hypotheses.

There also exists a natural class of nonzero `C-infinity`, compactly supported,
finite-energy divergence-free initial data

    u_0(x)=f(|x-D|) n cross (x-D),  f(s)=g(s^2),

for which the decaying physical pressure has exactly the exterior quadrupole
(9).  A separated initially resting spherical material tag has exactly the
centroid acceleration (10) and trace-free shape acceleration (12).  Comparing
with the zero Euler solution gives identical position and local velocity data,
including all local spatial velocity derivatives, but different exact initial
centroid and shape accelerations.  Thus no autonomous closure using only those
local tag position/velocity/spin/shape moments can govern the unrestricted
smooth finite-energy local Euler class.

For two such radial swirls with disjoint velocity supports, the pressure source
is additive at `t=0`.  Spherical material tags centered at the respective swirl
centers and disjoint from the other pressure source have the exact external
centroid acceleration (13), including the collinear and equatorial cases in
(14).  This is an exact initial, anisotropic, circulation-even pressure
interaction.  Ambient momentum remains part of the full Euler balance.

This result does not establish an invariant finite-dimensional closure, an
enduring carrier, an all-time two-object interaction, a reciprocal Newtonian
pair force, a Coulomb law, charge sign, particle mass, electron or neutrino
identification, stability, or global three-dimensional regularity.

## Independent centered-balance derivation

Let `dm=rho dx`, `r=x-X`, `c=u-V`, with the definitions in the source.  Reynolds
transport on a material tag and incompressibility give

    D_t r=c,
    D_t c=-grad(p)/rho-Xddot,
    integral r dm=integral c dm=0.

Therefore `Mdot=0`, `Xdot=V`, and

    M Xddot=-integral grad(p) dx=F.

Differentiating before using the two zero means gives

    Idot=B+B^T,
    Bdot=T-integral grad(p) r^T dx=T-P,
    Iddot=2T-P-P^T.

There is no missing mean-acceleration term: it multiplies `integral r dm`.
Similarly the `c cross c` and mean-acceleration contributions vanish, so

    Ldot=-integral r cross grad(p) dx.

For `T=integral cc^T dm`, the two mean-acceleration terms multiply
`integral c dm`; hence

    Tdot=-integral [grad(p)c^T+c grad(p)^T] dx.

Taking half the trace proves the internal kinetic-energy rate.  Adding
`V dot F` yields

    E_tag_dot=-integral u dot grad(p) dx
             =-integral_boundary p u dot n dA.

All factors of `rho` are correct: physical acceleration is `-grad(p)/rho`,
whereas `dm=rho dx`, so the pressure moments in the continuum balances carry
no extra density.  The finite API correctly represents this as
`sum m grad(p) r^T/rho` and `sum -m grad(p)/rho`.

Integration by parts gives, with outward tag normal,

    P_ij=integral_boundary p n_i r_j dA
         -delta_ij integral p dx,
    Ldot=-integral_boundary p r cross n dA.

These signs reproduce (2), (5), the scalar virial, and invariance under an
additive pressure constant.  The dimensions also close:
`[P]=[T]=mass*length^2/time^2`, `[F]=mass*length/time^2`, and every reported
rate has the dimension of its differentiated moment.

## Independent Cartesian source and Green calculation

Rotate to `n=e_z`, write `y=(x,y,z)`, `s=|y|`, and
`u=(-yf,xf,0)`.  Direct differentiation gives `div(u)=0`.  In
`tr[(Du)^2]=sum_ij partial_j u_i partial_i u_j`, the terms quadratic in
`f'` cancel, leaving

    tr[(Du)^2]=-2f^2-2 f f' (x^2+y^2)/s
              =-2f^2-2s f f' sin^2(theta).

Since `sin^2(theta)=2[1-P_2(cos theta)]/3`, the only spherical harmonics are

    S_0=-2f^2-(4/3)s f f',
    S_2=(4/3)s f f'.

Flat compact support (indeed, only the vanishing endpoint is needed for these
two integrations) gives

    integral_0^a s^2 S_0 ds=0,
    integral_0^a s^4 S_2 ds=-(10/3)J,
    J=integral_0^a s^4 f^2 ds>0.

For physical pressure, divergence of Euler gives
`-Delta p=rho tr[(Du)^2]`.  The exterior Green coefficient is
`rho/(2l+1) r^(-l-1) integral s^(l+2)S_l ds`.  Thus the monopole vanishes,
the `l=2` coefficient is `-2rho J/3`, and there are no higher exterior
multipoles.  Rotational covariance gives exactly (9).

The profile is a legitimate local-Euler datum when `g` is a smooth bump of
`s^2`: it is smooth at the origin, flat at its support edge, compactly
supported, divergence free, and belongs to every Sobolev space.  Standard
local well-posedness on `R^3` for divergence-free `H-infinity` data supplies a
smooth solution on a nonzero interval.  The cited Tao notes state the
normalized (kinematic-pressure) formula and local-existence corollary; multiplying
that pressure by the constant physical density gives the convention used here.
No global-existence premise enters the initial-acceleration argument.

## Harmonic-ball and closure calculation

Inside a ball separated from the swirl support, `p` and each `partial_i p` are
harmonic.  The mean-value identity therefore gives the centroid acceleration
as the center value of `-grad(p)/rho`, proving (10).

For a harmonic `h`, decompose its Taylor series into homogeneous harmonic
polynomials.  In the ball average of `x_j partial_i h`, orthogonality removes
every degree except the quadratic part of `h`; isotropy
`average(x_j x_k)=b^2 delta_jk/5` then gives

    average_B[x_j partial_i h]=(b^2/5) partial_ij h(0).

Direct differentiation of (9) for `D=d e_z` yields

    grad p(0)=-2rho J e_z/d^4,
    Hess p(0)=rho J diag(4,4,-8)/d^5.

At `t=0`, the tag has `T=0`, so substituting the harmonic-ball identity into
`Iddot=-P-P^T` gives exactly (12).  Its trace vanishes as required by
harmonicity, while its nonzero trace-free part distinguishes the two Euler
solutions even when centroid data alone are omitted from a proposed closure.

## Two-disjoint-swirl interaction calculation

For smooth compactly supported `u_1,u_2` with disjoint supports, the derivative
supports are also disjoint.  Equivalently the smooth cross tensors
`u_1 tensor u_2` and `u_2 tensor u_1` vanish identically.  Thus the quadratic
pressure source of `u_1+u_2` is exactly the sum of the two sources at `t=0`,
and decay fixes `p=p_1+p_2`.

Let `r=X_1-X_2`, `h=n_2 dot r`, and `s=r dot r`.  From (9),

    p_2(X_1)=-(rho J_2/3)(3h^2-s)s^(-5/2).

Taking the gradient before any axial specialization gives

    -grad p_2/rho
      =J_2[2h n_2 s+(s-5h^2)r]/s^(7/2),

which is (13).  With `r=-d n` and parallel axes this is
`+2J_2 n/d^4`; exchanging the labels gives `-2J_1 n/d^4`.  With `h=0`
it is `J_2 r/|r|^5`.  These checks independently fix the sign, normalization,
anisotropy, and `d^-4` scaling.

The self-pressure gradient is odd about its own swirl center and averages to
zero on a centered spherical tag.  The other pressure is harmonic on that tag,
so its volume average is exactly the center gradient used above.  The kinetic
energy is independently

    (rho/2) integral f^2 |n cross y|^2 dx=(4pi rho/3)J,

using `integral_S2 sin^2(theta)dOmega=8pi/3`.

In general `M_1 X_1ddot+M_2 X_2ddot` is not zero.  This is consistent, not a
momentum defect: the complement of the two material tags is a material ambient
region whose momentum rate supplies the negative tag-force sum (together with
the equivalent pressure flux across the tag boundaries; the infinity flux
vanishes under the stated decay).  Since pressure-driven velocity outside the
initial compact supports is generated immediately, later pressure sources have
cross terms.  Formula (13) is therefore exactly an initial interaction and no
pair-potential continuation is licensed.

## API and evidence audit

The additive module is faithful to the theorem at its declared finite-mass
algebra scope.  It imports the existing simultaneous moment calculation,
centers both position and velocity before rates are formed, divides supplied
physical pressure gradients by `rho`, and exposes rather than suppresses force,
pressure moment, torque, covariance work, and total kinetic-energy work.  It
does not pretend that arbitrary samples arise from one Euler field or solve the
Poisson equation.  Known nonpositive density, radial moment, support radius,
interior points, and nonunit axes are rejected; undecidable symbolic sign,
reality, and exterior-domain facts remain caller hypotheses, as is appropriate
for this conditional symbolic API.

The focused tests have meaningful but bounded roles:

- independent trajectory differentiation checks every returned centered rate,
  including nonzero force and torque;
- uniform acceleration exposes missed centering terms;
- the Cartesian polynomial profile exposes the source, harmonic projection,
  Green normalization, and pressure sign;
- exact pressure differentiation exposes density cancellation, the Hessian,
  harmonic trace, and circulation-even behavior;
- contract tests exercise representative invalid inputs.

The polynomial radial fixture is not itself a `C-infinity` zero extension at
the unit sphere.  The test and receipt state this correctly: it is an exact
coefficient regression, while smooth-bump admissibility is proved analytically
for the theorem.  Section 5 has no dedicated additional test, but the general
gradient above is exact and stronger evidence for (13); adding a regression
would improve maintenance only and is not missing scientific support.

The recorded `9 passed`/exit-zero receipt is consistent with the inspected
test sources.  I did not rerun it merely to reproduce a count at the unchanged
boundary.  A targeted trailing-whitespace scan is clean for the reviewed files
and this attempt.
No small-ratio numerical prescription binds because every reviewed quantity is
an exact symbolic identity and no numerical floor or fitted force is used.

## Findings and four-axis decision

No blocking scientific finding or source correction is required at this
boundary.  The strongest statement is already scoped to the smooth existence
interval, exact initial pressure response, and unrestricted local-moment
closure counterexample; its exclusions prevent particle, Coulomb, reciprocal
pair-law, stability, and all-time overclaims.

- Verification: `symbolic_verified` by independent exact calculus; focused
  tests are regression evidence.
- Review: `audited` in this non-author pass.
- Compatibility: `native` to constant-density incompressible Euler with
  physical pressure and the campaign's smooth finite-energy initial-data
  invariant.
- Epistemic: `proposed` campaign evidence; no claim is accepted or promoted.
- Relationship: exact P1 balance/closure evidence and bounded exact-initial P3
  interaction input, without upward inheritance to either obligation.

## Result and remaining dependency

The fixed theorem route is established at the precise scope above.  The next
P1/R3 achievement is an actual projected or retained-ambient dynamics whose
state includes the pressure/velocity information exposed by (1)--(5), followed
by a closure or controlled evaluated limit on an actual carrier family.  The
next P3 interaction achievement is persistence and a controlled same-family
interaction law for an actual stable R1 carrier.  Those constructions, plus
the independent P2 and P4--P7 obligations, remain open; this review does not
complete the electron/neutrino campaign.

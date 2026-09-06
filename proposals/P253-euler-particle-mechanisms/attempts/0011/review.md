# Independent separate-unit review of P253/0005 and P253/0009

Reviewer: `particle-balance-review`, a non-author/non-implementer reviewer
coordinated through Herdr pane `w3:p1`. I did not author or implement either
reviewed derivation, the shear-memory module, or its tests. The completed
`0004` and `0008` reviews are unrelated review work and confer no authorship or
automatic evidence on these units.

## Frozen transaction

The centrally activated `0011/README.md` was read before either source body.
Its SHA-256 is
`7d1d78412c9611f780fce6c352dd471781eb76768e56d799f6c145587cdf201e`;
the activation exit token is zero. This is one substantive equation-level pass
over two separately adjudicated units at the following working-tree boundary.

Unit A:

- `0005/derivation.md`:
  `0cde70fc4fe089654d42af7095d15380bab26a5f8b4cb8c3d1b81846913ff294`;
- `0005/source-audit.md`:
  `c91ae139b6d0f2db547c1d39422797c504145c1af7e075666e8f8a40d15f3d1a`;
- `0005/result.yaml`:
  `edf5b398f5c9ba281233db7ec50c1653d95f1af68606e5975c06b394d786e3cb`;
- exact verifier, captured output, and exit token:
  `36cc06849954d284905f058ca95431b2599ed2459e8d8c43dd99ea38b6fd5f64`,
  `378927d90259d43e2b95edfc18ef4bb31694c6595ce49917abc188af7a24d2ae`,
  and
  `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

Unit B:

- `0009/retained-euler-memory.md`, including the explicitly frozen initial-sine
  forcing:
  `ee560e1b12fbd94609f0867e9385ae9a2252f6a24663bd6e8606622f6e2b81be`;
- `src/substrate_framework/euler_shear_memory.py`:
  `821778ccd55e3f01ba55f769c33b6de7e5f5df745369fd51fe8d7a53e0b4d17f`;
- `tests/test_euler_shear_memory.py`:
  `ff8183a92d69b367b6c426e37e2245f9917657949572fd772cefbcbf359ca698`;
- captured four-test output and exit token:
  `9e4592019d9488a1b7ba853a93f4aa7e844f2f0b02615724b9c863ec28699691`
  and
  `9a271f2a916b0b6ee6cecb2426f0b3206ef074578be55d9bc94f6f3fe3ab86aa`.

The accepted base is `v0.183.0` at
`b6fc902a0942d07996f12a81028fbd3f7c909a43`; the active campaign branch was at
`a7676b7` during review. There is no accepted-claim, release, migration, or
existing-API contract change in this review. The recorded six-check and
four-test receipts are bounded regression evidence, not the basis of either
scientific verdict. I did not rerun them merely to reproduce their counts.

## Separate adjudications

### Unit A disposition

The full Euler--Poincare action, representation-A coadjoint-orbit action,
Hodge/Biot--Savart state completion, and compact-swirl rotation-sphere KKS
calculation are **established** with the stated density and sign conventions.
The prequantum and spinorial statements are established only as conditional
geometry on that finite sphere.

The representation-B route is **established locally and on regular Clebsch
coverage when the global circulation sector is already fixed**, but the claim
that the displayed labels/cocycle give the complete Euler orbit action on an
arbitrary manifold with nonzero `H^1(M)` is **refuted as written, with the
missing harmonic momentum row named**. The minimum correction is to scope
(13)--(18) to `H^1(M)=0`/decaying Euclidean data or to a fixed admissible
global-cocycle sector. A genuinely general statement must add and vary the
harmonic/circulation variables, or explicitly derive their evolution from the
global cocycle. The full Euler--Poincare and representation-A actions already
supply that row and are not affected.

The compact-cutoff rotation fields are valid tangent/orbit lifts and give the
reported KKS pairing, but they do not themselves form an `so(3)` Lie algebra in
the transition annulus. The global `SO(3)` action is the ordinary rigid action
in the larger volume-preserving group. Alternatively, compactly supported
time-dependent lifts equal to rigid rotations near the swirl support put every
rotated field on the same compact-fluid orbit. With that wording repair, the
genuine `S^2`, its stabilizer, form, period, and conditional `SU(2)` lift all
survive unchanged.

### Unit B disposition

The full nonlinear projected Euler split and energy exchange, the exact
finite-amplitude triangular Euler family, the bounded block elimination, the
Bessel kernel, and the deterministic unresolved-initial-state forcing are
**established as stated**. The nearby plus/minus profiles are a natural exact
counterexample only to prediction from the chosen resolved amplitude together
with its energy; the source does not turn this into a global no-go against all
finite or nonlinear closures.

One non-load-bearing operator wording correction is required. For
`L_k=-i beta sin(y)`, the skew-adjoint generator has spectrum
`i[-|beta|,|beta|]`. The real interval and arcsine measure in section 3 are the
spectral support/observation measure of the self-adjoint frequency operator
`i L_k=beta sin(y)`. The branch-point and finite-dimensional
constant-coefficient linear-system conclusion is unchanged.

Neither adjudication establishes a particle, force law, quantum theory,
relativistic dispersion, localized carrier, generic three-dimensional
stability, or parent-campaign completion.

## Unit A: coadjoint sign, density, and Hodge state

With

    <ad*_xi m,eta>=-<m,[xi,eta]>,

integration by parts for divergence-free tangent/decaying fields gives
`ad*_xi[m]=[L_xi m]`. The physical right-reduced tangent is therefore
`X_xi=-ad*_xi m`, and

    Omega(X_xi,X_eta)=-<m,[xi,eta]>.

For `m=[rho u^flat]` and `delta H/delta m=u`,

    dH(X_eta)=-<m,[u,eta]>=Omega(X_u,X_eta),

so `m_dot=-ad*_u m`. The representative equation is

    partial_t u^flat+L_u u^flat=-d Pi,

and differentiating gives material vorticity-form transport. All appearances
of physical density are consistent: `rho` belongs in `m`, `H`, and the KKS
pairing, while the kinematic vorticity form is `B=d u^flat`.

On a closed manifold, exact `B` satisfies

    u^flat=delta G_2 B+h,

because `d delta G_2 B=B`; the remaining co-closed, closed one-form is harmonic.
On decaying `R^3`, the harmonic term is absent and the displayed Biot--Savart
normalization is the standard one. On a bounded multiply connected domain,
the boundary-conditioned Hodge--Morrey inverse plus global circulation data is
indeed required.

The precise global interpretation matters. The harmonic component is
additional **instantaneous state data** needed to reconstruct velocity from
vorticity. Its coefficients in a fixed metric harmonic basis are not thereby
proved time-invariant. Kelvin invariants are circulations around material
loops; for nonzero vorticity they constrain the combination of the particular
inverse and harmonic part. Thus `(B,h)` is a state parametrization, while
membership in one coadjoint orbit imposes the corresponding global circulation
relations. Unit A should not promote `h` itself as a vector of fixed-basis
constants of motion.

## Unit A: the two action signs

For the material action, the constrained variation

    delta u=partial_t xi+[u,xi]

and endpoint/boundary conditions give

    delta S_EP=-integral <partial_t m+ad*_u m,xi> dt.

This independently fixes the sign in the Euler--Poincare equation. For the
inverse reconstruction `m=Ad*_g m0`, `delta g g^-1=-xi`, the Maurer--Cartan
calculation gives

    d<m0,g^-1 dg>(X_xi,X_eta)=Omega_KKS(X_xi,X_eta).

Since the paper's first-order convention is `Omega=-dTheta`, its local
potential is `Theta=-<m0,g^-1 dg>` and (11) has the correct overall sign.
Nonzero KKS periods forbid a single global potential but do not obstruct the
globally defined closed two-form or its local action charts.

For the Clebsch action

    S_B=-rho integral (theta_t+alpha beta_t+|u|^2/2),
    u^flat=dtheta+alpha dbeta,

the three variations give incompressibility and material advection of `alpha`
and `beta` with the displayed signs. Its canonical form is

    Omega_B=rho integral delta alpha wedge_field delta beta,

and substituting displacement variations
`delta_xi alpha=-xi(alpha)`, `delta_xi beta=-xi(beta)` yields exactly
`Omega_B=-<m,[xi,eta]>`. The local pullback and the single-chart helicity
obstruction are therefore correct.

### The missing harmonic row on nontrivial `H^1`

Taking `d` of a one-form momentum equation forgets every closed component, not
only a pressure exact form. The exposing case is a flat torus with

    u(t,x)=h(t),  B=0.

The local equations in (14) and `partial_t B+L_u B=0` allow arbitrary `h(t)`;
periodic Euler instead requires

    h_dot=-average(grad p)/rho=0.

A local potential `theta=h(t) dot x` is multivalued. Its time-dependent overlap
constants are precisely global variables, and the displayed bulk action does
not provide their variation or prove the required zero mode. Merely saying
that triple-overlap constants encode a cocycle does not close this Euler row.
The `H^1=0`/decaying case is unaffected; so is a fixed global sector whose
allowed variations have already been proved to preserve its circulation data.
For arbitrary `H^1`, representation A or the material action remains the
complete action until an explicit harmonic/cocycle term is supplied.

## Unit A: compact-swirl sphere, stabilizer, and period

For `u_n=f(r)n cross y`, direct angular integration gives

    H=(4 pi rho/3) integral r^4 f(r)^2 dr,
    L=(8 pi rho/3) integral r^4 f(r) dr n=j n.

There is no missing density or factor of two. For a nonzero profile, the
stabilizer under genuine rigid rotations is exactly the rotations fixing the
directed axis `n`, hence `SO(2)`, and the set of rotated fields is
`SO(3)/SO(2)=S^2`. This is only the rigid stabilizer; the full fluid
stabilizer is larger.

The ordinary rigid generators satisfy

    [r_a,r_b]=-r_(a cross b).

The physical tangent changes `L` by `delta L=a cross L`, so the KKS restriction
has the outward orientation

    Omega_rot(X_a,X_b)=(a cross b) dot L,
    Omega_rot=j sin(theta) dtheta wedge dphi.

For cutoff lifts `xi_a=chi(r)r_a`, radiality removes derivative cross terms but

    [xi_a,xi_b]=-chi(r)^2 r_(a cross b),

which is not `xi_(a cross b)` in the annulus. Thus the cutoffs are not by
themselves an `SO(3)` subgroup. They are nevertheless sufficient for the KKS
calculation because `chi=1` wherever `m` is supported, so their bracket pairing
is exactly the rigid pairing. A time-dependent cutoff lift of a rigid rotation
path is likewise rigid on a ball containing the support and proves compact-orbit
membership of each rotated field. The global homogeneous-space description
uses the actual rigid `SO(3)` action in the larger volume-preserving group.

The sphere integral is consequently

    integral_S2 Omega_rot=4 pi j
      =(32 pi^2 rho/3) integral r^4 f(r) dr.

This is a KKS/symplectic **period**, not a temporal period of an Euler orbit.
The Hamiltonian is constant on the rotation sphere, and the source correctly
does not claim this sphere is an invariant nontrivial Euler trajectory.

The north/south potentials differ by `2j dphi`, so prequantization with an
imported `hbar` requires `N=2j/hbar` integral. Equatorial curvature holonomy is
`exp(i2pi j/hbar)=(-1)^N`. An axial physical rotation is a constant base loop
and has connection holonomy one; separately, the standard equivariant `SU(2)`
lift has central element `-1` acting on a degree-`N` fibre by `(-1)^N`. Odd `N`
therefore permits a half-integer representation only after the declared line
bundle/polarization/quantum imports. Continuous Euler amplitude scaling changes
`j` continuously, so it selects neither `hbar` nor an odd class.

## Unit B: exact nonlinear split and energy exchange

Let `B(a,b)=-L[(a dot grad)b]`, with `L` the Leray projector, and let `Pi` be
orthogonal on divergence-free fields. Since `v=Pi u`,

    d(rho ||v||^2/2)/dt
      =rho <v,B(u,u)>
      =rho integral u_i u_j partial_j v_i.

Expanding `u=v+w`, the `vv` term and
`w_j v_i partial_j v_i` are divergences. The two surviving terms are

    rho integral w dot (v dot grad)v
      +rho integral w_i w_j partial_j v_i,

with exactly the positive signs in (2). Orthogonality and total Euler-energy
conservation give the negative of this expression for unresolved energy. The
law neither drops pressure physics nor asserts that either sector is closed.

## Unit B: exact finite-amplitude Euler residual

For

    u=(U(y),0,w0(x-tU(y),y)),

one has `div u=0`. The first two momentum components vanish for constant
pressure, while the last is exactly

    w_t+U(y)w_x=0.

No linearization or small-amplitude assumption enters. Smooth periodic input
remains smooth for every finite time, although transverse derivatives may grow
polynomially. The solution has finite energy per periodic cell, not finite total
energy as an isolated field on `R^3`, and it supplies no generic 3D stability.

The module implements exactly this composition and rejects known dependence
that leaves the invariant class. Reality, smoothness, and periodicity remain
caller hypotheses when symbolic assumptions cannot decide them, as its
documentation states.

## Unit B: bounded block elimination and conjugation

On `H=L2(S1,dy/(2pi))`, multiplication by real bounded `kU` makes
`L_k=-ikU` bounded and skew-adjoint. With `P=|1><1|`, `Q=1-P`, mean `U=0`,
`b=QL_k1`, and `A=QL_kQ` on `QH`, one has `A*=-A` and
`||A||<=||kU||_infinity`. There is no unbounded-domain or PDE-inverse issue.

Writing `g=a1+q`, the conjugate-first inner product fixes the block signs:

    a'=<1,L_k q>=-<b,q>,
    q'=a b+Aq.

Duhamel substitution therefore gives

    a'(t)=-<b,e^(tA)q0>
           -integral_0^t <b,e^((t-s)A)b>a(s) ds.

This is precisely (5), with `F=-<b,e^(tA)q0>` and a minus memory convolution.
The so-called noise term is deterministic unresolved-initial-state forcing;
no stochastic law or independence assumption is present.

## Unit B: Bessel kernel, forcing, and spectral scope

For `U_s sin y`, `beta=kU_s`, direct averaging gives, on the branch positive
for real `s>0`,

    a_0(t)=average exp(-i beta t sin y)=J0(beta t),
    Ahat(s)=average 1/(s+i beta sin y)=1/sqrt(s^2+beta^2).

The scalar Volterra resolvent then fixes

    Khat(s)=sqrt(s^2+beta^2)-s,
    K_beta(t)=beta J1(beta t)/t,
    K_beta(0)=beta^2/2.

Indeed

    Laplace[J1(beta t)]=(1-s/sqrt(s^2+beta^2))/beta

and differentiation with respect to `s`, together with decay at infinity,
gives `Laplace[J1(beta t)/t]=Khat/beta`. This proves the kernel rather than
fitting it.

For the newly explicit profiles

    g0_plus/minus=1+/-i epsilon sin y,

differentiating the same transported integral gives

    a_plus/minus(t)=J0(beta t)+/-epsilon J1(beta t).

Here `b=-i beta sin y` and `q0_plus/minus=+/-i epsilon sin y`, so

    F_plus/minus(0)=-<b,q0_plus/minus>
                   =+/-epsilon beta/2,

and the all-time block identity gives

    F_plus/minus(t)=+/-epsilon J1(beta t)/t,
    Laplace[F_plus/minus]=+/-epsilon Khat/beta

for nonzero `beta`, with the removable zero-frequency limit treated by the API.
The amplitude, kernel, and forcing returned by
`sinusoidal_shear_memory` match these formulas and limits.

For nonzero `beta`, the self-adjoint frequency multiplier
`iL_k=beta sin y` has absolutely continuous arcsine observation measure on
`[-|beta|,|beta|]` and no `L2` eigenvector at a single frequency; the
skew-adjoint generator `L_k` itself has the corresponding imaginary spectrum.
The observed resolvent has square-root branch points and is non-rational.
Since every finite-dimensional constant-coefficient linear realization has a
rational resolvent, no such finite linear system reproduces this exact all-time
transfer. Time-dependent, nonlinear, infinite-dimensional, and controlled
approximate realizations remain outside that negative result.

## Unit B: natural perturbation and energy normalization

The physical vertical velocity is real after taking
`Re[e^(ikx)g]`. The plus/minus initial profiles approach the reference field in
every fixed Sobolev norm as `epsilon` tends to zero, have the same resolved
amplitude `a(0)=1`, and satisfy

    ||g0_plus||^2=||g0_minus||^2=1+epsilon^2/2.

For the real Euler field, the vertical kinetic energy is this norm times
`rho*cell_volume/4`; the background x-shear energy is the same for both and is
separately constant. Their opposite initial resolved slopes therefore prove
that this resolved amplitude and the stated energy do not determine the
future. They do not rule out a richer resolved state that retains the sine
component, nor all nonlinear closure representations.

## Evidence and API boundary

The six symbolic checks recorded for unit A correctly corroborate the angular
integrals, KKS area, chart jump, equatorial phase, and action scaling. They do
not prove global Clebsch descent, existence of a compact `SO(3)` subgroup,
carrier stability, or quantum realization.

The four direct tests recorded for unit B have meaningful exposing roles: a
nontrivial two-harmonic profile checks the complete Euler residual; direct sine
moment integration checks amplitude and initial-force parity; the full Laplace
block equation checks the kernel and nonzero forcing; and the equal-energy
plus/minus fixture checks the natural unresolved-state distinction. They do not
prove localization, generic stability, or a particle reduction.

No small-ratio numerical prescription binds. Every reviewed quantity is an
exact algebraic, variational, integral, or bounded-operator identity; no soft
eigenvalue, numerical force, small energy splitting, or discretization floor is
used.

## Four-axis decisions and required corrections

Unit A:

- Verification: `symbolic_verified` for the full Euler/orbit signs, Hodge
  reconstruction, local Clebsch pullback, and compact-swirl KKS sphere.
- Review: `audited`, with the arbitrary-`H^1` Clebsch completion refuted by the
  constant-harmonic torus case.
- Compatibility: `native` for the Euler--Poincare/orbit actions and decaying
  `R^3` swirl; general-topology Clebsch is `required_refactor`.
- Epistemic: `proposed`; no accepted claim is promoted.
- Minimum repair: scope representation B to `H^1=0`/decay or an explicitly
  fixed admissible cocycle sector, or add the harmonic/circulation variational
  row. Describe cutoff rotations as tangent/path lifts rather than an `SO(3)`
  subgroup of `Diff_c`.

Unit B:

- Verification: `symbolic_verified` for the nonlinear split, exact Euler
  family, bounded elimination, Bessel kernel, and initial-state forcing.
- Review: `audited` and established at the stated periodic-cell/finite-linear
  scope.
- Compatibility: `native_periodic_euler`; it is not an isolated Euclidean
  finite-energy carrier.
- Epistemic: `proposed`; no accepted claim is promoted.
- Minimum repair: attribute the real arcsine interval to `iL_k`; state the
  spectrum of `L_k` on the imaginary axis. No equation or API change follows.

## Result and remaining dependencies

Unit A preserves a genuine exact full Euler orbit action and compact-swirl KKS
sphere after the bounded global-topology and cutoff-lift corrections above.
The next P4 achievement remains a stable localized intrinsic-swirl carrier with
a controlled invariant/adiabatic symplectic sector, a physical scale-selection
mechanism, an extension of the line bundle and polarization to that controlled
phase space, and an LP2 exchange configuration before statistics can be tested.

Unit B is an exact executed retained-state mechanism on a periodic triangular
Euler class. The next P1/R3 achievement is an analogous controlled projection
of an actual localized carrier/background system, retaining its complement,
pressure/current, and energy exchange. Neither unit inherits upward to the
electron or neutrino objective, and P1--P7 remain active at their existing
central scopes.


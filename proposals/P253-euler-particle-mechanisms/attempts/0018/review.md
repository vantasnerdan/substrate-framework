# P253/0018 independent equation and source-boundary review

## Boundary, activation, and independence

This review adjudicates two fixed units separately. Unit A is attempt 0013's
classical Lamb/Madelung/spinor field-and-action calculus and flat-interface
counterexample. Unit B is attempt 0017's global twisted initial carrier and
radial stationary/translation obstruction. Neither verdict inherits to the
other or to the parent particle objective.

The centrally registered activation receipt is zero. The frozen README hash is
`88729f438dacf7983c88330c4bf547b06a6e1230816ead93fec34f266f7d09f3`.
The accepted base is `v0.183.0` at
`b6fc902a0942d07996f12a81028fbd3f7c909a43`. The reviewer authored or
implemented neither target unit, module, tests, nor receipts. The
preregistration disclosure about isolated search-result snippets remains part
of the frozen independence record. Author verdicts were not used to establish
the decisions below.

## Primary-source applicability

The external bodies were used only to check the boundaries claimed by Unit A.
Their locally cached hashes are:

- Marmanis author PDF:
  `3f9d317c0f4e0617a30a0937c2e79d1363c679c4508ba98a9d812bae33bf12f0`;
- Fusca `arXiv:1512.04611v2` PDF:
  `e9c13297c2d406e3416034db50c0fa87abdd02f5368472b3d91aea1fe77c6ba7`;
- Meng--Yang `arXiv:2302.09741v1` PDF:
  `052d97dfc7a8becfd0771ee6e70db7d64353767ac5ad60561b7bf0a1b91f2d5f`;
- Zareei `arXiv:2105.12253v1` PDF:
  `c703983aef60861edd1559781b09042d88cee6662770f8b98dc0c011834b8660`.

Marmanis equations (3)--(13) use vorticity and the Lamb vector as formally
electromagnetic fields, but the radiative averaged system subsequently adds
filtering, modeled sources, and constitutive closure. It does not derive a
universal source-free speed for an arbitrary Euler flow. Fusca Theorem 3.5
makes the scalar Madelung map a momentum map into the semidirect-product phase
space for compressible quantum hydrodynamics; its own remarks identify the
gradient/potential-flow subspace. It does not identify fixed-density vortical
Euler with Schrödinger dynamics.

Meng--Yang Table I calls the fluid parameter corresponding to `hbar`
arbitrary. Its equations (24)--(31) distinguish a direct Euler encoding from a
Hermitian incompressible Schrödinger flow; the latter carries an extra
spin-gradient/Landau--Lifshitz-type force. Zareei's Heaviside construction
assumes a supplied free surface and surface tension. In deriving its real
interface equation, it restricts the test function by
`partial_n g|interface=0`, exactly the restriction tested below. Chern and
Slobodeanu supply context only; no load-bearing theorem from either is needed
for these exact results.

## Unit A: Lamb, Madelung, and spinor Euler action

### Lamb identities and their physical boundary

For constant material density `rho`, define

    B=curl u,
    E=B cross u,
    H=p/rho+|u|^2/2.

The vector identity
`(u.grad)u=grad(|u|^2/2)+B cross u` gives incompressible
Euler in Lamb form

    u_t=-E-grad H.

Taking curl and divergence therefore yields

    div B=0,
    B_t=-curl E,
    q=div E=-Delta H.                                   (A1)

For an arbitrarily supplied positive constant `c0`, setting

    J_c=c0^2 curl B-E_t                                (A2)

makes `E_t=c0^2 curl B-J_c` and
`q_t+div J_c=0` exact identities. The continuity row follows from
`div curl B=0`; it does not certify that the supplied field satisfies Euler.
Changing `c0` changes the divergence-free part of `J_c` without changing the
flow, so no universal propagation speed has been selected.

The periodic shear `u=(sin y,0,0)` gives

    B=(0,0,-cos y),
    E=(0,-sin y cos y,0),
    q=-cos(2y),
    J_c=(c0^2 sin y,0,0),

with zero Euler, Faraday, and charge-continuity residuals. Its nonzero defined
current exposes why (A1)--(A2) are not a source-free wave equation.

For smooth compact vorticity and the decaying Hodge velocity, `E=B cross u`
vanishes outside the vorticity support even though `u` has a dipole tail.
Gauss's theorem on any enclosing surface in the vorticity-free region gives

    integral q dx=surface_integral E.n dA=0.             (A3)

The same holds componentwise for disjoint compact-vorticity regions separated
by a vorticity-free surface. Thus the literal Lamb charge of such a localized
component cannot be a nonzero electron charge. This refutes that direct map,
not every structured-background or emergent electromagnetic construction.
Euler kinetic energy is also not the Maxwell quadratic in `E` and `B`:
differentiation and nonlinearity prevent equality without an additional map
and dimensional conversion.

### Scalar Madelung action

For `n>0`, `psi=sqrt(n) exp(iS/hbar)`, and positive supplied `m,hbar`, direct
differentiation gives

    hbar^2 |grad psi|^2/(2m)
      =n|grad S|^2/(2m)+hbar^2|grad n|^2/(8mn).           (A4)

The standard Schrödinger canonical term is `-n S_t`; varying the complete
action gives the continuity equation and

    S_t+|grad S|^2/(2m)+V
      -hbar^2 Delta(sqrt(n))/(2m sqrt(n))=0.              (A5)

The last term is exactly the Euler--Lagrange derivative of the positive
density-gradient energy in (A4). Here `n` is number/probability density, not
Euler's fixed mass density `rho`. Fixing `n` removes the gradient term but
leaves `v=grad S/m` locally curl-free. On all of `R3`, a regular single-valued
finite-energy incompressible scalar phase has an `L2` harmonic gradient and is
therefore trivial; torus harmonic/circulation sectors are separate. The scalar
map consequently does not cover a regular localized nonzero-vorticity Euler
carrier. Nodes, singular phase, multiple charts, or spinor texture are
different representations rather than exceptions silently contained in the
scalar theorem.

### Local normalized-spinor action

For a unit complex two-vector `z`, supplied circulation scale `kappa>0`, and
Pauli texture `s=z^dagger sigma z`, define

    a=-i z^dagger dz,
    u^flat=kappa a.

In the regular chart

    z=exp(i theta/kappa)
      (sqrt(1-f),sqrt(f) exp(i beta/kappa)), 0<f<1,

direct differentiation gives

    u^flat=d theta+f d beta,
    curl u=grad f cross grad beta.                        (A6)

The pointwise horizontal decomposition of `dz` gives the exact identity

    |grad z|^2=|u|^2/kappa^2+|grad s|^2/4.                (A7)

Hence the Euler Hamiltonian in these coordinates is

    H_E=(rho kappa^2/2) integral |grad z|^2
        -(rho kappa^2/8) integral |grad s|^2.             (A8)

The texture term must be subtracted. Keeping only positive spinor Dirichlet
energy changes the dynamics to a different spin-texture system, consistently
with the primary-source boundary.

The local first-order action

    integral rho[i kappa z^dagger z_t-|u|^2/2] dx dt
      =-rho integral[theta_t+f beta_t
                     +|grad theta+f grad beta|^2/2] dx dt (A9)

is the corrected 0005 local Clebsch action. Its variations give
`div u=0`, `D_t f=0`, and `D_t beta=0`; pressure is reconstructed by

    D_t theta=|u|^2/2-p/rho.

Therefore

    (partial_t+L_u)u^flat=d(|u|^2/2-p/rho),              (A10)
    i kappa D_t z=(p/rho-|u|^2/2)z.

Equation (A10) is exactly Euler's one-form equation. The displayed shear
labels reduce to `u=(kappa f(y),0,0)` and satisfy (A10), so this is an executed
vortical Euler chart rather than a prescribed transport field. The
construction is local/regular: it does not prove one global chart for nonzero
helicity or replace 0005's harmonic row on nontrivial `H1`.

### Flat-interface counterexample

For Zareei's flat resting zero-curvature surface at `z=0`, with zero velocity
potential, gravity, and pressure, the proposed wavefunction reduces to
`psi=H(-z)` and `V=0`. Distributionally,

    partial_z psi=-delta(z),
    partial_z^2 psi=-delta'(z),
    R=-(kappa^2/2) delta'(z).                             (A11)

For `phi(z)=z exp(-z^2)` times a compact tangential factor of unit integral,
`<R,phi>=kappa^2/2`, while a test with zero normal derivative gives zero. No
product of singular distributions is needed. The source's
`partial_n g=0` restriction deletes precisely this nonzero double layer and
cannot establish the unrestricted distributional Schrödinger equation.
Equation (A11) refutes the literal unregularized Heaviside transfer, while
leaving smooth-interface models or explicitly defined singular operator
domains open.

### Unit A verdict

Unit A is **established at its exact classical field/action scope, with the two
declared direct-map routes refuted by named mechanisms**. The Lamb identities,
neutrality theorem, scalar Madelung action split, local normalized-spinor Euler
action, and Heaviside counterexample are exact. The literal nonzero localized
Lamb-charge map is refuted by (A3); the unrestricted Heaviside Schrödinger map
is refuted by (A11). No scientific wording correction is required.

This unit licenses neither a source-free Maxwell sector nor a quantum theory.
`kappa`/`hbar` is supplied, spinor coordinates do not select physical spin,
and no Hilbert-space superposition, Born measure, noncommuting measurement,
exchange statistics, or relativistic dispersion follows.

## Unit B: global twisted Hodge carrier

### Smooth field and Hodge projection

Choose a smooth nonincreasing radial profile `F` equal to `pi` for `r<=a`,
zero for `r>=b`, and strictly decreasing through `(0,pi)` on the selected
transition shell. The displayed spinor

    z0=(cos F+i n_z sin F,(-n_y+i n_x)sin F)

has unit norm and is constant near the origin and infinity. Flat splicing
removes every apparent `1/r` singularity. Its Berry one-form is

    v^flat/kappa=
      cos(theta) F' dr
      -sin F cos F sin(theta) dtheta
      +sin^2 F sin^2(theta) dphi.                        (B1)

The corresponding compact smooth vector field is not divergence-free:

    div v=kappa cos(theta) S(r),
    S=F''+2F'/r-sin(2F)/r^2.                             (B2)

For `phi=kappa h(r)cos(theta)`, the regular-at-zero, decaying solution of
`Delta phi=div v` is

    h=-(1/3)[r^-2 integral_0^r s^3 S(s)ds
             +r integral_r^infinity S(s)ds].             (B3)

Direct differentiation returns
`h''+2h'/r-2h/r^2=S`. Since `S` vanishes near zero, `h` is proportional to
`r` there and `h cos(theta)` is smooth Cartesian linear data. Define

    u0=v-grad phi,
    omega0=curl v,
    z=exp(-i phi/kappa)z0.                               (B4)

The phase sign in (B4) gives
`-i kappa z^dagger dz=v^flat-dphi=u0^flat`. Thus `u0` is smooth,
divergence-free, finite-energy, and in `H^s(R3)` for every finite `s`, while
`omega0` is smooth and compactly supported. The velocity is correctly allowed
to have a noncompact tail.

### Tail and same-field invariants

Integration by parts in (B2) gives

    C=integral r^3 S dr
      =integral r[2F-sin(2F)]dr>0.                       (B5)

Outside `r=b`, `h=-C/(3r^2)`, so the field is exactly

    u0=(kappa C/(3r^3))[e_z-3n_z n].                    (B6)

Comparison with the compact-vorticity Biot--Savart expansion gives the
density-free impulse

    I=(1/2) integral x cross omega0 dx
      =-(4pi kappa C/3)e_z,                              (B7)

with the same sign as (B6). This is a dipole `r^-3` tail, not a Coulomb field
or mechanical force.

Using the same Hodge-corrected velocity and its curl,

    v^flat wedge dv^flat
      =2 kappa^2 F' sin^2(F) sin(theta)
        dr wedge dtheta wedge dphi.

The phase correction changes the helicity density by the compactly supported
exact form `-d(phi dv^flat)`, whose boundary integral vanishes. Since `F`
decreases from `pi` to zero,

    Hel(u0)=8pi kappa^2 integral F' sin^2(F)dr
           =-4pi^2 kappa^2.                             (B8)

The Hodge phase has no azimuthal component, so

    u_phi=kappa sin^2(F)sin(theta)/r,
    j_z=rho integral (x cross u0)_z dx
       =(8pi rho kappa/3) integral r^2 sin^2(F)dr>0.     (B9)

Equation (B9) is absolutely convergent because the swirl is compact. It does
not claim absolute convergence of every component of the full vector moment;
the meridional dipole tail is excluded from that stronger statement. The map
`z0:S3->S3` has degree `-1` for the stated Cartesian orientation: the angular
map is orientation preserving and `F` decreases from `pi` to zero. The global
real decaying phase in (B4) is null-homotopic and preserves that degree. This
topology does not select charge, `hbar`, or exchange statistics.

### Local-time Euler and global spinor lift

Because `u0` is smooth, divergence-free, finite-energy, and belongs in
`H^s(R3)` for every finite `s`, the standard local classical Euler theorem
applies for any explicitly chosen `s>5/2`. Let `X(t,a)` be that solution's
volume-preserving flow and set

    z(t,X(t,a))=exp[(i/kappa) integral_0^t
      (|u|^2/2-p/rho)(s,X(s,a))ds] z(0,a).               (B10)

Then `i kappa D_t z=(p/rho-|u|^2/2)z`. Its Berry one-form
`alpha=-i kappa z^dagger dz` satisfies

    (partial_t+L_u)alpha=d(|u|^2/2-p/rho),               (B11)

the same equation and initial data as `u^flat`. The transported difference is
therefore zero on the smooth interval, proving a global-in-space spinor lift
of the actual Euler solution. The compact vorticity support remains the
material image of the original shell and helicity is conserved on this
interval. This does not establish global-in-time regularity, stationary shape,
or particle persistence.

At fixed dimensionless shape, scaling gives

    H=rho kappa^2 L E_*,
    I=kappa L^2 I_*,
    j_z=rho kappa L^3 J_*,
    Hel=-4pi^2 kappa^2.                                 (B12)

Shrinking `L` at fixed degree and `kappa` sends the kinetic energy to zero,
but changes other transported data and generally moves between Euler leaves;
it is not a same-leaf instability. Fixing input helicity and `j_z` fixes
`kappa,L` only within this chosen shape family, not universally.

### Stationary and translation obstruction

The actual poloidal streamfunction and swirl label are

    Psi=kappa G(r)sin^2(theta),
    G=r^2(F'-h')/2,
    G'=sin F cos F-h,
    xi=r sin(theta)u_phi=kappa sin^2(F)sin^2(theta).      (B13)

The identity for `G'` uses the full radial Poisson equation, not a truncated
raw Berry field. Steady axisymmetric angular momentum transport requires the
Jacobian of `Psi` and `xi` to vanish. On the transition shell this implies

    (sin^2 F)'G-sin^2(F)G'=0,

and hence `G=c sin^2 F`. Flatness at `b` would require `G(b)=0`, whereas the
exterior Hodge solution gives `G(b)=-C/(3b)`, a contradiction.

For uniform axial translation at speed `U`, the relative streamfunction
coefficient is

    Q=kappa G-U r^2/2.

Flat-edge matching to `Q=c sin^2 F` requires both `Q(b)=0` and `Q'(b)=0`.
Using `G(b)=-C/(3b)` and `G'(b)=C/(3b^2)`, the first condition fixes

    U=-2 kappa C/(3b^3),

but then `Q'(b)=kappa C/b^2>0`. No common axial translation repairs this
radial family. The contradiction is the angular transport row for this exact
ansatz; it does not exclude persistent unsteady motion or other two-coordinate
twisted Euler carriers.

### Unit B verdict

Unit B is **established as a global smooth compact-vorticity, finite-energy
initial carrier with exact local-time Euler/spinor lift**. Its same-field
helicity, finite nonzero axial angular momentum, dipole tail, scaling, and
degree are exact. The stationary and uniformly translating realization of the
specified radial Hopf--Hodge family is **refuted with the angular-transport
edge mechanism named**. No scientific correction is required.

The remaining positive dependency is a controlled same-leaf neighborhood and
a nonlinear persistence/modulation or recurrent-orbit theorem for this or a
different stationary/deforming twisted carrier. The local existence theorem
and nonzero invariants alone do not provide stability or a particle sector.

## Oracle coverage and statuses

Unit A's repaired focused receipt records six tests passing in 2.62 seconds;
the preserved earlier failure was a SymPy half-angle simplification defect,
not a failed identity. Its separate flat-interface exact script exits zero.
The `e38a8e9` full receipt records 2695 tests passing and includes the 0013
module/tests. These checks corroborate the formulas above; primary-source
scope and global/non-inheritance boundaries were audited independently.

Unit B's focused receipt records four tests passing in 3.46 seconds and exit
zero. They expose the direct Berry form and spherical metric factors, full
`l=1` Poisson inverse and divergence cancellation, oriented helicity, and the
translation-edge derivative. The 0017 module and tests are absent from the
`e38a8e9` tree, so the 2695-test full receipt does not cover them and is not
claimed as evidence for Unit B. The focused receipt is the applicable code
oracle. No tests were rerun for this review.

For both units, verification is `symbolic_verified`, review is `audited`, and
epistemic status is `proposed`. Compatibility is `native_or_scoped_extension`
for the exact Euler identities/action charts and global initial carrier;
Schrödinger/Maxwell/particle interpretations outside the stated classical maps
are not compatibility results. No numerical or small-ratio prescription
binds.

Neither unit licenses a universal Euler--quantum equivalence, source-free
electromagnetism, nonzero electric charge, a selected action or length scale,
physical spin one-half, Born probabilities, measurement, statistics,
relativity, nonlinear stability, all-time localization, an electron or
neutrino, parent completion, or a global no-go.

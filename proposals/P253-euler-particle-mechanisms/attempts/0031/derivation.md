# Exact relative moments, helicity, and finite action boundary of the solitary column wave

## 1. Activated supplier boundary and conventions

Central registration is present and `0031/activation-schema.exit` is zero. The
frozen `0031` README hash is
`ba377906ac1f4e8b75444815c7389f6b86a92afa78bb16c16fc296a26bd7defe`.
This derivation uses the exact `0027` existence theorem and exterior estimates
at the boundary independently established by `0028`; it does not use the
concurrent `0029` or `0030` bodies.

The current `0027` construction and exterior hashes are the reviewed hashes

    solitary-wave-construction.md  c18d394a757f9fd967975ba23a486e62087c285c8fda4f08d8aeaa0159c6f850
    exterior-construction.md       0c2cda195dd5bea9163a7a6c582bfcff746a608a2e48d0a9cd7a007c1b260fc3.

The current `result.yaml` and `validation.md` have hashes
`dcef9066af6c83596af6fcba77d4ed1b0e679fc2b044ddbcc8664308460c0b33`
and `66d3ddee32976066e0a7a708d00bb9bf927da6dcab52237645da77fb4eec5486`.
Their difference from the hashes quoted in the earlier review is the
append-only recording of the independent `0028` verdict and validation state;
the proof and API/test boundary is unchanged.

Use the nondimensional variables of `0027` until Section 8. The background and
traveling streamfunction are

    U=(L(r)/r)e_theta,       psi_0=-c r^2/2,
    psi=psi_0+f(r,z),        F_c(psi_0(r))=L(r),                 (1)

where `L'=r w(r^2)`, `w>=0` is smooth, positive and constant near the
axis, and zero for large `r`. The exact laboratory-frame wave is

    u_r=-f_z/r,
    u_theta=F_c(psi_0+f)/r,
    u_z=f_r/r.                                                   (2)

Put

    G(r,z)=F_c(psi_0+f)-L(r),
    v=u-U=(-f_z/r)e_r+(G/r)e_theta+(f_r/r)e_z.                  (3)

The label derivative on the background is fixed, not fitted:

    F_c'(psi_0)=L'/psi_0'=-w(r^2)/c.                            (4)

The exterior has `G=0` and `Delta_* f=0`. The weighted construction gives
`f` and its core radial derivatives in `L^1_z` as well as the exact exterior
energy. These are the estimates used below; finite `||v||_2` alone is never
used as an `L^1` moment claim.

## 2. Vorticity and exact direct field integrals

With `Delta_*=partial_rr-r^-1 partial_r+partial_zz`, direct cylindrical curl
of (2) gives

    omega_r=-F_c'(psi) f_z/r,
    omega_theta=-Delta_* f/r,
    omega_z=F_c'(psi)(psi_0'+f_r)/r.                            (5)

For `f=0`, this reduces to the background vorticity
`Omega=w(r^2)e_z`, and `U dot Omega=0` pointwise. In the exact exterior both
`F_c'` and `Delta_*f` vanish, so the full wave vorticity vanishes there even
though the background velocity retains its `1/r` circulation tail.

### 2.1 Literal axial angular momentum

The axial component of the direct relative angular momentum is

    J_z^rel
      =rho_m integral_R3 (x cross v)_z d^3x
      =2 pi rho_m integral_R integral_0^infinity r G(r,z) dr dz. (6)

Only the swirl perturbation contributes. It is radially compact, and the
Lipschitz label map together with the weighted `L^1_z` control of `f` makes the
last integral absolutely convergent. Thus (6), with physical coordinates, is
a literal finite angular momentum rather than cancellation of two infinite
background values.

If `delta omega_z=(1/r)partial_r G`, the exact radial integration is

    J_z^rel=pi rho_m integral [r^2G]_0^infinity dz
              -pi rho_m integral r^3 delta omega_z dr dz.       (6a)

Both the axis and infinity surface terms vanish because `G=O(r^2)` at the
axis and is radially compact. Equation (6a) is therefore a valid vorticity
moment here; omitting its displayed surface row would not be a general
identity on an infinite-background field.

The transverse components vanish after the azimuthal integral on every
axisymmetric cylinder. The reviewed exterior estimate, however, establishes
finite energy rather than the absolute first moment
`integral |x cross v| d^3x` of its poloidal tail. Therefore this attempt claims
the literal absolutely convergent **axial component** (6), not an absolutely
convergent three-vector obtained from conditional azimuthal cancellation.

Changing the origin changes a full angular momentum by the corresponding
linear-momentum cross term. Axisymmetry leaves only axial translation impulse,
so the axial value (6) is unchanged by that term. Nevertheless the same-density
velocity perturbation supplies no localized mass density and no material
centroid. Equation (6) is a finite axial field moment; it is not thereby an
intrinsic rigid-body or quantum spin.

### 2.2 Literal helicity and its surface term

From (2) and (5), before any integration by parts,

    u dot omega
      ={F_c'(psi)[f_z^2+f_r(psi_0'+f_r)]-F_c(psi)Delta_*f}/r^2. (7)

On a rectangle `D=[epsilon,R_c] x [-Z,Z]` in the meridional half-plane,
self-adjointness of `Delta_*` in measure `dr dz/r` gives

    integral_D F_c(psi) Delta_*f dr dz/r
      =B_D-integral_D F_c'(psi)
          [f_z^2+f_r(psi_0'+f_r)] dr dz/r,                     (8)

where the complete signed boundary term is

    B_D= integral_-Z^Z [F_c(psi)f_r/r]_(epsilon)^(R_c) dz
        +integral_epsilon^R_c [F_c(psi)f_z/r]_(-Z)^Z dr.        (9)

The axis contribution vanishes because `F_c(psi)=O(r^2)` and `f_r=O(r)`.
The axial caps vanish by the weighted `L^1_z` regularity. In the exterior
`F_c=F_infinity`; if `M(r)=integral_R f(r,z)dz`, the exact zero-Fourier
exterior equation makes `M'(r)=0`, so the radial contribution is also zero.
Thus no hidden circulation-tail surface term remains in helicity, and

    H_kin[u]
      :=integral_R3 u dot omega d^3x
       =4 pi integral_R integral_0^infinity F_c'(psi)
          [(f_r^2+f_z^2)/r-c f_r] dr dz.                       (10)

This is absolutely convergent: `F_c'` has compact radial support, the
quadratic term is controlled by the wave energy, and the linear term by the
weighted `L^1_z` estimate. Since `U dot Omega=0`, (10) is simultaneously the
literal full helicity and the background-relative helicity; no
infinity-minus-infinity subtraction occurs. The density-free kinematic
convention is (10). The density-weighted physical helicity is
`rho_m H_kin`; inserting an additional `1/2` is a separate convention and is
not used here.

Equivalently, expansion before reduction gives all three background rows

    H_kin=integral [U dot curl v+v dot Omega+v dot curl v]d^3x. (10a)

The identity
`div(U cross v)=v dot Omega-U dot curl v` shows that the first two
cross terms are equal only after the full boundary flux
`integral_boundary (U cross v) dot n` is retained. Its axis, axial-cap, and
radial-infinity limits are the corresponding specializations of (9) and
vanish by the same estimates. Equation (10), rather than an assumed doubling,
is the primary formula.

## 3. Translation impulse and the exterior zero mode

The vorticity-defined axial hydrodynamic impulse is

    I_z=(rho_m/2) integral_R3 (x cross (omega-Omega))_z d^3x
       =pi rho_m integral r^2 omega_theta dr dz
       =-pi rho_m integral r Delta_*f dr dz.                    (11)

It is absolutely convergent because `Delta_*f=0` in the exact exterior and
the interior source is smooth and axially localized. Integrate in `z` first
and again set `M(r)=integral f(r,z)dz`. Axis regularity gives `M(0)=0` and
`rM'(0)=0`; the exact exterior zero mode gives `M(r)=M_infinity` and
`M'(r)=0` beyond the matching radius. Therefore

    I_z=-pi rho_m [rM'(r)-2M(r)]_0^infinity
       =2 pi rho_m M_infinity
       =2 pi rho_m integral_R f(R,z)dz.                         (12)

This is the load-bearing surface row. Dropping the exterior zero mode would
incorrectly set the impulse to zero. The cylinder-ordered velocity integral

    2 pi rho_m integral integral f_r dr dz

has the same limit (12), but the available decay does not make `v_z`
absolutely integrable on `R^3`. Hence the literal physical observable is the
vorticity impulse (11); the velocity integral is only its specified cylindrical
surface representation.

## 4. Relative kinetic energy and on-shell material action

The independently established literal kinetic excess reduces exactly to

    E^rel
      =(rho_m/2) integral (|U+v|^2-|U|^2)d^3x
      =pi rho_m integral
          [f_r^2+f_z^2+2L G+G^2] dr dz/r.                       (13)

Every term on the right is absolutely integrable. It is a finite Hamiltonian
value on the affine background even though both total kinetic energies are
infinite. Excess does not mean positive: its sign is determined below.

For a finite physical time interval `[0,T]`, volume preservation converts the
material kinetic action to its Eulerian integral. The background-subtracted
on-shell value of the incompressible material action is therefore exactly

    S_mat^rel[0,T]=T E^rel.                                    (14)

The pressure-multiplier term vanishes on both volume-preserving histories.
Equation (14) is a finite physical action with units of energy times time for
every finite `T`. It is not a periodic-orbit action and is not finite on the
whole time line unless its nonzero rate vanishes.

## 5. Actual fixed-background coadjoint geometry

Let `m=[rho_m u_flat]` in the regular dual of divergence-free vector fields,
modulo exact one-forms. Begin with compactly supported volume-preserving
diffeomorphisms, which leave the background and its circulation tail unchanged
outside a compact set. For the bracket

    [xi,eta]=xi dot grad eta-eta dot grad xi,

define

    <ad*_xi m,eta>=<m,[eta,xi]>,
    X_xi(m)=-ad*_xi m,
    Omega_m(X_xi,X_eta)=-<m,[xi,eta]>
      =-rho_m integral u dot [xi,eta] d^3x.                     (15)

Both generators in (15) are compactly supported, so the KKS pairing is
finite despite the infinite background energy. With
`H=E^rel`, background subtraction changes the value but not its derivative,
and

    dH(X_eta)=-<m,[u,eta]>=Omega_m(X_u,X_eta).                  (16)

Equation (16) is the weak Hamiltonian identity against compact generators; it
gives the full Euler equation modulo the pressure-exact row with one consistent
sign and density. It does **not** assert that the generally noncompact Euler
generator `u` belongs to the Lie algebra of `Diff_vol,c`. The physical domain
is `R^3`, so there is no nonzero global `H^1` harmonic row; the asymptotic
circulation class is instead fixed as boundary data. A contractible compact-
generator orbit chart has a finite symplectic potential and first-order action.
Nonzero global KKS periods, if any, require a separate two-cycle calculation
and are not inferred from this local statement. Extending the chart to an
arbitrary nearby finite-excess field requires a separately specified
asymptotic diffeomorphism domain and remains conditional.

The symmetry group can be extended by axial translations and axial rotations,
which preserve `U` and normalize the compactly supported group. If
`Z=e_z` and `R_z=e_z cross x=r e_theta`, their brackets with compact generators
remain compact. Moreover

    integral U dot [Z,eta] d^3x=integral U dot partial_z eta d^3x=0, (17)

so the translation KKS row is finite. With the convention (15), the momentum
maps have positive pairings `J_xi=<m,xi>` after the background/surface
renormalization. Equations (6) and (11)--(12) are respectively the finite
rotation and translation momentum maps. Consequently the exact traveling wave,
whose tangent is `c X_Z`, is a relative equilibrium and obeys

    d(E^rel-c I_z)=0                                             (18)

on the declared dynamically accessible fixed-background orbit. The finite
number `E^rel-cI_z` is an augmented Hamiltonian, not a closed-loop action.

For completeness, the impulse-to-momentum-map step is not assumed and does not
pretend that the velocity tangent is compact. For compact divergence-free
generator `eta`, the coadjoint vorticity tangent is

    delta_eta omega=curl(eta cross omega),                       (18a)

so integration by parts acts on the genuinely compact field
`eta cross omega` and gives directly

    delta_eta I_z
      =rho_m integral (eta cross omega)_z d^3x
      =rho_m integral eta dot partial_z u d^3x
      =Omega_m(X_Z,X_eta).                                      (18b)

The middle equality uses
`omega cross e_z=partial_z u-grad u_z`; the gradient integral vanishes because
`eta` is compact and divergence free. The last equality follows from (15) by
integrating `-integral u dot partial_z eta`. Thus (11)--(12), including the
exterior surface row, is the true translation momentum map on these tangents
even though `delta u=P(eta cross omega)` generally has a Hodge pressure tail.

The rotation row retains that tail as well. For
`R_z=r e_theta`, the Leray part contributes on every complete azimuthal
cylinder

    integral R_z dot grad p d^3x=integral partial_theta p d^3x=0. (18c)

Hence the variation of (6) equals `Omega_m(X_Rz,X_eta)` for the actual
coadjoint tangent, not only for an artificially compact velocity variation.

This establishes a finite local KKS form and weak Hamiltonian identity on the
compact-generator orbit, together with the exact noncompact axial-translation
orbit of the solitary solution and its finite symmetry moment maps. It does
not claim a general noncompact Euler generator lies in `Diff_vol,c`, an
arbitrary-nearby asymptotic group, or a nontrivial compact spin orbit:

- axial `SO(2)` fixes both `U` and the axisymmetric wave, so it is in the
  stabilizer and its state-space tangent is zero;
- axial translation moves the solitary profile nontrivially, but its group is
  `R`, the localized profile has no nonzero translation period, and the orbit
  is not a circle;
- a rotation tilting the column axis is not in the fixed-background phase
  space. On an infinite cylinder around the original axis the tilted column
  tends to zero as `|z|` grows while `U(r)` remains nonzero, so the squared
  difference has an infinite axial integral;
- the branch parameter `mu` (or `c`) is an interval parameter between distinct
  solutions, not a proved dynamically accessible cyclic coordinate.

In particular, writing `2 pi J_z^rel` because `J_z^rel` has action units would
describe a lift through an isotropy group, not the action of a nonconstant
closed state-space orbit. No such action or phase is promoted here.

## 6. Exact small-amplitude observable transfer

The existence proof gives

    f(r,z)=mu F_mu(r,z/L_mu),
    F_mu -> f_0(r) A_*(X),
    A_*(X)=3 sech^2(X/2)/(2 beta),
    integral_R A_*(X)dX=6/beta,                                 (19)

in a weighted space strong enough for every linear functional below. Also
`c->c_0`, `h_0=f_0(R)>0`, and

    mu L_mu^2=h_0^2 log L_mu,
    mu L_mu=h_0 sqrt(mu log L_mu)->0.                            (20)

Define four strictly positive radial constants

    C_J=integral_0^R r w(r^2)f_0(r)dr,
    C_H=integral_0^R w(r^2)f_0'(r)dr,
    C_E=integral_0^R L(r)w(r^2)f_0(r)dr/r,
    h_0=f_0(R).                                                  (21)

Positivity of `C_J,C_E,h_0` is immediate. The critical-mode identity in
`0027` gives `f_0'>0` on an initial interval where `w>0`, hence `C_H>0`.
Taylor expansion of the smooth label, with the exact weighted convergence and
not merely formal substitution, gives

    J_z^rel
      =-[12 pi rho_m/(beta c_0)] C_J mu L_mu+o(mu L_mu),         (22)

    H_kin
      =[24 pi/beta] C_H mu L_mu+o(mu L_mu),                     (23)

    I_z
      =[12 pi rho_m/beta] h_0 mu L_mu+o(mu L_mu),               (24)

    E^rel
      =-[12 pi rho_m/(beta c_0)] C_E mu L_mu+o(mu L_mu).        (25)

For (23), the linear part of (10) is exactly
`-4 pi c integral F_c'(psi_0)f_r=4 pi integral w f_r`; all
gradient and nonlinear-label terms are `O(mu^2 L_mu)+O(mu^2/L_mu)`.
For (22) and (25), `G=F_c'(psi_0)f+O(f^2)` and (4) fixes the signs. For
(24), the exterior trace in (12) supplies `h_0`. Thus these are controlled
asymptotics of the exact integrals, not replacements for them.

It follows that for every sufficiently small nonzero member of the positive
branch, the exact quantities satisfy

    J_z^rel<0,      H_kin>0,      I_z>0,      E^rel<0.           (26)

The augmented Hamiltonian is also nonzero:

    E^rel-c I_z
      =-[12 pi rho_m/beta]
         [C_E/c_0+c_0 h_0] mu L_mu+o(mu L_mu)<0.                (27)

Equations (6), (10), (11), and (13) are the exact full observables. Equations
(22)--(27) give leading behavior and, because their coefficients are strictly
nonzero with controlled little-`o` remainders, eventual exact nonvanishing.
No closed form is claimed for the full finite-`mu` integrals.

## 7. Physical units and all normalization factors

Restore `x_phys=ell_0 x`, `u_phys=U_0 u`, and
`t_phys=(ell_0/U_0)t`. If the radial constants and `c_0` in
(21)--(27) remain nondimensional, the prefactors are

| Quantity | Physical factor multiplying its nondimensional integral |
| --- | --- |
| `J_z^rel` | `rho_m U_0 ell_0^4` |
| `H_kin` | `U_0^2 ell_0^2` |
| `rho_m H_kin` | `rho_m U_0^2 ell_0^2` |
| `I_z` | `rho_m U_0 ell_0^3` |
| `E^rel` | `rho_m U_0^2 ell_0^3` |
| `S_mat^rel` and a KKS/action integral | `rho_m U_0 ell_0^4` |

Thus the `rho_m` displayed in (6), (11)--(13), and (22), (24), (25), (27)
is understood either with physical coordinates or replaced by the complete
factor in this table when using the nondimensional `0027` variables. The
`2 pi` and `4 pi` factors come from the actual azimuthal integral and are not
absorbed into `beta` or the radial constants.

## 8. Route verdicts and next construction

**Route A is established.** The exact solitary family carries a finite literal
axial relative angular momentum, finite literal helicity, finite absolute
hydrodynamic impulse, finite literal kinetic excess, finite finite-time
background-subtracted material action, and a fixed-background finite local KKS
form with the weak Hamiltonian identity (16). Their density, azimuthal, sign,
zero-mode, Hodge, and surface normalizations are explicit. The positive branch
transfers to nonzero exact observables through (22)--(27). The translation
orbit is an actual noncompact relative-equilibrium orbit; the axial rotation is
a stabilizer. A general nearby asymptotic orbit remains conditional on its
function-space group construction.

**Route B is blocked by one named missing construction.** No parameter already
present in `0027` is a nontrivial compact internal coordinate. What remains is
an actual same-background non-axisymmetric or otherwise internal `S^1` family
of exact Euler states, with a nonzero tangent modulo stabilizer, a derived least
period, finite KKS pullback, and a nonzero symplectic period. The present
calculation neither proves nor disproves existence of such a family in Euler.

The next executable construction is therefore to continue a non-axisymmetric
azimuthal mode of the exact solitary wave to a same-background relative-periodic
Euler family and compute its KKS pullback before assigning any internal action.
Stability, `0029` carrier reducibility, `0030` restoring dynamics, particle
identity, quantum spin/statistics, and any value involving `hbar` remain
outside this attempt.

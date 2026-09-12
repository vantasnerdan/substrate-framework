# Dynamically accessible swirl Hessian and Hamiltonian continuation

## 1. Conventions and the exact orbit tangent

Let `Pi={(r,z):r>0}`, `dnu=r dr dz`, and use the 0010 convention

    u=(-Psi_z e_r+xi e_theta+Psi_r e_z)/r,
    zeta=omega_theta/r=L Psi.

For meridional scalars define the weighted bracket

    {f,g}=(f_r g_z-f_z g_r)/r.                                      (1)

An axisymmetric divergence-free displacement has two independent scalar
generators `(g,a)`:

    v_r=-g_z/r,       v_z=g_r/r,       v^theta=a,
    v_theta=r a.                                                     (2)

Initially take `g,a` smooth and compactly supported.  Global regularity at the
axis is obtained with `g=O(r^2)` and finite smooth `a`; both actual carriers used
below lie in a torus separated from the axis, so all exposing generators may be
supported in that torus.  The dynamically accessible space is the `X`-closure
of the image found below.

The velocity one-form and vorticity two-form are

    u^flat=u_r dr+xi dtheta+u_z dz,
    B=d u^flat=-r zeta dr wedge dz
                +xi_r dr wedge dtheta+xi_z dz wedge dtheta.          (3)

For a coadjoint displacement, `delta B=-L_v B=-d(i_v B)`.  Direct contraction
gives

    i_v B=(zeta g_r-a xi_r)dr+(zeta g_z-a xi_z)dz-{xi,g}dtheta.

Comparing the `dr wedge dz` and angular coefficients in (3) proves the exact
axisymmetric coadjoint-orbit tangent map

    chi=delta xi={xi,g},
    eta=delta zeta={zeta,g}+{xi,a}.                                  (4)

This is an image statement derived from `delta m=-ad^*_v m`, not the converse
claim that all variations preserving a displayed list of Casimirs are
dynamically accessible.  In particular, `a` is the independent azimuthal
displacement, not a freely varied Eulerian swirl field.

Equation (4) preserves the complete smooth axisymmetric Casimir families used
in 0010.  If `G'=F`, integration of a compactly supported Jacobian gives

    delta int C(xi)dnu
      =int {C(xi),g}dnu=0,

    delta int zeta F(xi)dnu
      =int {zeta F(xi),g}dnu+int {G(xi),a}dnu=0.                      (5)

It also gives the reduced Lie--Poisson operator.  With a covector `(p,q)`,

    J_(zeta,xi)(p,q)=({zeta,p}+{xi,q},{xi,p}).                        (6)

Thus (4) is `J(g,a)`.  As an independent sign check, for
`delta E=(Psi,xi/r^2)`, equation (6) yields

    xi_t={xi,Psi},
    zeta_t={zeta,Psi}+{xi,xi/r^2}
          ={zeta,Psi}+r^-4 partial_z(xi^2),                           (7)

which is precisely the axisymmetric Euler evolution in 0010.

The axial-translation tangent is genuinely in this closure.  Choose a compact
cutoff generator equal to `g_Z=r^2/2`, `a_Z=0` on a neighbourhood of the
carrier.  Then

    q_Z=(-partial_z zeta,-partial_z xi).                              (8)

It is a Hessian kernel because the augmented functionals are translation
invariant and the carrier is critical.  Axial rotation acts trivially on the
axisymmetric fields; tilt rotations are non-axisymmetric and are not silently
inserted into this sector.  Stabilizer relabellings with zero coadjoint image
are quotient directions in generator space, not additional nonzero `(eta,chi)`
vectors.

## 2. Cao--Zhan: the physical leaf still has both Hessian signs

Write `e` for the fixed positive thickness parameter.  Strictly inside the
active core, 0010 gives

    psi=K zeta-U_e r^2/2-mu_e=e xi,
    zeta=alpha/e^2+xi/(e r^2),                                       (9)

and the Hessian of the conserved augmented Hamiltonian is

    Q_e(eta,chi)=int [eta K eta+chi^2/r^2-2e eta chi]dnu,             (10)
    ||(eta,chi)||_X^2=int [eta K eta+chi^2/r^2]dnu.                  (11)

Both signs occur on the actual image (4).

**Positive accessible directions.**  Set `g=0` and choose `a` supported in an
interior patch on which `grad xi` is nonzero.  Then

    (eta,chi)=({xi,a},0),
    Q_e(eta,0)=||(eta,0)||_X^2=int eta K eta dnu>0.                  (12)

The space of such `a` is infinite-dimensional, so finitely many smooth
symmetry-orthogonality conditions can be imposed without removing all of these
directions.

**Negative accessible sequence.**  Cao--Zhan's `partial_z psi<0` for `z>0`
supplies an interior rectangle `R` away from the axis and free boundary on
which `xi_z` does not vanish.  Take `phi in C_c^infinity(R)`,

    g_k=phi(r,z) sin(k r),       a_k=0.                               (13)

Then (4) and (9) give

    chi_k=-(k xi_z/r)phi cos(k r)+O(1),
    eta_k=-(k zeta_z/r)phi cos(k r)+O(1)
         =chi_k/(e r^2)+O(1).                                       (14)

Local ellipticity away from `r=0` makes `K` order minus two, hence
`int eta_k K eta_k dnu=O(1)` while `int chi_k^2/r^2 dnu=Theta(k^2)`.
Consequently

    Q_e(q_k)=-int chi_k^2/r^2 dnu+O(k),
    Q_e(q_k)/||q_k||_X^2 -> -1.                                     (15)

For any finite-dimensional smooth symmetry kernel `Z`, subtracting the
`X`-projection onto `Z` changes neither limit.  In particular `Q_e(q_Z,w)=0`
for the translation vector (8), and oscillatory pairings with every fixed
smooth mode are lower order.

Combining (12)--(15) with the frozen definition of `gamma_e`,

    inf Q_e/||.||_X^2 <= -1,
    inf (-Q_e)/||.||_X^2 <= -1,
    gamma_e <= -1 < 0.                                               (16)

Thus route A has an exact route-scoped verdict: the full axisymmetric
coadjoint leaf admits no globally consistent energy--Casimir sign and no
coercive gap in the frozen norm.  This does not imply temporal instability.
It instead activates the preregistered Hamiltonian continuation in section 4.

Independent review supplies a useful strengthening without changing that
verdict.  In (13), take `a_k=c g_k` with fixed real `c`.  Then

    eta_k=chi_k/(e r^2)+c chi_k+O(1),
    Q_e(q_k)/||q_k||_X^2 -> -1-2 e c R_phi,
    R_phi=[int w dnu]/[int (w/r^2)dnu]>0,
    w=(xi_z phi/r)^2.                                                (16a)

First taking `k` to infinity at fixed `c`, and then taking `c` to either
sign infinity, proves that `Q_e` is unbounded in both signs on the `X`-unit
accessible tangent sphere.  Thus `gamma_e=-infinity` in the extended-real
interpretation, strengthening the reported `gamma_e<=-1`.  This also shows
why the dynamical generator needs an explicit graph domain; it is not a
temporal instability argument.

## 3. Gavrilov: a local leaf saddle without touching the flat boundary

For a regular steady swirl patch with `xi=F(psi)`, `F'` finite and nonzero, the
general augmented Hessian from 0010 is

    Q(eta,chi)=int [eta K eta+2 A'(xi)eta chi+D(r,z)chi^2]dnu,
    A'=-1/F',
    D=1/r^2+C''(xi)+zeta A''(xi).                                    (17)

The compact Gavrilov carrier has such a patch.  Indeed,
`psi=-int^p omega_c(s)ds` and `xi=omega_c(p)b(p)` in the open cutoff shell.
The positive function `xi` vanishes at both flat edges and is nonconstant, so
there is an interior regular pressure patch with `F'!=0`, `grad psi!=0`, and
therefore `grad xi!=0`.  No statement is made at the flat support boundary.

The pure azimuthal-displacement direction in (12) again gives `Q>0`.  For an
interior high-frequency phase with covector `k`, put

    b=(xi_r k_z-xi_z k_r)/r,
    a_zeta=(zeta_r k_z-zeta_z k_r)/r,
    g_k=phi exp(i k dot x),       a_k=c g_k.                           (18)

Taking real or imaginary parts gives real compactly supported displacements.
At principal order,

    chi_k=i b phi exp(i k dot x),
    eta_k=i(a_zeta+c b)phi exp(i k dot x).                            (19)

The `K` term is two differential orders lower, whereas the local leading
coefficient is

    2 A' b(a_zeta+c b)+D b^2.                                       (20)

On a smaller patch with `b!=0`, the coefficient of `c` is
`2 A' b^2!=0`; choosing the sign and magnitude of the fixed real `c` produces
a negative sequence.  Hence the actual Gavrilov coadjoint leaf also has both
Hessian signs on a smooth interior patch.  The flat cutoff creates an
additional global-Casimir regularity problem, but it cannot restore a sign
already lost to compactly supported interior accessible perturbations.

This closes route B only as a sign-definite energy--Casimir strategy.  It is
not a no-go for Gavrilov persistence or for a smaller invariant perturbation
class justified by additional physics.

## 4. Same-carrier regular-core Hamiltonian continuation

The Hessian operator corresponding to (10) is

    L_e(eta,chi)=(K eta-e chi, chi/r^2-e eta).                         (21)

On the KKS quotient of the orbit by its stabilizer, set
`Omega_e^flat=J_e^{-1}` and use the 0005 convention
`i_X Omega_e=dH`, `dot q=J_e grad H`.  A canonical `dq wedge dp` test fixes
the quadratic orbit action in the moving frame as

    S_e^(2)[q]=int dt [-(1/2) Omega_e(q,dot q)-(1/2)Q_e(q,q) ].       (22)

The opposite kinetic sign gives `dot q=-J_e L_e q`.  With (22), the
Euler--Lagrange equation on a smooth regular active-core patch, where
`dA_e=0`, is

    dot eta={zeta,K eta-e chi}+{xi,chi/r^2-e eta},
    dot chi={xi,K eta-e chi}.                                        (23)

This is not an ambient arbitrary-variation equation: `J_e` maps covectors to
the orbit tangent (4).  Direct linearization of (7) in the translating frame
gives

    dot chi={chi,psi}+{xi,K eta},
    dot eta={eta,psi}+{zeta,K eta}
             +{chi,xi/r^2}+{xi,chi/r^2}.                             (24)

Using (9), equations (23) and (24) agree term by term on the regular active
core.  This independently checks the density and Poisson sign in (6).

It is not a global free-boundary derivation.  For a Poisson evolution
`F(q)=J(q)dA(q)`, the complete linearization is

    DF(q_*) delta q
      =J_* d^2A_* delta q+(dJ_*[delta q])dA_*.                       (24a)

The second term vanishes at an interior critical point, but `psi=e xi` is not
a global identity outside the active support.  Normal displacements of the
Cao vorticity jump can carry boundary distributions and nonzero traces.  The
global operator must therefore be derived from (24), with the physical
interface/transport conditions, or from both terms in (24a).  The compact
interior sign witnesses in sections 2--3 are unaffected.

The exact frozen-coefficient principal symbol exposes what the mixed Hessian
does dynamically.  For a local covector `k=(k_r,k_z)`, define

    b=(xi_r k_z-xi_z k_r)/r,
    a_zeta=(zeta_r k_z-zeta_z k_r)/r.                                (25)

Since `K` has order minus two, the first-order symbols of `J`, the local part
of `L_e`, and their product are

    J#=i [[a_zeta,b],[b,0]],
    L_0=[[0,-e],[-e,1/r^2]],
    (J L)#=i [[-e b,-e a_zeta+b/r^2],[0,-e b]].                      (26)

Relation (9) implies

    a_zeta=b/(e r^2)-2 xi k_z/(e r^4),

so

    (J L)#=i [[-e b,2 xi k_z/r^4],[0,-e b]].                         (27)

Its repeated characteristic value is

    lambda_pr=-i e b=-i u_mer dot k,                                 (28)

which is purely convective.  When `xi k_z!=0`, (27) has a genuine nilpotent
part.  The frozen system is solved exactly by

    chi(t)=exp(-i e b t) chi(0),
    eta(t)=exp(-i e b t)
             [eta(0)+i(2 xi k_z/r^4)t chi(0)].                       (29)

Therefore the raw first-differential-order characteristic has no
`O(|k|)` exponential rate.  It has a double convective eigenvalue and a raw
Jordan shear.  This does **not** exclude a finite `O(1)` short-wave rate in the
physical energy norm: there `eta` has inverse-elliptic weight, so the apparently
lower-order `K` coupling enters the energy-normalised amplitude system at order
zero.  Coefficient gradients, evolving wave covectors, free boundary, and
global orbit geometry enter at the same required continuation level; no
finite-rate or global spectral conclusion is inferred from (27).

The same conclusion holds on every regular interior Gavrilov patch.  For (17),

    L_0=[[0,A'],[A',D]],
    (J L)#=i [[A'b,A'a_zeta+D b],[0,A'b]].                            (30)

Because `b=F' a_psi` and `A'=-1/F'`, its double eigenvalue is again
`-i a_psi=-i u_mer dot k`; a nilpotent off-diagonal term is generic.  Thus the
Hamiltonian continuation retains the same compact carrier and distinguishes
indefinite energy from an order-`|k|` exponential characteristic; finite
energy-normalised rates remain open.

## 5. Strongest physical verdict and next construction

The positive bridge is now exact: the two scalar displacement generators map
to the physical leaf by (4), the conserved Hessian pulls back to that leaf,
and its KKS action generates (23).  For both actual intrinsic-swirl carriers,
sign-definite energy--Casimir restoration fails already on smooth interior
leaf tangents.  For Cao, the frozen target has the quantitative verdict
`gamma_e=-infinity` in the extended-real interpretation, and in particular
the original `gamma_e<=-1` bound remains true.  Hamiltonian continuation then
shows pure convection plus a raw Jordan shear at order `|k|`, excluding only
an order-`|k|` exponential characteristic.  It does not decide finite-rate
short-wave instability and supplies no restoring oscillation.

The concrete missing LP2 mechanism is consequently sharper than “prove
stability”: one must construct the global linearization (24a), including its
interface term, and control the energy-weighted Kelvin/trajectory system and
nonlocal `K` coupling by a resolvent/semigroup or phase-mixing estimate, or
exhibit its global growing mode.  A physical restoring frequency would require a
semisimple isolated imaginary pair and its KKS-normalised action; neither the
mixed Hessian nor the local Jordan rate in (29) supplies one.  This is a route
gap, not a physical impossibility and not a quantum-mechanical conclusion.

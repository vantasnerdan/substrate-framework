# Massive cylindrical fiber, common DA graph, and Fredholm gap

## 1. Physical reduction and normalization

Fix a sufficiently thin `p>=6` Cao ring and write

    r=R_delta+a_delta y_1,
    z=z_delta+a_delta y_2,
    delta=a_delta/R_delta,
    q_delta(y_1)=1+delta y_1,
    k_delta=delta n_delta.                               (1)

The ring remains in its translating, nonrotating frame. For a physical
orthonormal cylindrical vector, use

    v(r,z,theta)=(2*pi)^(-1/2)
       v_n(r,z) exp(i n theta).                          (2)

Then

    ||v||_(L2(R3))^2=integral |v_n|^2 r dr dz
      =R_delta a_delta^2
         integral |v_n(y)|^2 q_delta(y_1)dy.             (3)

Thus the exact unitary meridional reduction is multiplication by
`[R_delta a_delta^2 q_delta]^(1/2)`, after rotating the cylindrical frame to
the fixed reference frame. The straight column uses the normalized Fourier
transform `(2*pi)^(-1/2) integral exp(-ikZ)dZ` and transverse measure `dy`.
Only these fibers are compared. No map sends the ring's `S^1` to the column's
`R`.

All following differential rows are written before the harmless constant
physical length factors are restored. This keeps the exact kinetic/KKS
normalization in (3) and prevents a coordinate norm from being called physical
energy.

## 2. Exact scalar and full-vector cylindrical operators

Put `q=1+delta y_1` and `k=delta n`. For a scalar `f`, the exact reduced
gradient, divergence and scalar Laplacian, multiplied by the appropriate
powers of `a_delta`, are

    G_(delta,k)f
      =(partial_1 f, i k q^(-1)f, partial_2 f),           (4)

    D_(delta,k)v
      =partial_1 v_r+partial_2 v_z
        +delta q^(-1)v_r+i k q^(-1)v_theta,              (5)

    a_delta^2 Delta_n f
      =D_(delta,k)G_(delta,k)f
      =partial_1^2f+partial_2^2f
        +delta q^(-1)partial_1f-k^2q^(-2)f.              (6)

The exact physical curl is

    a_delta curl_n v = C_(delta,k)v,

    C_(delta,k)v=
      (i k q^(-1)v_z-partial_2v_theta,
       partial_2v_r-partial_1v_z,
       partial_1v_theta+delta q^(-1)v_theta
          -i k q^(-1)v_r).                              (7)

Direct calculation gives the differential-complex identities

    D_(delta,k)C_(delta,k)=0,
    C_(delta,k)G_(delta,k)=0,
    D_(delta,k)G_(delta,k)=a_delta^2 Delta_n.            (8)

The `delta q^(-1)v_theta` connection term in (7) is load bearing: reversing
its sign violates `D C=0`.

Let `v_+=v_r+i v_theta` and `v_-=v_r-i v_theta`. The cylindrical vector
Laplacian diagonalizes exactly:

    -a_delta^2(Delta_vec v)_+
      =H_(delta,k)^+ v_+,
    -a_delta^2(Delta_vec v)_-
      =H_(delta,k)^- v_-,
    -a_delta^2(Delta_vec v)_z
      =H_(delta,k)^0 v_z,                               (9)

    H_(delta,k)^s
      =-partial_1^2-partial_2^2-delta q^(-1)partial_1
        +(k+s delta)^2q^(-2),      s=+1,-1,

    H_(delta,k)^0
      =-partial_1^2-partial_2^2-delta q^(-1)partial_1
        +k^2q^(-2).                                     (10)

Thus the three exact axial masses are `delta^2(n+1)^2`,
`delta^2(n-1)^2`, and `delta^2n^2`; replacing all of them by `k^2` before
expansion loses genuine first-order blocks.

For later comparison,

    q^(-1)=1-delta y_1+delta^2 y_1^2+O(delta^3),
    q^(-2)=1-2delta y_1+3delta^2 y_1^2+O(delta^3).       (11)

Consequently the exact scalar row is

    a_delta^2 Delta_n
      =Delta_y-k^2
       +delta(partial_1+2k^2y_1)
       +delta^2(-y_1partial_1-3k^2y_1^2)+O(delta^3),    (12)

and the helical rows are obtained from (12) by replacing `k` by
`k+delta` or `k-delta` before expanding. Equations (4)--(12) display all
`k`, `delta*k`, `delta^2*k`, `k^2`, and `delta*k^2` terms. Since
`k_delta` remains in `[k_-,k_+]`, none is estimated by pretending that
`n_delta` is fixed.

With the physical weight `q dy`, `D_(delta,k)=-G_(delta,k)^*`. The scalar
pressure and Leray rows are therefore

    L_(delta,k)=-D_(delta,k)G_(delta,k),
    P_(delta,k)=1-G_(delta,k)L_(delta,k)^(-1)
                         [-D_(delta,k)],                 (13)

with the decaying/axis solution and finite harmonic rows fixed. Equivalently,
the velocity is the unique decaying solution of

    C_(delta,k)v=eta,
    D_(delta,k)v=0.                                     (14)

The two formulations agree by (8). The support edge of `eta` is a coefficient
interface in this whole-space problem, not an elliptic boundary. No artificial
core DtN condition is imposed there.

## 3. What positive mass does and does not control

For `s=0,+1,-1`, the scalar quadratic form associated with (10) is

    h_(delta,k)^s[f]
      =integral_(q>0)
         {q|grad f|^2+(k+s delta)^2q^(-1)|f|^2}dy.      (15)

On every fixed ring-centered collar `|y_1|<=M` and sufficiently small
`delta`, (15) obeys

    h_(delta,k)^s[f]
      >=c_M(||grad f||_2^2+k_-^2||f||_2^2).             (16)

There is no corresponding unrestricted global mass bound: as physical
`r->infinity`, `q->infinity` and the angular mass `(k+s delta)^2/q^2`
vanishes. Thus the valid comparison is from compactly supported vorticity (or
compact pressure source) to its global decaying field, not norm-resolvent
convergence of the free scalar Laplacians on arbitrary exterior data.

This restricted comparison is sufficient for the Euler generator. Let
`S subset {|y|<=M}` contain the common vorticity support after the meridional
Hanzawa map. In helical component `s=0,+1,-1`, define the scaled global norm

    ||v||_(V_(delta,k)^4)^2
      =sum_(j+ell<=4) integral_(q>0) q
         |partial_y^j[(k+s delta)q^(-1)]^ell v_s|^2dy,  (16a)

with the equivalent covariant-derivative terms included. For a fixed collar
`B_M`, let `nabla_(delta,k)` denote the exact covariant derivative in the
orthonormal cylindrical frame after (1): its angular rows are `ik/q` and its
two connection entries are `delta/q`. The energy-tail seminorm is

    E_tail,M(v)^2
      =sum_(j=0)^4 integral_((q>0)\B_M) q
            |nabla_(delta,k)^j v|^2dy.                  (16b)

This is a physical scaled `H^4` tail, not a norm containing only ordinary
`y` derivatives. On a fixed collar it is uniformly equivalent to (16a),
while near the axis it retains the genuine `n`, `n+1`, and `n-1` polar rows.

The scalar Newton fiber has the exact angular integral kernel

    G_(delta,n)(y,y')
      =[sqrt(q q')/(4*pi)] integral_(-pi/delta)^(pi/delta)
        exp(-i k t)
        [|y-y'|^2+4q q' delta^(-2)sin^2(delta t/2)]^(-1/2)dt, (17)

for the inverse of `-a_delta^2 Delta_n` after the half-density unitary in
(3). To pin the constant independently, projection of the physical Newton
kernel with (2) leaves one angular integral and no residual Fourier factor:

    ((-Delta_n)^(-1)f_n)(r,z)
      =(4*pi)^(-1) integral exp(-in alpha)
       f_n(r',z')r'dr'dz'dalpha/|x-x'|.                 (17a)

Conjugating `a_delta^(-2)(-Delta_n)^(-1)` by
`[R_delta a_delta^2q]^(1/2)` contributes
`(R_delta/a_delta)sqrt(q q')=delta^(-1)sqrt(q q')` after
the factor `a_delta` is extracted from `|x-x'|`. Finally
`alpha=delta t` contributes `dalpha=delta dt`. The two powers of `delta`
cancel and leave exactly `sqrt(q q')/(4*pi)` in (17). For the
helical components replace `k` in the phase by `k+delta` or `k-delta`. On
compact `y,y'`, splitting the integral into
`|t|<=T` and its complement gives

    2delta^(-1)sin(delta t/2)->t,
    G_(delta,n)->(2*pi)^(-1)K_0(|k||y-y'|).              (18)

The operator-norm statement needs more than pointwise convergence. Choose
angular cutoffs `chi_0+chi_far=1`, with `chi_0` supported before the endpoints.
For the near kernel, subtract the common logarithmic fundamental solution.
The inverse is order minus two; no extra smoothing is inserted. A more robust
operator proof uses nested collars. If `R_delta^s=(H_delta^s)^(-1)`,
localization and the exact resolvent identity give

    chi(R_delta^s-R_0^s)chi
      =chi R_delta^s(H_0^s-H_delta^s)R_0^s chi
        +S_(delta,chi),                                  (18a)

where `H_delta^s-H_0^s:H^5->H^3` is `o(1)` by (10)--(11). On a slightly
larger collar, (16), the exact global decaying Green response, and interior
parameter-elliptic regularity give uniform localized `H^3->H^5` bounds.
The commutator source is separated from the observation collar, so its exact
kernel is smooth. The chord square may there be differentiated through order
five, giving `||S_(delta,chi)||_(H^3->H^5)=o(1)`. Thus (18a) is an
order-minus-two operator-norm argument, not differentiation of the pointwise
limit (18).

For the far angular kernel integrate in the original periodic angle
`alpha=delta t`. At `alpha=+/-pi`, the amplitude and each differentiated
amplitude agree, while
`exp(-i n pi)=exp(+i n pi)=(-1)^n`; hence every endpoint term cancels for the
integer toroidal harmonic. Repeated integration by parts supplies uniform
smooth-kernel seminorms tending to zero after first increasing `T` and then
letting `delta->0`. Therefore

    ||chi(R_delta^s-R_0^s)chi||_(H^3->H^5)->0.           (18b)

Applying the exact curl row (7) upgrades (18b) to the order-minus-one
`H^3->H^4` Hodge comparison. This proves differentiated singular-kernel
operator convergence rather than merely pointwise `K_0` convergence.

The pieces outside a fixed collar do not create a missing boundary term.
Toward the remote axis `y_1=-delta^(-1)`, regularity of the helical components
is

    v_+=O(r^|n+1|),
    v_-=O(r^|n-1|),
    v_z=O(r^|n|).                                       (19)

The angular barrier and (16) give exponential Agmon decay across the distance
`Theta(delta^(-1))`. More explicitly, a cutoff energy identity for (10) with
weight `exp(c*dist(S,y))`, stopped before `q` changes by a fixed factor, gives
`C exp(-c/delta)` at the remote axis. On the outward side the decaying
harmonic solutions have
the reciprocal powers `r^(-|n+s|)` before the true far-field multipole region;
matching this exact exterior solution at the stopped collar gives the same
`C exp(-c/delta)` bound beyond it. The axial far field follows directly from
(17). Derivatives through order four obey the same estimate by the homogeneous
equation and local elliptic regularity. Thus

    lim_(M->infinity) sup_(small delta)
       E_tail,M(Q_delta B_(R3,n_delta)Q_delta^(-1)eta)=0 (19a)

uniformly for `||eta||_(H^3)<=1`, with the analogous column estimate. Hence,
for compact sources,

    ||Q_delta B_(R3,n_delta)Q_delta^(-1)-B_col(k_delta)||
       _(H0^3(S)->V_(delta,k)^4)=epsilon_delta,
    epsilon_delta->0.                                  (20)

The same kernel argument applied to (13) gives the pressure/Leray comparison.
Equation (20) includes the axis, infinity and exterior energy tail. It does
not assert a global inverse on arbitrary exterior forcing.

There is one low-frequency row that must not be hidden in this statement. In
the `m=+/-1` column sectors, one helical transverse component has shifted
angular number zero, so its scalar exterior response is `K_0(kr)` and an
uncancelled monopole would grow like `log(1/k)`.  The uniform estimate is
proved directly on the fixed-support `q`-space, without choosing a curl
potential.

For every smooth compact DA curl, Stokes' theorem gives

    integral_(R2) q(x_perp,0)dx_perp=0.                 (20a)

The integral functional is continuous on `H^3` for support in the fixed core,
so (20a) survives closure and the finite-row projection, whose correcting
profiles are themselves chosen from the same compact DA-curl space.  For a
general axial fiber `k`, transverse divergence reads

    i xi dot qhat_perp(xi,k)+ik qhat_z(xi,k)=0.          (20b)

The axial component still has zero transverse mean by the curl/Stokes
identity.  Differentiating (20b) at `xi=0` gives

    qhat_j(0,k)=-k partial_(xi_j)qhat_z(0,k),  j=1,2.    (20c)

Since the support lies in one fixed compact set `S`, Cauchy--Schwarz gives

    |partial_xi qhat(0,k)|
       <=integral_S |x_perp||q|dx_perp
       <=C_S||q||_(L2),

and the mean-value theorem gives the same bound for
`qhat(xi,k)-qhat(0,k)`.  Hence, with no gauge or potential ambiguity,

    |qhat(xi,k)|
       <=C_S(|xi|+|k|)||q||_(H3).                       (20d)

The zero-shift helical scalar is exactly the zero angular Fourier coefficient
of `q_x+i q_y` or `q_x-i q_y` (the choice depends on whether `m=-1` or
`m=+1`): `q_x+/- i q_y=exp(+/- i beta)(q_r+/- i q_beta)`.
Its radial monopole is therefore the corresponding Cartesian component of
`integral q dx_perp`, so (20a)--(20d) cancel precisely the scalar `K_0`
monopole, not merely the total vector after summing channels.

At low transverse frequency the Hodge multiplier obeys
`|vhat|<=|qhat|/sqrt(|xi|^2+k^2)`.  Equation (20d) makes it uniformly bounded;
fixed support controls the low-frequency integral and the standard order-minus-
one estimate controls high frequency.  This proves the uniform low-frequency
part of the `H^3->H^4` estimate.  The cancellation is distinct from the radial
`m=1` translation pencil, whose exterior velocity has the genuine `K_1` DtN
and whose `k^2 log k` term is retained in (47).

At the physical support edge, the fixed-carrier jet
`zeta_delta~d^p`, `p>=6`, makes the zero-extended generator coefficients and
DA curls `H^3`; their traces through order two vanish. The whole-space fields
from (17) pass through that coefficient interface with the usual global
Sobolev traces. Splitting the kernel for estimates and recombining it before
the final norm preserves the core/exterior cancellation.

## 4. Ambient DA closure and the genuinely common column domain

Let `K` be the fixed meridional vorticity support. For fixed `n!=0`, define

    X_n^3=closure_(H^3(R3)) Range(C_n),
    C_n xi=-[xi,zeta partial_theta],                     (21)

with zero extension and the physical finite rows. The boundary jet gives

    Range(C_n) subset H_0^3(K;C^3) intersect ker div.   (22)

The reverse density does not invert `C_n` at the edge. If
`q in C_c^infinity(K interior)` is divergence-free, action--angle components
give

    q^I=i n zeta xi^I,
    q^beta=i n zeta xi^beta,
    q^theta=i n zeta xi^theta-zeta'(I)xi^I.             (23)

Since `zeta>0` on the compact support of `q`, (23) defines a compact smooth
`xi`. Moreover

    div q=i n zeta div xi,                              (24)

so `xi` is divergence-free. Harmonic-preserving Bogovskii/profile corrections
give density with the finite rows. Therefore

    X_n^3=H_0^3(K;C^3) intersect ker div
           intersect {finite rows},                    (25)

for every fixed nonzero `n`. This proves that the support edge adds no passive
exterior-vorticity sector to the ambient DA closure.

For the limiting column define `T_col q=-[W_col,q]` and

    D_col={q in X_col^3:T_col q in X_col^3},
    ||q||_(D_col)=||q||_(X_col^3)+||T_colq||_(X_col^3). (26)

For `k in [k_-,k_+]`, write `A_col(k)=T_col+K_col(k)`. The compact-source
Hodge estimate gives

    sup_(k in [k_-,k_+])||K_col(k)||_(X_col^3->X_col^3)<infinity. (27)

Thus bounded perturbation gives the exact domain identity and two-sided graph
equivalence

    D(A_col(k))=D_col,
    c||q||_(D_col)
      <=||q||_X+||A_col(k)q||_X
      <=C||q||_(D_col),                                 (28)

uniformly in `k`. This discharges the noncircular column-domain obligation in
the 0074 README.

The ordinary Hanzawa map is not by itself an admissible transport-domain
identification: a small transverse component differentiates `q` in `I`, which
the anisotropic graph norm does not control.  The exact exposing sequence is,
in a regular cell with `T_0=omega(I)partial_beta`,

    q_N=N^(-3)chi(I)exp(iNI)e_theta.                     (28a)

After a same-cell divergence/finite-row correction of lower order,
`||q_N||_(H^3) asymp 1` and `T_0q_N=0`, whereas on any subcell where
`a(I,beta)!=0`,

    ||delta a partial_Iq_N||_(H^3) asymp delta N.        (28b)

Thus a transverse first-order remainder is not even bounded
`D(T_0)->H^3`, however small its `C^4` coefficient is.  The needed map comes from the
actual flow.  With physical meridional area form

    mu_delta=r dr wedge dz,

define `I_delta` to be enclosed `mu_delta` area divided by `2*pi`, increasing
from the elliptic center to the support edge, and let `beta_delta` be elapsed
flow time multiplied by `2*pi/T_delta(I)`.  The Hamiltonian identity
`i_(W_delta)mu_delta=-dP_delta` and the coarea formula give

    mu_delta=d beta_delta wedge dI_delta,
    W_delta=omega_delta(I_delta)partial_(beta_delta).    (29a)

This orientation makes `i_(W_delta)mu_delta=-dP_delta` when `I_delta`
increases from the maximum at the center toward the boundary.  Reversing the
angle reverses both displayed signs and changes no spectrum.

The symplectic Morse lemma supplies the regular-center extension; the nonzero
boundary normal derivative supplies the regular-edge extension.  Normalize
`x=I_delta/I_(b,delta)`.  The map matching equal `(x,beta)` between the column
and ring cores then satisfies

    phi_delta^*mu_delta=c_delta mu_0,
    c_delta=I_(b,delta)/I_(b,0),                         (29b)

and is `C^5` on the closed core.  This constant-Jacobian statement is not a
false assertion that differently sized physical carriers have equal volume.

The center statement is constructive.  A parameter-dependent Morse map first
writes `P_c-P=h_delta(rho^2)` and
`mu_delta=f_delta(rho,theta)rho drho wedge dtheta`.  At each `rho`, replace
`theta` by the normalized cumulative integral of `f_delta` around the circle;
the `drho` part of this angle change drops out when wedged with `drho`.  Then

    I_delta(rho)=(2*pi)^(-1)integral_(disk rho)mu_delta

gives `mu_delta=dbeta_delta wedge dI_delta`.  Smooth polar compatibility at
`rho=0`, the uniformly nondegenerate Hessian, and the `C^5` coefficient bounds
give a `C^5` Cartesian extension with uniform inverse bounds.  On the regular
annulus the same `beta_delta` is the normalized flow time; equality on an
overlap fixes the phase.  This supplies the isochore-Morse map rather than
assuming one by name.
The divergence-preserving identification is

    U_delta^P q=c_delta^(-1)(phi_delta)_*q,              (29c)

followed by the exact cylindrical component rotation and a convergent finite-
row correction.  It is distinct from the half-density unitary used to derive
(17).  In particular, the transition includes composition by a diffeomorphism
and is not merely a matrix multiplier.

Here is the required common extension.  The nested-ball Schauder iteration
below gives `C^(7,alpha)` convergence of the carrier on a fixed ball.  The
parameter-dependent isochore Morse construction at the center, the normalized
flow-time construction on the regular annulus, and the Hopf boundary chart
therefore give a core map `phi_delta` with uniform `C^5` bounds and
`||phi_delta-id||_(C^4)->0` after the fixed centering/scaling.  On overlaps,
two canonical flow angles differ by a function of `I`; changing one phase by
that function glues them without changing `dbeta wedge dI` or (29a).

Extend `phi_delta` across a two-sided normal collar by matching its boundary
value and first normal jet and applying a one-dimensional positive-Jacobian
interpolation.  An isotopy extension then makes it the identity in the
outward far chart.  On the inward side choose a monotone map

    h_delta:R->(-delta^(-1),infinity)                    (29c.1)

which equals the ring-centered coordinate on every fixed collar and tends to
`-delta^(-1)` as its argument tends to `-infinity`; join it to the collar by
the same isotopy.  This gives a global orientation-preserving meridional
diffeomorphism `Phi_delta` from the reference plane to the open physical
half-plane.  It is uniformly `C^5` and converges in `C^4` on every fixed
collar.  Its derivatives need not be uniformly invertible at the negative
end; that end is the physical axis and is controlled by (19), not by a false
global Euclidean norm equivalence.

Pulling the exact scalar/vector system (4)--(14) back by `Phi_delta`, including
its Jacobian, cofactor and cylindrical component rotation, gives an elliptic
operator `H_delta^Phi`.  On nested fixed collars,

    H_delta^Phi-H_0:H^5->H^3=o(1),                       (29c.2)

because all metric coefficients converge in `C^4`.  The localized resolvent
identity (18a), now written for `H_delta^Phi`, proves the same
`H^3_comp->H^5_loc` scalar estimate.  Exact curl (7) gives
`H^3_comp->H^4_loc` for the Hodge map.  The collar commutators are supported
where `Phi_delta` is still uniformly regular.  Beyond that collar the exact
integer-`n` kernel cancellation and the regular-axis/exterior estimate
(19)--(19a) give the common `E_tail,M` bound. Consequently the version of
(20) used below is, with `B_delta^#` denoting the full Piola/metric pullback
under this same `Phi_delta`,

    ||B_delta^#-B_col(k_delta)||_(H0^3(S)->V^4)
       =epsilon_delta^Phi ->0.                           (29c.3)

No composition operator has been omitted. The metric Leray formula (13)
gives the identical pressure conclusion.

For completeness, before applying (29a), the exact scaled cylindrical Lie
brackets for an axisymmetric meridional field `W=(W_r,0,W_z)` and axial
vorticity `zeta e_theta` are

    a[W,q]_r=W_r partial_1q_r+W_z partial_2q_r
               -q_r partial_1W_r-q_z partial_2W_r,
    a[W,q]_theta=W_r partial_1q_theta+W_z partial_2q_theta
               -delta q_geom^(-1)W_r q_theta,
    a[W,q]_z=W_r partial_1q_z+W_z partial_2q_z
               -q_r partial_1W_z-q_z partial_2W_z,       (29d)

and, for `v=Bq`,

    a[v,zeta e_theta]_r=-i k q_geom^(-1)zeta v_r,
    a[v,zeta e_theta]_theta
       =v_r partial_1zeta+v_z partial_2zeta
          -delta q_geom^(-1)zeta v_r
          -i k q_geom^(-1)zeta v_theta,
    a[v,zeta e_theta]_z=-i k q_geom^(-1)zeta v_z.        (29e)

Here `q_geom=1+delta y_1`; it is not the perturbation `q`.  Equations
(29d)--(29e), together with (4)--(13), are the full vector transport,
stretching, curl and pressure ledger.  In particular, the stretching block is
zeroth order after the `H^3->H^4` Hodge gain.

The nested-ball bootstrap in corrected 0066 can be iterated without crossing
a nonsmooth boundary equation: `(t_+)^p` is `C^5` for `p>=6` and the Cao
equation is posed on the ambient meridional plane. Starting from the proved
`C^(2,alpha)` convergence, the sequence

    C^(2,alpha) profile -> C^(1,alpha) source
      -> C^(3,alpha) profile -> C^(2,alpha) source
      -> C^(4,alpha) profile -> C^(3,alpha) source
      -> C^(5,alpha) profile -> C^(4,alpha) source
      -> C^(6,alpha) profile -> C^(5,alpha) source
      -> C^(7,alpha) profile.                            (29e.1)

On successively smaller balls this follows from composition by the `C^5`
nonlinearity and interior Schauder estimates. Thus the carrier velocity
converges in `C^(6,alpha)`.
The quantitative isochore-Morse/flow/edge construction above then gives
parameter-uniform `C^5` maps and `C^4` convergence of (29b), of
`omega_delta`, and of the edge chart. Thus, in the common action variables,

    (U_delta^P)^(-1)[W_delta,U_delta^Pq]
      =[omega_delta(x)partial_beta,q],                  (29f)

and the principal difference is exactly

    (omega_delta-omega_0)partial_beta q
      -q^x partial_x(omega_delta-omega_0)partial_beta.   (29g)

There is no transverse `partial_x q` remainder.  Equation (29g) maps
`D_col->X` with norm `o(1)` by definition of the transport graph.  All
coordinate/frame terms left in (29e) are zeroth order or contain only the
bounded multipliers `k`, `k+delta`, and `k-delta` already displayed in
(4)--(12).

Define

    epsilon_prof(delta)
      =||omega_delta-omega_col||_(C^4([0,1]))
        +||zeta_delta^#-zeta_col||_(C^4),
    epsilon_H(delta)
      =||B_delta^#-B_col(k_delta)||_(H0^3->V^4),
    epsilon_row(delta)
      =||P_delta^fin-P_col^fin||.                        (29h)

The profile bootstrap and flow chart give `epsilon_prof->0`; (20) gives
`epsilon_H->0`. The finite row matrices converge because their physical
moment, harmonic and KKS functionals are continuous under (29c) and the
normalization (3). For `q in D_col`, the full difference is

    (Ahat_delta-A_col(k_delta))q
      =-[ (omega_delta-omega_col)partial_beta,q]
       -[(B_delta^#-B_col)q,omega_col]_0
       -[B_delta^#q,omega_delta^#-omega_col]_0
       +R_(curv,delta)(B_delta^#q,omega_delta^#)
       +R_(row,delta)q.                                  (29i)

The first row maps `D_col->X` with norm `C epsilon_prof` by (29g).  The next
three rows map `X->X` with norm
`C(epsilon_H+epsilon_prof+delta)`, using (29e), the `H^3->H^4` Hodge gain and
`C^4` carrier multiplication.  Every `n_delta` derivative is one of the exact
bounded factors `k_delta`, `k_delta+delta`, or `k_delta-delta`; no hidden
`delta partial_theta` is placed in the remainder. The last row has norm
`epsilon_row`. Thus the
bounded isomorphisms

    U_delta:X_col^3->X_(delta,n_delta)^3,
    U_delta:D_col->D(A_(delta,n_delta)),                 (29)

obey uniform two-sided graph bounds, and the leading massive graph convergence

    ||(U_delta^(-1)A_(delta,n_delta)U_delta
          -A_col(k_delta))q||_(X_col^3)
      <=rho_delta||q||_(D_col),
    rho_delta
      =C[epsilon_prof(delta)+epsilon_H(delta)
          +epsilon_row(delta)+delta+exp(-c/delta)]
      ->0.                                               (30)

Every `delta*n_delta` contribution is retained inside (4)--(12) and
`A_col(k_delta)` before (30) is estimated. Equation (30) is a leading
graph/norm-resolvent theorem, not the unproved second-order `HJ2` expansion of
0054 and not a product two-index commutator estimate.

## 5. Route B: converse transport spectrum on the ambient q-space

The 0066 construction already gives the Weyl inclusion. Equality follows
directly on the ambient space (25), without dividing by `zeta` at the edge.
The physical streamline flow is periodic on every level, with one elliptic
center, no separatrix, and

    0<omega_min<=omega(I)<=omega_max.                    (31)

In the volume action--angle frame, the `m`-th streamline block of
`T_n=-[W,·]` is

    T_(n,m)q=-i m omega(I)q+omega'(I)q^I partial_beta.   (32)

The second term is nilpotent of square zero. Equivalently, along a physical
closed orbit of period `tau(I)`, the monodromy is the unipotent derivative of
the flow. For

    lambda notin Sigma_T
      :={0} union closure{-i m omega(I):m in Z\{0}},     (33)

the periodic characteristic inverse exists because
`exp(lambda tau)I-Dphi_tau` is invertible. In (32) it is simply a scalar
denominator plus its squared-denominator nilpotent correction. On compact
subsets of the complement of `Sigma_T`, these denominators and their first
three `I` derivatives are uniformly bounded over all `m`; for large `|m|`
the bound improves like `|m|^(-1)`.

At the center this is not an appeal to singular polar coordinates. On a center
disk, write the characteristic inverse directly using the smooth physical
flow `phi_t` and its derivative. The nondegenerate elliptic fixed point gives
`phi_t(y)=exp(t J D^2P(0))y+O(|y|^2)` with derivatives through order four;
the denominator stays uniformly invertible away from `Sigma_T`. Composition
and multiplication therefore preserve the Cartesian `H^3` regular-center
rows. At the outer edge, `W` is `C^4` and tangent. Its flow preserves a smooth
defining function up to a positive `C^3` factor, so composition/pushforward
preserves the vanishing normal traces of orders zero, one and two that
characterize `H_0^3(K)`. The time integral and monodromy inverse preserve the
same traces. Hence the characteristic inverse maps `X_n^3` into `D(T_n)` and
proves the converse Fredholm inclusion

    sigma_ess(T_n)=Sigma_T.                             (34)

The `m=0` block supplies the infinite-dimensional point `0`; every nonzero
band has the singular sequences constructed in corrected 0066. Finite rows
are finite-codimensional and do not alter (34), although they can alter the
ordinary point spectrum. No equality for the full `sigma(T_n)` is claimed.

For each fixed `n`, `K_n=-[B_n·,omega]` is compact on `X_n^3`: common compact
curl support gives zero integral and controlled first moments, hence removes
the low-frequency Biot--Savart singularity; global `H^3->H^4` gain, restriction
to the coefficient support, fixed-`n` toroidal multiplication and Rellich give
compactness. Fredholm essential spectrum is invariant under this compact
perturbation, so

    sigma_ess(A_n)
      ={0} union closure{-i m omega(I):m!=0}.             (35)

This proves, for every fixed sufficiently thin Cao carrier and every fixed
nonzero toroidal harmonic, the genuine Fredholm essential gap

    lambda=i nu,       0<|nu|<omega_min.                 (36)

It does not say that (36) is free of discrete eigenvalues; those are precisely
the possible Kelvin modes.

## 6. Intermediate full-cluster contour before sector isolation

The limiting compact-column `m=0` Sturm--Liouville problem has simple
positive-Krein frequencies `sigma_j(k)->0`. Choose fixed `J` so

    0<sigma_J(k_0)<omega_min/3.                          (37)

By (35) for the limiting compact DA column, this point belongs to a finite-
rank isolated full-column spectral cluster, even if a nonzero-`m` eigenvalue
happens to coincide with it. Indeed, in the gap `T_col-lambda` is invertible
and

    A_col(k)-lambda
      =[1+K_col(k)(T_col-lambda)^(-1)](T_col-lambda)     (37a)

is Fredholm of index zero. The compact factor is analytic in `k` for
`k in [k_-,k_+]`; analytic Fredholm theory makes the gap spectrum discrete
with finite algebraic multiplicity. Fix `Gamma_J` around the **full**
`k_0` cluster, not merely its `m=0` line, and no other spectral point.
Common-domain analyticity and compactness of `Gamma_J` give

    sup_(|k-k_0|<=eta) sup_(z in Gamma_J)
      ||(z-A_col(k))^(-1)||_(X->D_col)<=C_Gamma          (37b)

after reducing `eta` if necessary. Combining (30) with (37b) gives the actual
relative bound

    sup_(z in Gamma_J)
      ||(Ahat_delta-A_col(k_delta))
          (z-A_col(k_delta))^(-1)||
      <=rho_delta C_Gamma<1.                            (38)

for small `delta`. The exact resolvent identity then gives graph-norm
resolvent convergence and

    ||U_delta^(-1)Pi_(delta,J)U_delta-Pi_(col,J)||
       _(X->D_col)->0,
    rank Pi_(delta,J)=rank Pi_(col,J)>0.                 (39)

At this intermediate rung the combined routes construct only a nonzero actual
Cao Riesz cluster in the massive harmonic `n_delta`.  The real and reversible
conjugation rows transfer, but an unknown coincident streamline sector would
leave its KKS inertia unresolved.  Section 7 performs the missing isolation;
this paragraph is not the final Route-A verdict.

## 7. Whole-column isolation and two definite modes

The full-cluster statement in (39) is not the Route-A target.  We now isolate
the `m=0` eigenlines rather than declaring an unknown-sector cluster a
success.

### 7.1 Uniform removal of large streamline harmonics

Decompose the column operator by the integer streamline harmonic `m`, keeping
the axial Fourier number `k in [0,k_1]`.  In the action frame,

    T_m-im omega(I)

is a triangular multiplication operator whose only off-diagonal entry is
`-omega'(I)q^I partial_beta`; in the `m` block this entry is multiplication,
not an additional derivative.  For `|lambda|<=omega_min/4`, direct inversion
of this triangular matrix and three `I` derivatives gives

    ||(T_m-lambda)^(-1)||_(X_m^3->X_m^3)
       <=C/|m|,             |m|>=1.                      (40)

The parameter-elliptic estimate for the full column div--curl system gives,
uniformly in `m` and `0<=k<=k_1`,

    ||B_(m,k)q||_(H^4)<=C||q||_(H^3),
    ||K_m(k)q||_(H^3)<=C||q||_(H^3).                    (41)

To see the uniformity rather than assume it, diagonalize the transverse
velocity into `v_r+/-iv_beta`.  The scalar div--curl forms then contain
`(m+/-1)^2/r^2+k^2` and the axial component contains `m^2/r^2+k^2`.
Testing the equations, differentiating through three physical covariant rows,
and using the regular-axis polar conditions gives the first estimate in (41).
The apparent factors `m` in the cylindrical stretching row consume precisely
the one derivative gained by `B_(m,k)` and do not make the second constant
grow.  The exterior is homogeneous and its decaying Bessel response obeys the
same form estimate.  For `m=+/-1`, the only zero-shift helical channel is
uniform at `k=0` by the exact compact-curl/zero-moment calculation (20a), not
by a false uniform bound for arbitrary `K_0` sources.  Thus (41) is the
physical core-plus-exterior Hodge bound.
Therefore

    ||K_m(k)(T_m-lambda)^(-1)||<=C/|m|.                 (42)

Choose `M` so the right side is below one for `|m|>=M`.  The Neumann
factorization of `A_m-lambda` proves a uniform punctured neighborhood of zero
free of the entire large-`|m|` spectrum.  This is stronger than excluding only
large-harmonic eigenvectors.

### 7.2 The finite harmonics and the genuine translation row

At `k=0`, the zero-frequency equation in any nonzero `m` sector reduces
exactly to the planar steady equation.  If `phi` is the perturbation
streamfunction, `q=-Delta phi`, and `W=f(U)=U_+^p`, then

    {U,q}+{phi,W}=0
      implies {U,q-f'(U)phi}=0.                          (43)

The last bracket has no nonzero angular harmonic, hence

    (-Delta-pU_+^(p-1))phi=0.                            (44)

Cao's quoted Lane--Emden nondegeneracy says the bounded kernel of (44) is
exactly `span{partial_1U,partial_2U}`.  Thus every `|m|>=2` block is invertible
at zero.  The `m=+/-1` kernel is precisely physical translation; rotations
and Casimir/stabilizer directions lie in `m=0` and are not silently removed
as generic finite rows.

For each of the finitely many `2<=|m|<M`, the regular-axis row and decaying
exterior Bessel row have a common form domain as `k->0`; their DtN maps
converge to the `r^(-|m|)` row.  The Hodge block is norm-continuous on that
domain.  Hence (44) and bounded inversion leave these blocks zero-free on one
common interval `0<k<k_f`.

The translation block needs its actual response rather than quotienting away
the `k!=0` bending mode.  In the `m=1` scalar convention

    T(b)=-partial_r^2-r^(-1)partial_r+r^(-2)
          +W'(r)/[r(Omega(r)-b)],                        (45)

with exterior row `psi'(a)=-psi(a)/a`, the exact kernel is
`psi_0=r Omega`.  The ground-state transform gives a simple zero and a
positive complement, while the sign is

    <psi_0,T'(0)psi_0>_(rdr)
      =integral_0^a r^2W'(r)dr=-2F(a)<0.                (46)

For `k!=0`, the exact `K_1` exterior row
`-1/D=1+(ka)^2[log(ka/2)+gamma_E]+O(k^4log^2k)` gives the
translation Feshbach determinant

    d_1(k,0)=F(a)^2k^2log(ka)+O(k^2).                   (47)

It is nonzero for every sufficiently small positive `k`; equivalently the
translation eigenvalue moves to

    b_1(k)=-(F(a)/2)k^2log(1/(ka))+O(k^2) !=0.          (48)

The conjugate `m=-1` block has the same conclusion.  Equations (43)--(48)
are the nonidentity witness: no blanket zero-kernel assertion was made before
the true translation row was separated.

Choose one `k_0 in (0,min(k_f,k_1))`.  Equations (42)--(48) prove
`0 in rho(A_m(k_0))` for every `m!=0`.  The essential bands of each such block
are at distance at least `|m|omega_min` from zero; analytic Fredholm
discreteness then gives a number `epsilon_*>0` such that

    sigma(A_m(k_0)) intersect {|lambda|<epsilon_*}=empty,
                 every m!=0.                            (49)

Only a finite minimum is taken after (42), so `epsilon_*` is genuine even
though the `m=1` bending eigenvalue can be very small.

### 7.3 Two `m=0` modes and the exact KKS scope

The simple positive axisymmetric Sturm--Liouville frequencies satisfy
`sigma_J(k_0)->0`.  Choose two fixed, distinct, sufficiently large indices
`J_1,J_2` such that

    0<sigma_(J_2)<sigma_(J_1)<epsilon_*/3.              (50)

They are now isolated from every `m!=0` point and essential spectral value,
and simplicity of the Sturm problem separates them from each other.  Draw
disjoint contours `Gamma_1,Gamma_2` around the full-column eigenvalues
`-i sigma_(J_1)` and `-i sigma_(J_2)`.  The projectors in (37)--(39) on these
contours have complex rank one.

The Fourier convention (2) and (3) fixes the kinetic norm
`rho_0 integral |v_J|^2rdrdz`; it is not identified with the orbit Hessian.
The reviewed 0048/0053 Sturm wave calculation establishes positive
constrained-energy/Krein **sign** for every nonzero `m=0` mode, but it did not
carry the density, axial Fourier volume, real/complex factor and second
orbit-variation term far enough to identify its displayed Sturm normalization
with the absolute physical relative Hessian `E_2`.

Let `h_J=E_2(e_J,conjugate(e_J))>0` denote that actual, convention-dependent
value.  Rescaling the eigenvector by `h_J^(-1/2)` is legitimate but is not a
calculation of `h_J` from the Sturm amplitude.  In this energy-unit convention,

    E_2(e_J,conjugate(e_J))=1,
    A e_J=-i sigma_J e_J,                               (51)

and the right-reduced Hamiltonian identity
`Omega_KKS(Ae,conjugate(e))=E_2(e,conjugate(e))` gives

    Omega_KKS(e_J,conjugate(e_J))=i/sigma_J,
    |Omega_KKS|=1/|sigma_J|.                            (52)

Equation (52) is an exact energy-unit identity; it is not an independently
derived absolute physical normalization of the original Sturm eigenfunction.

Only the sign is needed for the contour transfer.  For each fixed selected
`J`, the first-order column eigen-equations reconstruct a smooth DA
displacement on the compact core because `sigma_J!=0`; the apparent edge
division cancels against the `W,W'` factors in those equations.  The action-
flow Piola map and graph-Riesz convergence carry these representatives to the
ring.  In the physical KKS integral

    Omega_KKS(Cxi,Ceta)=rho_0 integral omega dot (xi cross eta)dx, (52a)

the carrier, representatives, positive Jacobian and normalized fiber factors
then converge on each finite-dimensional range.  Hence its nonzero sign is
locally constant and both actual Cao modes retain the positive Krein sign.
This proves sign transfer, but leaves the absolute conversion between the
Sturm normalization and (52a) open.  Equations (30), (38) applied separately
to `Gamma_1,Gamma_2` preserve rank; Hamiltonian/reversible symmetry keeps each
simple definite eigenpair on the imaginary axis.

Therefore Route A constructs two isolated positive-Krein modes on every
sufficiently thin actual Cao ring at the massive harmonic
`n_delta=nearest(k_0/delta)`, together with their real `+/-n_delta`
representatives.  What remains for the 0051 supplier is not spectral
existence or sign: it is the absolute physical KKS/Hessian conversion needed
for a common analyzer normalization, frequency matching, and evaluation of
two noncommuting physical deformation couplings on these particular modes.
No nonlinear branch, restoring
neighborhood, P2/LP2, particle identity or quantum inference follows here.

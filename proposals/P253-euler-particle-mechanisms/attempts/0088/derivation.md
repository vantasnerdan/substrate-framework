# Physical Hessian-row compression and the first response obstruction

## 1. Conventions and the physical orbit form

Let `[u,v]=u dot grad(v)-v dot grad(u)`.  In the right-reduced convention
used by P253/0005, a displacement `eta` generates

    C eta=-[eta,bar_omega]=curl(eta cross bar_omega),       (1)

and

    Omega(C eta,C xi)
      =rho_0 integral bar_omega dot (eta cross xi) dx,      (2)

with `i_(X_H) Omega=dH`.  The frozen 0088 README writes a control as
`h=[xi,bar_omega]`; this is (1) with generator `eta=-xi`.  Both generators
change sign, so (2), expressed in the README representatives, is still

    Omega(h_1,h_2)
      =rho_0 integral bar_omega dot (xi_1 cross xi_2) dx.   (3)

Thus no sign is changed silently when passing between the canonical orbit
map and the README control convention.

Use physical orthonormal cylindrical components and write the Cao
streamfunction convention as

    bar_omega=omega_theta e_theta=r zeta(r,z)e_theta.     (3a)

Use
`phi_n=(2 pi)^(-1/2)exp(i n theta)`.  If
`xi=hat_xi(r,z)phi_n` and `eta=hat_eta(r,z)phi_n`, then

    Omega([xi,bar_omega],[conjugate(eta),bar_omega])
      =rho_0 integral zeta r^2
        (xi_z conjugate(eta_r)-xi_r conjugate(eta_z)) dr dz. (4)

There is no extra `2 pi` in (4): it has been consumed by the normalized
characters.  In the fixed meridional variables, put

    tilde_xi=(R a^2 q)^(1/2) hat_xi,
    q=1+delta y_1.                                        (5)

Since `r dr dz=R a^2q dy`, the two half-density factors cancel the volume
Jacobian but not the additional `r` in `omega_theta=r zeta`.  Thus

    Omega=rho_0 integral r zeta
      (tilde_xi_z conjugate(tilde_eta_r)
       -tilde_xi_r conjugate(tilde_eta_z))dy.             (5a)

The residual `r=Rq` is load bearing.  Equations (4)--(5a) pin the density,
Fourier, core-volume, and physical-component conventions and agree with the
unnormalized-character row `2 pi rho_0 integral r^2 zeta(...)dr dz`.

## 2. The actual relative Hessian and absolute normalization

Set

    F(omega)=H(omega)-c_g P_z(omega),
    H=(rho_0/2) integral omega dot (-Delta)^(-1)omega dx,
    P_z=(rho_0/2) integral (e_z cross x) dot omega dx.     (6)

For compact vorticity its ambient first derivative is the one-form

    mu_F=rho_0(-Delta)^(-1)bar_omega
         -(rho_0 c_g/2)(e_z cross x).                     (7)

If `q=[xi,bar_omega]` and `p=[eta,bar_omega]`, differentiating two
push-forwards and symmetrizing gives the physical orbit Hessian

    E2_phys(q,p)
      =rho_0 integral Bq dot Bp dx
       +(1/2)<mu_F,
          [xi,[eta,bar_omega]]+[eta,[xi,bar_omega]]>.     (8)

The second line is the orbit-variation/energy--impulse term.  It is the term
lost by identifying `E2_phys` with kinetic energy.  At the relative
equilibrium it is independent of the displacement representatives modulo
the stabilizer and is the relative Hessian on the reduced leaf.

Let `tilde e_a` be either of the two reviewed 0074/0078 positive-Krein
eigenlines, reconstructed with the full physical Piola/Hodge map and a
displacement `tilde xi_a`.  Define, without a fitted constant,

    H_a=E2_phys(tilde e_a,conjugate(tilde e_a))>0,
    e_a=tilde e_a/sqrt(H_a),
    xi_a=tilde xi_a/sqrt(H_a).                            (9)

Equations (4)--(9), with the residual radius in (5a), are the absolute
Sturm-to-physical normalization
algorithm: `rho_0`, normalized `theta`, the half density, and the second
orbit variation are all inside the displayed integrals.  For a real field
`q=z e_a+conjugate(z e_a)`, Hessian invariance first gives
`E2(e_a,e_a)=0`: indeed `E2(Ae_a,e_a)+E2(e_a,Ae_a)=0` and
`nu_a!=0`.  Together with Fourier orthogonality this gives

    (1/2)E2_phys(q,q)=|z|^2.                              (10)

No real/complex factor remains implicit.

The reviewed fixed-fiber gap also supplies the regularity needed here.  On a
fixed Cao member, the simple mode satisfies

    e=-(T+i nu)^(-1)K e.                                 (11)

The characteristic inverse is bounded on the regular polar/edge rows, `K`
gains one Hodge derivative, and the `p>=6` compact-edge traces allow zero
extension.  One bootstrap therefore puts these fixed-fiber modes in
`D(A_phys) intersect H^4`; their Hodge velocities are locally `H^5`.
Likewise, for `xi in H^5_c(C)`, `h=[xi,bar_omega]` is `H^4`, and the two
brackets in `A_phys h` are `H^3`.  Hence the modes and compact controls used
below lie in the frozen smooth graph core `D_4`, without identifying the
full anisotropic generator domain with `H^4`.

## 3. Hessian-row left covector: exact sign and clock

On the reduced physical domain,

    E2_phys(q,p)=Omega(A_phys q,p),
    Omega(A_phys q,p)+Omega(q,A_phys p)=0.                (12)

At the double frequency choose the KKS frame continued from the two simple
incoming branches and normalize it by (9):

    A_phys e_a=-i nu_* e_a,
    A_phys conjugate(e_a)=+i nu_* conjugate(e_a),
    E2_phys(e_a,conjugate(e_b))=delta_ab,                 (13)

where `nu_*=Omega_a sigma_*`, not `sigma_*`.  From (12)--(13),

    E2_phys(q,conjugate(e_a))
      =Omega(A_phys q,conjugate(e_a))
      =-Omega(q,A_phys conjugate(e_a))
      =-i nu_* Omega(q,conjugate(e_a)).                  (14)

Therefore the exact physical left covector is

    ell_a(q):=E2_phys(q,conjugate(e_a))
             =-i nu_* Omega(q,conjugate(e_a)),           (15)

and `ell_a(e_b)=delta_ab`.  In particular

    Omega(e_a,conjugate(e_b))=i delta_ab/nu_*.            (16)

This proves the sign and core-clock factor requested in the frozen contract.
It also removes the need to solve a global `C_0^*` potential merely to define
the physical dual.  If a local `L2` density `a_a` is desired, it is any
solution, modulo gradients and the actual finite-row annihilator, of

    C_0^* a_a
      =-i nu_* rho_0 P_slice(conjugate(xi_a) cross bar_omega)
        +row covectors.                                  (17)

Equation (15), rather than an arbitrary smooth `a_a`, fixes which solution
of (17) is physical.

## 4. Covariant Feynman--Hellmann and invariant two-by-two statement

Let `nabla_h` include the fixed-slice derivative, or the full recentering
connection and Riesz derivative when the slice moves.  Differentiating
`(A_phys+i nu_a)e_a=0` and applying (15) gives

    ell_a((nabla_h A_phys)e_a)=-i D_h nu_a,              (18)

because `ell_a((A+i nu_*)nabla_h e_a)=0` and
`ell_a(e_a)=1`.  Equation (18) is an eigenfrequency derivative only along a
differentiable steady operator/eigenbranch.  For an arbitrary same-leaf DA
seed at a displaced, nonsteady background, its left-hand side remains the
exact instantaneous compression, but the notation `D_h nu_a` is not used.

For any real admissible seed define

    M(h)_ab=ell_a((nabla_h A_phys)e_b).                   (19)

The physical positive-sector metric is `G_s=E2_s`, and
`A_s^*G_s+G_s A_s=0`.  Differentiating at
`A_0|P=-i nu_* I` gives on `P`

    P(A'_0^*G_0+G_0A'_0)P
      +i nu_* P G'_0 P-i nu_* P G'_0 P=0.               (20)

Thus, for the full covariant/fixed-slice derivative in the same frame,

    M(h)^*=-M(h),                                        (21)

so its diagonal entries are imaginary and `iM(h)` is Hermitian.  Omitting a
speed, projector, or slice-connection row destroys the hypothesis behind
(20), and away from exact degeneracy the `G'_0` cancellation is unavailable.

The basis-invariant instantaneous control condition is consequently

    [M(h^(1)),M(h^(2))] !=0.                             (22)

A common unitary change of the continued KKS frame conjugates both matrices
and does not change (22).  In the same frame, exact toroidal selection makes
an `ell=n_1-n_2` seed off diagonal, while a toroidal-zero seed is diagonal.
If their nonzero entries are `c` and `-i(d_1,d_2)`, then

    [M_off,M_diag]_12=-i(d_2-d_1)c.                      (23)

Thus (22) is exactly the nonparallel-Pauli condition, without treating
`G_12` and `G_3` as basis-independent scalars.

There is also a useful Lagrange form of an off-diagonal row.  Differentiate
`E2=Omega(A .,.)` with a symplectic Riesz connection.  Since
`nabla Omega=0`,

    M(h)_ab=-i nu_* (nabla_h E2_phys)
                       (e_b,conjugate(e_a)).             (24)

For a non-symplectic coordinate connection, (24) acquires the explicit
`nu_*^2 (nabla_h Omega)(e_b,conjugate(e_a))` correction.  Equation (24) is
the exact third-relative-Hessian/Lagrange representation and automatically
retains the orbit-constraint term in (8).

## 5. Direct full-Hodge density and its exact scope

Suppose temporarily that (15) has an `L2` representative `a`, and put
`v=B e`.  With `V_h e=-[Bh,e]-[Be,h]`, divergence freedom and decay give

    ell(V_h e)=integral g_(e,a) dot h dx,                (25)

    g_(e,a)
      =B[-(grad e)^T a-(e dot grad)a]
        +(v dot grad)a+(grad v)^T a.                     (26)

The first term contains both the induced whole-space Hodge velocity and its
exterior tail; no artificial core boundary is introduced.  For the README
control convention `h=C_0 xi=[xi,bar_omega]`, a second integration by parts
gives

    integral g dot C_0 xi=integral xi dot C_0^*g,
    C_0^*g=(grad bar_omega)^Tg+(bar_omega dot grad)g.     (27)

Equations (25)--(27) independently reproduce the frozen dual-density
formula and its sign.  Gradients, stabilizer/Casimir rows, and slice rows are
then removed by the actual constrained projection, not by an unconstrained
Helmholtz slogan.

Let the finite displacement rows be `L_j`, and choose fixed smooth regular-
cell correction profiles `chi_j` with `L_i(chi_j)=delta_ij`.  The bounded
primal oblique projection is

    P_row xi=xi-sum_j L_j(xi)chi_j,                       (27a)

so `Range(P_row)` is the constrained displacement slice and
`||P_row||<=1+C_rows`, where `C_rows` is the norm of the row map times the
inverse moment matrix and the profile norms.  The correct dual witness is
`P_row^* C_0^*g`, not `P_row C_0^*g`.

If a selected component of
`curl(P_row^* C_0^*g)` has magnitude at least `gamma` on a meridional ball of
radius `r_p`, choose `psi=A phi((y-y_*)/r_p)`, a physical cylindrical curl
test `xi_0`, and set `xi=P_row xi_0`.  Adjointness gives the exact numerator

    <xi,C_0^*g>=<xi_0,P_row^*C_0^*g>,                    (27b)

so the row restoration cannot cancel the local response.  Then

    |integral xi_0 dot P_row^*C_0^*g|
       >=c_0 A gamma r_p^2,
    ||P_row xi_0||_(H5)
       <=(1+C_rows)C_0 A r_p^(-5).                       (28)

With local `C_0` and Hodge bounds `C_C0,C_Hodge`, (28) yields

    |G(C_0 xi)|/||C_0 xi||_Y
      >= c_geom gamma r_p^7/
          [C_C0 C_Hodge(1+C_rows)].                      (29)

The exponent seven is therefore the exact two-meridional-dimensional `H5`
count frozen in the README.  Equation (29) is a fixed-carrier physical
estimate.  In a core-rescaled family, the statement that angular derivatives
cost `|ell|/R=O(1/a)` avoids a coordinate `N^5` loss only after the physical
Piola and Sobolev scaling factors have been included in
`C_C0,C_Hodge,C_rows`; (29) by itself is not uniform in `N`.  It remains
conditional on the source-specific adjoint-projected curl witness.

## 6. An exact refutation inside the frozen diagonal slice

On a positive regular cell use volume action--angle coordinates
`dV=dI d beta d theta` and

    bar_omega=zeta(I)partial_theta.                       (30)

An `(ell,d)=(0,0)` displacement has coefficients depending only on `I`.
Divergence freedom says `partial_I xi^I=0`; compact support forces
`xi^I=0`.  Hence

    [xi,bar_omega]
      =xi^I zeta'(I)partial_theta-zeta partial_theta xi=0. (31)

After quotienting by `ker C_0`, the frozen compact-cell space
`Y_(0,0)^4` is zero.  Consequently

    G_3|_(Y_(0,0)^4)=0.                                  (32)

This is a route-scoped refutation, not a statement that the actual curved
ring has no diagonal response.  It identifies why an IVT frequency crossing
or the value-level Sturm derivative cannot supply the required compact
same-orbit seed: the only source-compatible zero/zero displacement is a
stabilizer.

The failure-derived physical route is toroidal order `ell=0` with a genuine
meridional harmonic `d!=0`.  In the straight axisymmetric column, a pure
`m=0` eigenmode makes its diagonal matrix element vanish by meridional
selection.  The exact cylindrical ledger in P253/0074 contains at first
order the factors `q^(-1)=1-delta y_1+...`, the connection `delta q^(-1)`,
and shifted helical masses `(k+/-delta)^2q^(-2)`.  Since `y_1` has
meridional character `+/-1`, every first-order row has character `+/-1`.
Thus the smallest possible diagonal supplier is exactly `d=+/-1`, and its
first permitted response order is `O(delta)`.

On the selected smooth source core, write the source-derived expansion

    A_delta=A_0+delta A_1+o(delta),
    V_delta[h]=V_0[h]+delta V_1[h]+o(delta).              (33a)

where `A_1,V_1` contain the differentiated full Hodge/Leray block,
cylindrical connection, carrier first jet, finite rows and Riesz connection.
In Kato gauge,

    e_a(delta)=e_a^0+delta e_a^1+O(delta^2),
    e_a^1=-R_a^Q(A_1+i nu_a^1)e_a^0,                    (33b)

where `R_a^Q=Q_a(A_0+i nu_a^0)^(-1)Q_a`.  Meridional selection kills the
order-zero diagonal row.  The exact first coefficient is

    D_curv=(d/d delta)|_(0)
      {M_11(h_(0,1))-M_22(h_(0,1))},                    (33)

    D_curv
      =[ell_1^1(V_0[h]e_1^0)+ell_1^0(V_1[h]e_1^0)
         +ell_1^0(V_0[h]e_1^1)]
       -[ell_2^1(V_0[h]e_2^0)+ell_2^0(V_1[h]e_2^0)
         +ell_2^0(V_0[h]e_2^1)].                        (33c)

where the derivative includes the full Piola/Hodge, carrier, slice, and
Riesz connection.  The displayed coefficient is a source-derived smooth-core
coefficient.  Promoting (33a) to a `C1` statement on the full action-flow
graph additionally requires a uniform differentiated Hodge/edge, finite-row,
and left/right Riesz remainder, which reviewed 0074/0078 did not supply.
Equations (33a)--(33c) identify the first permitted power `delta^1`; it is
the actual leading power only if `D_curv!=0`.  They do not prove that the six
source terms fail to cancel.
Neither the source paper nor reviewed 0074/0078 evaluates (33c).
Establishing `D_curv!=0` (with row-corrected `h_(0,1)`) is the precise
missing diagonal construction; absence of that calculation is not a
physical no-go.

### 6.1 Fixed-`k` massive curvature coefficient

P253/0044 equations (13)--(22) were expanded with the toroidal integer held
fixed.  They cannot be inserted verbatim when `n=Theta(delta^(-1))`.
Introduce the two-parameter cell operator `A(delta,k)` by replacing the
integer Fourier row `delta n` in the exact fiber equations by an independent
continuous `k` in a fixed compact massive interval.  Differentiate this
symbolic/analytic extension at fixed `k`; this is not a derivative of one
physical fixed-integer ring family.  For every physical integer sequence,

    A_(delta_N,n_N)
      =A_0(k_N)+delta_N A_1(k_N)+o(delta_N),
    k_N=delta_N n_N,                                    (33c.1)

with the remainder required uniformly for `k_N` in that compact interval.
If `k_N-k_*` is not `o(delta_N)`, the separate term
`A_0(k_N)-A_0(k_*)` is retained and is not folded into `D_curv`.

Hold `k` fixed in this two-parameter sense and re-expand the exact
P253/0074 operators.
Put `c=cos(alpha)`, `s_a=sin(alpha)`, and `h_g=1+delta s c`.  Directly from
the orthogonal tube metric,

    G_0^k f=(partial_s f,s^(-1)partial_alpha f,ik f),
    G_1^k f=-ik s c e_theta f,                            (33d)

    D_0^k V=partial_sV_s+s^(-1)V_s
             +s^(-1)partial_alpha V_alpha+ikV_theta,
    D_1^k V=cV_s-s_aV_alpha-ik s c V_theta.              (33e)

Consequently, if `L_delta=D_delta G_delta`,

    L_1^k=C+2k^2 s c,
    C=c partial_s-s_a s^(-1)partial_alpha.               (33f)

The `2k^2 s c` term is leading first curvature at fixed massive `k`; the
fixed-integer bookkeeping would incorrectly place it later.  With
`S_0^k=(L_0^k)^(-1)` (including the decaying whole-space row),

    P_1^k=-G_1^kS_0^kD_0^k
       +G_0^kS_0^kL_1^kS_0^kD_0^k-G_0^kS_0^kD_1^k.     (33g)

For a field `Y` with axial wave number `k_Y`, define the cylindrical frame
rotation

    C_theta Y=(-cY_theta,s_aY_theta,cY_s-s_aY_alpha).

Regrading the exact covariant derivative gives

    B_0^k(V,Y)=V_s partial_sY+(V_alpha/s)D_alpha Y
                  +i k_Y V_theta Y,
    B_1^k(V,Y)=V_theta C_theta Y
                  -i k_Y s c V_theta Y.                 (33h)

Thus the first full velocity-generator block is

    N_0^k v=B_0^k(U_0,v)+B_0^k(v,U_0),
    N_1^k v=B_1^k(U_0,v)+B_1^k(v,U_0)
       +B_0^k(U_1^odd,v)+B_0^k(v,U_1^odd),               (33i)

    A_1^k=-(P_1^kN_0^k+P_0^kN_1^k).                    (33j)

Here `U_1^odd` is not a fitted profile: it is the odd solution of the
linearized Lane--Emden cell problem in Cao equation (3.29), with the odd
forcing constructed immediately before that equation.  Equations
(33d)--(33j) include pressure, frame rotation, fixed-`k` transport, and the
actual first carrier cell.

For a toroidal-zero, `d=1` physical seed velocity `w=B h`, set

    V_0^k[w]v=-P_0^k{B_0^k(w,v)+B_0^k(v,w)},             (33k)

    V_1^k[w]v=-P_1^k{B_0^k(w,v)+B_0^k(v,w)}
       -P_0^k{B_1^k(w,v)+B_1^k(v,w)
               +B_0^k(w_1,v)+B_0^k(v,w_1)},             (33l)

where `w_1` is the first Piola/Hodge and finite-row correction of the same
physical seed.  Let `v_(a,0)(s)` be the physical-Hessian-normalized `m=0`
column velocity, and let `R_(a,1)^Q` be the full `m=1` reduced resolvent at
`-i nu_a`.  The right-mode first jet is

    v_(a,1)=-R_(a,1)^Q Pi_1 A_1^(k_a)v_(a,0).            (33m)

Let `lambda_(a,0)` denote the Hessian row (15), expressed in the same
physical velocity coordinates.  Differentiating
`lambda_a(A+i nu_a)=0` with `lambda_a(e_a)=1` gives on the `m=1` complement

    lambda_(a,1)|_(m=1)
      =-lambda_(a,0) A_1^(k_a)R_(a,1)^Q.                (33m.1)

The remaining multiple of `lambda_(a,0)` is fixed by normalization and does
not pair with the `m=1` row below.  Angular integration now reduces (33c) to
the explicit one-dimensional coefficient

    D_a[w]=integral_0^infinity s ds {
       lambda_(a,0) dot Pi_0 V_1^(k_a)[w]v_(a,0)
      -lambda_(a,0) dot Pi_0 V_0^(k_a)[w]
          R_(a,1)^Q Pi_1 A_1^(k_a)v_(a,0)
      +lambda_(a,1) dot Pi_1 V_0^(k_a)[w]v_(a,0)},       (33n)

    D_curv[w]=D_1[w]-D_2[w].                             (33o)

The measure and the `lambda` rows in (33n) are the physical ones obtained
from (4)--(9); in fixed fibers they retain the residual `r=Rq` from (5a).
For the rational ray the two copies of (33n) use their actual massive
sequence values `k_(1,N),k_(2,N)`.  Any mismatch from their limiting values
is kept through `A_0(k_(a,N))`.  Variation along the connected carrier path
adds its bordered slice derivative to `U_1^odd,w_1,lambda_(a,1)`, but remains
a thin-ring asymptotic derivative.  By contrast, `nabla_h A` in (18)--(24)
is the instantaneous background derivative along a physical same-leaf seed
at one fixed carrier.  The two derivatives are composed in (33n), never
identified.  A derivative with respect to `kappa` or `R` is not substituted
for a same-orbit seed.

All possible first-order cancellations are visible in (33n): the direct
pressure/frame row, the right Riesz-vector response, and the physical left
Hessian response.  The scalar part of (33f) satisfies

    C[psi(s)exp(i alpha)]
      =(psi'+psi/s)/2+(psi'-psi/s)exp(2i alpha)/2.       (33p)

For `psi=s Omega`, the two forcings are exactly `W/2` and
`s Omega'/2`, reproducing the 0044 geometric check.  Hence the allowed
`m=1 -> m=0` constituent is genuinely nonzero.  But the nodal modes make the
three signed radial terms in (33n) indefinite, and `W/2!=0` alone cannot
exclude their cancellation or prove `D_1!=D_2`.  The route has therefore
advanced from an abstract derivative to a replayable source-mode integral.
At full graph scope it also needs the `C1` remainder named after (33c); at
smooth-source scope its remaining analytic datum is the sign/nonidentity of
(33o), not an unspecified operator theorem.

### 6.2 The `d=1` functional and its finite cokernel

The seed in (33n) need not be guessed.  On a compact regular annulus
`B subset {W'!=0}`, write a divergence-free `d=1` displacement as

    xi_s=i psi(s)e^(i alpha)/s,
    xi_alpha=-psi'(s)e^(i alpha).                        (33q)

Its column coadjoint tangent is

    h_z=xi_s W'=q(s)e^(i alpha),
    psi=-i s q/W'.                                       (33r)

Thus every compact smooth radial coefficient `q` on `B`, after its finite
moment rows are removed, is an actual local DA seed.

The source-bearing principal coefficient can be read directly from (33f):

    Pi_0[2 k^2 s cos(alpha)q(s)e^(-i alpha)]
       =k^2 s q(s).                                      (33s)

After the `m=0` elimination this is the explicit multiplication contribution

    delta Phi_eff,principal=s q.                         (33t)

The full coefficient map for branch `(J,k)` is instead

    Q_(J,k)=s I+R_(J,k),                                 (33t.1)

where direct substitution of (33g)--(33n) shows that `R_(J,k)` contains the
`C,G_1,D_1,B_1` frame/Hodge rows, the Cao odd cell, and the branch-dependent
right and left reduced resolvents.  Those substitutions do **not** yet prove
that every order-zero symbol in `R_(J,k)` cancels, nor that the remaining
family is collectively compact as `J` grows.  In particular, pointwise
compactness for one fixed branch would not control reduced-resolvent norms
when adjacent Sturm gaps shrink.  Equation (33t) is therefore the earned
principal constituent, not an asserted full symbol theorem.

For the normalized generalized Sturm problem

    H_k u_J=lambda_J Phi u_J,
    integral Phi |u_J|^2s ds=1,
    sigma_J=k/sqrt(lambda_J),                            (33u)

Feynman--Hellmann gives

    delta sigma_J
      =(sigma_J/2)integral delta Phi |u_J|^2s ds
       -(sigma_J/(2lambda_J))<u_J,delta H u_J>_(s ds).   (33v)

At an equal-frequency pair, (33s)--(33v) make the direct row density
proportional to

    s^2{|u_(J_1)(k_1)|^2-|u_(J_2)(k_2)|^2}.             (33w)

For a fixed pair, Sturm nodal count proves that this density is not zero on
every interval: distinct Sturm modes cannot share all zeros, so at a zero of
one which is not a zero of the other the difference has a strict sign.  But
that fact alone does not exclude cancellation by the two distinct full
remainders `R_(J_1,k_1)^*` and `R_(J_2,k_2)^*`, or membership in the finite
row annihilator.  This is the exact fixed-pair cokernel condition:

    P_row^*{(nu_*/2)s^2(|u_1|^2-|u_2|^2)
      +R_(J_1,k_1)^*d_1-R_(J_2,k_2)^*d_2}=0.            (33x)

The curvature nonidentity route is therefore reduced, but not completed.
One can finish it either by deriving the full order-zero principal symbol and
a collective high-`J` resolvent-gap estimate, or by contradicting (33x)
directly for one fixed pair with its adjoint Sturm equations.  No shrinking-
gap uniform compactness is assumed.  The independent near-axis
off-diagonal route below does not use (33x) and supplies the positive response
result of this attempt.  The full action-flow graph-`C1` remainder after
(33c) remains separately open for any actual-Cao curvature transfer.

## 7. Off-diagonal response and the remaining noncancellation

For the rational-ray crossing, `ell=n_1-n_2!=0`.  On a positive regular cell
the orbit map has the exact principal row

    [xi,zeta(I)partial_theta]
      =xi^I zeta'(I)partial_theta-i ell zeta(I)xi.        (34)

Thus, modulo the genuine stabilizer, the nonzero-`ell` compact DA map is
order zero and has a bounded local inverse; the off-diagonal test space is
nonempty.  Toroidal integration in (25) is nonzero only for
`ell=n_1-n_2`; the real conjugate supplies `-ell`.  A partner-frequency
matrix element instead requires `ell=+/-(n_1+n_2)` and belongs to the
squeezing block, not (19).

The Hessian shortcut reduces the full nonvanishing problem to

    gamma_12=sup_(||h||_Y=1)|M(h)_12|.                  (35)

The following calculation evaluates it in the fixed-column limit without
choosing an arbitrary `C_0^*` potential.  First construct the representation
which is actually needed.  Use velocity variables for the two `m=0` modes.
The column base is

    U=s Omega(s)e_alpha,       bar_omega=W(s)e_z,
    W=2 Omega+s Omega'.                                  (35a)

Fix a positive-core disk `B_0={s<s_1}` on which `W>=w_0>0` and a nonzero
axial character `k`.  If `q=curl v=[xi,W e_z]`, the local orbit inverse is

    xi_s=i q_s/(kW),       xi_alpha=i q_alpha/(kW),
    xi_z={xi_sW'-q_z}/(ikW).                            (35a.1)

Denote this regular-axis operator by `T_k v=C_(0,k)^(-1)curl v`.  It is a
first-order differential operator with smooth coefficients on `B_0`, acting
on regular cylindrical fields; it makes no division by a vanishing edge
profile.

Put `p=conjugate(e_1)=[eta_1,bar_omega]` and retain the full energy--impulse
one-form `mu_F` from (7).  The two orbit terms in (8) obey

    integral mu_F dot [xi,p]
      =integral xi dot f_p,
    f_p=(grad p)^Tmu_F+(p dot grad)mu_F,

    integral mu_F dot [eta_1,q]
      =integral q dot g_p,
    g_p=-(eta_1 dot grad)mu_F-(grad eta_1)^Tmu_F.        (35a.2)

For every smooth divergence-free velocity `v` compactly supported in `B_0`,
substitution of (35a.1), `q=curl v`, and one integration by parts gives the
actual Hessian row

    ell_1(curl v)=integral lambda_1 dot v dx,

    lambda_1=P_div{rho_0 Bp
       +(1/2)T_k^*f_p+(1/2)curl g_p}.                    (35a.3)

All three terms are locally `L2` (indeed smooth) by the fixed-mode `H4/H5`
regularity.  Equation (35a.3) contains the kinetic and both second
push-forward pieces of (8), so it represents the actual orbit-Hessian row on
compact velocity tests in `B_0`.  It is not yet the whole-space velocity
representer used in (35c).  The symbol `P_div` in (35a.3) means the local
solenoidal class modulo gradients: applying the whole-space Leray projector
would first require an extension of the expression in braces outside `B_0`,
and different extensions need not restrict to the same field.

The ambiguity is exposed exactly.  For the `-k_1` character,

    phi_h(s,z)=I_0(|k_1|s)exp(-i k_1 z),
    Delta phi_h=0.                                       (35a.4)

Adding `c grad phi_h` changes no compact divergence-free pairing in `B_0`,
but changes the radial Taylor coefficient by `c k_1^2/2`.  Therefore the
local formula cannot by itself fix the coefficient `A_1^sharp` below.

There is also a graph-domain distinction.  If `q=curl v` with `v` compact in
`B_0`, then the vorticity `Aq` is compact, but its Biot--Savart velocity is a
decaying whole-space field, not a compact velocity.  The Hamiltonian identity

    ell_1(Aq)=-i nu_* ell_1(q)                            (35a.5)

implies the local adjoint equation only after the intrinsic expression for
`ell_1(Aq)` is integrated through the edge and exterior with no boundary
row, or equivalently after a cutoff-independent whole-space representer has
been constructed.  Compactness of `Aq` does not supply that step.

If that global representer theorem holds, the left field obeys

    A^*lambda_1
      =P_div{(U dot grad)lambda_1-(grad U)^Tlambda_1}
      =-i nu_* lambda_1,                                (35a.6)

and for the `-k_1` character its component equations are

    -i nu_* lambda_r=-W lambda_alpha-partial_s pi,
    -i nu_* lambda_alpha=2 Omega lambda_r,
    -i nu_* lambda_z=+i k_1 pi,
    partial_s^*lambda_r-i k_1lambda_z=0.                (35a.7)

Eliminating `lambda_alpha,lambda_z,pi` gives

    -partial_s partial_s^*lambda_r
      +k_1^2(1-Phi/nu_*^2)lambda_r=0.                   (35a.8)

This conditional singular-axis statement has an explicit recurrence.  Write

    1-Phi(s)/nu_*^2=sum_(j>=0) psi_j s^(2j),
    lambda_r=sum_(j>=0) a_j s^(2j+1).

Then

    a_(j+1)={k_1^2/[4(j+1)(j+2)]}
             sum_(l=0)^j psi_l a_(j-l).                 (35a.9)

The second indicial branch is `s^(-1)`.  Hence regularity excludes that
branch, and `a_0=0` kills every regular coefficient.  This proves that a
nonzero global left eigenfield which reaches the axis has
`A_1^sharp=a_0!=0`; it does not prove that the intrinsic Hessian row has such
a global field rather than an edge/dual component.

The physical representative is nevertheless fixed without choosing a local
`P_div` extension.  Let `q=[xi,W e_z]` have character `+k_1`, and let
`p=conjugate(e_1)=[eta_1,W e_z]` have character `-k_1`.  In the column core,

    q_s=-i k_1W xi_s,       q_alpha=-i k_1W xi_alpha.

Substitution in the KKS formula and then (15) gives the exact cancellation

    ell_1(q)=-i nu_* Omega(q,p)
      =(rho_0 nu_*/k_1)integral
         (q_s eta_(1alpha)-q_alpha eta_(1s))dx
      =integral a_1 dot q dx,                            (35a.10)

    a_1=(rho_0 nu_*/k_1)
          (eta_(1alpha)e_s-eta_(1s)e_alpha)              (35a.11)

inside the compact core.  Thus the apparently singular factors `1/W` cancel
before any edge limit is taken.

This row has the required global domain.  For the conjugate eigenmode, whose
velocity is denoted `v_1`, the transverse vorticity equations give

    eta_(1s)=i v_(1s)/nu_*,
    eta_(1alpha)=i v_(1alpha)/nu_*
                   +(s Omega'/nu_*^2)v_(1s).             (35a.12)

The fixed-mode regularity and the `p>=6` regular core boundary therefore give
an `H4` transverse coefficient up to the edge.  A bounded Sobolev extension
of (35a.11) through a fixed exterior collar, followed by a compact cutoff,
produces `a_1^E in H1(R3)`.  Put

    lambda_1=curl a_1^E in L2_sigma(R3).                 (35a.13)

It is decaying and divergence free, and integration by parts gives
`ell_1(curl v)=integral lambda_1 dot v` for every physical velocity in the
fixed fiber.  Two exterior extensions differ only where the physical
vorticity is zero.  More precisely, if `delta a=0` on the closed vorticity
core and `v=Bq`, then the decaying whole-space/interface rows give

    integral curl(delta a) dot v dx
      =integral delta a dot curl(v) dx
      =integral delta a dot q dx=0.                      (35a.13a)

Thus the represented functional is extension independent, although this
first construction does not select a unique field in the larger ambient
`L2_sigma` space.  This resolves the functional ambiguity in (35a.4); the
dual-Riesz construction below selects the physical adjoint line.

The dual-Riesz route gives an independent uniqueness check.  On the exact
fixed-fiber space `X_(0,k_1)^3`, reviewed 0074/0078 algebraic simplicity makes
the dual Riesz range one dimensional.  Equation (35a.10) is continuous on
that space, is normalized by `ell_1(e_1)=1`, and (12) gives

    ell_1(Aq)=-i nu_*ell_1(q).                            (35a.14)

Hence it is the normalized element of the rank-one dual Riesz range.  No
continuity of a raw `C_k^(-1)` is assumed: continuity follows from the
explicit `L2` density (35a.11).  Rank-one dual Riesz makes this functional
unique; it does not make an ambient `L2` representative unique.

The cutoff-independent decaying adjoint representative can also be built
directly.  Solve on `(0,infinity)`

    [-partial_s partial_s^*+k_1^2]lambda_s
       =(k_1^2/nu_*^2)Phi lambda_s,                      (35a.14a)

with the regular `s` branch at zero, continuity of `lambda_s` and
`partial_s^*lambda_s` at the compact coefficient edge, and the decaying
exterior `K_1(|k_1|s)` branch.  This is the same scalar Sturm pencil that
defines the simple right mode, so its simple eigenvalue supplies a unique
nonzero solution up to scale.  Reconstruct

    lambda_alpha=2i Omega lambda_s/nu_*,
    lambda_z=-k_1 pi/nu_*,                               (35a.14b)

with `pi` from the radial row and divergence.  The regular-axis, interface,
and exterior Bessel rows put this field in the whole-space velocity `L2`
domain.  Its radial component is fixed by matching the explicit core
representative: `lambda_s=rho_0 v_(1s)` in the core, and the Sturm interface
rows continue that value uniquely to the exterior.  Direct integration of
the velocity equations, with the interface
terms cancelling by the two matching rows, proves that

    L_St(q)=integral lambda_St dot Bq dx                 (35a.14c)

is a continuous nonzero left eigenfunctional on `X_col^3`.  Here
`B A_q=A_v B` on the graph, and the global pressure drops against the
divergence-free decaying `lambda_St`.  Both `ell_1` and `L_St` therefore lie
in the one-dimensional dual Riesz range supplied by algebraic simplicity, so
they are proportional.  Start the Sturm solution with the positive-core data
from `curl(a_1^E)`.  Because its explicit radial component
`rho_0 v_(1s)` is nonzero, a compact solenoidal positive-core test `v` can be
chosen with nonzero pairing; `q=curl v` is DA by (35a.1).  The two functionals
agree on that test, fixing the proportionality.  Hence

    ell_1=L_St,             L_St(e_1)=1.                (35a.14d)

If the Sturm field were instead started with an arbitrary nonzero scale, its
pairing with `e_1` could not vanish: orthogonality of left and right vectors
in an algebraically simple Riesz block would create a generalized vector.
Thus (35a.14d) identifies the intrinsic orbit-Hessian row with the global
adjoint Sturm field and selects, from the extension class (35a.13), the
cutoff-independent representative needed below.  No density of compact
positive-core tests in the full edge-degenerate DA closure is asserted.

The graph test is now legitimate.  Although `Aq` has a noncompact
Biot--Savart velocity, (35a.13) is global and
`curl(A_vv)=A(curl v)`.  Therefore

    integral lambda_1 dot A_vv
      =ell_1(A curl v)=-i nu_*ell_1(curl v)

for compact regular-axis tests.  The pressure pairs to zero with
`lambda_1`, and (35a.6)--(35a.9) follow distributionally on the positive-core
component.  More directly, (35a.12)--(35a.13) give one representative there

    curl(a_1^E)_s=rho_0 v_(1s),
    curl(a_1^E)_alpha=rho_0 v_(1alpha)
       -i rho_0(s Omega'/nu_*)v_(1s).                    (35a.15)

The conjugate right-mode rows
`v_(1alpha)=iWv_(1s)/nu_*` and
`partial_s^*v_(1s)-ik_1v_(1z)=0`, together with
`W-s Omega'=2 Omega`, turn (35a.15) into

    curl(a_1^E)_alpha=2i rho_0 Omega v_(1s)/nu_*,
    curl(a_1^E)_z=rho_0 v_(1z),                          (35a.15a)

which verifies every core component of the adjoint reconstruction
(35a.14b), not only its radial row.

The two representatives give the same functional by (35a.13a) and
(35a.14d).  On the positive core, every compact divergence-free test velocity
`v` has `q=curl v` in the actual DA tangent space by (35a.1), and its
whole-space Hodge velocity is exactly `Bq=v`.  Hence

    integral [lambda_St-curl(a_1^E)] dot v dx=0

for every such test.  The difference is therefore a local gradient; because
both representatives are divergence free, it is a harmonic gradient.  For
the `m=0,-k_1` character its alpha component is zero.  Both fields satisfy the
local adjoint system, whose alpha row gives
`delta lambda_alpha=2i Omega delta lambda_s/nu_*`; `Omega>0` forces the radial
difference to vanish, and `k_1!=0` plus divergence forces the axial difference
to vanish.  Consequently `lambda_St=curl(a_1^E)` pointwise on the positive
core, although their ambient exterior representatives are not identified.
Their difference also annihilates every physical response velocity, including
`B(delta A[w]e_2)`.

Since `s Omega'=O(s^2)` and the regular right Sturm mode has
`v_(1s)=A_1s+O(s^3)` with `A_1!=0`,

    A_1^sharp=rho_0 A_1!=0.                              (35a.16)

Equation (35a.15) fixes the core representative and its interface data, while
(35a.14a)--(35a.14d) gives its unique decaying adjoint continuation.  This
identifies the nonzero regular-axis radial datum directly; it is not inferred
from positive Krein sign.  The recurrence (35a.9) is the independent
singular-axis check.

Let `v_2` be the right velocity at axial wave number `k_2`.  If `h=curl w`
and the divergence-free seed
velocity `w=Bh` has wave number `k_h=k_1-k_2`, the full fixed-slice derivative
is

    delta A[w]v_2
      =-P{(w dot grad)v_2+(v_2 dot grad)w}.              (35b)

Because `lambda_1` is divergence free and has the regular-axis and decaying
whole-space rows, self-adjointness of the physical Leray projector and one
integration by parts remove the pressure exactly and give

    M_12(w)=integral w dot F_12^v dx,
    F_12^v=(v_2 dot grad)lambda_1-(grad v_2)^Tlambda_1.  (35c)

Both Euler convection terms and the full orbit-Hessian left row are already
inside (35c); there is no omitted Hodge term with which its leading
coefficient could later cancel.

Write regular-axis expansions

    v_(2r)=A_2s+O(s^3),       v_(2alpha)=B_2s+O(s^3),
    v_(2z)=C_2+O(s^2),
    lambda_(1r)=A_1^sharp s+O(s^3),
    lambda_(1alpha)=B_1^sharp s+O(s^3),
    lambda_(1z)=C_1^sharp+O(s^2).                       (35d)

The divergence and right/left eigen-equations at the common nonzero physical
frequency `nu_*` give

    C_2=2i A_2/k_2,
    B_2=-i W(0)A_2/nu_*,
    B_1^sharp=2i Omega(0)A_1^sharp/nu_*.                (35e)

Here `W(0)=2 Omega(0)`.  Directly from (35b), the velocity adjoint of the
linearized Euler operator is

    A^*lambda=P{(U dot grad)lambda-(grad U)^Tlambda}.    (35e.1)

Its `m=0` alpha row gives the second identity in (35e), while divergence gives
`C_1^sharp=-2iA_1^sharp/k_1`.  The dual-Riesz/Sturm identification
(35a.14d), the recurrence (35a.9), and continuation to the axis give
`A_1^sharp!=0`.  The same
regular right-mode recurrence gives `A_2!=0`; no
pointwise proportionality between the nonlocal Hessian row and the right
mode is assumed.

The exact orthonormal-cylindrical alpha component of (35c) is

    F_(12,alpha)^v
      =v_(2r)(lambda_(1alpha)'-lambda_(1alpha)/s)
       -i k_1 v_(2z)lambda_(1alpha)
       +2 v_(2alpha)lambda_(1r)/s.                      (35f)

Substitution of (35d)--(35e) yields the full leading coefficient

    F_(12,alpha)^v
      ={4i Omega(0)/nu_*}
        A_2 A_1^sharp(k_1/k_2-1)s+O(s^3).               (35g)

It is nonzero because the rational-ray pair has `k_1!=k_2`.  The first term
in (35f) is `O(s^3)` by regular-axis compatibility; the remaining two terms
are exactly the two possible order-`s` convection contributions and have
already combined into (35g).  No pressure, Hodge, or Hessian row remains at
that order to cancel the displayed coefficient.

Choose a sufficiently small punctured annulus
`B={s_0<s<2s_0}`.  Radial monotonicity of the Lane--Emden/Cao core gives
`W'(s)!=0` there, equivalently `zeta'(I)!=0` in its regular action coordinate,
and (35g) has fixed nonzero phase on `B`.  Let

    w=w_alpha(s)e^(i k_h z)e_alpha,
    w_alpha in C_c^infinity(B).                          (35h)

This field is divergence free, and `h=curl w` is compact with

    h_s=-i k_h w_alpha,
    h_z=s^(-1)(s w_alpha)'.                              (35i)

The exact coadjoint inverse on `B` is

    xi_s=w_alpha/W,
    xi_alpha=0,
    xi_z={xi_s W'-h_z}/(i k_h W),                       (35j)

and `div xi=0` follows from `div h=0` and
`div[xi,W e_z]=-i k_h W div xi`.  Hence
`h=[xi,bar_omega]` is an actual compact smooth DA tangent.  Taking the phase
and sign of `w_alpha` from (35g) gives `M_12(h)!=0`.

For the rational-ray sequence, choose the dense rational value
`r=k_2/k_1=Q/P` outside the finite exceptional set `{1/2,1,2}`.  Then
`ell=n_1-n_2` is eventually of absolute value at least two and is distinct
from `+/-n_1,+/-n_2`.  Circulation, impulse, recentering and translation rows
have toroidal characters zero or one, while the modal gauge rows have
characters `+/-n_1,+/-n_2`.  Taking their correction profiles in the same
characters, Fourier orthogonality makes every declared finite row vanish on
(35j).  Thus `P_row xi=xi`, `h=C_0xi` is already an admissible DA tangent,
and `w=Bh` is its velocity representative.  The correctly typed identity
`M_12(h)=integral w dot F_12^v` is nonzero.  Hence (35g), evaluated on the
punctured annulus, is a nonzero obstruction for the
**final constrained functional**, not a raw convection component.  In
particular

    gamma_12^col>0.                                     (35k)

The fixed-`k` massive operators and the seed (35h) are analytic on simple-mode
cells, so this nonzero response persists on a neighborhood of the selectable
column crossing.  The local construction alone would not have proved this;
the KKS cancellation, dual-Riesz uniqueness, and global Sturm representative
in (35a.10)--(35a.16) are the required bridge.

### 7.1 Source-specific fixed-`k` transfer

For the one compact seed above, the common-domain remainder can be derived
without claiming an operator-norm derivative for every control.  Let
`U_N` be the reviewed 0074/0078 action-flow/Piola intertwiner at an exact Cao
crossing on the independently reviewed 0083/0089 path.  In the fixed massive fibers,

    ||U_N^(-1)A_(a,N)U_N-A_(a,col)||_(D_col->X)->0,
    ||U_N^(-1)P_(a,N)U_N-P_(a,col)||_(X->D_col)->0.      (35l)

The second convergence uses the fixed simple contours.  Choose the incoming
KKS frame by these Riesz maps.  Formula (8), pulled back by the same unitary
half-density/Piola chart, has only two terms: the whole-space Hodge kinetic
pair and the compact-core second-orbit pair.  The reviewed Hodge convergence,
carrier `C4` convergence on the support, and exterior energy-tail estimate
give, on the two finite Riesz ranges,

    E2_N(U_N q,U_N p)=E2_col(q,p)+o(1)||q||_D||p||_D.   (35m)

No kinetic-energy identification is made: the second term of (8) is included
in (35m).  Positivity and finite rank make the Gram matrices uniformly
invertible, so physical Hessian normalization commutes with the limit and
the left rows `lambda_(a,N)` converge to `lambda_(a,col)` in the dual graph
norm.

Now push the compact displacement (35j), not the velocity, through the
volume action-flow chart and denote it by `xi_N`.  Then

    h_N=[xi_N,bar_omega_N]                               (35n)

is exactly DA, has the required integer difference character, and remains
inside the corresponding punctured regular torus.  The character exclusions
above make every speed, circulation, impulse, center, and modal phase row
zero, so no moving-slice or `dc_g` term is hidden.  The reviewed whole-space
Hodge intertwiner gives

    U_N^(-1)B_N h_N -> w                                (35o)

in local `H5`, with the collar/exterior energy tail tending to zero.

The bilinear Euler map obeys on the smooth graph core

    ||[v,q]||_(H3)<=C||v||_(H5_loc)||q||_(H4),           (35p)

and the full pulled-back covariant derivatives contain only the already
retained fixed masses `k_1,k_2,k_h` plus `o(1)` metric/frame coefficients.
Equations (35l),(35o)--(35p) therefore give the source-specific remainder

    U_N^(-1)V_N[h_N]e_(2,N)
       =V_col[h_col]e_(2,col)+o_(H3)(1).                 (35q)

Pairing (35q) with the convergent physical Hessian row in (35m) yields

    M_(12,N)(h_N)=M_12^col(h_col)+o(1).                 (35r)

The right side is nonzero by (35g)--(35k), so every sufficiently thin exact
crossing on that connected path has an actual compact DA seed with

    M_(12,N)(h_N)!=0.                                   (35s)

Thus the off-diagonal nonvanishing transfers at source-specific fixed-`k`
scope.  The exact-crossing supplier is the reviewed 0083/0089 path.  The
transfer still does not provide a two-sided
bound after normalizing `h_N` in the physical `Y^4` norm: the residual
`r=R_Nq`, Piola factors, and possible `N` dependence of that norm remain the
separate high-index ledger in Section 8.

There is a second, non-curvature control continuation.  Two real seeds in
the same difference character can have

    M(h^(r))=[[0,c_r],[-conjugate(c_r),0]],   r=1,2,     (36)

and direct multiplication gives

    [M(h^(1)),M(h^(2))]
      =2i Im(c_1 conjugate(c_2)) diag(1,-1).             (37)

Thus `Im(c_1 conjugate(c_2))!=0` supplies two complex-linear,
action-preserving noncommuting axes without the `delta`-suppressed diagonal
coefficient.  This is not earned by calling the cosine and sine phases of
one static seed two controls.  It requires two independently switchable,
phase-locked same-carrier Euler histories (or a retained material
orientation which physically breaks the axial phase orbit), and the two
response functionals must remain independent after all slice rows.  Equation
(37) is an exact algebraic alternative for R2; at R1 it is only a candidate
because neither autonomous history has yet been constructed.

The first missing autonomous equation is explicit.  If `f_r`, `r=1,2`, are
retained physical Euler sectors in the difference harmonic with amplitudes
`a_r`, phase-locked mixing requires

    dot a_r=-i(nu_1-nu_2)a_r+N_r(a_r,z,q_Q),
    M_r(t)_(12)=c_r a_r(t),                              (37a)

with `Im(c_1 conjugate(c_2))!=0`, while the carrier and complement equations
keep both sectors coherent through the gate.  A single simple retained mode
has only one axial phase orbit.  Two chosen initial phases of it, a rigid
rotation, or free waiting do not make two switches.  The first construction
needed by this route is therefore two distinct invariant retained sectors,
or one material orientation with two internally generated switching states,
which solves (37a) on the same nonlinear interval.  This precisely separates
the established compression algebra from the still-open autonomous-history
mechanism.

## 8. High-harmonic and gate boundary

For the rational ray `ell_N=Theta(N)` and `R_N=Theta(N)a`, so the physical
wave number `ell_N/R_N` is `Theta(1/a)`.  The local estimate (29) has no
coordinate `N^5` loss.  However the absolute normalization (9) alone does
not prove uniform upper and lower bounds for the raw mode Hessians `H_a(N)`
or the dual source coefficient `gamma_12(N)`: in particular the residual
factor `r=R_Nq` in (5a) must be carried through the raw displacement
amplitudes before any cancellation or growth in `R_N` can be asserted.
Hence no two-sided scale for
`g_(12,N)^4` is claimed.

For the curvature-diagonal route, (33c) gives, if its normalized coefficient
has a uniform nonzero limit,

    g_(3,N)^4=delta_N |D_curv,N|+O(delta_N^2)
              =Theta(N^(-1)).                           (38)

This is an additional `N^(-1)` response loss beyond whatever physical
normalization appears in the off-diagonal row.  A diagonal rotation of fixed
angle therefore needs `A_N T_gate=Omega(N)` even before complement or
nonlinear errors.  If `D_curv,N` itself decays, that decay multiplies (38).
No finite-`delta` nonzero value is substituted for this uniform ledger.

The three frozen verdicts remain separate.  The Hessian-row identity and
matrix-skew oracle are exact, but R1 needs `gamma_12>0` and `D_curv!=0` with
uniform physical normalization.  R2 still needs two autonomous histories.
R3 still needs both outbound and return kernels, squeezing, and a nonlinear
remainder through the same gate interval.  No instantaneous identity here
licenses those later conclusions.

## 9. Route verdicts and strongest result

- **Route A (global local-potential inversion): blocked as a construction,**
  because a globally regular representative of (17) on the constrained
  whole-space graph has not been built.  The exact density (25)--(27) and
  patch theorem (29) remain valid whenever that representative exists.
- **Route B (Hessian/Lagrange representation): established.**  Equations
  (14)--(24) construct the physical left covector with the correct clock,
  density, Fourier and orbit-Hessian normalization, prove the covariant
  Feynman--Hellmann row, and give the basis-invariant skew-Hermitian
  compression criterion.
- **Frozen diagonal `(0,0)` subroute: refuted** by the stabilizer identity
  (31)--(32).  The generated `(0,+/-1)` curvature route remains active at
  coefficient (33n); `delta^1` is the first permitted loss and becomes the
  leading loss only when `D_curv!=0` and the full graph `C1` remainder closes.
- **Off-diagonal nonzero-response subroute: established at fixed-column and
  source-specific transfer scope.**  Equations (35a.1)--(35a.4) expose the
  local projection ambiguity; (35a.10)--(35a.16) remove it through the exact
  KKS cancellation, rank-one dual Riesz identification, global adjoint Sturm
  row, and singular-axis recurrence.  Equations (35b)--(35k) then give a
  nonzero constrained DA seed.  The reviewed 0083/0089 path then transfers
  it to sufficiently thin exact same-family Cao crossings.  No uniform
  normalized high-index response is claimed.
- **Route C (uniform high-index useful response): blocked** until (33), (35),
  and the physical mode/dual bounds are uniform in the rational-ray sequence.
- **Two independent off-diagonal histories: blocked at physical
  realization,** while the exact noncommutation criterion (37) is
  established.  It avoids the curvature loss only if two histories, not two
  hand-selected phases, are constructed.

The strongest exact verdict is therefore a useful physical bridge: the
positive-mode KKS covector is the actual relative-Hessian row, its full
degenerate compression is skew-Hermitian, and noncommuting analyzer axes
reduce invariantly to two source-specific Hessian derivatives.  The exact KKS
cancellation and rank-one dual Riesz/Sturm identification supply the global
physical left row; its near-axis coefficient and an actual constrained DA
seed prove fixed-column nonvanishing.  One frozen
candidate diagonal slice is exactly a stabilizer, exposing the materially
different curvature-mixed seed required next.  R1, P2, P4, particle,
quantum, and relativistic conclusions remain open.

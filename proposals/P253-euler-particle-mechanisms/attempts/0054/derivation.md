# Actual Cao moving-boundary Green/Leray graph jet

## 1. Conventions and theorem statement

Write `delta=a_delta/R_delta`, `L=log(1/delta)`, fix a physical toroidal
Fourier number `ell>=2`, and retain `k_delta=ell delta` in the straight-column
operator. Vorticity is a two-form (or its Euclidean divergence-free vector
density), and the exact linearization is

    partial_t eta=A_delta eta,
    A_delta eta=curl(u_delta cross eta
                     +B eta cross Omega_delta),           (1)
    B=curl(-Delta_R3)^-1.                                 (2)

The sign in (1) follows by differentiating
`partial_t Omega=curl(u cross Omega)`. The decaying whole-space operator (2),
not an interior inverse, defines velocity and pressure closure.

For an integer Cao exponent `p>=6`, after centering, fixing circulation and
factoring the exact physical `R_delta,a_delta`, the following objects are
constructed below:

1. a polyhomogeneous `C^2` carrier/free-boundary chart with independent
   variables `(delta,tau)`, restricted at the end to
   `tau=delta^2 log(delta)`;
2. a fixed Hodge representation of the complete DA quotient, including its
   harmonic circulation row;
3. the twice-differentiated three-dimensional whole-space Hodge and metric
   Leray operators on one graph domain; and
4. the expansion

       Ahat_delta=A_col(k_delta)+delta C_1
          +delta^2 L C_(2,log)+delta^2 C_2+R_delta,        (3)

   with

       ||R_delta^loc||_(D_* -> X_*)<=C delta^3 L,
       ||R_delta^H||_(Psi^(-1))<=C delta^3 L.             (4)

The order-minus-one mode estimates and their order-minus-two Schur products
then imply the contour-relative estimate required by 0052. This is an HJ2
theorem for the frozen `p>=6` subfamily. It is not a nonlinear rotating-wave
existence theorem.

## 2. The exact carrier jet rather than differentiability by assertion

### 2.1 Fixed-core equation

Let `x_epsilon=(x_epsilon,1,0)` and `s_epsilon` be Cao's center and core scale,
and put

    q=s_epsilon/x_epsilon,1.                              (5)

Equation (3.12) of the source becomes, exactly on the rescaled meridional
half-plane,

    -Delta w+q/(1+q y_1) partial_1 w
       =(1+q y_1)^2 (w_+)^p.                             (6)

Indeed (6) follows by multiplying
`-(1+qy_1)^-1 div((1+qy_1)^-1 grad w)=(w_+)^p` by
`(1+qy_1)^2`. Expanding (6) gives the already source-fixed cells

    L_U V=2 y_1 U_+^p-partial_1 U,                       (7)

    L_U Z={p(p-1)/2}U_+^(p-2)V^2
       +2p y_1 U_+^(p-1)V+y_1^2 U_+^p
       -partial_1 V+y_1 partial_1 U,                     (8)

where `L_U=-Delta-pU_+^(p-1)`.

Equation (6) alone does not impose the axis boundary or the far logarithmic
coefficient. To retain them, rescale the exact source kernel (2.2). For
physical points

    r=R(1+q y_1),  z=Rq y_2,
    r'=R(1+q y'_1), z'=Rq y'_2,                          (9)

its small parameter is exactly

    rho=q^2 |y-y'|^2/[(1+q y_1)(1+q y'_1)].              (10)

The elliptic-integral expansion evaluated in 0052 therefore makes the fixed
core integral operator a smooth function of the two independent variables

    (q,tau),       tau=q^2 log q,                         (11)

after its circulation-dependent rank-one logarithm is put into the Bernoulli
constant. Its remainder through the displayed second cells is
`O(q^3 |log q|)` in `C^(3,alpha)` from compactly supported
`C^(1,alpha)` data. The `rho^2 log rho` kernel remainder is actually
`O(q^4|log q|)`; the larger order in this statement comes from the geometric
prefactors and the first cell.

### 2.2 Augmented implicit map

Work on `B_2` in the affine logarithmic-potential space

    W={w in C^(3,alpha): w-c_log log|y| is harmonic off B_2},   (12)

with fixed `c_log` (circulation), even symmetry in `y_2`, and center gauge
`partial_1 w(0)=0`. The exact core integral equation, circulation condition,
center condition, and the two physical scale/speed equations obtained from
source system (3.36) define

    F(q,tau;w,mu,c_log)=0.                                (13)

Because `p>=6`, `w -> (w_+)^p` is `C^4` between the spaces needed for two
parameter derivatives and one graph derivative. The kernel calculation
following (10) makes `F` `C^3` in `(q,tau,w,mu,c_log)` near
`(0,0,U,mu_0,c_0)`.

The profile part of `D F` is `L_U`. Source Lemma 3.8 says its bounded kernel
is exactly

    ker L_U=span{partial_1 U,partial_2 U}.                (14)

Evenness removes `partial_2 U`, and the center gauge removes
`partial_1 U`. The far-log/circulation row is nonzero because

    -2 pi U'(1)=integral_(B_1) U^p dy>0,                 (15)

while the boundary/scale row is nonzero by the Hopf sign `U'(1)<0`. The
derivative of the first source (3.36) row with respect to the centered major
radius has leading coefficient `W log(1/epsilon)`, while the second row has a
nonzero scale derivative by (15); after using `q=s/R`, these rows are uniformly
triangular in the factored physical variables. Consequently the augmented
derivative is triangular modulo a compact
potential perturbation, has index zero, and has trivial kernel by (14)--(15);
it is an isomorphism. This is the missing surjectivity step that uniqueness
alone would not supply.

The Banach implicit-function theorem applied with `(q,tau)` independent gives
a unique `C^3` solution

    w(q,tau)=U+qV+tau Z_log+q^2 Z+Rem_w,                 (16)
    ||Rem_w||_(C^(3,alpha))
       <=C(|q|^3+|q tau|+tau^2).                         (17)

Restricting `tau=q^2 log q` turns (17) into

    ||Rem_w||_(C^(3,alpha))<=C q^3 |log q|.              (18)

The source's centered estimate in Lemma 3.9 supplies uniformity of the
inverse as `q->0`; (14)--(15) identify its limit. The constructed solution has
the concentration, circulation, symmetry, and speed class of Theorem 1.6, so
the source uniqueness theorem identifies it, up to the removed axial
translation, with the actual Cao family. Thus (16)--(18) are a newly derived
carrier jet, not a differentiability claim imported from the source's
`O(epsilon^2|log epsilon|)` estimate.

The physical thinness `delta=a_delta/R_delta` and `q` differ by the smooth
free-boundary radius gauge `delta=q(1+O(q))`; its inverse has the same
`(q,q^2 log q)` regularity. Hence (16)--(18) hold with `delta` after a change
of the finite coefficients. Slowly varying exact major-radius and clock
factors are retained in the physical normalization rather than Taylor-expanded
as if they were constants.

### 2.3 Interface regularity and absence of a generated sheet

In the fixed Hanzawa coordinate let `d=1-s` be inward normal distance. Hopf's
lemma and (16) give

    Omega_delta=d_+^p a_delta(s,beta) tau_delta,
    0<c<=a_delta<=C,                                     (19)

with `a_delta` uniformly `C^(3,alpha)`. Parameter differentiation yields

    partial_delta Omega_delta
       =p d_+^(p-1)(partial_delta d)a_delta
          +d_+^p partial_delta a_delta,                  (20)

and a second derivative whose lowest power is `d_+^(p-2)`. For `p>=6`, (20)
and its second derivative have all traces needed by the graph calculation;
distributional differentiation creates no boundary delta measure. For a DA
tangent

    eta=-L_xi Omega_delta=d(i_xi Omega_delta),            (21)

the lowest boundary power is `d_+^(p-1)`, again with zero vortex-sheet trace.
This closes the physical core/interface support statement before applying the
global Hodge operator.

## 3. A genuine common DA quotient

Let `mu_delta` and `varpi_delta=i_(Omega_delta)mu_delta` denote the pulled-back
volume and vorticity two-forms. A no-normal-flow displacement satisfies

    d(i_xi mu_delta)=0.                                  (22)

On a solid torus, the relative Hodge decomposition of the closed two-form in
(22) is

    i_xi mu_delta=d alpha+h xi_harm,                     (23)

where the single harmonic row records meridional flux/circulation. Retaining
`h` is essential; (23) is not a local Clebsch completion.

In the action--angle chart,

    varpi_delta=f_delta(I)dI wedge d beta,                (24)

and `f_delta/f_0` extends as a positive `C^2` function to the boundary after
the common factor `d^p` in (19) is removed. Equations (21)--(24) identify the
stabilizer as the closed-one-form kernel of `d(i_xi varpi_delta)`. Fix its
orthogonal complement together with the impulse and centering rows.

The one-dimensional Hardy estimate

    ||d^(p-1)g||_(H^-1(0,1))
       <=C ||partial_d(d^p g)||_(H^-1(0,1))               (25)

and its tangential Fourier version show that multiplication by
`f_delta/f_0`, the Hodge projections in (23), and their inverses are uniformly
bounded on the quotient completion. The finite harmonic row is controlled
directly by circulation. Therefore

    T_delta:[Y_*/Stab_*] -> [Y_delta/Stab_delta]          (26)

is an isomorphism with uniform norms, and so is the induced map of closed DA
energy ranges. This constructs the fixed `X_*`; it does not assume all
divergence-free vorticities are dynamically accessible.

Cross-`delta` maps in (26) are operator identifications. Only the physical
same-carrier displacement generating (21) is volume preserving and belongs to
one coadjoint orbit.

## 4. Exact moving-coordinate Biot--Savart derivatives

### 4.1 General shape formula and the singularity cancellation

For a compactly supported orientation-preserving map `Phi_t`, write

    F_t=D Phi_t,  J_t=det F_t,
    P_t eta(Phi_t(y))=J_t(y)^-1 F_t(y)eta(y).             (27)

Pull velocity back as a one-form,

    b_t(y)=F_t(y)^T B(P_t eta)(Phi_t(y)).                 (28)

The volume Jacobian cancels exactly inside the Biot--Savart integral:

    b_t(y)=F_t(y)^T integral K(Phi_t(y)-Phi_t(y'))
                         F_t(y')eta(y')dy',              (29)
    K(r)a=(4 pi)^-1 a cross r/|r|^3.                     (30)

Let `Phi_dot=h`, `Phi_ddot=k`, `H=Dh`, `Q=Dk`, and
`d_h=h(y)-h(y')`, `d_k=k(y)-k(y')`. Direct differentiation of (29) gives

    B'[h]eta=H(y)^T B eta
       +integral {K(r)H(y')eta+DK(r)[d_h]eta}dy',         (31)

and

    B''[h,k]eta=Q(y)^T B eta
       +2H(y)^T integral {K(r)H(y')eta+DK(r)[d_h]eta}dy'
       +integral {K(r)Q(y')eta
          +2DK(r)[d_h]H(y')eta+DK(r)[d_k]eta
          +D^2K(r)[d_h,d_h]eta}dy'.                      (32)

Every term in (31)--(32) has the original order. Indeed

    |K(r)|<=C|r|^-2,
    |D K(r)| |d_h|<=C||Dh||_infty |r|^-2,
    |D^2 K(r)| |d_h|^2<=C||Dh||_infty^2 |r|^-2.          (33)

The differences `h(y)-h(y')` in (31)--(33) are the load-bearing cancellation;
omitting either produces a false order-zero shape derivative.

### 4.2 Thin torus rather than an axisymmetric surrogate

For dimensionless cross-section coordinates `y=(x,z)` and angular separation
`vartheta=theta-theta'`, the exact circular-tube distance is

    |X_delta(y,theta)-X_delta(y',theta')|^2/R_delta^2
      =delta^2 |y-y'|^2
       +4(1+delta x)(1+delta x')sin^2(vartheta/2).         (34)

In the near region put `vartheta=delta t`. Then

    delta^-2 R_delta^-2 |X-X'|^2
      =|y-y'|^2+t^2
       +delta(x+x')t^2
       +delta^2{x x't^2-t^4/12}+O(delta^3 |t|^4
                                      +delta^4 |t|^6).    (35)

The toroidal phase is retained exactly as

    exp(i ell vartheta)=exp(i k_delta t).                 (36)

Thus the leading near kernel is precisely the straight-column Hodge kernel at
`k_delta`, not its `k=0` value. Expanding the tensor numerator and (35), then
matching to the region `delta<<|vartheta|<<1`, produces the finite first
curvature tensor, the `delta^2 L` local-induction tensor, and the finite second
tensor. The far angular region is smooth. Three Taylor remainders integrated
over the matching annulus obey

    integral_0^1 delta^3/(vartheta+delta) d vartheta
       <=C delta^3 L.                                    (36a)

More precisely, after the homogeneous order-minus-one kernel is factored,
the front-face remainder (`vartheta=delta t`) and the side-face remainder
(`vartheta` fixed) are both bounded by the integrand in (36a), with the same
bound after the cross-sectional derivatives defining the order-minus-one
symbol seminorms. Their overlap coefficients agree because both are Taylor
expansions of the exact expression (34). A partition of unity on these two
regions therefore leaves no unmatched intermediate scale and gives

    ||R_delta^H||_(H^s -> H^(s+1))
       <=C_s delta^3 L,                                  (37)

uniformly for the finite Sobolev interval used below. The axisymmetric Cao
kernel (2.2) gives the carrier coefficients; (34)--(37) give the actual
nonaxisymmetric perturbation operator.

Combining (31)--(37), the collar Hanzawa map, and its identity exterior
extension proves

    ||Lambda^(s+1) partial_(delta,tau)^j
          B_delta Lambda^(-s)||<=C_s,  |j|<=2.            (38)

The global construction has no artificial interface boundary condition:
compact vorticity is extended by zero, (2) is applied in all of `R^3`, and
the resulting velocity has the unique decaying interior/exterior traces.
Equations (19)--(21) ensure that this zero extension contains no sheet.

## 5. Metric Leray jet and equality of the two closures

Let `g_delta=F_delta^T F_delta`. On pulled-back vector fields define

    D_delta v=J_delta^-1 partial_i(J_delta v^i),
    G_delta phi=g_delta^-1 d phi,
    Delta_delta=D_delta G_delta,                          (39)

and on the decaying/mean-zero range

    P_delta^L=I-G_delta Delta_delta^-1 D_delta.           (40)

If `R_delta=Delta_delta^-1`, direct inverse differentiation gives

    R'=-R Delta' R,
    R''=2R Delta' R Delta' R-R Delta'' R.                (41)

Substitution of (41) into (40) gives, for example,

    (P^L)'=-G' R D+G R Delta' R D-G R D',                (42)

with the second derivative obtained by differentiating these three displayed
products. All `J`, cofactor, gradient, divergence, and inverse-Laplacian terms
are therefore present. Uniform ellipticity and the identity exterior extension
give the same `H^s` bounds as (38), at Leray order zero.

For divergence-free compact vorticity, `v=B_delta eta` satisfies
`curl_g v=eta`, `D_delta v=0`, and decay. Uniqueness of this Hodge system gives

    P_delta^L v=v.                                       (43)

Thus (29) and (40) are the vorticity and velocity representations of one
global closure; pressure has not been discarded merely because (1) is written
in vorticity form.

## 6. Common graph domain and the full jet

Differential-form covariance rewrites (1) on the fixed domain as

    Ahat_delta eta=-L_(uhat_delta)eta
                    -L_(B_delta eta)varpi_delta.          (44)

The action coordinate makes

    uhat_delta=omega_delta(I) partial_beta.               (45)

For `m!=0`, `omega_delta>=omega_min>0` makes the graph norm of the first term
uniformly equivalent to one `beta` derivative. For `m=0`, it vanishes and the
compact Kelvin graph from the DA quotient is used. The remaining stretching,
connection, Hodge, collar, and exterior terms are relatively graph bounded.
Equations (16)--(18), (31)--(38), and (44)--(45) consequently give one domain

    D(Ahat_delta)=D_*                                    (46)

and the expansion (3).

Separate the local graph and Hodge remainders. Taylor's formula in the smooth
variables `(delta,tau)` and `tau=delta^2 log delta` gives

    ||R_delta^loc||_(D_* -> X_*)<=C delta^3 L.            (47)

The matched kernel calculation gives (37). Hence (4) follows. The local
transport coefficient variations are block diagonal in exact streamline
Fourier number; longitudinal `i ell delta` terms remain inside
`A_col(k_delta)`. Curvature factors proportional to `cos beta,sin beta` shift
`m` by one, while connection terms without those factors preserve `m`. This
derivation replaces the false rule that all of `C_1` shifts `m`.

The pulled-back KKS form is a coefficient integral of (19)--(24). It is `C^2`
on the same quotient; the harmonic and impulse duals remain finite dimensional.
Thus the KKS-biorthogonal projections used by 0052 share the domain (46).

## 7. Modewise order and the resonant estimate

Let `Lambda` be the positive first-order elliptic weight obtained from the
fixed Kelvin Sturm--Liouville operator and the streamline Fourier weight, with
eigenvectors `e_(m,n)` and eigenvalues comparable to `1+|m|+n`. For an
order-minus-one block `B_j`, repeated commutation gives

    (lambda_a-lambda_b)^N <e_a,B_j e_b>
       =<e_a,ad_Lambda^N(B_j)e_b>.                        (48)

The kernel bounds (33), (37), and `p>=6` control the three commutators needed
for summation. Therefore, for `N=0,1,2,3`,

    |<e_(m,n),B_j e_(m',n')>_E|
      <=C_N <m-m'>^-N <n-n'>^-N
       /<|m|+|m'|+n+n'>.                                 (49)

First curvature maps between the appropriate neighboring streamline sectors.
The inverse nonzero-`m` transport block is order zero in the Kelvin variable,
so the product

    P_0 C_1 Q_perp[Q_perp(A_0-z)Q_perp]^-1
        Q_perp C_1 P_0                                  (50)

is order minus two. Equation (49) with the last denominator squared follows
by convolution; `N=3` makes the double mode sum finite.

The remainder obeys (49) with coefficient `C delta^3 L`. In particular, at

    n_*=Theta(1/(delta L)),                               (51)

the Hodge remainder on a Kelvin vector is

    O(delta^3 L/n_*)=O(delta^4 L^2),                     (52)

while the entire local graph remainder is `O(delta^3 L)`. Both are smaller
than the 0052 contour distance

    r_delta=gamma delta^3 L^2                            (53)

by respectively `O(delta)` and `O(1/L)`. An order-minus-two second Schur block
is `O(delta^2/n_*^2)=O(delta^4L^2)`, also `O(delta)` relative to (53).

This proves the precise resonant bridge. No unweighted `o(delta^2)` is compared
directly with (53).

## 8. Contour-relative resolvent conclusion

Let `A_app,delta` be (3) without `R_delta`, with the exact diagonal Kelvin
frequencies and every mode in the 0052 four-coupling window retained. The
all-sector Grushin result of 0052 gives

    sup_(z in Gamma_delta)
      ||(A_app,delta-z)^-1||_(X_* -> D_*)<=C/r_delta.     (54)

Equations (47), (52), and the far-from-resonance spectral split yield

    sup_(z in Gamma_delta)
      ||R_delta(A_app,delta-z)^-1||_(X_* -> X_*)
       <=C/L+o(1)<1/2                                   (55)

for small `delta`. The resolvent identity then gives

    (Ahat_delta-z)^-1=(A_app,delta-z)^-1
       [I+R_delta(A_app,delta-z)^-1]^-1,                 (56)

and

    sup_(z in Gamma_delta)
      ||(Ahat_delta-z)^-1-(A_app,delta-z)^-1||_(X_* -> D_*)
       <=C/[L r_delta].                                  (57)

Since the contour length is `Theta(r_delta)`, its Riesz projection differs in
graph norm by `O(1/L)`. Rank, real conjugation, and KKS nondegeneracy persist.
Transport by (26)--(29) reconstructs compact DA vorticity and its exact
finite-energy exterior velocity.

Equations (54)--(57) also execute the frozen norm-resolvent fallback. Here it
is a consequence of the strong jet, not a replacement invoked without the
modewise estimate.

## 9. Exact and conditional boundaries

Established in 0054 for the frozen `p>=6` Cao subfamily:

- the centered polyhomogeneous carrier/free-boundary jet (16)--(18), identified
  with the actual source family by uniqueness;
- the full Hodge DA quotient with its harmonic circulation row;
- exact Piola shape derivatives (31)--(32), the genuine nonaxisymmetric
  toroidal distance expansion (34)--(37), and the metric Leray derivatives;
- one common graph domain and the order-minus-one/order-minus-two estimates;
- the contour-relative norm-resolvent and Riesz transfer (55)--(57).

Consequently the `HJ2` condition isolated in 0052 is discharged for this
subfamily. The retained crossing, coupling, spacing, and `O(1/L)` bad-window
scales are unchanged.

Still conditional on the earlier finite-core spectral inputs are the existence
and positivity of the particular bending/Kelvin cluster selected by the 0052
Feshbach theorem. This attempt supplies the graph transfer, not an independent
reproof of those spectral inputs. Even with that cluster, the fixed-order
nonlinear residual is not an exact solution: a nonlinear solid-torus
relative-periodic branch still requires the separate range/transparency and
convergent bifurcation construction. No quantum or particle inference follows.

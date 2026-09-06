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

For integer Cao exponent `p>=6`, this attempt establishes the exact local Cao
cells, the pulled-back Biot--Savart first and second shape identities, the
exact toroidal distance, and a limited no-vortex-sheet trace fact for the
vorticity and its DA tangent. It also derives the candidate expansion

    Ahat_delta=A_col(k_delta)+delta C_1
       +delta^2 L C_(2,log)+delta^2 C_2+R_delta.           (3)

The claimed estimates

    ||R_delta^loc||_(D_* -> X_*)<=C delta^3 L,
    ||R_delta^H||_(Psi^(-1))<=C delta^3 L                 (4)

remain conditional on four named constructions: an exact augmented Cao
carrier map, a uniform common DA graph domain, genuine two-index symbol
bounds, and the full nonnormal all-sector graph resolvent. Consequently 0054
does not discharge `HJ2` and is not a nonlinear rotating-wave existence
theorem.

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

### 2.2 Candidate augmented implicit map and the source gap

Work formally on `B_2` in the affine logarithmic-potential space

    W={w in C^(3,alpha): w-c_log log|y| is harmonic off B_2},   (12)

with fixed `c_log` (circulation), even symmetry in `y_2`, and center gauge
`partial_1 w(0)=0`. The exact core integral equation, circulation condition,
and center condition are available. Cao's system (3.36), however, is an
auxiliary two-parameter asymptotic system defining comparison parameters. It
is not an exact scale/speed row of the steady problem: Corollary 3.12 supplies
the true equations only with `O(epsilon^2|log epsilon|)` errors, and
Proposition 3.13 records nonzero differences between the actual and auxiliary
parameters.

A successful continuation must construct from the exact steady integral
equation and exact parameter relation an augmented map

    F(q,tau;w,mu,c_log)=0,                                (13)

with an explicit Banach domain/codomain and exact restriction relating
`epsilon`, `q`, `tau`, circulation, and prescribed speed. The formal kernel
expansion following (10) identifies the candidate `(q,tau)` coefficients, but
does not prove that `F` is a `C^3` map on an open neighborhood.

The limiting profile part would be `L_U`. Source Lemma 3.8 gives

    ker L_U=span{partial_1 U,partial_2 U}.                (14)

Evenness removes `partial_2 U`, the center gauge removes `partial_1 U`, and

    -2 pi U'(1)=integral_(B_1) U^p dy>0                  (15)

provides a nonzero limiting boundary coefficient. Lemma 3.9 is an a-priori
estimate for the source operator `L_epsilon`; it does not prove Fredholm index
zero or surjectivity of the newly lifted augmented derivative. The
low-frequency/far-log compatibility, index, and parameter-row Schur
determinant remain to be constructed.

If those exact mapping and invertibility statements were proved, the implicit
function theorem would produce the candidate expansion

    w(q,tau)=U+qV+tau Z_log+q^2 Z+Rem_w,                 (16)
    ||Rem_w||_(C^(3,alpha))
       <=C(|q|^3+|q tau|+tau^2),                         (17)

and restriction `tau=q^2 log q` would give

    ||Rem_w||_(C^(3,alpha))<=C q^3 |log q|.              (18)

Equations (16)--(18) are therefore a formal polyhomogeneous candidate, not a
constructed Cao branch. Theorem 1.6 can identify an independently constructed
exact family only after its concentration and parameter hypotheses have been
verified; it cannot create the missing branch or supply augmented
surjectivity. The exact cells (7)--(8) remain established. The claimed smooth
conversion from `q` to the physical thinness `delta=a_delta/R_delta` also
belongs to the missing exact augmented-map construction.

### 2.3 Limited interface regularity and absence of a generated sheet

In a fixed Hanzawa coordinate let `d=1-s` be inward normal distance. For the
source's actual smooth polynomial profile, Hopf behavior gives

    Omega_delta=d_+^p a_delta(s,beta) tau_delta,
    0<c<=a_delta<=C.                                     (19)

Formal parameter differentiation yields

    partial_delta Omega_delta
       =p d_+^(p-1)(partial_delta d)a_delta
          +d_+^p partial_delta a_delta,                  (20)

and a second derivative whose lowest power is `d_+^(p-2)`. For `p>=6`, these
coefficients have zero-extension traces and create no interface delta measure.
For a DA tangent

    eta=-L_xi Omega_delta=d(i_xi Omega_delta),            (21)

the lowest boundary power is `d_+^(p-1)`, again with zero vortex-sheet trace.
This establishes the limited no-sheet statement for `Omega_delta` and `eta`.
It does not establish the required traces of `A_delta eta` after two parameter
derivatives or freeze the full Sobolev interval of the graph calculation.

## 3. Candidate common DA quotient

Let `mu_delta` and `varpi_delta=i_(Omega_delta)mu_delta` denote the pulled-back
volume and vorticity two-forms. A no-normal-flow displacement satisfies

    d(i_xi mu_delta)=0.                                  (22)

On a solid torus, the relative Hodge decomposition would take the form

    i_xi mu_delta=d alpha+h xi_harm,                     (23)

where the harmonic row records meridional flux/circulation. In action--angle
coordinates the candidate density is

    varpi_delta=f_delta(I)dI wedge d beta.                (24)

The one-dimensional estimate

    ||d^(p-1)g||_(H^-1(0,1))
       <=C ||partial_d(d^p g)||_(H^-1(0,1))               (25)

is compatible with this weighted representation. It does not by itself prove
the reverse estimate, closedness of the DA energy range, uniform boundedness
of multiplication by `f_delta/f_0` and its inverse on the `H^-1`/graph
completion, or boundary/gauge control of the harmonic row. Therefore

    T_delta:[Y_*/Stab_*] -> [Y_delta/Stab_delta]          (26)

is the required common-domain isomorphism, not an established conclusion.
Its proof must freeze a Sobolev interval, retain tangential zero modes and the
harmonic row, and estimate the composite whole-space Hodge projectors on the
actual zero-moment/DA spaces. Cross-`delta` maps remain operator
identifications; only physical same-carrier displacements are coadjoint.

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
symbol seminorms. Their overlap coefficients agree formally because both are Taylor
expansions of the exact expression (34). The scalar bound (36a) does not by
itself prove the cross-sectional derivative seminorms required for

    ||R_delta^H||_(H^s -> H^(s+1))
       <=C_s delta^3 L.                                  (37)

That estimate remains conditional on explicit front, side, and overlap symbol
bounds with every anisotropic `delta` factor retained on a frozen Sobolev
interval. The axisymmetric Cao
kernel (2.2) gives the carrier coefficients; (34)--(37) give the actual
nonaxisymmetric perturbation operator.

The desired consequence of (31)--(37), the collar map, and the exterior
extension is

    ||Lambda^(s+1) partial_(delta,tau)^j
          B_delta Lambda^(-s)||<=C_s,  |j|<=2.            (38)

but (38) is not established until the missing symbol seminorm estimates and
interface graph domains are supplied.

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
are therefore present. The algebraic identities are exact, but `Delta^-1` is not separately
bounded on general inhomogeneous whole-space Sobolev spaces. Bounds must be
proved for the composite Hodge/Leray operators on the actual zero-moment or
homogeneous DA domain. Uniform ellipticity alone does not give the asserted
global estimates.

For divergence-free compact vorticity, `v=B_delta eta` satisfies
`curl_g v=eta`, `D_delta v=0`, and decay. Uniqueness of this Hodge system gives

    P_delta^L v=v.                                       (43)

Thus (29) and (40) are the vorticity and velocity representations of one
global closure; pressure has not been discarded merely because (1) is written
in vorticity form.

## 6. Conditional common graph domain and full jet

Differential-form covariance formally rewrites (1) as

    Ahat_delta eta=-L_(uhat_delta)eta
                    -L_(B_delta eta)varpi_delta.          (44)

The action coordinate would make

    uhat_delta=omega_delta(I) partial_beta.               (45)

A proof of one common closed graph domain

    D(Ahat_delta)=D_*                                    (46)

requires the missing two-sided DA isomorphism (26), traces for
`A_delta eta`, composite whole-space Hodge bounds, and uniform control of the
metric pressure operator. The exact identities (31)--(35) are inputs to that
proof, not a substitute for it. Accordingly the remainders

    ||R_delta^loc||_(D_* -> X_*)<=C delta^3 L,            (47)

and (37), as well as twice differentiable KKS dual projections on the same
domain, remain conditional graph estimates. The exact streamline and Fourier
selection rules remain valid formal bookkeeping: longitudinal `i ell delta`
stays in `A_col(k_delta)`, curvature factors proportional to
`cos beta,sin beta` shift `m` by one, and connection terms without those
factors preserve `m`.

## 7. Required two-index mode bounds

Let `Lambda` be a candidate positive first-order elliptic weight with
`lambda_(m,n)` comparable to `1+|m|+n`. The exact commutator identity is

    (lambda_a-lambda_b)^N <e_a,B_j e_b>
       =<e_a,ad_Lambda^N(B_j)e_b>.                        (48)

A single `Lambda` controls only the combined eigenvalue difference. Large
changes in `m` and `n` can cancel, so (48) does not imply the product bound in
the former equation (49). The required repair is either separate bounded
commutators with `M=-i partial_beta` and a radial Kelvin weight, including
mixed commutators, or one directly proved summable two-dimensional bound in
`|m-m'|+|n-n'|`, together with uniform eigenbasis/Riesz-basis constants.

Conditional on those bounds and the common graph domain, two first-curvature
factors would make the nonzero-`m` Schur product order minus two. At

    n_*=Theta(1/(delta L)),                               (51)

its target scale is `O(delta^4 L^2)`, while the contour distance is

    r_delta=gamma delta^3 L^2.                            (53)

These scale comparisons are preserved. The actual product-decay convolution,
modewise Hodge remainder, and full order-minus-two Schur estimate remain
unproved.

## 8. Conditional contour algebra and missing nonnormal inverse

Let `A_app,delta` denote the formal truncation of (3), retaining every mode in
the finite coupling window. The necessary estimate is

    sup_(z in Gamma_delta)
      ||(A_app,delta-z)^-1||_(X_* -> D_*)<=C/r_delta.     (54)

Scalar form gaps and spectral separation do not prove (54) for this nonnormal
Hamiltonian operator. It requires an all-sector graph-resolvent/Grushin
construction including every low `m`, the full pressure/Hodge block,
axis/interface/collar/exterior pieces, and uniform KKS projection constants.
This is an independent condition in addition to `HJ2`.

If (37), (47), the two-index estimates, and (54) are established, the standard
Neumann-resolvent algebra would give

    (Ahat_delta-z)^-1=(A_app,delta-z)^-1
       [I+R_delta(A_app,delta-z)^-1]^-1,                 (56)

and the advertised graph-resolvent and Riesz differences. Those implications
are exact conditional algebra. Their hypotheses have not been proved here,
so no exact Cao Riesz projection is constructed in 0054.

## 9. Exact and conditional boundaries

Established in 0054:

- the exact Cao local cells (7)--(8) and the source's actual
  `O(epsilon^2|log epsilon|)` parameter remainder;
- the limited `p>=6` zero-extension/no-sheet fact for `Omega_delta` and its DA
  tangent `eta`, without the stronger `A_delta eta` graph traces;
- the exact Piola-conjugated Biot--Savart shape identities (29)--(33), including
  the point-difference cancellation;
- the exact nonaxisymmetric toroidal distance and retained phase (34)--(36);
- the exact algebraic metric-Leray differentiation identities (39)--(43).

Blocked by named constructions are the exact augmented `C^3` Cao carrier map,
the common closed DA graph domain, the full symbol seminorm estimate, genuine
two-index mode decay and order-minus-two Schur bound, and the nonnormal
all-sector graph resolvent. Thus `HJ2`, the norm-resolvent/Riesz transfer, and
the exact nonlinear branch remain open. The corrected crossing, coupling,
spacing, and `O(1/L)` scale algebra from 0052 is unchanged. This route verdict
does not imply instability, quantization, an electron, or a neutrino.

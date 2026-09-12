# Boundary topology, physical seeds, and a massive one-contour Kelvin transfer

## 1. Conventions and the source boundary jet

Write the fixed Cao ring in its translating, nonrotating frame as

    omega_0=zeta(r,z) partial_theta,
    W_0=B_R3 omega_0-c_0 e_z,
    A_0 eta=-[W_0,eta]-[B_R3 eta,omega_0],
    C_0 xi=-[xi,omega_0]=curl(xi cross omega_0).          (1)

All brackets and KKS signs are those of 0005 and 0062.  Put

    P_epsilon=G_1 zeta_epsilon
       -(W/2)r^2 log(1/epsilon)-mu_epsilon,
    zeta_epsilon=epsilon^(-2)(P_epsilon)_+^p.             (2)

The reviewed input was only `zeta=O(d^p)`.  The stronger boundary coefficient
is, however, an output of the exact source equations.  Cao et al. Appendix A,
Lemma A.2 parameterizes the scaled free boundary and proves in (A.8)

    partial_t P_tilde_epsilon
      =C_epsilon U'(1)+o(1)<0,                            (3)

uniformly on that boundary; `C_epsilon>0` and the limiting radial ground
state has `U'(1)<0`.  The local equation

    -Delta P_epsilon+r^(-1) partial_r P_epsilon
       =r^2 epsilon^(-2)(P_epsilon)_+^p                  (4)

is uniformly elliptic because the thin torus stays away from the axis.
For each fixed sufficiently thin ring, the bounded right-hand side in (4)
first gives local `W^(2,q)`/`C^(1,alpha)` regularity. Since `(t_+)^p` is at
least `C^5` for `p>=6`, local Schauder bootstrapping gives `P_epsilon` at
least `C^2` near its compact regular level set. Equation (3), the implicit
function theorem and normal Taylor expansion therefore give

    P_epsilon=kappa_P(s)d+O(d^2),       kappa_P(s)>0,
    zeta_epsilon=kappa_zeta(s)d^p+O(d^(p+1)),
    0<c<=kappa_zeta(s)<=C.                               (5)

The constants are uniform along that fixed compact boundary. This derives the
nonzero jet and two-sided `d^p` comparability for the fixed carrier. It does
not provide uniform-in-`epsilon` fixed-domain `C^2` convergence or global
critical-point/period control. The upper estimates used by the 0062 collar
counterexample require only the exact positive-part form and do not depend on
the lower bound in (5).

The relative velocity satisfies

    W_0=r^(-1) grad^perp P_epsilon.                       (6)

The missing `C^2` control follows from the exact transformed source equation,
not from Proposition 3.2 alone. On a fixed scaled ball containing the full
cross-section, the Cao profile obeys

    -Delta w_delta+[delta/(1+delta y_1)]partial_1 w_delta
       =(1+delta y_1)^2(w_delta)_+^p.                    (6a)

Proposition 3.2 gives `w_delta->U` in `C^1` there. Since `(t_+)^p` is `C^5`
for `p>=6`, the required Schauder upgrade is obtained on nested balls rather
than directly from `C^1` convergence. Choose
`B_L compactly contained in B_(L+1) compactly contained in B_(L+2)`. The
uniform `L^infinity` bound on `w_delta` makes the semilinear right-hand side
uniformly `L^infinity`; the coefficient
`delta/(1+delta y_1)` and its derivatives are uniformly bounded on
`B_(L+2)`. Interior two-dimensional elliptic estimates first give, for any
fixed `q>2`,

    ||w_delta||_(W^(2,q)(B_(L+1)))<=C_(L,q),
    ||w_delta||_(C^(1,gamma)(B_(L+1)))<=C_L,
    gamma=1-2/q>0.                                      (6b)

Fix `0<alpha<gamma`. Interpolation of the uniform `C^(1,gamma)` bound with
the source `C^1` convergence shows convergence of the nonlinear term in
`C^(0,alpha)`. The full lower-order term is retained:

    [delta/(1+delta y_1)]partial_1 w_delta
       ->0 in C^(0,alpha)(B_(L+1)),                       (6c)

because its coefficient is `O(delta)` in `C^1` and
`partial_1w_delta` is uniformly `C^(0,gamma)`. Thus the complete right-hand
side of the difference of (6a) and `-Delta U=U_+^p` tends to zero in
`C^(0,alpha)`. Interior Schauder on the smaller ball gives

    ||w_delta-U||_(C^(2,alpha)(B_L))->0.                 (6d)

The support is contained in `B_L` for fixed large `L`, and the outer regular
level is covered by (3), so this fixed-domain estimate covers the center,
every positive level, and its boundary. The radial ground state has one
nondegenerate maximum at the origin. On a small center ball, (6d) and the
quantitative implicit-function theorem give one and only one nearby critical
point `y_delta`, with

    D^2w_delta(y_delta)->D^2U(0)=-lambda_c I,
    lambda_c>0.                                         (6e)

On its compact complement, `|grad U|` has a positive lower bound, so (6d)
excludes every secondary critical point. The boundary is a single
regular closed curve by Lemma A.2. Hence every intermediate positive level is
a regular closed curve and there is no separatrix.

The period integral is continuous on regular levels. At the center, (6e) and
the fixed positive cylindrical weight give two-sided bounds on the
linearized elliptic frequency, hence on the limiting center period. At the
outer level, (3) gives two-sided bounds on `|grad P_epsilon|`; the ring stays
away from the axis and the regular boundary lengths have positive lower and
finite upper bounds, so

    0<T_c^-<=lim_(I->0)T(I)<=T_c^+<infinity,
    0<T_b^-<=lim_(I->I_b)T(I)<=T_b^+<infinity.           (6f)

Continuity on the intervening compact interval supplies the same two-sided
bound there. Thus, with `omega=2*pi/T`, in the fixed carrier's action--angle
coordinates,

    W_0=omega(I) partial_beta,
    0<omega_min<=omega(I)<=omega_max<infinity.            (7)

The constants are for each sufficiently thin fixed carrier; no uniform raw
physical-time normalization across the family is asserted. If a different
Cao member falls outside the `C^2` perturbative regime, it must instead be
partitioned into regular cells with every separatrix frequency retained.

## 2. Route A: the three topologies really separate

The 0063-corrected 0062 construction supplies, at `s=4`, `p>=6` and fixed
`|n|>=2`, divergence-free collar displacements with

    ||[xi_h]||_(H^5/ker C_0)>=c,
    ||C_0 xi_h||_(H^3)=O(h^(p+1)),
    ||C_0 xi_h||_(H^4)=O(h^p),
    ||B_R3 C_0 xi_h||_(H^4)=O(h^(p+1)),
    ||A_0 C_0 xi_h||_(H^3)=O(h^p).                       (8)

Thus both ordinary ambient and ordinary generator-graph lower bounds are
refuted.  The exact Hodge identity used here is

    B_R3 C_0 xi=P_L(xi cross omega_0),                    (9)

with the compact-curl zero moment and decaying low-frequency row; (9) is not
an unconditional inhomogeneous zero-frequency Biot--Savart gain.

There is nevertheless a natural candidate orbit tangent topology. Let
`G_vol,c^5` be the local group of orientation- and volume-preserving `H^5`
diffeomorphisms equal to the identity outside a fixed compact neighborhood of
the ring, with the translation, rotation, phase, center, circulation and
stabilizer rows fixed.  Define

    O_w^4={g_*omega_0:g in G_vol,c^5},
    ||C_0 xi||_(X_orb,w^5)
       =inf_(k in ker C_0)||xi+k||_(H^5).                 (10)

The completion of the tangent image in (10) is separated by definition of the
quotient; smooth compact displacements are dense. This is not a disguised
ambient `H^3` norm, as (8) proves. Equation (10) is only a tangent norm until
reverse Hardy/interface control and descent through distinct displacement
representatives are proved.

For `H^5` maps in three dimensions, the chain rule, the `H^4` algebra
estimate and the inverse-function formula give, on a sufficiently small
bounded neighborhood,

    ||g o h-id||_(H^5)
       <=C(||g-id||_(H^5)+||h-id||_(H^5)),
    ||g^(-1)-id||_(H^5)<=C||g-id||_(H^5),                (11)

and the corresponding difference estimates with one displayed high norm and
lower norms on the other factors.  Volume is preserved under both operations.
These are operations on the ambient volume-preserving diffeomorphism group;
they do not by themselves descend to the quotient completion (10). Formula
(5) identifies candidate two-sided `d^p` coordinate weights, but equivalence
to (10) remains the reverse Hardy/interface theorem.

For any `H^5` volume-preserving path `g_t`, put
`u_t=dot g_t o g_t^(-1)` and `omega_t=(g_t)_*omega_0`.  With the convention
in (1),

    dot omega_t=-[u_t,omega_t]=curl(u_t cross omega_t),
    omega_1-omega_0=curl F_g,
    F_g=integral_0^1 u_t cross omega_t dt.                (12)

Consequently the complete velocity, including its exterior tail, is

    B_R3(omega_1-omega_0)=P_L F_g.                        (13)

For every individual smooth path, the whole-space Leray projection gives an
`H^4` velocity difference and an `H^3` pressure gradient. Compact support of
the curl, zero circulation change and decay fix the harmonic row. Uniformity
and descent of these maps to the completed quotient, including the moving
interface, are not proved. Equations (10)--(13) establish the tangent norm and
exact smooth-path Hodge identity, not the frozen nonlinear weighted chart.

## 3. Route B: an exact physical DA test seed and the symmetry audit

Fix any `|n_0|>=2` and a positive-core coordinate patch. Components in (15)
are physical orthonormal cylindrical components in `(e_r,e_theta,e_z)`, not
coordinate coefficients in `(partial_r,partial_theta,partial_z)`. Choose real
`f,g in C_c^infinity` in that patch with

    integral r^2 zeta f g dr dz !=0.                      (14)

Set

    xi^r=f,       xi^z=i g,
    xi^theta=(i r/n_0)(partial_r f+f/r)
                 -(r/n_0)partial_z g,                    (15)

times `exp(i n_0 theta)`.  Direct substitution gives `div xi=0`.  After the
finite-dimensional row correction of 0062, which can be chosen without
destroying (14),

    q_test=C_0 xi in C_c^infinity(K_0) intersect D(A_0)  (16)

is a real-conjugate, finite-energy, whole-space-Hodge physical DA seed.  Its
KKS pairing is explicitly

    Omega_KKS(q_test,conjugate(q_test))
      =4 pi i rho_0 integral r^2 zeta f g dr dz !=0.      (17)

Thus absence of a DA seed or of a physical symplectic normalization is not the
obstruction.

The eigenstate condition is much stronger.  The exact intertwiner gives

    (A_0-i nu)C_0 xi
       =C_0{P_L(xi cross omega_0)+[xi,W_0]-i nu xi}.       (18)

On the positive nonzero-`n` core, `C_0` is injective.  Equation (18) is the
actual full-pressure eigenproblem; (15) does not solve it for any asserted
`nu`.

All exact symmetry candidates can be classified.  For a constant translation
`a`, Euclidean equivariance gives

    B_R3 C_0 a=-[a,B_R3 omega_0],
    A_0 C_0 a=0.                                         (19)

The transverse translations occupy `n=+/-1`, the axial translation occupies
`n=0`, and

    Omega_KKS(C_0 a,C_0 b)
       =rho_0 integral omega_0 dot(a cross b) dx=0        (20)

by angular integration of the purely azimuthal vorticity.  They are isotropic
zero modes removed by the center/impulse quotient.  For a rotation generator
`R`, equivariance and (18) give

    A_0 C_0 R=C_0{-c_0[R,e_z]},                           (21)

so tilted rotations are generalized zero chains ending in translations;
axial rotation stabilizes the ring.  Hence the symmetry-seed route is
refuted for the required nonzero-frequency, nonzero-Krein internal state.  It
does not refute internal Kelvin modes.

The toroidal nonlinear selection rule is exact even before a mode is found.
Because the carrier has toroidal number zero, a real seed in `+/-n_0` produces
at homogeneous order `k` only

    n=(k-2j)n_0,        j=0,...,k.                        (22)

In particular the quadratic rows are `0,+/-2n_0`.  There is no corresponding
finite cross-sectional `m` lattice: the exact toroidal Green kernel and the
curved action--angle coefficients have generally infinite Fourier matrices in
`m` already at linear order.  The honest support recursion is the convolution
of the current source support with those source-defined coefficient supports.
Thus a finite nearest-neighbor `m` truncation is refuted as an exact selection
rule, while the toroidal rule (22) survives.

On smooth compact DA tangents, the relative Hessian

    E_2=d^2(H-c_0P_z)|_(omega_0,orbit)                    (23)

and `A_0=J_0E_2` imply

    E_2(u,v)=Omega_KKS(A_0u,v),
    Omega_KKS(A_0u,v)+Omega_KKS(u,A_0v)=0.                (24)

There is no physical free-boundary flux at this smooth-core scope: the
`d^p`, `p>=6`, traces
vanish, the whole-space Green operator is self-adjoint in the energy pairing,
and its core/exterior interface contributions cancel.  A nonzero flux appears
only if the core is artificially truncated without its exterior Hodge row.
Finite circulation and symmetry rows are fixed before applying (24).
Extension of (24) to the completed weighted/interface graph is part of the
open descent theorem. On any domain where that extension is proved, an
isolated mode `A_0q=i nu q` has left functional
`v -> Omega_KKS(conjugate(q),v)` up to residue normalization.  This constructs
the Green identity and adjoint trace formula, but no `q`, `V_*`, or first
source-bearing resonance is named by (24).

## 4. Route C: the fixed-`n` essential gap

For a fixed integer `n!=0`, define the ambient physical DA space by

    X_n^3=closure_(H^3(R3)){
      C_0 xi: xi smooth, divergence-free, n-harmonic,
      C_0xi compactly supported in K_0, finite rows fixed}. (25a)

All elements are zero-extended physical vorticities; no displacement norm is
put on `X_n^3`. Define

    T_n eta=-[W_0,eta],       K_n eta=-[B_n eta,omega_0]. (25)

The closed transport domain and full generator domain are

    D(T_n)={eta in X_n^3:T_n eta in X_n^3},
    D(A_n)=D(T_n),             A_n=T_n+K_n,               (25b)

where equality in (25b) follows from the bounded compactness result below.
Use the upper-semi-Fredholm essential spectrum

    sigma_ess^+(A_n)={lambda:A_n-lambda is not left Fredholm}. (25c)

Failure of left Fredholmness implies failure of Fredholmness, so (25c) is
contained in every standard non-Fredholm essential spectrum used here. In the
action--angle frame, the local flow is

    (I,beta,theta) -> (I,beta+omega(I)t,theta),
    D phi_t=[[1,0,0],[t omega'(I),1,0],[0,0,1]].          (26)

Thus the vector transport monodromy is unipotent.  In streamline harmonic
`m`, `T_n` is multiplication by `-i m omega(I)` plus a nilpotent shear; the
shear changes Jordan growth but not the spectral set.

For fixed `n`, `K_n:X_n^3->X_n^3` is compact. Every bounded sequence in
`X_n^3` has curl support in the same compact `K_0`. Divergence freedom and
compact support give

    integral_R3 eta dx=0,
    eta_hat(k)=O(|k|),                                  (26-)

with the first moments controlled uniformly by the fixed-support `H^3` norm.
Thus the whole-space multiplier `i k cross/|k|^2` has no low-frequency
singularity on this sequence. The decaying harmonic/circulation row is fixed,
so

    B_n:H^3_(K_0)->H^4(R3)

is bounded, including the low-frequency part. Its restriction to a fixed
neighborhood of `K_0` is compact into `H^3` by Rellich. In

    [B_n eta,zeta partial_theta]

the only derivative falling on `B_n eta` in the second bracket leg is the
fixed multiplier `i n` plus the cylindrical connection, hence is zero order
for this fixed fiber. Multiplication by the `C^4` carrier coefficients and
the local Rellich embedding prove compactness of `K_n`; the exterior velocity
and regular-axis limit are retained by the global `B_n` construction. This is
a fixed-`n` statement only and gives no bound uniform in a varying `n_0`.

The action coordinate is not an arbitrary streamline label. Define

    I(s)=(2*pi)^(-1) integral_{P_0>s} r dr dz,            (26-0)

with orientation chosen so `I` increases outward, and choose periodic
`beta in R/(2*pi Z)` by normalized travel time. Then

    r dr dz=dI d beta,
    dV=dI d beta d theta,

and both streamline and toroidal characters `m,n` are integers.
The constrained band inclusion can be proved on every regular cell. Use
volume action--angle coordinates, normalized so that the physical volume form
on the cell is `dI d beta d theta`. Fix a regular `I_0`, choose
`epsilon_j->0`, `N_j->infinity` with `N_j epsilon_j->infinity`, and set

    a_(m,n)=partial_beta-(m/n)partial_theta,
    xi_j=A_j a_(m,n) exp(i m beta+i n theta)
          exp(i N_j(I-I_0)) chi((I-I_0)/epsilon_j).       (26a)

Here `m,n` are integers and `n!=0`; the cutoff is supported in a compact
subannulus of the regular cell. This polarization is exactly tangent to
`I=constant`, lies in
the kernel of the nilpotent shear in (26), and is exactly divergence-free:
`partial_beta xi_j^beta+partial_theta xi_j^theta=0`, while
`xi_j^I=0`. Thus no merely formal high-frequency divergence projection is
being used. In these coordinates

    C_0 xi=-[xi,zeta(I)partial_theta]
          =i n zeta xi-xi^I zeta'(I)partial_theta.       (26b)

Consequently `C_0 xi_j=i n zeta xi_j` exactly. Positive-core `C_0` is an
order-zero multiplier in the high-`I` oscillation, preserves the
tangent-to-`I` polarization, and is bounded above and below on the packet.
Equivalence of coordinate and physical Sobolev norms on the compact cell and
`N_j epsilon_j->infinity` give, for unnormalized point amplitude `A_j`,

    ||C_0xi_j||_(H^3)
       =Theta(A_j N_j^3 epsilon_j^(1/2)).                 (26c)

Choose `A_j=(N_j^3 epsilon_j^(1/2))^(-1)` and divide by the remaining fixed
Theta-factor, so the resulting raw packet `q_j^0=C_0xi_j` is smooth, compact,
belongs to `D(A_n)`, and has `H^3` norm one.
Because `q_j^{0,I}=0`, the shear term also vanishes exactly and the product
rule gives the exposing estimate

    ||(T_n+i m omega(I_0))q_j^0||_(H^3)
       <=C(epsilon_j+N_j^(-1)).                           (26d)

Indeed the undifferentiated factor `omega(I)-omega(I_0)` costs
`O(epsilon_j)`, while each derivative falling on that factor removes at least
one power of `N_j`. Derivatives of the cutoff are subordinate because
`N_j epsilon_j->infinity`. The normalized packets are weakly null in `H^3`
by oscillation (and shrinking support), first for smooth tests and then by
density.

All divergence correction is therefore exactly zero in this volume
action--angle gauge. The fixed toroidal character already annihilates the
axisymmetric symmetry, impulse, circulation and center rows. If the chosen DA
closure contains additional continuous same-character finite rows, their
values on `q_j^0` are `o(1)` by weak convergence. Choose the correcting
profiles in the same fixed-`n` smooth compact DA sector and in `D(A_n)`, with
an invertible row matrix. Since
`A_nq_j^0=-i*m*omega(I_0)q_j^0+o_(H^3)(1)`, the packets are weakly null also
in the graph topology. Subtracting the row coefficients times those profiles
therefore produces normalized packets `q_j=q_j^0+o_(D(A_n))(1)` and changes
the residual by `o_(H^3)(1)`. Hence

    ||q_j||_(H^3)=1,
    q_j weakly ->0,
    ||(T_n+i m omega(I_0))q_j||_(H^3)
       <=C(epsilon_j+N_j^(-1))+o(1).                     (26e)

Compactness of `K_n` gives `K_nq_j->0`. Hence, with the sign absorbed into
`m`,

    i m omega(I_0) in sigma_ess(A_0|X_n^3)               (27)

for every regular cell, regular `I_0`, and integer `m`. To expose the
essential-spectrum implication, suppose `A_n-lambda` had a bounded left
regularizer modulo compact operators on its graph: then

    R(A_n-lambda)=1-C,       C compact.                  (27a)

Applying (27a) to the normalized graph-weak packet gives a right-hand side
tending to zero, contradicting `||q_j||_(H^3)=1`. Thus `A_n-lambda` is not
left Fredholm, so `lambda` lies in `sigma_ess^+(A_n)` and hence in the
standard non-Fredholm essential spectrum. This is the required constrained DA
Weyl inclusion. It does not prove the converse
direct-integral/domain inclusion, nor classify boundary or separatrix
spectrum. Conditional on the missing converse constrained-domain theorem, the
positive essential-gap candidate from the established period bound (7) is

    G_ess=(0,omega_min).                                  (28)

The `m=0` compact Kelvin block can have discrete eigenvalues in (28)
accumulating only at zero. Thus (28) is not yet a proved essential gap,
resolvent interval, or mode-existence assertion. If the foliation theorem
fails, every regular cell and separatrix frequency must instead enter (27).

## 5. The positive massive high-toroidal transfer target

The small-`k` bending scaling is not the only same-carrier limit.  Choose
`k_0>0` and, for every sufficiently thin Cao member with core-to-major-radius
ratio `delta`, choose the integer

    n_delta=nearest(k_0/delta),
    k_delta=n_delta delta,
    |k_delta-k_0|<=delta/2.                               (29)

This is one fixed physical toroidal harmonic on each fixed ring.  The
positive mass `k_delta>=k_0/2` removes the low-frequency logarithmic Hodge
singularity.

For the limiting compact column, the `m=0` velocity equation is the positive
Sturm--Liouville pencil

    H_k u_j=lambda_j(k) Phi u_j,
    sigma_j(k)=k/sqrt(lambda_j(k)),
    H_k=-partial_r partial_r^*+k^2,
    Phi=2 Omega W>0,                                     (30)

with regular axis and decaying `K_1(kr)` exterior matching. Multiplication by
`r` gives the generalized Sturm--Liouville problem

    -(r u')'+(r^(-1)+k^2 r)u=lambda r Phi u,
    L_Phi=integral_0^a sqrt(Phi(r))dr.                   (30a)

Here `0<L_Phi<infinity`; the compact-edge weight vanishes like `(a-r)^p`.
Dirichlet--Robin bracketing away from the regular-axis and vanishing-weight
endpoints, followed by shrinking the endpoint collars, gives the generalized
Sturm--Liouville Weyl law

    sqrt(lambda_j(k))=pi*j/L_Phi+O(1).                  (30b)

This is the 0048 equation (13)--(16) calculation retained by the independent
0053 review, now with `L_Phi` explicit. Endpoint singularities alter the
bounded phase, not the coefficient. The eigenvalues are simple and therefore

    sigma_j(k_0)=k_0 L_Phi/(pi j)+O(j^(-2)) ->0.          (31)

The associated wave equation

    H_k partial_t^2u+k^2 Phi u=0                         (32)

gives positive constrained-energy/Krein sign. Using (7), one can choose
consecutive indices `J,J+1` with both
frequencies below `omega_min/3`. Their simplicity holds inside the `m=0`
Sturm--Liouville sector only. Isolation from the **entire** column spectrum is
a separate obligation: put every `m` sector on an analytic common domain,
prove uniform large-`|m|` resolvent control, and show that the finitely many
remaining nonzero-`m` branches neither coincide identically with nor enter the
two proposed circles. Merely choosing generic `k_0` outside a named discrete
crossing set does not prove this.

Here a one-contour curved transfer may avoid the conditional small-`k`
all-sector `HJ2/GR` theorem. Rescale the torus cross-section and use
the orientation-preserving Cao free-boundary chart supplied by (3)--(5), with
the full Piola Jacobian/cofactor factors.  On a fixed core/collar and for
`k_delta in [k_0/2,2k_0]`, the Fourier-reduced Hodge form is uniformly
coercive:

    a_delta(h,h)>=c(||grad h||_2^2+k_0^2||h||_2^2).      (33)

The required vector Hodge estimate is

    ||B_(delta,n_delta)-B_(col,k_0)||_(H^3_comp->H^4_loc)
       <=epsilon_delta,          epsilon_delta->0.        (34)

To prove (34), one must write the full vector Piola transform of curl,
divergence and Leray pressure, retain every `delta*n_delta=O(1)` contribution,
and establish the moving-interface, remote-axis and exterior maps together
with the exact `D(A_col)->H^3` coefficient ledger. Scalar coercivity (33),
Cao's `C^1` profile convergence and formal elliptic bootstrapping do not by
themselves prove (34).

The remaining curved first-order coefficients have the form
`delta a(I,beta)partial_beta` plus coefficients tending to zero.  Since
`omega>=omega_min`, uniformly over every nonzero `m`,

    ||delta a partial_beta u||_(H^3)
       <=C delta(||T_col u||_(H^3)+||u||_(H^3));          (35)

the `m=0` row contains no such derivative. Conditional on the complete vector
ledger for (34), the desired
common-domain graph-relative estimate is

    ||(Ahat_(delta,n_delta)-A_(col,k_0))u||_(H^3)
      <=rho_delta(||A_(col,k_0)u||_(H^3)+||u||_(H^3)),
    rho_delta->0.                                        (36)

Raw physical `H^3` norms and the rescaled graph norms are equivalent for each
fixed `delta`; their constants are not asserted uniform as
`n_delta->infinity`.

Conditional also on whole-column isolation, (36) and the exact resolvent
identity would yield

    (z-Ahat_delta)^(-1)
      =(z-A_col)^(-1)
       [1-(Ahat_delta-A_col)(z-A_col)^(-1)]^(-1),         (37)

once `rho_delta sup_Gamma||(z-A_col)^(-1)||_(X->D)<1`. Only after that estimate
would the Riesz ranks be preserved and the two column eigenpairs transfer to
Cao modes

    A_0 q_(delta,a)=i nu_(delta,a)q_(delta,a),
    nu_(delta,a)->sigma_(J+a-1)(k_0),       a=1,2.        (38)

Preservation of the DA range by the transformed Piola graph is part of the
unproved interface ledger. Literal `C^infinity` regularity across the support
edge is also not claimed: for finite polynomial `p`, the nonzero jet (5)
makes a generic boundary-moving tangent only finitely regular. The physically
classical frozen `H^3` graph scope is the intended target. If the transfer is
proved, continuity of (23)--(24) and positivity in (32) would preserve the
positive Krein sign and permit the exact normalization

    E_2(q_a,conjugate(q_a))=1,
    Omega_KKS(q_a,conjugate(q_a))=-i/nu_a.               (39)

Equations (38)--(39) are a conditional positive construction, not established
objects in 0066. No production numerics or conditional 0052/0054 theorem is
used to promote them.

## 6. Nonlinear and P2 boundary

For any pair eventually produced by (38), the toroidal selection rule (22)
is exact. The candidate modes share the same `n_delta` and differ in radial
Kelvin index. A first-order
axisymmetric area-preserving cross-sectional generator has only tangential
`m=0` action on radial column modes; the obvious quadrupole generators have
`m=+/-2` and hence zero direct `m=0` matrix element.  Consequently two
noncommuting projected deformation matrices required by the 0051 doublet do
not follow from mode existence.  Their first possible contribution is the
second-order Schur path

    m=0 -> m=+/-2 -> m=0,                                (40)

with the full exterior Hodge resolvent between the two vertices.  Its two
ordered quadrupole phases need not commute, but the physical integrals have
not been evaluated here.

Even a proved (38) would supply only a linear center pair, not a nonlinear
rotating wave. The immediate work remains the reverse Hardy/interface chart,
the converse constrained direct-integral/domain theorem, full vector high-`n`
Piola ledger and entire-column isolation. Only after those may a
two-mode source-specific Lyapunov--Schmidt map use (24), (37), (22), and (40).

## 7. Exact input boundary for a later small-charge continuation

Let a compact scalar charge label be chosen from the actual carrier foliation,

    chi_0=F(I),             W_0 dot grad chi_0=0.          (41)

For `|c_0|<c_EM`, a gauge-fixed comoving Maxwell solve has principal operator

    L_EM=-Delta+(c_0^2/c_EM^2)partial_z^2,                (42)

whose Fourier symbol is
`|k_perp|^2+(1-c_0^2/c_EM^2)k_z^2>0`.  Standard decaying elliptic estimates
therefore give `E,B=O(g)` for charge `g chi_0`, and its Lorentz self-force is
`O(g^2)`.  Equation (41) proves advection only.  It does not make
`rho_q E+J cross B` a pressure gradient; an arbitrary advected charge is not
a stationary charged Euler solution.

There is a shorter **steady axisymmetric** route than the full dynamical graph.
In Lorenz gauge, first solve the decaying comoving potentials

    L_(EM,c) phi=g chi/epsilon_EM,
    L_(EM,c) A=mu_EM g chi u,
    E=c partial_z A-grad phi,       B=curl A,             (43)

where `L_(EM,c)` is (42), together with its gauge row. Put
`f_EM=g chi(E+u cross B)`. With
`P=psi-(c/2)r^2-mu` and
`{a,b}=a_r b_z-a_z b_r`, the genuine no-swirl steady/free-boundary map is

    L_axis psi-zeta=0,
    {P,zeta}=(curl f_EM)_theta,
    {P,chi}=0,                                            (44)

plus (43), the support/positive-part constitutive row at `g=0`, and center,
impulse, circulation, charge and translation gauges. This is not Cao's
auxiliary parameter equation (3.36).

At order `g^2`, write `zeta_0=f(P_0)` and
`f_EM=g^2 f_2+O(g^4)`. The curl row of (44) becomes exactly

    {P_0,h_2}=(curl f_2)_theta,
    h_2=zeta_2-f'(P_0)P_2.                               (45)

On every closed level `P_0=s`, (45) has a single-valued solution only if

    integral_(P_0=s) (curl f_2)_theta/|grad P_0| dl=0.   (46)

This is a continuum of streamline compatibility rows; total Maxwell
self-force or the one global translation Noether identity does not imply
(46). If (46) holds, the zero-mean solution `h_2` is obtained by integration
along the streamline, up to a profile function of `P_0`, and the remaining
equation is the Cao steady elliptic linearization

    [L_axis-f'(P_0)]P_2=h_2+F_2(P_0),                    (47)

with its actual translation kernel. Cao's steady elliptic nondegeneracy and a
translation Lyapunov--Schmidt row can then solve (47), after the center,
impulse and charge constraints determine the finite parameters. Thus a full
time-linearized Euler inverse is not required for this route.

Equations (43)--(47) define the immediate separate charged-continuation
attempt. It must evaluate (46), construct the nonlinear free-boundary map and
prove its tame remainder before claiming existence. Route A's tangent norm,
Route B's smooth-core Green identity and Route C's conditional high-`n`
contours do not supply those statements. Even a successful steady branch
would not prove tag persistence: the passive label has rearrangement and
filamentation directions. Nor does Cao's constrained kinetic maximizer plus
positive Maxwell field energy automatically give a definite joint Lyapunov
functional; the coupled second variation may be a saddle and needs its own
signature/coercivity audit.

No nonlinear branch, stability neighborhood, P2/LP2, particle identity,
stationary charged continuation, quantum rule, or relativistic bridge is
claimed in this attempt.

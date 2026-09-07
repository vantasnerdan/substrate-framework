# Exact Green bordered map and the charged Cao branch

## 1. Source normalization and physical rows

Write `x=(r,z)` and

    dnu(x)=r dr dz.                                               (2)

Cao's Green operator includes this measure:

    psi(x)=G_1 zeta(x)=integral G_1(x,y) zeta(y)dnu(y)
          =integral G_1(x,y)zeta(y)y_1dy.                         (3)

Let `L=log(1/epsilon_core)` and let `c=W L` be physical translation speed.
The uncharged equation is exactly

    zeta=epsilon_core^(-2)(G_1zeta-c r^2/2-mu)_+^p.              (4)

The circulation and physical axial impulse are

    kappa(zeta)=integral zeta dnu,
    I_z(zeta)=pi rho_m integral r^2 zeta dnu
             =pi rho_m integral r^3 zeta dr dz.                  (5)

A circulation-preserving radial shift of a core concentrated at `R` has

    delta I_z=2 pi rho_m R kappa delta R+o(delta R).             (6)

These conventions fix every radial and `2*pi` factor below.

## 2. A fixed compact-source Green space

Fix `p>=6`, `0<alpha<1`, and one small `epsilon_core`.  Choose compact
annuli `K_0 compactly contained in K compactly contained in {r>0}`, both even
in `z`, such that `K_0` contains the Cao support and a collar.  The base
threshold is strictly negative off `K_0`, and the same holds in a small
neighborhood.  Use

    X={zeta in C_even^(2,alpha)(R^2_+): supp zeta subset K_0}.    (7)

The norm is the global zero-extension norm.  Thus traces through order two
vanish at `partial K_0`; arbitrary `C^(2,alpha)(K)` functions are not extended
by zero.

Cao (2.2)--(2.6) and elliptic regularity away from the axis give

    G_1:X -> C_even^(4,alpha)(K)                                (8)

continuously and compactly back into `X` after multiplication by a supported
coefficient.  Its global representative has the exact regular-axis and
vanishing far-field rows.  No free inhomogeneous `Delta^(-1)` bound is used.

Set

    S(zeta,mu,c)=G_1zeta-c r^2/2-mu,
    F_0(zeta,mu,c)=zeta-epsilon_core^(-2)S_+^p.                 (9)

For `p>=6`, zero extension across `S=0` preserves the two source derivatives
in (7), so `F_0` is `C^2`.  On the tagged band `|grad P|` is bounded below.
A fixed flow box and the coarea formula then make

    zeta -> P -> I_P -> chi_P                                  (10)

`C^2`: the second level-set variation uses curvature controlled by the
`C^(4,alpha)` image in (8).

The P253/0077 subluminal Maxwell inverse maps the compact current to affine
Coulomb potentials and finite-energy fields.  On `K_0`, its eliminated
Lorentz correction is a `C^1` map `M:X times R -> X`.  Hence

    F(zeta,mu,c,tau)=F_0(zeta,mu,c)-tau M(zeta,c),
    tau=g^2,                                                   (11)

is the exact `C^1` charged map needed below.

## 3. Exact linear signs and the soft cell

At the uncharged base define

    a_epsilon=p epsilon_core^(-2)S_0,+^(p-1),
    A_epsilon=I-a_epsilon G_1.                                (12)

Direct differentiation gives

    partial_mu F_0=+a_epsilon,
    partial_c F_0=+a_epsilon r^2/2.                            (13)

The Green kernel is symmetric in `L^2(dnu)`.  Let `s_epsilon` and `R_epsilon`
be the core and mean radii and choose a cutoff, equal to one on the active
core, for

    phi_t,epsilon=cutoff partial_1U((x-(R_epsilon,0))/s_epsilon),
    q_t,epsilon=a_epsilon phi_t,epsilon,
    ell_t,epsilon(h)=<phi_t,epsilon,h>_nu
                     /<phi_t,epsilon,q_t,epsilon>_nu.          (14)

Then `P_epsilon h=q_t,epsilon ell_t,epsilon(h)` and `Q_epsilon=I-P_epsilon`
are an explicitly normalized right-cell/adjoint-functional split in the
physical `r dr dz` pairing.  This is a soft finite-core cell, not an exact
radial symmetry.  The exact Green expansion gives

    ||A_epsilon q_t,epsilon||_Xepsilon
       =o(1)||q_t,epsilon||_Xepsilon,                           (15)

where the scaled density norm is

    ||h||_Xepsilon=epsilon_core^2||h||_infinity
      +epsilon_core^3||grad h||_infinity
      +epsilon_core^4||D^2h||_infinity
      +epsilon_core^(4+alpha)[D^2h]_alpha.                     (16)

On `ell_t,epsilon(h)=0`,

    ||Q_epsilon h||_Xepsilon
      <=C||Q_epsilon A_epsilon Q_epsilon h||_Xepsilon.         (17)

Here is the complete contradiction mechanism for (17).  Normalize the sum of
(16) and the analogous scaled `C^(4,alpha)` potential norm.  The linear
equation expresses the source as `a_epsilon` times its potential plus a
vanishing residual.  Compact Green estimates exclude mass escape to the
axis, exterior, or outside the rescaled core, so the pulled-back potentials
have a nonzero local limit.  That limit solves

    (-Delta-pU_+^(p-1))phi=0.                                 (18)

Cao Lemma 3.8 leaves the two translations.  Evenness removes `partial_2U`.
The approximate right and left cells in (14) converge to
`pU_+^(p-1)partial_1U` and `partial_1U`, so the normalized orthogonality
removes `partial_1U`.  The Green decay and maximum principle remove an
exterior remainder, contradicting normalization.  Since `A_epsilon` is
identity minus compact, `Q A Q` is Fredholm of index zero on the two stated
complements, and (17) gives its onto inverse.

The triangular speed column follows exactly from (13):

    partial_cF_0-(R_epsilon^2/2)partial_muF_0
      =a_epsilon[R_epsilon(r-R_epsilon)+(r-R_epsilon)^2/2].   (19)

## 4. Direct bordered isomorphism

Fix the uncharged base `zeta_0` and the target numbers

    kappa_*=kappa(zeta_0),
    I_*=I_z(zeta_0),
    M_1(zeta)=integral r zeta dnu,
    R_*=M_1(zeta_0)/kappa_*.

In the blow-up notation of this section, `R=R_*`; every centered finite row
uses this frozen target number.

Define the circulation--impulse border

    B_I,epsilon(h,delta_mu,delta_c)
      =(A_epsilon h+a_epsilon delta_mu+a_epsilon r^2 delta_c/2,
        integral h dnu,
        pi rho_m integral r^2h dnu).                         (20I)

It maps `X times R^2` to `X times R^2` and is Fredholm of index zero: the
infinite block is identity minus compact and equally many finite columns and
rows are added.

The operator (20I) is injective for every sufficiently small core.  Suppose a
kernel sequence contradicted a uniform inverse.  Let `H_epsilon` be the
positive threshold-amplitude scale in the Cao blow-up and set

    D_epsilon=epsilon_core^(-2)H_epsilon^p,
    h_tilde(y)=h(R+s y_1,s y_2)/D_epsilon,
    phi_tilde(y)=G_1h(R+s y_1,s y_2)/H_epsilon.

Use the affine source coordinates

    eta_0=delta_mu+R^2 delta_c/2,
    eta_1=R s delta_c.                                       (21)

Normalize in the product norm

    ||h||_Xepsilon/H_epsilon^p
      +|eta_0|/H_epsilon+|eta_1|/H_epsilon=1.                (21a)

Indeed, because `s/epsilon_core` stays bounded in the Cao blow-up, (16) gives

    ||h||_Xepsilon
      asymp H_epsilon^p||h_tilde||_(C^(2,alpha)_y),           (21a')

whereas `D_epsilon` is the raw pointwise density scale.  Dividing the scaled
source norm by `D_epsilon` would lose an extra `epsilon_core^2` and would not
exclude a source-dominated sequence.

The exact first row, divided by `D_epsilon`, is

    h_tilde=p U_+^(p-1)
      (phi_tilde-hat_eta_0-hat_eta_1 y_1)+o(1),              (21b)

where `hat_eta_j=eta_j/H_epsilon`; the remainder is in the fixed local source
norm and includes the quadratic `s^2 delta_c y_1^2/2` term.  Thus a sequence
dominated by the parameter columns is included in the contradiction.

This limit has no unlisted cell.  Put

    w=phi_tilde-hat_eta_0-hat_eta_1 y_1.

Equations (18) and (21b) give

    (-Delta-pU_+^(p-1))w=0.                                 (21c)

The normalized whole-space Cao Green row allows logarithmic `m=0` behavior
and affine plus decaying `m=1` behavior outside the unit positive set.  In
the regular `m=0` class, matching value and normal derivative through
`|y|=1` leaves the one-dimensional scaling solution

    v_a=(2/(p-1))U+y dot grad U.

In the regular `m=1` class it leaves the translations `partial_1U` and
`partial_2U`; evenness in the axial coordinate removes `partial_2U`.
For `m>=2`, the bounded/decaying class is removed by the Lane--Emden kernel
classification in Cao Lemma 3.8.  The `p>=6` positive-part density has the
traces required for this `C^1` interior/exterior matching and introduces no
interface distribution.  Thus `v_a` and `partial_1U` exhaust the allowed
limit.  Applying `-Delta` gives exactly the two density cells below.

The two limiting density cells are explicit.  With

    U_lambda(y)=lambda^(2/(p-1))U(lambda y),
    Lambda_p=integral U_+^p dy=-2 pi U'(1),

set

    q_a=d/dlambda[(U_lambda)_+^p]|_(lambda=1)
       =(2p/(p-1))U_+^p+y dot grad(U_+^p),
    q_t=partial_(y_1)(U_+^p).                                (22)

Integration by parts on the compact positive set gives

    integral q_a dy=2 Lambda_p/(p-1)=:M_a !=0,
    integral q_t dy=0,
    integral y_1 q_t dy=-Lambda_p.                           (23)

Retain the physical Jacobian by using the normalized cells

    Q_a(x)=(R s^2)^(-1)q_a((x-(R,0))/s),
    Q_t(x)=(R s^2)^(-1)q_t((x-(R,0))/s),                     (24)

with the source cutoffs and `o(1)` corrections furnished by (17).  Since
`dnu=(R+s y_1)s^2dy`, first replace `Q_t` by

    Q_t^0=Q_t-[Dkappa(Q_t)/Dkappa(Q_a)]Q_a.                  (25)

This is only a triangular normalization of the same two-cell space and makes
its circulation row exactly zero.  Radial symmetry and (23) then give the
physical moment matrix

    [ Dkappa(Q_a)                         Dkappa(Q_t^0) ]
    [ (D I_z-pi rho_m R_*^2 Dkappa)(Q_a)  (D I_z-pi rho_m R_*^2 Dkappa)(Q_t^0) ]

      = [ M_a+O(s^2)                         0 ]
        [ O(s^2)   -2 pi rho_m R s Lambda_p+O(s^2) ].        (26)

Thus circulation removes the profile cell.  It does not remove radial
translation at leading order.  After circulation is enforced, the centered
physical impulse row removes radial translation; in physical shift
coordinates it is (6).  The determinant in (26) is nonzero for every
sufficiently small fixed core.  Evenness has already removed axial
translation.

It remains to remove the affine source columns, rather than naming an
unspecified matching step.  In the rescaled first row, the affine part is

    eta_0+eta_1 y_1+O(s^2|delta_c||y|^2).

After (26) removes the two density-cell coefficients, (17) and compact Green
convergence give `h_tilde,phi_tilde ->0`.  Projecting (21b) on `1` and `y_1`
gives

    [ A_0  0   ] [hat_eta_0] =o(1)(|hat_eta_0|+|hat_eta_1|),
    [ 0    A_2 ] [hat_eta_1]

    A_0=integral pU_+^(p-1)dy>0,
    A_2=integral pU_+^(p-1)y_1^2dy>0.                        (27)

This is the leading constant/linear affine matching matrix; its determinant
is `A_0 A_2>0`.  It forces `hat_eta_0=hat_eta_1=0`, contradicting the
normalization (21a) after the density part has vanished.  For each fixed core,
the triangular matrix in (21), whose determinant is `R s`, then forces
`delta_mu=delta_c=0`.  Hence (20I) is injective and, by index zero, onto.

The centered impulse expression in (26) is a fixed-base row operation on the
nonlinear target equations `kappa(zeta)=kappa_*` and `I_z(zeta)=I_*`; `R_*`
is not varied inside that derivative.

For the spectral-continuity application, use instead the nonlinear target
rows

    kappa(zeta)=kappa_*,
    M_1(zeta)-R_* kappa(zeta)=0,

and the circulation--mean-radius border

    B_R,epsilon(h,delta_mu,delta_c)
      =(A_epsilon h+a_epsilon delta_mu+a_epsilon r^2 delta_c/2,
        Dkappa(h),(D M_1-R_* Dkappa)(h)).                   (20R)

On the exactly mass-zero cell `Q_t^0` from (25),

    (D M_1-R_* Dkappa)(Q_t^0)
      =-s Lambda_p+O(s^2),                                  (27a)

whereas the same row on `Q_a` is `O(s^2)`.  Replacing the second row of
(26) by (27a), the same product normalization and the same positive affine
matrix (27) prove that (20R) is an index-zero isomorphism.  Thus `B_I` fixes
the nonlinear targets `(kappa_*,I_*)` for the charged branch, while `B_R`
fixes `(kappa_*,kappa_*R_*)` and therefore supplies the exact target-mean-radius
external-core path used below.  The two are distinct bordered realizations;
only their leading radial-cell matrices are triangularly related.

Equation (17) is the projected-complement HSE of P253/0077: `Q A Q` is an
index-zero isomorphism on the declared right and left complements.  The direct
augmented theorems prove branch existence without assuming invertibility of
the full fixed-`(mu,c)` operator `A_epsilon` or a sharp two-by-two Schur
coefficient.

## 5. Sharp two-by-two Schur jet and remaining theorem

Use circulation and mean radius as finite physical coordinates,

    kappa=integral zeta dnu,
    R=kappa^(-1)integral r zeta dnu.                          (28)

The value-level exact Green and local Pohozaev calculations give

    log s_epsilon=-L+O(1),
    A_LE=kappa R/Lambda_p+O(epsilon_core L),
    a_log=kappa R/(2 pi)+O(epsilon_core L),                  (29)

where `A_LE` multiplies `U` and `a_log` is the matched logarithmic
coefficient.  They are different parameters.  The exact radial balance and
Green value are

    c_epsilon R-kappa/(4 pi)
      [log(8R/s_epsilon)+(p-1)/4]=O(epsilon_core^2L),        (30)

    mu_epsilon=a_log L-c_epsilon R^2/2
      +kappa R/(2 pi)[log(8R)-2]+O(epsilon_core^2L^2).       (31)

Hence, at value level,

    c_epsilon=kappa L/(4 pi R)+O(1),
    mu_epsilon=3 kappa R L/(8 pi)+O(1),
    I_z=pi rho_m kappa R^2+O(epsilon_core^2).                (32)

The last remainder is the centered second radial moment; its `C^1` version
must include the derivative of that variance.

If (29)--(32) are proved with the corresponding uniform `C^1` remainders,
then

    det D_(kappa,R)(mu,c)
      =-3 kappa L^2/(16 pi^2 R)+O(L),                        (33)

and, ordering columns `(mu,c)` and rows `(kappa,I_z)`,

    (partial_mu kappa)|_c=4 pi/(3RL)+O(L^-2),
    (partial_c I_z)|_kappa=-8 pi^2 rho_m R^3/L+O(L^-2),     (34)

    S_epsilon=-32 pi^3 rho_m R^2/(3L^2)+O(L^-3).            (35)

This is the exact leading algebra checked by the oracle.  Its finite-core
sign remains conditional on a differentiated-remainder theorem: differentiate
the exact mass and local Pohozaev identities; control derivatives of the
scale, mean center, profile, and centered variance; prove the `p>=6`
interface cancellation; and bound differentiated near/far Green terms.
The distinct LS kernel-to-sharp-Schur bridge needs the exact tangent identity
taking every full `A_epsilon` kernel into the nullspace of (33).  Cao
Proposition 3.13 supplies value estimates, not that bridge or the differentiated
remainder theorem.  This does not reopen the projected-complement HSE proved
in (17).

## 6. The subluminal window and charged branch

At fixed positive circulation and bounded radius, (32) makes `c_epsilon`
diverge with `L`.  A fixed Maxwell speed therefore cannot support a theorem
for every sufficiently small core.  Fix `eta in (0,1)` and positive compact
parameter intervals.  If `C_c` bounds the speed remainder, define

    L_EM=(4 pi R_- / kappa_+)[(1-eta)c_EM-C_c].              (36)

The exact result uses the nonempty window

    L_border<=L<=L_EM,                                      (37)

with the declared foundation-scale inequality `L_border<L_EM`.  Here
`L_border` is the threshold of the bordered proof.  Throughout this window
`c<=(1-eta)c_EM`, and the Maxwell ellipticity margin is uniformly
`1-c^2/c_EM^2>=1-(1-eta)^2>0`.

On (37), the derivative of (11) augmented by (5) is (20).  The `C^1` IFT
therefore gives, for both signs of `g`, an exact subluminal axisymmetric
no-swirl relative equilibrium with

    zeta_g-zeta_0=O(g^2),   c_g-c_0=O(g^2),
    chi_g-chi_0=O(g^2),    E_g,B_g=O(g).                    (38)

The tag integral is fixed while `Q_g=g C_chi`; fixed-charge language applies
to variations at each `g`, not across the branch.  The nonisochronous support
condition from P253/0077 remains in force.

With `epsilon_core` treated as an external parameter, the uncharged bordered
IFT supplies a local continuous exact path at fixed circulation and fixed mean
radius, solving for `mu` and `c`.  The value asymptotic is

    s_epsilon/R=C(kappa,R)epsilon_core(1+o(1)).              (39)

This keeps the limiting column profile fixed.  To turn (39) into P253/0079's
exact geometric-`delta` bracket one still needs either a `C^1` scale remainder
giving a nonzero leading derivative, or a direct uniform two-sided coverage
estimate whose error is `o(h_N)`.  Varying `R` or `c` instead would require
including the induced column-frequency derivative in the transversality and
proving it does not cancel the fixed-column derivative.

The charged continuation adds a finite hierarchy that the uncharged
`N->infinity` crossing does not automatically satisfy.  Since
`delta_N=Theta(N^(-1))` and `c_epsilon=Theta(log N)`, one must exhibit an
integer in the nonempty interval

    max(N_graph,N_IVT,N_response)<N<N_EM,                    (40)

where `N_EM` is the finite ceiling corresponding to (36)--(37).  The lower
thresholds come respectively from graph--Riesz transfer, the endpoint-sign
bracket, and the normalized mixer/leakage bounds.  Neither (37) nor (39) by
itself proves that (40) is nonempty.

Equation (38) proves existence, not joint energy--Casimir coercivity or
all-time tag localization.  It selects neither charge nor action scale and
does not provide analyzer dynamics, event statistics, a full-field Lorentz
cone, electron, or neutrino identity.

This branch belongs to the declared Euler--Maxwell extension and uses the new
`U(1)` field and the explicit coefficients `g`, `epsilon_EM`, and `mu_EM`.
It neither derives those objects from bare incompressible Euler nor adopts
them as a replacement for the Issue 203 Euler objective.  The bare-Euler
carrier and asymptotic-tail routes remain active.  Euler similarity still
leaves the overall action scale continuous; the Maxwell coefficients do not
select a universal nonzero action without a separate compact gauge/topological
matter mechanism or another declared selection law.

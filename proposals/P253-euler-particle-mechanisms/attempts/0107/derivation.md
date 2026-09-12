# Exact fixed-leaf completion and the surviving P2 routes

## 1. Authority, conventions, and corrected invariant row

This body opens only after the central momentum-domain schema replay with exit
`0` on README SHA-256
`59cfe8d5ec379fdded8f34ecdf951ae9a66c6b910603f8dbef1fc1e94f91f345`.
The conserved translation row is

    J_ren=Delta I_h+Delta P_EM,
    I_h(omega)=(rho_m/2) integral x cross omega dx,
    P_EM=epsilon_EM integral E cross B dx.                (1)

Hydrodynamic impulse is retained below only as an initial carrier/modulation
slice.  It is not a second conserved quantity of the coupled dynamics.

Use

    qhat(k)=integral exp(-i k dot x)q(x)dx,
    uhat=i k cross qhat/|k|^2,
    C_g xi=curl(xi cross omega_g)=-[xi,omega_g].           (2)

The last equality uses `[a,b]=a dot grad b-b dot grad a` and divergence-free
fields.

## 2. Exact forced-impulse and Maxwell-momentum balance

Let `f=rho_q E+j cross B` be the Lorentz-force density.  The forced vorticity
equation is

    partial_t omega=curl(u cross omega)+(1/rho_m)curl f.  (3)

For a decaying vector field `A`, integration first on a ball gives

    (1/2) integral x cross curl A dx=integral A dx.       (4)

The surface term vanishes in the localized class.  Also

    u cross omega=grad(|u|^2/2)-(u dot grad)u,
    integral u cross omega dx=0.                          (5)

Consequently

    d I_h/dt=integral f dx.                               (6)

For the Maxwell convention

    div(epsilon_EM E)=rho_q,      div B=0,
    epsilon_EM partial_t E=(1/mu_EM)curl B-j,
    partial_t B=-curl E,                                  (7)

define

    sigma_EM=epsilon_EM(E tensor E-|E|^2 I/2)
             +(1/mu_EM)(B tensor B-|B|^2 I/2).

Substitution of (7), with no time averaging, gives

    partial_t(epsilon_EM E cross B)=div sigma_EM-f.       (8)

Thus `-sigma_EM` is the field-momentum flux and

    d P_EM/dt=-integral f dx.                             (9)

For the difference from the traveling carrier, the initial charge difference
has zero monopole, its Coulomb row is `O(|x|^-3)`, and `E_g=O(|x|^-2)`.
The compact transverse correction is initially supported in a fixed ball.
On a fixed finite interval `[0,T]`, finite Maxwell propagation confines newly
generated transverse radiation from the compact correction/current to a
finite expanding wave zone; outside its backward cone the initial dipole tail
retains its `O(|x|^-3)` order.  A radiative `O(|x|^-1)` profile inside that
finite wave zone therefore never reaches the spheres `R->infinity` at fixed
`t`.

There is also a flux proof that does not require pointwise asymptotics.  Use
smooth cutoffs `vartheta_R` with `|grad vartheta_R|<=C/R`.  Finite field energy
makes the Maxwell stress `L1`; hence

    |integral sigma_EM grad vartheta_R dx|
       <=(C/R)||sigma_EM||_L1 ->0.                       (9a)

The same estimate applies after subtracting the carrier, including every
base--difference cross term by `L2 times L2 -> L1`.  The charge/current and
Lorentz force are supported with the compact transported tag on `[0,T]`.
The smooth finite-energy local Euler--Maxwell solution used by the reviewed
same-target bridge preserves these properties for each fixed `T`; no
time-uniform scattering estimate is claimed.  Equations (6)--(9) therefore
hold on every fixed finite time interval by the cutoff limit, and

    d J_ren/dt=0.                                        (10)

This is also the translation Noether balance of the material-plus-Maxwell
system.  It explains why `I_h` and `J_ren` cannot both be called independent
conserved rows at nonzero charge.

## 3. The velocity-integral domain diagnostic

Let compact smooth `q` be divergence-free.  Then `integral q=0`, since
`q_i=partial_j(x_i q_j)`.  Put

    M_(ji)=integral x_j q_i dx.

Expanding `qhat` at zero and using `k dot qhat(k)=0` gives

    qhat_i=-i k_j M_(ji)+O(|k|^2),
    M_(ji)+M_(ij)=0.                                     (11)

The adopted impulse normalization gives

    M_(mn)=rho_m^(-1) epsilon_(lmn) I_h,l(q),
    qhat=i rho_m^(-1) k cross I_h(q)+O(|k|^2).            (12)

If `I_h(q)=0`, the exact Biot--Savart multipole expansion has its monopole and
first moment removed, hence

    u(x)=O(|x|^-4),     u in L1(R3).                      (13)

Conversely, if a divergence-free `u` belongs to `L1`, its Fourier transform is
continuous.  From `k dot uhat(k)=0`, set `k=t e` and send `t->0`; every unit
vector `e` obeys `e dot uhat(0)=0`, so

    integral u dx=uhat(0)=0.                             (14)

Equations (13)--(14) show that `rho_m integral(u-u_g)` is identically zero on
the zero-impulse `L1` slice.  Equations (6) and (10) show why that slice need
not be invariant under Lorentz exchange.  It is not the coupled momentum row.

## 4. A quantitative compact witness lemma

Let `F` be a nonzero, divergence-free, localized `L2` vector field.  For a
constant vector `y`, define

    S_F y=P_L(F cross y),
    Q_F,ij=<S_F e_j,S_F e_i>_L2.                          (15)

If `S_F y=0`, then `F cross y` is a gradient and

    0=curl(F cross y)=(y dot grad)F.                      (16)

For `y!=0`, the Fourier transform of (16) is supported on the measure-zero
plane `{k:k dot y=0}`.  Since `Fhat in L2`, this forces `F=0`, a contradiction.
Thus `Q_F` is positive definite.  Write

    lambda_F=lambda_min(Q_F)>0,
    s_F=||S_F||_(R3->L2).                                (17)

This injectivity proof uses the actual localized field and excludes a constant
translation stabilizer; it is not a generic-density assertion.

The ideal columns in (15) need not be compact.  Set

    A_F,j=(-Delta)^(-1)curl(F cross e_j),
    curl A_F,j=S_F e_j.                                  (18)

Choose radial `chi_R` equal to one on `|x|<=R`, zero on `|x|>=2R`, and with
`|grad chi_R|<=C/R`.  Then

    w_F,j,R=curl(chi_R A_F,j)                            (19)

is compact and divergence-free, but for a merely `L2` field `F` it is not yet
smooth.  The exact cutoff error is

    w_F,j,R-S_F e_j
       =(chi_R-1)S_F e_j+grad chi_R cross A_F,j.          (20)

Here `A_F,j in dot H1`; Hardy's inequality gives
`||A_F,j/|x|||_2<=2||grad A_F,j||_2`.  Therefore

    eta_F(R):=||[w_F,1,R-S_Fe_1,...,w_F,3,R-S_Fe_3]||
      <=(sum_j ||1_(|x|>R)S_Fe_j||_2^2)^(1/2)
        +C(sum_j||1_(R<|x|<2R)A_F,j/|x|||_2^2)^(1/2)
      ->0.                                                (21)

Choose first `R_F` so that the right-hand side of (21) is at most
`lambda_F/(4s_F)`.  Let `varphi_epsilon` be a compactly supported standard
mollifier and define

    w_F,j,R,epsilon=varphi_epsilon * w_F,j,R.            (21a)

Convolution preserves divergence-free character and compact support, and
places every column in `C_c^infinity`, hence arbitrarily more than the two
derivatives above the frozen `H^s`, `s>=4`, required for the packet path.  For
this fixed `R_F`, choose `epsilon_F>0` so that

    kappa_F(R_F,epsilon_F):=
      ||[w_F,j,R_F,epsilon_F-w_F,j,R_F]_(j=1)^3||
      <=lambda_F/(4s_F).                                 (21b)

This is the second, independent regularity/error budget.  The total column
error is at most `lambda_F/(2s_F)`.

For the response `T_(c,F)w=c integral w cross F dx`, the matrix formed from
the compact smooth columns (21a) is

    M_F(R,epsilon)_(ij)=e_i dot T_(c,F)w_F,j,R,epsilon
       =c <w_F,j,R,epsilon,S_F e_i>.                     (22)

Consequently

    ||M_F(R,epsilon)-c Q_F||_op
       <=|c|s_F[eta_F(R)+kappa_F(R,epsilon)].            (23)

Choose the finite pair `(R_F,epsilon_F)` satisfying

    eta_F(R_F)<=lambda_F/(4s_F),
    kappa_F(R_F,epsilon_F)<=lambda_F/(4s_F).             (24)

Then

    sigma_min(M_F(R_F,epsilon_F))>=|c|lambda_F/2,
    ||M_F(R_F,epsilon_F)^(-1)||<=2/(|c|lambda_F).        (25)

This supplies an explicit cutoff-error threshold.  It also makes transparent
that a fixed-member bound can degenerate when either the physical coefficient
or the field tends to zero.

## 5. The electromagnetic matrix G

For the fixed nonzero charged Cao member, `B_g` is a nonzero divergence-free
finite-energy field.  Here is the required transverse-current argument.  The
reviewed axisymmetric solution has `B_g=rH_g e_theta` and the gauge-fixed
traveling vector-potential equation has transverse source
`P_L(g chi_g u_g)`.  If this source vanished, the current would be a gradient.
On any closed strict-band streamline with `chi_g(P)!=0`, however,

    integral_orbit (g chi_g u_g) dot dx
      =g chi_g(P) integral_0^T |W_g|^2 dt !=0,           (25a)

because the translation contribution integrates to
`c_g integral_orbit dz=0`.  A gradient has zero closed circulation, so the
transverse current is nonzero.  Subluminal invertibility of the reviewed
vector-potential block then gives a nonzero transverse potential and
`B_g=curl A_g!=0`; equivalently this is the nontrivial `H_g` row in the
reviewed toroidal Ampere primitive.  Apply Section 4 with

    F=B_g,       c=epsilon_EM,
    e_T,j=w_B,j,R_B,epsilon_B.                           (26)

The fields `e_T,j` are compact and divergence-free.  Thus adding them to the
electric difference changes neither Gauss nor magnetic divergence, and with
`b=0` changes field momentum exactly linearly.  Its physical response is

    G_(ij)=epsilon_EM integral(e_T,j cross B_g)_i dx,
    sigma_min(G)>=epsilon_EM lambda_B/2,
    ||G^(-1)||<=2/(epsilon_EM lambda_B).                  (27)

Because `B_g=O(g)`, `lambda_B=O(g^2)` under this normalization; no uniform
`g->0` inverse is claimed.  At `g=0` the electromagnetic correction is the
separate automatic-zero row.

## 6. The hydrodynamic matrix H and the center slice C

For a displacement `xi`, (2) and (4) give

    D I_h[C_g xi]=rho_m integral xi cross omega_g dx.    (28)

Apply Section 4 with

    F=omega_g,       c=rho_m,
    xi_I,j=w_omega,j,R_omega,epsilon_omega.              (29)

Then

    H_(ij)=rho_m integral(xi_I,j cross omega_g)_i dx,
    sigma_min(H)>=rho_m lambda_omega/2,
    ||H^(-1)||<=2/(rho_m lambda_omega).                  (30)

For the center slice use the compact transported tag, with mass
`Q_chi=integral chi_g dx!=0`, and

    X_chi=Q_chi^(-1)integral x chi dx.                   (31)

Choose `psi=1` on a neighborhood of
`supp omega_g union supp chi_g`, with compact support, and define

    xi_C,j=curl[psi(x)(e_j cross x)/2].                  (32)

Since `curl[(e_j cross x)/2]=e_j`, this field is compact, divergence-free, and
equals `e_j` on both supports.  The transported-density variation
`delta chi=-div(chi_g xi)` gives

    D X_chi[xi_C,j]=e_j,                                 (33)

so `C=I_3`.  Also `integral omega_g=0` for compact divergence-free vorticity,
and hence

    D I_h[C_g xi_C,j]
       =rho_m e_j cross integral omega_g dx=0.           (34)

If `c_(aj)=D X_chi,a[xi_I,j]`, replace

    xi_I,j by xi_I,j-sum_a c_(aj)xi_C,a.                 (35)

This makes the impulse columns center-orthogonal without changing `H`.

The complete frozen row inventory inherited from P253/0095/0104 is as
follows; it also fixes which rows are automatic and which enter the finite
correction.

1. `div omega=0`, `div B=0`, the difference Gauss law, and the zero difference
   charge monopole are automatic exact rows of the curl, Maxwell, and
   Gauss-reconstruction formulas.
2. The coadjoint/Kelvin-circulation data and the complete transported-tag
   distribution are automatic exact rows of the common material push-forward.
3. Scalar axisymmetric branch rows, including circulation and mean-radius
   profile pairings used to select the Cao member, carry toroidal character
   zero.
4. The three components of initial hydrodynamic impulse are auxiliary
   carrier/tail-slice rows, not conserved coupled rows.
5. The three components of the tag barycenter (31) are the translation slice.
6. The three components of `J_ren` are the conserved translation row.  Its
   full linear field part is displayed in (39a).
7. Electromagnetic gauge is absent in the physical `(E,B)` variables.  Axial
   rotation/phase is a stabilizer and quotient parameter, not an additional
   scalar conservation row.

The reviewed fixed circulation--impulse bordered uniqueness selects the
equilibrium member; it does not turn hydrodynamic impulse into a conserved
coupled row.  The chosen tag barycenter replaces any general theta-dependent
center covector.  If a later modulation convention deliberately uses a
non-`SO(2)`-equivariant center, phase, or carrier covector, its packet row is
only `o_N(1)` by nonstationary phase and must be removed through (38); it is
not silently called exact.  No anonymous finite row remains after fixing
`J_ren` and the auxiliary rows (4)--(5).

## 7. Lower-triangular joint-momentum derivative

Under a displacement `xi`, the tag and its Gauss-reconstructed electric field
vary by

    delta chi_xi=-xi dot grad chi_g,
    delta E_L,xi=-(g/epsilon_EM)grad(-Delta)^(-1)delta chi_xi
       =(g/epsilon_EM)grad(-Delta)^(-1)(xi dot grad chi_g). (36)

With `delta B=0`, the exact displacement derivative of (1) is

    K(xi)=rho_m integral xi cross omega_g dx
          +g integral [grad(-Delta)^(-1)(xi dot grad chi_g)]
                    cross B_g dx.                        (37)

Let `K_I` and `K_C` be (37) on the orthogonalized columns (35) and (32).
With rows `(Delta I_h,X_chi,J_ren)` and columns `(xi_I,xi_C,e_T)`, the actual
matrix is

    M_leaf=[[H,0,0],
            [0,I,0],
            [K_I,K_C,G]],
    det M_leaf=det H det G!=0.                           (38)

The lower-left blocks are retained; no block-diagonal claim is made.

## 8. Exact packet rows and nonlinear fixed-J_ren curve

The reviewed growing packets use the exact phase

    exp[i N(S(I)+ell theta)],       n=N ell.              (39)

Choose `|n|>=2` and theta-independent tube cutoffs.  Inventory items 1--2 are
exact nonlinear constraints.  The scalar axisymmetric branch rows in item 3
see only character zero.  The vector rows in items 4--6 are physical
`SO(2)`-covariant rows.  In Cartesian components an axisymmetric covariant
vector integrand has only the cylindrical packet character shifted by the
frame characters `0,+1,-1`.  In particular, the full linear field part of
`J_ren` is

    epsilon_EM integral(e_N cross B_g+E_g cross b_N)dx, (39a)

while the quadratic `e_N cross b_N` is absent at the linear tangent.  The
Gauss-slaved longitudinal `e_N` has the same covariant toroidal character as
the tag tangent; multiplying by the axisymmetric `E_g,B_g` produces no zero
Cartesian vector component for `|n|>=2`.  Therefore Fourier orthogonality
gives, packet by packet,

    Delta I_h'(q_N)=0,
    D X_chi(q_N)=0,
    D J_ren(q_N,eta_N,e_N,b_N)=0                         (40)

exactly, including the Cartesian frame shifts.  The longitudinal Gauss field
has the same nonzero toroidal character, so it does not reintroduce a vector
zero mode.  Thus the frozen finite-rank tangent projection is the identity on
this actual high-character sequence.  If a different modulation convention
adds a smooth nonaxisymmetric covector, its row is not called exactly zero:
direct nonstationary-phase integration makes it `o_N(1)`, and (38) removes it
with a vanishing smooth correction.  The corrected sequence remains weak null
and its fixed-time observed correction is `o_N(1)`.  No continuity of total
velocity momentum on `X_in^s` is used.

Let `Phi_tau` be the volume-preserving flow of the smooth packet generator.
Compose it with the six compact flows (35),(32), and transport

    omega_tau=(Phi_tau)_*omega_g,
    chi_tau=chi_g composed Phi_tau^(-1).                 (41)

The first-stage map `(Delta I_h,X_chi)` has derivative `diag(H,I)`.  By (40),
its IFT coefficients are `O(tau^2)` and it sets the auxiliary initial impulse
and center rows exactly.  Define

    E_tau=E_g-(g/epsilon_EM)grad(-Delta)^(-1)(chi_tau-chi_g)
              +sum_j gamma_j(tau)e_T,j,
    B_tau=B_g.                                           (42)

This satisfies Gauss, fixed charge, and magnetic divergence exactly.  After
the first stage, solve the conserved row explicitly:

    gamma(tau)=-G^(-1)J_ren(omega_tau,chi_tau,E_tau|_(gamma=0),B_g).
                                                               (43)

Equation (40) gives `gamma=O(tau^2)`.  Thus (41)--(43) are an exact curve in
the invariant fixed-`J_ren` leaf, with the declared initial carrier/center
slice, and with precisely the original unstable tangent.

The fixed-member solution of the coupled equations conserves (43) by (10).
The auxiliary impulse and center normalizations are allowed to evolve and are
handled by modulation; they are not falsely promoted to conserved rows.

## 9. Fixed-observation nonlinear contradiction

Retain the quantifier order

    fixed circuit j -> fixed observation tube -> sufficiently large N
      -> construct (41)--(43) -> tau->0 -> N->infinity. (44)

The reviewed same-`H^s` pathwise differentiability applies because the packet
generator is chosen two derivatives smoother and (41)--(43) have the same
linear tangent.  The `O(tau^2)` finite rows vanish after division by `tau`.
The modulation derivative removes only finitely many symmetry characters;
the sequence `|n|->infinity` is weak null and is detected by the fixed local
observation.

If the frozen all-time Lipschitz bound held with constant `C`, differentiation
along the **same normalized corrected packet curve** (41)--(43), for fixed
`j,N`, would give only the packetwise estimate

    ||J_loc S_g(jT_*) Q_(j,N)||_(X_DA,loc^s/sym)
       <=C ||Q_(j,N)||_(X_in^s).                          (45)

No full operator-norm estimate is inferred.  The reviewed Unit A--F
construction supplies, for each fixed `j`, a normalized weak-null sequence
with

    liminf_(N->infinity)
      ||J_loc S_g(jT_*) Q_(j,N)||_(X_DA,loc^s/sym)
       >=c_obs lambda_+^j,
    ||Q_(j,N)||_(X_in^s)=1,       lambda_+>1.             (46)

Choose `j` with `c_obs lambda_+^j>C`, differentiate first at each fixed `N`,
and only then send `N->infinity`.  Equations (45)--(46) contradict each other
on that sequence.  Route A therefore establishes the fixed-member nonlinear
uniform-Lipschitz obstruction after consuming the reviewed same-target smooth
pathwise solution map; Section 2 supplies its finite-time momentum-flux
closure.  It is not a full semigroup-norm claim, a one-datum nonlinear escape
theorem, a verdict on larger charge, or a refutation of P2 for another carrier.

## 10. Route B: whole-core charged-column classification

The exact slow block retained from P253/0095 is

    C_col=x[[0,-Z_tau],[R_tau,0]],
    K(s,tau)=4 Z_tau R_tau/Omega_tau^2.                  (47)

At `tau=0`,

    Z_0=zeta=2Omega+s Omega_s,      R_0=2Omega,
    K(s,0)=8 zeta/Omega.                                  (48)

At a smooth positive center `s=0`, regularity gives
`zeta(0)=2Omega(0)`, and therefore

    K(0,0)=16.                                           (49)

This is the first Hattori--Fukumoto resonance value.  The reviewed invariant
first-curvature coefficient satisfies

    |b(0)|=15|U_theta|/256>0.                             (49a)

On a fixed regular center annulus, continuity of the charged column/profile
and returned coefficient gives `|b(tau)|>=|b(0)|/2` for sufficiently small
`|tau|`.  Hence `K>=1` there opens an actual hyperbolic tongue, not merely a
kinematic frequency equality.  Charge must make an order-one change in (47),
not a perturbative `O(tau)` change, before the center can enter a uniform
`K<1` window.

At the compact Cao support edge, `zeta->0` with the `p>=6` free-boundary jet,
and the strict-band tag and its Lorentz correction vanish no later than that
edge.  Enclosed circulation generally leaves `Omega_edge!=0`.  Hence

    Z_tau->0,       R_tau->2Omega_edge,
    K(s,tau)->0                                      (50)

for the direct strict-band charged construction.  No positive lower bound
`kappa_0<=K` can cover the closed core through its support edge.  The interior
regular-band classifier remains useful, but the edge requires its own
weighted-accessibility theorem.

For `R_tau,Z_tau>0`, the order-zero column symmetrizer

    H_0=diag(R_tau,Z_tau)                                (51)

satisfies

    H_0 C_col=x[[0,-R_tau Z_tau],[R_tau Z_tau,0]].        (51a)

Thus if `R_tau,Z_tau` have common positive upper and lower bounds, the frozen
column evolution is uniformly symmetrizable even through `x=0`.  The physical
KKS Hessian, however, is proportional to

    |x| H_0                                               (51b)

up to the already fixed positive physical factors.  It degenerates as the
accessible wave angle `x->0`.  This blocks direct KKS-Hessian coercivity in the
frozen physical topology but, by (51a), does **not** itself imply linear growth
or loss of order-zero symmetrizability.  The open question is whether the
global action-flow graph norm realizes `H_0`, or whether only the singular KKS
form descends on the physical DA closure.

The compact edge is sharper.  For `R,Z>0`, put
`nu=|x|sqrt(RZ)`.  Exactly

    exp(t C_col)=cos(nu t)I+[sin(nu t)/nu]C_col.          (51c)

At `t_s=pi/(2|x|sqrt(RZ))`, its lower off-diagonal entry has magnitude
`sqrt(R/Z)`.  Hence any frozen-column accessible sequence approaching the edge
with `R->R_edge>0`, `Z->0` has no uniform-in-core all-time bound in the
unweighted velocity-amplitude norm.  At `Z=0` exactly,

    exp(t C_col)=I+t x[[0,0],[R,0]],                     (51d)

the Jordan shear.  Transferring (51c)--(51d) to the finite Cao DA graph still
requires an accessible localized edge packet and control of the
interface/free-boundary operator.  The accessibility weight might remove or
renormalize this channel; that is the active B2 construction.

The missing charged constitutive row is not free notation.  In the local
column scaling it must be the reviewed modified Grad--Shafranov law

    zeta_tau=h(P)+tau/rho_m
       [chi_P'(P) Phi_hat-chi_P(P) H_hat],                (52)

together with `P_s=R_core s Omega` in the fixed physical convention and the
radial Maxwell equations already recorded in P253/0095.  Equations (47)--(52)
give Route B a precise split: its perturbative center has an actual hyperbolic
first tongue, the frozen-column edge has the exact conditioning mechanism
(51c), and the weighted whole-closure/interface transfer remains active.

## 11. Route C: Hill remains a full bulk/interface problem

Write Hill's interior vorticity and its zero extension in the translating
frame as

    omega_in=C_H(e_z cross x),
    omega_H=omega_in 1_B,                                (53)

where `B` is the spherical core.  For a smooth compact divergence-free
displacement `xi`, the exact DA tangent is

    delta omega=-[xi,omega_in] 1_B
                 +(xi dot n)omega_in delta_(partial B).  (54)

Indeed

    delta omega=curl(xi cross omega_in 1_B),
    delta u=P_L(xi cross omega_in 1_B).                  (55)

The second identity is the whole-space Hodge reconstruction; it contains the
interior, interface, and exterior harmonic velocity in one operator and does
not impose an artificial boundary condition at `partial B`.

Let

    F_H(omega)=rho_m/2 integral omega dot B_R3 omega dx
                -c_H I_z(omega),
    mu_H=rho_m B_R3 omega_H
          -(rho_m c_H/2)(e_z cross x).                   (56)

Differentiating two material push-forwards gives the exact weak relative-orbit
Hessian

    Q_H(xi,eta)
      =rho_m integral P_L(xi cross omega_in 1_B) dot
                       P_L(eta cross omega_in 1_B) dx
       +(1/2)<mu_H,
          [xi,[eta,omega_H 1_B]]
          +[eta,[xi,omega_H 1_B]]>.                      (57)

Equation (57), rather than a boundary-only VSH energy, is the frozen Route-C
object.  For the diagonal, distributional integration by parts gives a useful
single-core representation.  Define

    D_xi mu_H=(xi dot grad)mu_H+(grad xi)^T mu_H.

Then

    Q_H(xi,xi)
      =rho_m ||P_L(xi cross omega_in 1_B)||_2^2
       +integral_B (xi cross omega_in) dot
          curl(D_xi mu_H) dx.                            (58)

The interface distribution in (54) has not disappeared: it is encoded by the
zero-extension in the first term and by the trace of the global `mu_H` in the
integration leading to (58).  Formula (58) is valid first for smooth compact
`xi`; its closure norm must include both the interior displacement and the
surface trace.  Exact translations are null directions of (57) and must be
removed before any sign statement.

The surface term alone is not a closed phase space.  A complete constrained
second variation must use (57)--(58), including the interior stretching term,
the interface trace, the exterior harmonic velocity, and the translation
quotient.  Choi's rearrangement controls the axisymmetric no-swirl scalar class
only and does not decide the sign of (58) on arbitrary nonaxisymmetric interior
displacements.

Route C therefore remains **active** at an improved analytic rung: the full
weak quadratic form is now explicit, and the next calculation is its coupled
vector-spherical-harmonic radial matrix after removal of the `l=1` translation
block.  A materially different candidate is a smooth swirling ring, but it
must first pass the full BAS filter: Lifschitz--Hameiri derive localized
exponential/algebraic criteria for swirl, Lifschitz--Suters--Beale report such
instabilities on exact computed rings, and Turkington supplies variational
existence rather than full-three-dimensional coercivity.  Swirl is therefore
not assumed to cure Hill's interface problem.

## 12. Route verdicts

- **Route A: active with its constructive core established.**  The finite
  translation moment map, compact smooth H/C/G witnesses, lower-triangular row
  inverse, and exact fixed-`J_ren` initial curve are derived.  The route-level
  nonlinear obstruction remains prospective until the finite-time weighted
  propagation hypothesis and the packetwise same-target solution-map bridge
  are pinned together at the final evidence boundary.
- **Route B1: refuted for perturbative whole-core standard-norm coercivity.**
  The mechanisms are the nonzero center tongue (49a) and the frozen-column
  edge amplification (51c)--(51d).
- **Route B2: active.**  The exact missing construction is the finite-Cao
  accessible edge-packet/interface transfer and a decision whether the global
  action-flow graph realizes the regular symmetrizer `H_0` through `x=0` and
  `Z=0`.
- **Route C: active.**  Equations (57)--(58) construct the full weak Hessian;
  its nonaxisymmetric VSH radial sign after the translation quotient remains
  to be decided.

No route yet establishes P2, a stable particle, a quantum bridge, or a parent
campaign verdict.  The prospective Route-A obstruction, the active weighted
charged-column construction, and the materially different Hill/swirl carrier
tests remain separately scoped.

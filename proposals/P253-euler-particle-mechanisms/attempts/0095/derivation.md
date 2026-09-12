# Complete-state P2 competition: exact principal systems and carrier tests

## 1. Charged Cao linearization on the affine fixed-charge leaf

Use the reviewed 0068 bracket and the translating-frame convention whose
free Maxwell roots are those in 0085.  Put

    W=u_g-c_g e_z,                 F_g=E_g+u_g cross B_g,
    rho_q=g chi_g,                 c_EM^2=(epsilon_EM mu_EM)^(-1).

For a tangent `(v,eta,e,b)` with `div v=div b=0` and
`div(epsilon_EM e)=g eta`, the exact linearized equations are

    v_t=-P_L[(W dot grad)v+(v dot grad)u_g]
        +(g/rho_m)P_L[eta F_g
          +chi_g(e+v cross B_g+u_g cross b)]+R_slice,       (8)

    eta_t=-W dot grad eta-v dot grad chi_g+R_tag,slice,     (9)

    e_t=-c_g partial_z e+c_EM^2 curl b
        -(g/epsilon_EM)(chi_g v+eta u_g)+R_EM,slice,       (10)

    b_t=-c_g partial_z b-curl e+R_gauge,slice.             (11)

On the smooth compact packet core the declared linear impulse, center, gauge
and carrier-chart functionals are defined, have axisymmetric coefficients,
and vanish on a nonzero toroidal character.  Their formal finite-rank
microlocal rows are included in `R_slice` and do not change the principal
symbol.  No continuous extension or submersion on the completed
`mathcal X_in^s` is asserted.  The coadjoint/Kelvin-circulation data, tag
distribution, charge and Gauss rows are built into the constrained tangent.
The physical total-momentum functional is **not** included on the topology
used here: a generic compact-vorticity perturbation has a whole-space Hodge
velocity tail `O(|x|^-3)`, which belongs to the declared weighted `L^2` rows
but need not be in `L^1`, so `rho_m integral u` is not a continuous
functional there.  Equations (8)--(11), rather than the sign of the relative
Hessian, are the starting point for stability; this finite-row scope does not
alter their principal symbol.

For a phase `exp(i N phi)`, write `xi=grad phi` and restrict the velocity and
Maxwell amplitudes to `xi`-transverse components.  Gauss gives
`xi dot e=O(N^(-1))`.  The homogeneous degree-one symbol is block diagonal:

    p_1(x,xi)=diag(-i W dot xi I_2, -i W dot xi,
                   M_EM(xi)),                              (12)

where the four-dimensional transverse Maxwell block has roots

    lambda_EM,+/-=-i c_g xi_z +/- i c_EM |xi|             (13)

with the two physical polarizations.  Thus the characteristic polynomial,
up to a nonzero convention factor, is

    (lambda+i W dot xi)^3
    [(lambda+i c_g xi_z)^2+c_EM^2|xi|^2]^2.              (14)

Every current, Lorentz, material-tag gradient and moving-slice coupling is
degree zero.  This does not mean it can be deleted from the amplitude
cocycle.  It means that away from a coincidence

    -W dot xi=-c_g xi_z +/- c_EM|xi|,                     (15)

the Maxwell amplitude forced by a fluid WKB packet is one order smaller.
At (15) the full mode-conversion block must instead be retained.

## 2. The Stokes-streamfunction-normal ray is exactly neutral but not a
## leading DA packet

For the charged axisymmetric no-swirl carrier, reviewed 0077 gives

    W=r^(-1) grad P cross e_theta,
    B_g=r H e_theta,
    F_g=-grad Phi-H grad P,
    chi_g=chi(P).                                         (16)

Here `P` is the Stokes streamfunction, not physical pressure.  Set

    n=grad P/|grad P|,  t=W/U,  U=|W|=|grad P|/r,
    xi=grad P.                                            (17)

Then `W dot xi=0`.  Strict subluminality implies that (15) has no nonzero
solution on this ray, because

    |-c_g xi_z|<c_EM|xi|.                                 (18)

The same-leaf tag amplitude is zero at principal order: a solenoidal
displacement amplitude is tangent to the `P` level and hence
`eta=-d dot grad chi=0`.

Let `L=grad W` in the physical cylindrical orthonormal frame.  The ambient
Euler BAS equation, with the whole-space Leray pressure symbol retained, is

    dot a=-L a+2 n(n dot L a),       n dot a=0.            (19)

Write `a=A t+C e_theta`.  Since the carrier has no swirl, `L` preserves the
meridional plane and

    L e_theta=(W_r/r)e_theta,       dot theta=0.           (20)

There is therefore no cylindrical-frame rotation of `e_theta` and no
off-diagonal entry between `t` and `e_theta`.  Projecting (19) gives

    dot A=-(t dot L t)A=-(dot U/U)A,
    dot C=-(W_r/r)C=-(dot r/r)C.                           (21)

The possible speed integral is displayed rather than hidden:

    integral_0^T t dot L t dt
      =1/2 integral_0^T W dot grad log|W|^2 dt
      =log[U(T)/U(0)].                                    (22)

On a closed regular `P` contour, `U(T)=U(0)` and `r(T)=r(0)`.  In the
returned physical frame the exact two-polarization matrix is

    M_P(T)=diag(U(0)/U(T),r(0)/r(T))=I_2.                 (23)

No Bernoulli assertion about `|W|^2` was used; it may vary along the contour.
The off-diagonal shear is zero pointwise by (20), not inferred from the two
diagonal periods.

The base magnetic term also vanishes separately on both polarizations.  With
`t=n cross e_theta`,

    P_n[t cross(B_theta e_theta)]=P_n[-B_theta n]=0,
    P_n[e_theta cross(B_theta e_theta)]=0.                (24)

The `eta F_g` term vanishes because `eta=0`, while the nonresonant perturbed
Maxwell field is `O(N^(-1))`.  Hence (23) remains the degree-zero coupled
return.

This neutral matrix is not yet a physical DA stability sector.  The Cao
vorticity is `bar_omega=r zeta e_theta`, so

    xi dot bar_omega=0.

For a principal solenoidal displacement `d`,

    delta omega_pr=i[(xi dot bar_omega)d-(xi dot d)bar_omega]=0.  (25)

Thus the leading coadjoint symbol degenerates on the Stokes-streamfunction-
normal ray.
Equation (23) is an exact ambient BAS comparator; a DA realization would need
a next-order construction and cannot be used as an LP2 theorem.

## 3. The actual toroidal DA cocycle

The nondegenerate local DA ray is the integer toroidal phase

    xi=(n/r)e_theta,       n!=0.                           (26)

It is again separated from the Maxwell characteristics, now by
`|lambda_EM|=c_EM|n|/r`, since both `xi_z` and `W dot xi` vanish on this
special ray.  Strict subluminality alone does not exclude (15) for an
arbitrary covector; that would additionally require `sup|u_g|<c_EM` or a
direct exclusion of the exact characteristic equation.  Let `a=(a_r,a_z)` be the meridional velocity
amplitude and let `d` be the displacement amplitude.  Comparing
`curl(a exp(i n theta))` with `-[d exp(i n theta),bar_omega]` gives

    d=(e_theta cross a)/(r zeta),
    eta=-d dot grad chi=-(chi'(P)/zeta)W dot a.           (27)

This is the same-leaf constraint; `eta` is not an independent triangular
tag variable.  On the meridional plane let

    L_m=grad_m W,       J a=a cross e_theta.

Because the Leray symbol for (26) is the identity on meridional vectors, the
full nonresonant degree-zero coupled DA cocycle is

    dot a=C_g(t)a,                                        (28)

    C_g=-L_m+(g chi B_theta/rho_m)J
        -(g chi'/(rho_m zeta)) F_g tensor W.              (29)

All functions in (29) are evaluated on the same closed material contour.
The two electromagnetic coefficients are `O(g^2)` on the reviewed branch,
because `E_g,B_g=O(g)`.  The determinant is not guessed from the uncharged
problem.  Directly,

    tr C_g=W_r/r-(g chi'/(rho_m zeta))W dot F_g
          =d(log r)/dt+(g chi'/(rho_m zeta))d Phi/dt.     (30)

The last equality uses `W dot grad P=0` and (16).  Its period integral has
zero `O(g^2)` part: at `g=0`, `zeta=h(P)` is constant on the orbit and
`integral_0^T dot Phi dt=0`.  Thus

    det M_theta(T)=1+O(g^4),                              (31)

uniformly on a fixed regular contour.  An exact all-`g` unit determinant
would require the full KKS density and an invariant two-dimensional reduction
in these velocity coordinates and is not asserted by (31).  Consequently,
small charge preserves a real hyperbolic pair, and in particular an expanding
multiplier, when the gap estimate below holds; (31) does not establish exact
reciprocity for the charged pair.

Equation (29) is the first actual charged-Cao BAS classifier.  Its uncharged
thin-column limit can be evaluated completely.  Let

    W_col=s Omega(s)e_alpha,
    zeta_col(s)=2 Omega+s Omega'>0                        (32)

on a regular positive-vorticity annulus.  In the frame rotating with the
particle at angular rate `Omega`, a returned covector with zero azimuthal
component has `k=(k_s,0,k_z)` and `K=|k|`.  In the transverse basis

    e_1=e_alpha,       e_2=(k_z e_s-k_s e_z)/K,

the exact full-pressure matrix is constant:

    C_col=[[0,-zeta_col k_z/K],[2 Omega k_z/K,0]],
    C_col^2=-nu_K^2 I,
    nu_K^2=2 Omega zeta_col k_z^2/K^2.                    (33)

Therefore

    M_col(T)=
      [[cos(nu_K T),-(zeta_col k_z/(K nu_K))sin(nu_K T)],
       [(2 Omega k_z/(K nu_K))sin(nu_K T),cos(nu_K T)]],  (34)

with the continuous value `I_2` at `k_z=0`.  The positive metric

    G_col=diag(2 Omega |k_z|/K,zeta_col |k_z|/K)          (35)

makes `G_col C_col` skew when `k_z!=0`, under the explicit regular-cell signs
`Omega>0` and `zeta_col>0`.  The quantitative elliptic gap is uniform only
on a compact cell where both have positive lower bounds; it collapses as
`zeta_col` tends to zero at the positive-core interface.  Thus every returned
covector in this regular column subfamily is elliptic or scalar, not
hyperbolic.  Small-`g` continuity can preserve a nonresonant member only on
such a compact cell, not across the collar or over all covectors.  For a nonzero
azimuthal covector component,

    dot k_s=-s Omega' k_alpha,                            (36)

so the cotangent fiber is not returned and the actual object is a moving-
fiber skew product, not (34).

The finite-curvature charged ring has neutral/scalar limits in (23), (34)
and nonreturned fibers in (36).  Consequently a small coefficient
perturbation cannot give a uniform all-ray conclusion: near a scalar return,
an arbitrarily small periodic perturbation can be elliptic or hyperbolic.

### 3.1 First curvature-harmonic Melnikov comparator

The invariant returned column covector family used in (33)--(34) is
explicitly

    k_0=(p,0,q)                                           (36a)

in the rotating `(e_s,e_alpha,e_theta)` frame: `p` is radial, `q` is the
straight-column axial/toroidal component, and the cross-sectional angular
covector is zero.  It is this zero component, not a generic fixed wavevector,
that makes the unperturbed ray returned even when `Omega'!=0`.

To isolate geometry from the physical Cao profile jet, first bend a radial
column streamfunction into the toroidal metric

    h=1+delta s cos alpha,       W_alpha=s Omega/h.        (36b)

This is an exactly divergence-free metric-bent comparator.  Its physical
wavevector is `k=(p,0,q/h)`, so it also returns after one orbit.  Expanding
the full Leray BAS equation in the same moving transverse frame gives

    C(delta)=C_col+delta C_geom,1(alpha)+O(delta^2),       (36c)

    C_geom,1=
      [[-Omega s sin alpha,
        q s{5 Omega p^2+3 Omega q^2
             +2s Omega' p^2+s Omega' q^2}cos alpha/K^3],
       [-2 Omega q s(2p^2+q^2)cos alpha/K^3,
        Omega p^2 s sin alpha/K^2]],                     (36d)

where `K=sqrt(p^2+q^2)` and
`zeta_col=2 Omega+s Omega'`.  This formula includes the variation
`q/h`, cylindrical strain, frame rotation and the order-zero pressure
projection; deleting any one changes (36d).

Put

    a=2 Omega q/K,       b=zeta_col q/K,
    nu=sqrt(a b)>0                                       (36e)

for `q>0`.  With right modes
`v_+=(1,-i nu/b)`, `v_-=(1,i nu/b)` and the dual
`w_+^T=(1/2,i b/(2nu))`, the first Fourier coefficient driving `v_-` into
`v_+` is

    mu_geom=w_+^T[(C_cos-i C_sin)/2]v_-
      =i s Omega/4[
          sqrt(2 Omega/zeta_col) q/K+(2p^2+q^2)/K^2].    (36f)

The `d`th parametric resonance is

    2 nu(p,q,s)=d Omega(s),
    |q|/K=d sqrt[Omega/(8 zeta_col)]<=1.                 (36g)

For the first harmonic `d=1`, it exists on every regular cell satisfying
`Omega<=8 zeta_col`.  At that resonance,

    mu_geom,res
      =i s Omega/4[2+3 Omega/(8 zeta_col)]!=0.            (36h)

Thus toroidal metric bending alone opens the standard first-order
parametric coupling rather than cancelling by symmetry.  This is an exposing
comparator, not yet the physical Cao Floquet gap.  The actual first
coefficient is

    mu_Cao=mu_geom+mu_profile+mu_translation+mu_Leray-frame. (36i)

The reviewed 0044 formulas check the geometric `m -> m+/-1` character, but
their fixed-integer-harmonic grading is not the returned massive-ray grading
used here; they leave the odd Cao first cell as an unsolved linearized
Lane--Emden response.  They therefore do not determine the last three terms
in (36i), and cancellation against (36h) is still possible.  The next exact
Route-A construction is the Cao first curvature/profile/translation jet in
the action-flow common domain, inserted into (36f).  Independently, the
moving-fiber sectors in (36) and the DA realization at (25) remain part of
the all-covector Lyapunov problem.  Interface loss
`zeta_col ->0` is a separate mechanism from the interior resonance (36g).

### 3.2 Quarantined pre-recovery source-cell calculation

This subsection records the exact local operator expansion and bordered
translation moment that survived the recovery audit, but its stronger divided-
source and endpoint-remainder assertions are **not consumed below**.  In
particular, equations (36p.1)--(36u) were drafted during a Luna-low fallback.
They remain visible as provenance, not as premises for the physical growth
coefficient.  The independent compact-interior argument in Section 3.5 uses
only Cao's printed value estimates, subsequential elliptic compactness, and
the cancellation theorem in Hattori--Fukumoto.

The last paragraph can be sharpened without asserting the `C^1`
external-`epsilon` path that 0083 deliberately left open.  Use the exact
source equation and the reviewed `B_R` bordered realization.  In centered
core coordinates

    r=R+a y_1,       z=a y_2,       delta=a/R,

the streamfunction operator in Cao (2.1), after the Lane--Emden amplitude is
removed, is

    A_delta w=-(1+delta y_1)^(-1)
       div_y[(1+delta y_1)^(-1)grad_y w].                 (36j)

Consequently

    A_delta w=-Delta w
      +delta[partial_1 w+2y_1 Delta w]+O(delta^2)         (36k)

on every fixed core ball.  If `U` is the radial Lane--Emden cell,
`-Delta U=U_+^p`, the raw odd first source is therefore

    f_geom=(-partial_1U+2y_1U_+^p).                       (36l)

It is not orthogonal to the translation kernel.  Put

    E_p=integral_(B_1) U^(p+1)dy,       Lambda_p=integral_(B_1)U^pdy.

Radial symmetry and integration by parts give the exact moment

    <partial_1U,f_geom>
      =-integral(partial_1U)^2
        +2 integral y_1U^p partial_1U
      =-{1/2+2/(p+1)}E_p
      =-(p+5)E_p/[2(p+1)],                                (36m)

whereas

    <partial_1U,pU^(p-1)y_1>=-Lambda_p.                   (36n)

Thus the affine far-field/speed column is load-bearing.  With the convention

    L_p=-Delta-pU_+^(p-1),

its unique coefficient that makes the decaying odd equation solvable is

    beta_p=-(p+5)E_p/[2(p+1)Lambda_p],                    (36o)

and the complement solution `phi_1` is defined by

    L_p phi_1=f_geom+beta_p pU_+^(p-1)y_1,
    <phi_1,U_+^(p-1)partial_1U>=0,                        (36p)

with odd-`y_1`, even-`y_2`, regular-core and decaying exterior rows.  Equation
(36p) is Cao's Appendix-C Fredholm problem (C.1)--(C.4) applied to the first
**non-affine** coefficient of the explicit odd Green correction (A.1);
(36o) is the separate affine harmonic coefficient that makes the leading
Appendix-C source orthogonal to `partial_1U`.  In Cao's notation the remaining
`b_h partial_1U` centering multiplier is then zero at first order (Lemma C.3
is consistent with this through `b_epsilon=O(epsilon^2|log epsilon|)`).

The qualifier is essential.  Before the Kelvin--Hicks rows are imposed, the
displayed Green expression (A.1) is polyhomogeneous:

    mathcal F_delta=mathcal F_0
       +delta log(delta)[A_0+A_1y_1]
       +delta F_1^odd+O(delta^2|log delta|).              (36p.1)

The coefficient of its affine logarithmic part can be read off rather than
postulated.  The self-Green radial derivative contributes

    [psi_self]_(delta L,odd)
       ={kappa R/(4pi)}delta L y_1,       L=log(1/epsilon), (36p.2)

while the exact physical speed polynomial contributes

    -c_delta r^2/2
      =constant-W_delta R^2 delta L y_1
        -W_delta R^2 delta^2 L y_1^2/2.                  (36p.3)

Here `c_delta=W_delta L`.  Cao (3.34), or equivalently the leading part of
(3.36), is the exact coefficient identity

    W_delta R=kappa/(4pi)+O(1/L).                         (36p.4)

It cancels the `delta L y_1` coefficients in (36p.2)--(36p.3); their finite
difference belongs to the order-`delta` affine border below.  The constant
logarithmic coefficient is absorbed by the chemical-potential row in (3.16).
More explicitly, the first bracket in (A.1) contains
`a_delta L/(2R)-W_delta R L`; the mass/scaling relation (3.35) and the
chemical row of (3.16) give
`a_delta L/(2R)=kappa L/(4pi)+O(1)`, while (36p.4) gives
`W_delta R L=kappa L/(4pi)+O(1)`.  Multiplication by
`x_1-x_epsilon,1=R delta y_1` proves (36p.2)--(36p.4) directly from the
displayed A.1 terms.
This is also the cancellation used in Appendix A to improve the preliminary
`F=O(epsilon|log epsilon|)` bound to (A.3), `F=O(epsilon)`.

The dynamical cancellation must retain the quadratic term in (36p.3): for
the complete polynomial,

    r^(-1)grad(-c_delta r^2/2) cross e_theta
       =-c_delta e_z,       grad(-c_delta e_z)=0.         (36p.5)

A constant chemical row gives zero velocity.  Thus the comoving subtraction
removes the full speed velocity and no logarithmic strain is left; truncating
the threshold to its affine Taylor term would give a false residual through
the cylindrical `1/r`.  The reviewed mean-radius estimate
`R-x_epsilon,1=O(epsilon^2)` excludes a hidden first-order recentering term.
Thus no `delta log(delta)` coefficient survives in the relative-velocity BAS
cocycle; the first possible parametric term is the order-`delta` odd cell.

The divided odd source really converges, but this is a sequential source
calculation rather than endpoint differentiability.  Put

    H_0=(R C_s)^(-2/(p-1)),
    C_p=3R kappa/(4pi Lambda_p H_0).

Substitute `x=x_epsilon+s_epsilon y` in the displayed (A.1), use
`s_epsilon=C_s epsilon+O(epsilon^2)`, (3.34)--(3.36), and dominated
convergence for its compact logarithmic convolution.  After the affine
logarithmic rows above are removed,

    F(x_epsilon+s_epsilon y)/(delta H_0)
       -> F_1(y)
        =y_1[U(y)/2+gamma_p]
          +C_p integral eta_1 log(1/|y-eta|)U_+^p(eta)deta, (36p.6)

in `C^(1,alpha)` on compact core subsets.  The finite constant `gamma_p` is
the explicit finite affine remainder of (3.16)/(3.36); it is subsequently
fixed equivalently by the solvability number (36o), so it is not a fitted
datum.  The singular logarithm is locally integrable and the source is
compact; differentiating once and using the Poisson equation gives the stated
convergence.  Appendix-C solution continuity then sends this divided source
to the unique complement solution (36p).  Any sequence `delta_j->0` has this
same limit, and the reviewed bordered uniqueness fixes the same `a_t`; hence
the whole divided odd cell converges without selecting a differentiable
carrier branch.

On a compact regular subcell strictly inside `{U>0}`, (3.29)--(3.33) and
Proposition 3.11 give an `O(delta^2|log delta|)` scaled `C^1` remainder.
The exact local elliptic equation and interior Schauder estimates upgrade it
to

    P_delta/H_delta
       =U+delta T_1+O_C2,loc(delta^2|log delta|).          (36p.7)

This is exactly strong enough for `W` in local `C^1`, its strain in local
`C^0`, and the returned orbit/frame ODE.  It is not a global estimate through
the free boundary, collar, axis or exterior, and no ordinary derivative of
that entire map at `delta=0` is asserted.

The solution of (36p) is still defined modulo how the physical center is
represented.  Write

    T_1=phi_1+beta_p y_1+a_t partial_1U,
    Z_1=pU_+^(p-1)T_1.                                   (36q)

Circulation is automatic by oddness.  Differentiating the *normalized* exact
mean-radius row

    integral y_1(1+delta y_1)Z_delta dy=0

fixes the remaining kernel coefficient, rather than choosing it by hand:

    a_t=Lambda_p^(-1){
          integral y_1 pU_+^(p-1)(phi_1+beta_p y_1)dy
          +integral y_1^2U_+^pdy}.                        (36r)

Equations (36o)--(36r) are the blow-up form of the reviewed
circulation--mean-radius bordered inverse

    B_R(h,delta_mu,delta_c)
      =(Ah+a delta_mu+a r^2 delta_c/2,
        Dkappa h,(DM_1-R Dkappa)h).                       (36s)

They include the speed/affine column, radial centering, circulation and exact
mean-radius rows; axial centering is the even slice.  Let `Ren_1` mean:
subtract the two explicit affine `delta log(delta)` rows in (36p.1), divide by
`delta`, and take the local coefficient supplied by (36p.3).  If
`mathcal R_delta` denotes the two physical rows, the whole first tuple is

    (Z_1,mu_1,c_1)
      =-B_R,0^(-1) Ren_1(mathcal F_delta,mathcal R_delta). (36t)

This is not a use of the false unrestricted fixed-parameter inverse.  The
existence and uniqueness are exactly on the reviewed bordered product space.
Nor is it an import of 0083's open differentiated-map Route A: (36t) is the
source-derived polyhomogeneous coefficient on an interior regular cell, not a
Fréchet derivative of the global carrier path.  Here `B_R,0^(-1)` abbreviates
the unique limiting Lane--Emden bordered solution proved by the normalized
sequence argument; it is not an independently asserted endpoint Banach-space
inverse.  Equations (36p)--(36r) are its radial `m=1` ODE and finite rows.

To insert this cell without mixing clocks or normalizations, let `alpha` be
the returned orbit phase and put the full BAS equation in that phase and the
same returned orthonormal frame:

    d a/dalpha=K_delta(alpha)a,
    K_delta=K_0+delta K_1+O(delta^2|log delta|)            (36u)

in operator norm on this finite-dimensional returned amplitude bundle over
the chosen compact regular orbit.  Equation (36u) is deliberately not a
global action-flow graph expansion.

Normalize `w_+^*v_+=1`, `w_+^*v_-=0` as in (36f), and define the one common
resonant functional

    M_d[K_1]=(1/(2pi)) integral_0^(2pi)
       exp(-id alpha) w_+^* K_1(alpha)v_- dalpha.          (36v)

The physical-time coefficient in (36f) is `Omega M_1`; hence

    M_geom=mu_geom/Omega,
    M_phys=M_geom+Lambda_profile[T_1]+Lambda_frame,        (36w)

    Lambda_profile[T_1]=M_1[D_WK_0(B_0T_1)],

where `B_0T_1` is the column velocity reconstructed from the threshold
streamfunction.  `Lambda_frame` is the same functional applied to the
order-`delta` coefficient of the orbit-time, returned frame, covector and
whole-space Leray rows after the affine `delta log(delta)` translation has
been removed.  In coordinate-free form its integrand is the renormalized
coefficient of

    C(W,k)a=-(grad W)a
       +2k{k dot(grad W)a}/|k|^2                         (36x)

together with the returned-frame connection.  Thus every term in (36w) is a
specified integral of `U`, the unique radial solution (36p), the two numbers
(36o), (36r), and the fixed right/left column modes.  This is an exact
source-defined polyhomogeneous scalar functional.  It does not silently
promote the continuous 0083 path to a differentiable one.

There is, however, a sharp distinction between definition and noncancellation.
At leading core order, changing the physical major radius at fixed normalized
Lane--Emden cell only rescales the common clock.  It does **not** vary the
dimensionless number `M_phys`, so radius scaling is not a transversality
parameter.  For each fixed integer `p>=6`, (36p) reduces `M_phys` to an
analytic function of the regular orbit radius and returned wave angle on a
simple cell.  A proof that this function is not identically zero would make
its zeros isolated and would permit the rational-ray shell to avoid them.
The nonzero comparator (36h) alone does not prove that statement because the
profile and frame terms have the same `O(delta)` order.  The outstanding
identity is now the concrete one-dimensional assertion

    M_1[D_WK_0(B_0T_1)+K_frame,1] != -M_geom              (36y)

at one regular radius, with `T_1` fixed by (36p)--(36r).  If (36y) cancels at
an isolated radius, analyticity leaves neighboring shells; if it were an
identity, the next physical harmonic is `d=2` at order `delta^2`, requiring
the second bordered cell rather than reusing the fixed-`l` 0044 expansion.
Absence of such an identity in Cao's prose is not a proof.  The next
subsection evaluates the identity in a controlled endpoint regime rather than
assuming genericity.

### 3.3 Withdrawn raw two-entry computation

The remaining local noncancellation can be decided analytically near, but not
at, the regular core axis.  Let `A=U(0)>0`.  The Lane--Emden equation gives

    U(s)=A-A^p s^2/4+O(s^4),
    Omega(s)=A^p/2+O(s^2),       zeta_col(s)=A^p+O(s^2).  (36z)

Centering means that the total odd threshold cell has no linear term.
Writing `T_1=t(s)cos alpha`, (36l) therefore gives

    t(s)=t_3s^3+O(s^5),       -8t_3=5A^p/2,
    t_3=-5A^p/16.                                      (36aa)

The affine `beta_p s` has zero strain, and the translation-kernel piece is
removed covariantly by the center slice.  The remaining profile velocity is

    w_1,s=-(t/s)sin alpha,       w_1,alpha=-t'cos alpha.  (36ab)

Put

    A_t=t'/s-t/s^2,       B_t=t''.

In the same transverse basis as (33), direct differentiation of the
full-pressure BAS matrix gives the profile-strain block

    C_prof,1=
      [[-A_t sin alpha,       B_t(q/K)cos alpha],
       [-A_t(q/K)cos alpha, A_t(q^2/K^2)sin alpha]].      (36ac)

The perturbed streamline also changes the returned covector.  In the
solid-rotation limit, with `k_0=(p,0,q)`, its unique zero-mean periodic
correction is

    k_1,r=-(p A_t/Omega)cos alpha,
    k_1,alpha=-(p A_t/Omega)sin alpha,
    k_1,theta=0.                                         (36ad)

This follows directly from `dot k_1=-(grad W_0)^T k_1-
(grad w_1)^T k_0` in the rotating frame.  Differentiating the pressure symbol
`2kk^T(grad W)/|k|^2` and adding (36ac), while treating periodic transverse-
basis changes as a Floquet gauge, yields

    mu_prof+ray
      ={i/4}[A_t(-1+x+3x^2-2x^3)+B_t x]+O(s^3),
    x=q/K.                                                (36ae)

At the `d=1` resonance, (36g) and (36z) give `x=1/4+O(s^2)`.
Equations (36aa) and (36ae) then give

    mu_prof+ray=-25 i A^p s/1024+O(s^3),                 (36af)

whereas the already full-metric/frame comparator (36h) gives

    mu_geom=280 i A^p s/1024+O(s^3).                     (36ag)

The phase-speed correction multiplies `K_0`; its resonant off-diagonal row is
zero because `w_+^*K_0v_-=0`.  The base solid-rotation strain is spatially
constant, so the periodic orbit displacement contributes only `O(s^3)`.
The Luna-low scratch then combined only these selected entries as

    mu_scratch=255 i A^p s/1024+O(s^3).                  (36ah)

Equation (36ah) is **not** a physical exponent and is withdrawn.  It omits the
Hill-variable factor, the `d log(k_perp/k)/dt` and `tr L_perp` terms, and their
derivatives in the invariant scalar equation.  No reciprocal entry or
similarity normalization was derived for this truncated matrix, so it is not
retained even as a physical off-diagonal coefficient.  The recovery receipt
pins the pre-correction hash and the superseded 9/9 execution; the next
subsection replaces the claim from the primary equations.

### 3.4 Exact Hattori--Fukumoto reduction and smooth-center exponent

Use Hattori--Fukumoto's toroidal coordinates `(r,theta,s_parallel)` and put
`rho=1+delta y`, where their thinness `epsilon` is exactly the present
geometric ratio

    delta=a/R.                                             (36ai)

Write `k_perp` for the cross-sectional covector and `k_parallel` for its
ring-tangent component.  With

    A=d/dt log(k_perp/k)+tr L_perp,
    B=2 rho k_parallel^2
       (k_perp H k_perp^T)/(k_perp^2 k^2),
    C=-omega_s/rho,                                       (36aj)

their printed equation (7), including its parentheses, is

    d/dt [p;q]=[[A,B],[C,-A]][p;q],                       (36ak)

where

    p=(k/k_perp)rho k_perp dot a_perp,
    q_HF=(k/k_perp)(k_perp cross a_perp) dot e_parallel.

For a steady axisymmetric ring without swirl,
`D_t(omega_s/rho)=0`; hence `C` is constant along the particle orbit.  Direct
elimination, with no amplitude normalization suppressed, gives

    q_HF''=(A^2+BC-A')q_HF,
    q_HF''+Vq_HF=0,
    V=A'-A^2-BC.                                          (36al)

Substituting (36aj) yields the exact printed equation (9):

    V=2 omega_s k_parallel^2
        (k_perp H k_perp^T)/(k_perp^2 k^2)
      +[d^2/dt^2 log(k_perp/k)+d/dt(tr L_perp)]
      -[d/dt log(k_perp/k)+tr L_perp]^2.                 (36am)

The last two bracketed rows are precisely what the withdrawn calculation
missed.

For an arbitrary smooth leading circular-column profile, put

    Omega=U_theta^(0)(s)/s,
    zeta=omega^(0)(s),
    G(s)=integral_0^s tau[U_theta^(0)(tau)]^2 d tau,
    c_chi^2=cos^2 chi.                                    (36an)

Hattori--Fukumoto Appendix C carries arbitrary first-order fields
`U_r^(1)` and `U_theta^(1)` through the particle and covector equations.
Substitution into (36am) cancels them and gives their equation (22):

    V(t)=2 zeta Omega c_chi^2
      +delta U_theta^(0) F(s) sin(Omega t)+O(delta^2),    (36ao)

    F=Omega+2c_chi^2[(1-4sin^2 chi)zeta-Omega
       -G/(s^2[U_theta^(0)]^2)
          {(2sin^2 chi-1)zeta+Omega/2}].                 (36ap)

Thus the first-order invariant depends only on the leading profile; it does
not require an endpoint derivative of the first Cao correction.  The `d=1`
resonance is the printed equation (23)

    zeta/Omega=1/(8c_chi^2),                              (36aq)

and the standard resonant Hill calculation applied to (36ao) gives their
equation (24)

    gamma_hat=delta |U_theta^(0)| |
       (1/2)sin^2 chi-3/8
       -(3/8-(1/4)sin^2 chi)G/(s^2[U_theta^(0)]^2)|
       +O(delta^2).                                      (36ar)

For a smooth center,

    U_theta^(0)=Omega_c s+u_3s^3+O(s^5),
    zeta=2Omega_c+4u_3s^2+O(s^4),
    G/(s^2[U_theta^(0)]^2)=1/4+O(s^2).                   (36as)

Equation (36aq) gives `c_chi^2=1/16+O(s^2)`.  At the endpoint values,
`sin^2 chi=15/16`, and the full forcing coefficient is

    F/Omega_c
      =1+(1/8)[(-11/4)2-1-(1/4){(7/8)2+1/2}]
      =15/128.                                           (36at)

The unperturbed Hill frequency is `nu_0=Omega_c/2`, so the exponent is also
`|delta U_theta F|/(4nu_0)`.  Both evaluations give

    gamma_hat=delta |U_theta^(0)(s)| 15/256
       +O(delta s^3+delta^2).                            (36au)

This uses both resonant couplings, not one matrix entry.  If
`q=A_+exp(i nu_0t)+A_-exp(-i nu_0t)` and
`f=delta U_theta F`, projection of (36ao) gives

    d/dt [A_+;A_-]
       ={f/(4nu_0)}[[0,1],[1,0]][A_+;A_-]+O(delta^2),    (36au.1)

whose eigenvalues are `+/-f/(4nu_0)`.  Equivalently the returned Hill
discriminant is `2+{fT/(4nu_0)}^2+O(delta^3)>2`.  This paired invariant is the
hyperbolicity certificate.

For the Cao Lane--Emden cell `Omega_c=A^p/2`, hence

    gamma_hat=delta (15/512)A^p s
       +O(delta s^3+delta^2)>0                           (36av)

on a punctured sufficiently small regular orbit.  If
`t_hat=Omega_a t_phys` with the positive Cao core clock

    Omega_a=kappa/(2 pi a^2),                            (36aw)

then `gamma_phys=Omega_a gamma_hat`.  This pins the clock and shows why the
withdrawn `255/1024` could not be a physical exponent.

### 3.5 Cao compact-interior applicability without endpoint differentiation

The physical coefficient (36av) needs less than the quarantined claim
(36p.7).  Cao Proposition 3.10 and equation (A.1) give an odd correction of
size `O(epsilon_core)`; Proposition 3.11 leaves an
`O(epsilon_core^2|log epsilon_core|)` remainder.  The reviewed scale rows give

    x_epsilon,1=R+O(epsilon_core^2),
    s_epsilon=C_s epsilon_core+O(epsilon_core^2),
    delta=s_epsilon/R.                                   (36ax)

The apparent `delta log(delta)` rows in (A.1) are affine.  The self-Green
odd coefficient is `+(kappa R/(4pi))delta L y_1`, while the full speed
polynomial is `-c r^2/2`.  Equations (3.34)--(3.36) give
`cR/L=kappa/(4pi)+O(1/L)` in the notation where `c=W L`, so the odd logarithmic
coefficients cancel.  The constant logarithmic row is chemical potential.
The cancellation is dynamical only for the full polynomial:

    r^(-1)grad(-c r^2/2) cross e_theta=-c e_z,
    grad(-c e_z)=0.                                      (36ay)

An affine truncation would miss the quadratic term whose two physical
derivatives enforce the zero strain.

Here is the topology actually used.  After the affine rows are removed, the
rescaled A.1 source divided by `delta` is bounded in `W^(-1,q)` on its fixed
support.  Its nonlocal part has the form

    T f(y)=integral eta_1 log(1/|y-eta|)f(eta)deta.       (36az)

The scale and center expansions (36ax), the `C^1` normalized-profile
convergence, and `p>=6` give convergence of the zero-extended densities in
`C^(0,alpha)` on every nested positive-core set.  Split the integral into a
small disk and its complement.  On the complement the differentiated kernel
converges uniformly; on the disk `|grad log|=O(|y-eta|^-1)` is integrable in
two dimensions, while the Holder difference supplies the derivative modulus.
Therefore `Tf_epsilon->Tf_0` in `C^(1,alpha)` on nested compact sets.  This is
the promised explicit compact-log argument, not an appeal to pointwise
dominated convergence.

Appendix C maps the bounded odd source to a bounded weighted `W^(1,q)`
solution after its translation row is removed.  On a compact set
`K compactly contained in {U>0}`, the exact semilinear equation has uniformly
`C^(0,alpha)` coefficients.  Interior `W^(2,q)`, Sobolev embedding, and
Schauder estimates give the uniform bound

    sup_(0<epsilon<epsilon_K)
      ||(P_epsilon/H_epsilon-U)/delta||_(C^(2,alpha')(K))
       <=C_(K,alpha'),       0<alpha'<alpha.             (36az.1)

Indeed the divided odd source is uniformly bounded in `W^(-1,q)` by (A.1)
and Proposition 3.10.  On a nested positive-core pair `K compactly contained
in K'`, the compact logarithmic row maps the zero-extended, uniformly
`C^(0,alpha)` density to `C^(1,alpha')`; the local semilinear equation then
upgrades the `W^(1,q)` Appendix-C bound first to `W^(2,q)` on `K'` and then
to (36az.1) on `K`.  This is an estimate for normalized differences, not
endpoint differentiability of a selected carrier family.  Consequently every
sequence has a subsequence converging in `C^(2,alpha')` on `K`.  Bordered
uniqueness identifies the affine and translation rows.  More importantly,
Appendix C's explicit calculation shows that every
allowed first-order `U_r^(1),U_theta^(1)` cancels from (36ao).  Hence every
subsequence has the same scalar limit (36av).  To make uniformity explicit,
suppose the scalar remainder were not `o(1)` on a fixed small orbit tube.
Choose a bad `epsilon_j` and orbit point.  Compactness of the tube and the
uniform `C^(2,alpha')` bound give a convergent point and a first-cell
subsequence.  Hattori--Fukumoto Appendix C then gives the same coefficient
independently of that cell, contradicting the bad lower bound.  This proves
the uniform compact-interior asymptotic

    gamma_hat,epsilon/(delta s)=15 A^p/512+O(s^2)+o_epsilon(1)              (36ba)

on a fixed sufficiently small regular orbit and a small surrounding tube.
No global `C^1` endpoint carrier path, free-boundary estimate, or collar
uniformity is inferred.

### 3.6 Actual DA/tag/Maxwell packet and fixed observed-propagator global-input-to-local-output essential norm

Choose such a regular orbit and action--angle coordinates `(P,alpha,vartheta)`
on a tube compactly contained in the positive core.  The ring is axisymmetric,
so `vartheta` is a global toroidal angle and the comoving no-swirl flow has no
`vartheta` component.  For an integer `ell!=0` choose a smooth real `F(P)` so
that on the central orbit

    phi=F(P)+ell vartheta,
    k=dphi,
    W dot k=0,
    k_parallel/|k|=cos chi.                              (36bb)

Then `exp(iNphi)` is globally single valued for every integer `N`: the
toroidal harmonic is `Nell`, while `F(P)` is an ordinary single-valued
function.  This is the required tube-first/frequency-last phase, not a local
plane wave pasted around the torus.

This phase also repairs the higher-order periodicity caveat following
Hattori--Fukumoto equation (20).  It is a circle-valued exact first integral:

    (partial_t+W dot grad)phi=0,
    dot k=-(grad W)^T k,       k=dphi.                   (36bb.1)

Therefore, if `x(T_*)=x(0)` on the closed streamline, then

    k(T_*)=dphi(x(T_*))=dphi(x(0))=k(0)                 (36bb.2)

exactly, not merely through first order.  In the leading circular-column
coordinates, `dP` is radial and `d vartheta` is ring-parallel.  The compact-
interior `C^2` estimate (36az.1) consequently identifies (36bb) with the HF
initial direction `varphi_HF=0+O(delta)`; the `O(delta)` term is the physical
streamline-normal correction, not an unreturned arbitrary wavevector.

There remains one continuous tuning variable: `F'(P_0)`.  With physical
normal and toroidal components

    k_perp=F'(P_0)|grad P|,       k_parallel=ell/r,

it ranges the wave angle through an open interval.  At `delta=0`, define

    D(0,chi)=2 sqrt(2 zeta Omega) cos(chi)-Omega.

The HF resonance (36aq) is `D(0,chi_0)=0`, and

    partial_chi D(0,chi_0)
      =-2 sqrt(2 zeta Omega)sin(chi_0)!=0.               (36bb.3)

Use the two resonant solutions of the unperturbed Hill equation as a basis.
Its one-period multiplier is `-I` at (36bb.3).  After removing `-I`, putting
`chi=chi_0+delta kappa`, and dividing the logarithm of the returned matrix by
`delta`, periodic Floquet reduction gives

    H_1(kappa)=
      [[i{d_1+partial_chi D(0,chi_0)kappa}, b_HF],
       [conjugate(b_HF),-i{d_1+partial_chi D(0,chi_0)kappa}]],              (36bb.4)

up to a fixed harmless real/complex basis convention.  Here `d_1` is the
finite first orbit-period/direction detuning and the paired source calculation
(36au.1) gives `|b_HF|=15|U_theta|/256>0`.  Equation (36bb.3) uniquely chooses

    kappa=-d_1/partial_chi D(0,chi_0),
    chi_delta=chi_0+O(delta),

equivalently a value of `F'(P_0)`, that zeros the diagonal.  The two returned
multipliers are therefore

    -exp(+/-delta |b_HF|T_*+O(delta^2)),                 (36bb.5)

and are genuinely reciprocal and hyperbolic.  Equation (36bb.5) is the
**uncharged `g=0`** returned pair; the nonzero-charge conclusion below is only
continuity persistence of a real hyperbolic pair and expanding multiplier.
Smoothness of the off-diagonal
coefficient in `chi` gives

    delta A_HF(chi_delta)
      =delta A_HF(chi_0)+O(delta^2).                    (36bb.6)

The derivative in (36bb.3) stays bounded away from zero on a sufficiently
small compact regular `P`-interval.  Hence the same argument gives a smooth
choice `F'(P)` throughout that interval; integrating it gives one
single-valued `F(P)`.  Packet cutoffs are supported strictly inside this
interval, so their commutators are lower WKB order and do not change
(36bb.5).

The preservation of the leading `15/256` value can also be seen directly from
Appendix C.  Its wavevector starts as
`(sin chi,0,cos chi)+O(delta)` and its arbitrary first-order **flow** fields
cancel from the off-diagonal invariant in equation (22).  The exact
streamline-normal covector correction here is periodic; it contributes the
detuning in (36bb.4), rather than being declared to cancel.  Under the corresponding
periodic near-identity change `a=(I+delta S(t))b`, the first coefficient
changes by

    A_1 -> A_1+[A_0,S]-dot S.

At resonance, pairing its off-diagonal part between the resonant left/right
solutions and integrating over one period gives only a zero boundary term;
the first Melnikov coupling is invariant.  The diagonal covector correction
is removed by the choice of `kappa` in (36bb.4).  Thus the exact finite-Cao covector lies in
the Appendix-C `phi_HF=0+O(delta)` class and differs in exponent only at
`O(delta^2)`.  The saturation mechanism discussed after HF (20) applies to
their nonreturned approximate wavevector; it does not apply to
(36bb.1)--(36bb.2).  Without this exact-first-integral tuning the coefficient
(36av) would license only a finite transient, and none of the repeated-
circuit conclusions below would follow.

At the resonance, `k dot bar_omega=zeta k_parallel!=0`.  For either growing
Hattori--Fukumoto velocity polarization `a`, `k dot a=0`, set

    d=(k cross a)/(k dot bar_omega),
    xi_N=N^(-1)chi_tube Re[d exp(iNphi)],
    q_N^raw=curl(xi_N cross bar_omega).                   (36bc)

The leading *rescaled* Biot--Savart velocity is exactly `a`.  Indeed,

    q_pr=i(k dot bar_omega)d=i k cross a,
    N B q_N^raw=(i k cross q_pr)/|k|^2+O(N^-1)
               =a+O(N^-1).                              (36bd)

A compact Bogovskii correction restores `div xi_N=0` one order lower.  Choose
it equivariantly, so the corrected generator retains the same nonzero toroidal
character.  Every frozen *linear* moment and slice row has an axisymmetric
coefficient at the axisymmetric carrier.  Its value on the packet is therefore
an integral of an axisymmetric coefficient times `exp(iNell vartheta)` and
vanishes exactly for `Nell!=0`.  There is no same-character restoration
profile and hence no fixed low-frequency component in the weak-null sequence.
Exponentiating the corrected divergence-free field gives an exact volume-
preserving coadjoint path, not merely a formal tangent.

The transported tag tangent and initial Gauss correction are

    eta_N=-xi_N dot grad chi_g,
    e_N=-(g/epsilon_EM)grad(-Delta)^(-1)eta_N,
    b_N=0.                                                (36be)

They have zero charge variation.  In the same `H^s` vorticity topology used
below, `||q_N^raw||=Theta(N^s)`, while the tag row is
`O(N^(s-1))` and its Gauss field gains one additional inverse derivative.
Thus neither changes the principal normalization.  More explicitly,
`integral eta_N=0` because it is the derivative of a volume-preserving
transport of `chi_g`.  Its Coulomb field therefore starts at dipole order,
`e_N=O(|x|^-3)`, and

    integral_(|x|>1)<x>^(2sigma)|e_N|^2 dx<infinity
       exactly when sigma<3/2.                           (36be.a)

The weight mechanism is exact.  In three dimensions
`w(x)=<x>^(2sigma)` is a Muckenhoupt `A_2` weight precisely throughout the
power-law range `|sigma|<3/2`; Calderon--Zygmund Riesz, Leray, and Hodge
operators are consequently bounded on its weighted `L^2` Sobolev rows.  The
anisotropic action-flow graph row uses the same scalar weight on its transport
derivative, which is equivalent on the compact tube and has the same `A_2`
exterior Hodge row.  The lower choice `sigma>1/2` is not a singular-integral
condition: it separates the excluded Coulomb monopole (`|x|^-2`) from the
zero-charge dipole tangent (`|x|^-3`).  The compact-curl zero moment gives the
fluid Hodge tail the same admissible decay.  Hence `1/2<sigma<3/2` is an actual
interval for these packet paths, not a claim that the Coulomb base belongs to
the tangent space.

The coupled Maxwell assertion can be made explicitly on the constrained
field fiber.  Decompose `e=e_T+e_L`, `b=b_T`, with
`k dot e_T=k dot b_T=0`; the magnetic constraint removes `b_L`, while Gauss
fixes rather than evolves

    e_L={g eta/(i N epsilon_EM |k|^2)}k+O(N^-2).          (36be.0)

Thus `e_L` is a longitudinal order-minus-one slaved row.  After the common
transport phase is removed, the degree-one symbol on the four-dimensional
transverse field space and the slow Euler/tag space is

    D_1(k)=diag(0_slow,D_EM(k)),
    spec D_EM={-i c_g k_z+i c_EM|k|,
               -i c_g k_z-i c_EM|k|}.                   (36be.1)

The slow block contains the Euler polarization and its same-leaf tag; the tag
is not discarded as an independent triangular mode.  On the phase (36bb),
`W dot k=0`, and strict subluminality gives the uniform Sylvester denominator

    dist(0,spec D_EM)>=(c_EM-|c_g|)|k|.                  (36be.2)

If `C_sM,C_Ms` are the complete order-zero transverse current/Lorentz
couplings after inserting (36be.0), the order-minus-one change
`I+N^(-1)S` with

    S_sM D_EM=C_sM,       D_EM S_Ms=-C_Ms               (36be.3)

cancels both off-diagonal rows at order zero.  The longitudinal substitution
and induced slow Schur row are `O(N^-1)`; the direct change of the charged
base fluid/tag cocycle is `O(g^2)`.  This proves, rather than assumes, that
the Maxwell response is one WKB order lower away from (15).

More explicitly, let `T_N=I+N^(-1)Op(S)` on the compact tube and reconstruct
`e_L` by (36be.0).  In the common action-flow graph core its conjugated
generator has the fixed-frequency form

    T_N^(-1)G_(delta,g)T_N
      =N diag(0,D_EM)
       +diag(C_slow(delta,g),C_EM,0)+N^(-1)R_N,           (36be.3c)

where the estimate actually used is

    sup_N ||J_loc R_N||_(mathcal X_in^s -> X_DA,loc^s)<infinity.              (36be.3d)

This follows after a compact cutoff and the whole-space Hodge continuation;
no propagation of a positive spatial weight by the free Maxwell group is
claimed.  The commutator of
the order-minus-one `Op(S)` with the first-order diagonal block is order zero
and is exactly what solves (36be.3); the remaining commutators are order
minus one.  The longitudinal Gauss reconstruction maps the slow tag row to
the field space with the same order-minus-one gain.  This supplies the
common-domain Duhamel remainder used below rather than importing the
Euler-only 0039 estimate unchanged.

The charge quantifier is fixed here.  For each fixed sufficiently thin
`delta`, the reviewed fixed-parameter charged IFT and the finite-dimensional
returned cocycle give a finite constant `C_mon(delta)` and

    ||M_(delta,g)-M_(delta,0)||
       <=C_mon(delta)g^2+o_delta(g^2).                   (36be.3a)

If `c_HF delta` is a lower bound for the uncharged Floquet gap on the chosen
compact orbit tube, that gap persists whenever

    C_mon(delta)g^2<c_HF delta.                          (36be.3b)

Equivalently, for each such fixed `delta` there is a positive
`g_0(delta)`.  Neither reviewed 0080/0084 nor the compact-interior argument
above proves `sup_delta C_mon(delta)<infinity`; therefore no uniform
`g^2<=c_*delta` wedge is asserted.  The regime in which
`C_mon(delta)g^2` is comparable to or larger than `delta` remains an active
charged-cocycle problem.  No assertion is made at a Maxwell characteristic
coincidence.

For the same-norm statement, let `mathcal X_in^s` denote the weighted global
initial graph topology (1).  Fix once and for all a window `Pi_obs` supported
in a compact positive-core action chart and equal to one on a smaller master
tube around the selected closed orbit.  Define `X_DA,loc^s` to be the
`Pi_obs`-localized fluid vorticity action-flow/`H^s` observation, with the
weight removed but **without losing a derivative**.  Let the single fixed
operator `J_loc` collect this row and the corresponding instantaneous local
tag/field terms of (2).  Neither `Pi_obs` nor `J_loc` will depend on the
circuit count `j`; only the packet amplitude support may be narrowed inside
the master tube.
Normalize the complete exact state by

    Q_N=(q_N^raw,eta_N,e_N,0)/
        ||(q_N^raw,eta_N,e_N,0)||_(mathcal X_in^s).       (36be.4)

The packet vorticity is supported in one fixed compact tube, so the positive
weight is uniformly bounded there.  It and the equivariant Bogovskii term are
compact curls and hence have zero spatial mean.  No finite linear row
correction is present.  The zero Fourier row makes its whole-space Hodge velocity one order
smoother, while (36be) makes the tag and Gauss rows lower by at least one WKB
order.  Therefore, for the fixed carrier and tube,

    ||(q_N^raw,eta_N,e_N,0)||_(mathcal X_in^s)=Theta(N^s),
    ||J_loc(q_N^raw,eta_N,e_N,0)||=Theta(N^s).            (36be.5)

The global Hodge tail is harmless in the candidate weight interval because
the fixed-support zero row gives the exact estimate

    |qhat_N(xi)|
      =|integral (exp(-i x dot xi)-1)q_N(x)dx|
      <=|xi| integral |x q_N(x)|dx,                     (36be.6)

at low frequency and
the oscillatory toroidal character gives rapid far-field cancellation at
high `N`; the graph transport row has the same compact principal support.
Because the covector and physical frame return, the same `N^s|k|^s` factor
occurs at input and output; the Floquet ratio is therefore `lambda_+`, not
`lambda_+/N`.  All packets are smooth members of the full action-flow graph
domain.  On the transverse plane the principal map
`a -> i k cross a` is invertible and returns with `k`; therefore the HF
velocity multiplier is exactly the vorticity multiplier detected by the
local `H^s` term.

The full whole-space Leray symbol is already in (36am).  The `A_2` estimate
above bounds the global weighted Leray/Hodge entrance map, while restriction
to the compact tube is bounded from that space to the unweighted local row.
On the compact tube,
the order-zero Leray/Hodge and Sylvester symbols have uniformly bounded
derivatives; the carrier is smooth there by Section 3.5.  Standard
composition through `s+3` symbol derivatives supplies the first corrector.
Between the compact vorticity tube and the collar/exterior, Biot--Savart and
Leray kernels are smooth, and integration by parts in the nonstationary phase
is `O(N^-M)` for arbitrary fixed `M`.  Maxwell propagation is handled by
(36be.1)--(36be.3), not by deleting it.  Consequently, for each fixed `T`,

    sup_(0<=t<=T)||J_loc[S_g(t)Q_N-Q_N^BAS(t)]||
       <=C_(T,tube,s)N^(-1).                             (36bf)

This is the coupled fixed-time Euler--Maxwell Egorov estimate.  It invokes no
torus spectral theorem, no compact-velocity assumption, and no estimate
uniform in growing `T`.

Let `T_*` be the period of the exact closed orbit.  For one fixed sufficiently
thin charged carrier satisfying (36be.3b), let `lambda_+>1` denote the real
expanding multiplier that persists continuously from the uncharged
exact-return reciprocal pair (36bb.5); (36ba) supplies only its positive
leading gap, and (31) does not make the charged pair exactly reciprocal.
For every integer circuit count `j`, first fix `T=jT_*` and, if needed, shrink
only the packet tube inside the region where `Pi_obs=1`.  The observation
operator remains the same.  Now let `N` tend through the integers.  After the
global normalization (36be.4), the leading order-`s` derivatives are bounded
compactly supported amplitudes times `exp(iNphi)`.  Since `dphi` is nonzero on
the tube, nonstationary phase makes them weakly null; lower derivatives, tag,
Gauss and Hodge-tail rows tend strongly to zero at their relative WKB orders.
The same holds for the graph row because `W dot dphi=0`, so it introduces no
extra order-`N` factor.  Thus `Q_N` is weakly null in `mathcal X_in^s` and
every compact operator, as well as every finite-dimensional symmetry or
modulation projection, sends it strongly to zero.  Equation (36bf) and the
returned covector/frame then prove

    ||J_loc S_g(jT_*)||_essential,
         mathcal X_in^s -> X_DA,loc^s/sym
       >=c_obs lambda_+^j,       c_obs>0,                 (36bg)

for each fixed `j`, with one `c_obs>0` independent of `j`.  Indeed, at the
central returned ray the same principal `H^s` symbol norm occurs before and
after every circuit; choose `c_obs` below this fixed norm-equivalence ratio.
Shrinking the packet tube changes cutoff derivatives and the finite-time
Egorov constant, but these are lower WKB order.  They change only the required
threshold `N_0(j)`, not `J_loc` or `c_obs`.  This is an essential-norm lower
bound for the fixed observed propagator `J_loc S_g` from the **complete-state
weighted global input topology** to the local term appearing
in `D`.  It does not use one packet uniformly for growing time and does not by
itself assert nonlinear instability.

The exact coadjoint/tag/Gauss path needed to test Route A is explicit.  It is
not yet an exact fixed-conserved-row leaf path.  If `Phi_tau` is
the volume-preserving flow of the corrected smooth packet generator, set

    omega_tau=(Phi_tau)_*omega_g,
    chi_tau=chi_g composed Phi_tau^(-1),
    B_tau=B_g,
    E_tau=E_g-(g/epsilon_EM)grad(-Delta)^(-1)
                         (chi_tau-chi_g).                (36bh)

Then `div(epsilon_EM(E_tau-E_g))=g(chi_tau-chi_g)` exactly and total charge is
unchanged.  The coadjoint/Kelvin-circulation data and tag-distribution
Casimirs are exact under this push-forward, and `div omega`, `div B`, the
difference Gauss law and its zero monopole are exact.  These are the retained
automatic rows.

The declared nonautomatic finite rows on the smooth compact core are:

1. the hydrodynamic impulse components
   `I_omega=(rho_m/2)integral x cross omega dx`;
2. the center and modulation-slice rows represented by fixed smooth compact
   covectors; and
3. any explicitly declared carrier-parameter row represented by such a
   covector.

There is no additional Gauss-reconstructed electromagnetic finite row after
charge and zero monopole, which are already exact.  In particular the total
momentum

    P_tot=rho_m integral u dx
          +epsilon_EM integral E cross B dx              (36bh.1)

is not a continuous row on `mathcal X_in^s` because the generic
`O(|x|^-3)` Hodge velocity tail need not be `L^1`.  The frozen README fixes
`P_tot`; simply deleting it enlarges the leaf and cannot prove a statement on
that smaller fixed-`P_tot` leaf.

For the remaining declared rows `R_a`, an exact correction would first
require their typed continuous traces (or an augmented row norm) and then
carrier-specific compact smooth divergence-free generators `psi_j` for which

    M_(a j)=D R_a(Z_g)[
       curl(psi_j cross omega_g),
       -psi_j dot grad chi_g,
       -(g/epsilon_EM)grad(-Delta)^(-1)
             (-psi_j dot grad chi_g),
       0]                                                (36bh.2)

has a displayed nonzero determinant after the actual symmetry quotient.
No such witnesses or determinant are constructed in this attempt.  Density
of a smooth core does not establish independence or surjectivity of (36bh.2).
Accordingly the `O(tau^2)` impulse/slice defects of (36bh) are not corrected,
and no exact curve in the frozen fixed-row leaf has been proved.  A valid
repair is either to put (36bh.1) on a relative-momentum domain and extend
(36bh.2) by its electromagnetic row, or to prove `P_tot` automatic/redundant
for this particular curve and then construct the displayed witness matrix.

The same-target differentiability calculation survives as a conditional
bridge once such an exact fixed-leaf curve exists.  For a fixed packet and
finite interval, choose the initial generator two derivatives smoother.  If
`Y_tau=(Z_tau-Z_g)/tau` is formed from that corrected curve and `Y` solves the
full constrained linearized system, then `R_tau=Y_tau-Y` obeys the same linear
principal operator with a source of the form

    tau Q_2(Y_tau,grad Y_tau)+tau Q_EM(Y_tau,Y_tau),

including the transported-tag and longitudinal Gauss substitutions.  The
standard `H^s` product estimate, using the uniform `H^(s+1)` bound, gives

    sup_(0<=t<=jT_*)
      ||R_tau(t)||_(X_DA,loc^s plus field/tag terms)
       <=C_(j,N)|tau| ->0.                               (36bi)

Hence an exact corrected smooth path would be directionally differentiable in
the **same** local norm detected by `D`; no one-derivative-weaker solution-map
claim is needed.  This conditional statement is pathwise and does not assert
Frechet differentiability of the Euler solution map on the full weighted ball.

If the exact fixed-row curve were supplied, the derivative of distance to the
finite-dimensional modulated carrier manifold would be the quotient distance
to its tangent.  Choosing `Nell` outside its finite set of Fourier characters
would make the modulation derivative vanish while the instantaneous fluid
term of `D` detects the same-norm gain.  Dividing the hypothetical bound by
`tau` and sending `tau` to zero would then give
the correctly typed observation estimate

    ||J_loc S_g(jT_*)||_(mathcal X_in^s -> X_DA,loc^s/sym)<=C,               (36bi.1)

contradicting (36bg).  But (36bh.1)--(36bh.2) leave the antecedent open, so the
combined fixed-member nonlinear uniform-Lipschitz refutation is **blocked**.
No bound on the unobserved global operator `S_g` is inferred.  Units A--F,
including the fixed-`J_loc` linear observed-propagator theorem, do not depend
on this missing exact-leaf correction and remain established at author scope.

### 3.7 Complementary charged scaling: an exact rescaled Melnikov target

The absence of a uniform `C_mon` bound leaves the regime `g^2/delta=O(1)`
active.  At each fixed `delta`, use the algebraic parameter `tau=g^2` from
reviewed 0080.  If `B_(I,delta)` is its reviewed circulation--impulse
bordered derivative and `M_delta` the eliminated Maxwell forcing in the
steady map, the charged carrier derivative is exactly

    dot Z_delta:=partial_tau Z_(delta,tau)|_(tau=0)
       =B_(I,delta)^(-1)M_delta(Z_(delta,0)).             (36bj)

This is a fixed-`delta` derivative and makes no endpoint claim in `delta`.
Let `K_E(Z_delta)` denote the returned Euler BAS coefficient in the covariant
orbit/frame gauge.  The order-`tau` slow coefficient is

    K_ch,delta
      =D_Z K_E(Z_(delta,0))[dot Z_delta]
       +K_Lorentz[chi_0,F_1,B_1]
       +K_slice[dot c_delta,dot center_delta],            (36bk)

where `F_1=lim_(g->0)F_g/g` and `B_1=lim_(g->0)B_g/g` are supplied by the
fixed-parameter Maxwell inverse.  In displacement variables the direct row
is obtained from

    rho_m^(-1){-(d dot grad chi_0)F_1
                    +chi_0 a cross B_1},                 (36bl)

followed by the full pressure projection.  The transverse radiative Maxwell
field generated by the packet remains in the `N^(-1)` Schur row (36be.3c),
so it does not enter (36bk) in the frequency-last BAS limit.

Let `Phi_delta(t)` be the uncharged returned fundamental matrix and
`v_-,w_+` the covariantly normalized resonant right/left pair.  The exact
fixed-`delta` charged Melnikov row is

    m_ch(delta)=(1/T_delta) integral_0^Tdelta
       exp(-i Omega_delta t)
       w_+^* Phi_delta(t)^(-1)K_ch,delta(t)
                    Phi_delta(t)v_- dt.                  (36bm)

Together with the analogous diagonal detuning row, (36bm) determines the
first derivative of the returned discriminant in `tau`.  For the joint scale
`tau=lambda delta`, the candidate leading slow matrix is

    H_eff(lambda)=H_HF+lambda H_ch.                       (36bn)

Here `H_HF` has the paired nonzero coefficient (36au.1).  Equation (36bn) is
an earned fixed-parameter target, not yet a joint-limit theorem: completing
it requires a `delta`-uniform bound for `B_(I,delta)^(-1)M_delta`, convergence
of the orbit/frame rows in (36bk), and convergence of `delta^(-1)` times the
returned discriminant.  If those land, `det H_eff(lambda)` gives the exact
hyperbolic/elliptic tongues and tests whether charge closes the HF gap.  If
they do not, (36be.3b) is the strongest charged conclusion.  This is the
active continuation for the complementary charge/thinness regime.

A materially different way to construct that missing uniform limit is to
solve the charged straight column first.  Let its leading fields be

    W_tau=s Omega_tau(s)e_alpha,
    bar omega_tau=zeta_tau(s)e_parallel,
    B_tau=g mathcal B_tau(s)e_parallel,
    F_tau=g mathcal F_tau(s)e_s,       tau=g^2,

and put `k=p e_s+q e_parallel`, `x=q/sqrt(p^2+q^2)`.  In the transverse basis
`e_1=e_alpha`, `e_2=(q e_s-p e_parallel)/|k|`, the DA inverse is

    d={k cross(a_1e_1+a_2e_2)}/(q zeta_tau),
    eta=-d dot grad chi_tau=a_1 partial_s chi_tau/zeta_tau.

The electromagnetic and radial-force rows of this column are not placeholders.
Writing `E_tau=g mathcal E_tau e_s` and taking the parallel direction to be
`z`, the stationary Maxwell constraints and Euler radial balance reduce to

    (s mathcal E_tau)'=s chi_tau/epsilon_EM,
    mathcal B_tau'=-mu_EM chi_tau s Omega_tau,
    mathcal F_tau=mathcal E_tau+s Omega_tau mathcal B_tau,
    p_tau'=rho_m s Omega_tau^2+tau chi_tau mathcal F_tau.   (36bo.0)

The signs follow from `curl(mathcal B e_z)=-mathcal B'e_alpha` and
`e_alpha cross e_z=e_s`.  These equations determine the Maxwell fields and
pressure once the charged vorticity/streamfunction constitutive row and its
compact free boundary are supplied.  That remaining elliptic row, rather than
Maxwell algebra, is the existence problem for the charged-column-first route.

Projecting the direct Lorentz force onto `k`-transverse velocity gives the
exact constrained slow column block

    C_col(tau)=x[[0,-Z_tau],[R_tau,0]],                  (36bo)

    Z_tau=zeta_tau+tau chi_tau mathcal B_tau/rho_m,
    R_tau=2Omega_tau+tau{chi_tau mathcal B_tau
             +(partial_s chi_tau)mathcal F_tau/zeta_tau}/rho_m.            (36bp)

The two terms in `R_tau` come respectively from `chi a cross B` and
`eta F`; the upper correction comes from the other polarization of
`chi a cross B`.  Whole-space pressure is included by the transverse
projection.  The radiative Maxwell amplitude is still the order-`N^(-1)`
Sylvester row, so (36bo) is the complete frequency-last slow principal block.
When `Z_tau R_tau>0`, its inertial frequency is

    nu_tau=|x|sqrt(Z_tau R_tau),                         (36bq)

and the first curvature resonance can be tuned by

    2nu_tau=Omega_tau,
    |x_tau|=Omega_tau/{2sqrt(Z_tau R_tau)}<=1.           (36br)

The orbit-by-orbit coverage quantity is

    K(s,tau)=4Z_tau(s)R_tau(s)/Omega_tau(s)^2.           (36br.1)

If `Z_tau R_tau<0`, (36bo) is already hyperbolic at column order.  If
`Z_tau R_tau>0`, the first curvature resonance is admissible precisely when
`K>=1`, at `|x|=1/sqrt(K)`; more generally the `d`th resonance requires
`|x|=d/sqrt(K)<=1`.  Thus `0<K<1` excludes every integer `d>=1` curvature
resonance at that point, without excluding other channels or giving a uniform
gap.  A restoring candidate must instead prove the required common
positive-Krein signs of `Z_tau,R_tau` and constants `kappa_0,eta>0` with
`0<kappa_0<=K(s,tau)<=1-eta` uniformly on its declared core.  A carrier-level
P2 analysis must evaluate `K` throughout the whole regular core and evaluate
`b(s,tau)` on the entire resonant set `{K>=1}`, rather than select one orbit;
one resonant point with `b!=0` already supplies a refutation mechanism.
Even the uniform `K` bound classifies only curvature-parametric resonance:
`nu_tau=|x|sqrt(Z_tau R_tau)` tends to zero as the accessible angle `x` tends
to zero.  Complete-DA restoring control additionally requires either an
action-flow graph/accessibility weight controlling this sector or a
fixed-Casimir zero-frequency complement with a resolvent-modulation estimate.

For a charged-column branch with a `delta`-uniform first curvature jet, write
the periodic coefficient as

    C(delta,tau,alpha)=C_col(tau)
       +delta C_1(tau,alpha)+o(delta).

In the covariantly continued resonant frame the paired coefficient is the
explicit functional

    b(tau)=(1/(2pi))integral_0^(2pi)
       exp(-i alpha)w_+(tau)^*C_1(tau,alpha)v_-(tau)dalpha.                 (36bs)

It includes the geometric metric row, the charged profile/translation jet,
the tag terms in (36bp), and the returned Leray/frame rows in one
normalization.  Reality and the Hamiltonian positive-frequency pairing give
the conjugate paired leading off-diagonal entry; after detuning is tuned as in
(36bb.4), the first slow exponents are
`+/-delta|b(tau)|+o(delta)`.  This leading exponent pairing does not establish
exact charged monodromy reciprocity without an invariant symplectic
two-dimensional reduction.

At `tau=0`, (36au.1) gives
`|b(0)|=15|U_theta^(0)|/256!=0`.  Hence on any connected analytic
charged-column cell on which (36bo)--(36bs) are defined there is
`tau_resp>0` with `|b(tau)|>=|b(0)|/2` for `|tau|<tau_resp`; there is no stable
first-harmonic window in this neighborhood.  Globally on a simple-mode cell,
`b` is not identically zero and its zeros are isolated.  Every nonzero point
produces a hyperbolic tongue.  At an isolated zero `tau_0`, the next harmonic
is not merely named: with

    K_j(t)=Phi_0(t)^(-1)C_j(t)Phi_0(t),
    B_2=int_0^T K_2(t)dt
       +int_(0<t_2<t_1<T)K_1(t_1)K_1(t_2)dt_2dt_1
       +B_period/frame,                                   (36bs.1)

the `d=2` paired coefficient is the resonant off-diagonal entry of `B_2`
after its diagonal detuning is removed.  This is the exact ordered-integral
test at a zero of (36bs), including the second curvature cell, the double
first-cell interaction, and orbit-period/returned-frame correction.  A zero
of both coefficients activates the next Fourier harmonic; it is not evidence
for stability.
What is missing is now exact and constructive: either prove a
`delta`-uniform charged Cao-to-column expansion, or solve the charged-column
elliptic/free-boundary system and its first curvature jet directly.  The
fixed-parameter 0080 IFT alone supplies neither, so (36bo)--(36bs) are a
conditional positive route rather than a verdict on the complementary
charged regime.

## 4. What this decides for Route A

The exact streamfunction-normal comparator and straight-column sectors remain
neutral/elliptic, but the invariant Hattori--Fukumoto calculation finds a
different, genuinely DA, returned polarization with positive first-order
Floquet exponent on every sufficiently small punctured center annulus.  Its
fixed observed-propagator, global-input-to-local-output lift gives (36bg) for
every fixed sufficiently thin carrier and every charge satisfying
(36be.3b).  This is a linear complete-state theorem and is stronger than any
inference from the 0085 Hessian saddle.  The frozen nonlinear
uniform-Lipschitz conclusion would follow from the same-target calculation
(36bi)--(36bi.1) only after an exact path in the README's fixed conserved-row
leaf is constructed.  Equations (36bh.1)--(36bh.2) show that this Unit-G input
is absent, so the combined nonlinear verdict remains blocked.

The linear theorem remains carrier- and observable-scoped.  It proves neither
one nonlinear datum escaping every neighborhood nor a general nonlinear
instability theorem, and it decides no nonlinear orbital-persistence notion.
It also does not decide charged members for which
`C_mon(delta)g^2` is comparable with or larger than the `c_HF delta` gap.
Those members remain an active part of Route A.  A center-stable manifold for
the linearly hyperbolic small-charge members would at most describe a prepared
subset; without Unit G it neither proves nor refutes their fixed-leaf
open-ball estimate.

Even after the principal theorem, LP2 requires more: the full closed
Euler--Maxwell/tag generator, `P_int/P_rad` with cumulative decay only on
`P_rad`, bounded internal dynamics, modulation, a quasilinear invariant
neighborhood, and the all-time open-ball estimate.  No positive-codimension
graph can replace that estimate.

## 5. Route B: two exact refutations and the surviving Hill problem

### 5.1 Constant-lambda whole-space Beltrami is empty

Suppose `u in L2(R3)`, `div u=0`, and `curl u=lambda u` with constant real
`lambda`.  Fourier transformation gives

    i k cross uhat=lambda uhat.                           (37)

For `lambda!=0`, a nonzero value requires `|k|=|lambda|`, a sphere of
Lebesgue measure zero.  An `L2` Fourier function supported there is zero.
For `lambda=0`, curl-free, divergence-free `L2` fields also vanish.  Hence
there is no nonzero finite-energy constant-lambda whole-space carrier.

The formal energy-helicity Hessian fails independently.  On the two helical
polarizations `curl h_+/-=+/-|k|h_+/-`, the symbol of `I-alpha curl` is

    1-alpha|k|,        1+alpha|k|.                        (38)

For every nonzero `alpha` it has both signs at high frequency.  Equations
(37)--(38) refute only the constant-lambda/naive helicity subroute, not
generalized Beltrami or swirl carriers.

### 5.2 The smooth compact Gavrilov swirl carrier is not coercive

Gavrilov supplies a genuinely smooth compact finite-energy generalized-
Beltrami carrier with meridional flow and swirl.  Corrected 0032/0038 and
0039 establish on one member a dynamically accessible full-pressure
hyperbolic Kelvin return and the exact whole-space linear lower bound

    ||S(jT_*)||_ess>=lambda_+^j,       lambda_+>1.         (39)

If a conserved constrained Hessian controlled the kinetic/shape norm by
`E2(q,q)>=gamma||q||_E^2` modulo its finite symmetry kernel, conservation
would make the quotient linear group uniformly bounded.  Equation (39)
contradicts that consequence.  Therefore this concrete smooth swirl carrier
cannot satisfy Route-B coercivity.  This is a carrier-scoped refutation, not
a generalized-swirl no-go.

### 5.3 Hill requires a bulk/interface three-dimensional Hessian

Inside Hill's ball the vorticity has the Cartesian form

    omega_H=C e_z cross x,                                (40)

and jumps to zero outside.  For a smooth divergence-free displacement `xi`,
the physical coadjoint tangent has both pieces

    delta omega=-[xi,omega_H] 1_B
                 +(xi dot n)omega_H delta_(partial B),    (41)

up to the fixed bracket sign convention.  A compact `xi` supported strictly
inside `B` has `xi dot n=0` but generally
`-[xi,omega_H]!=0`.  Thus boundary vector-spherical-harmonic deformations do
not span the full three-dimensional DA leaf.

Choi's Theorems 1.1--1.2 give all-time stability only for nonnegative
axisymmetric no-swirl relative vorticity in the stated `L1+L2+impulse`
metric.  They neither evaluate (41) nor provide a coercive Hessian for its
interior nonaxisymmetric modes.  The surviving Hill construction is the
coupled bulk/interface/exterior Hodge second variation on (41), sector by
sector after translations.  Until that form is derived, Hill is neither a
full-3D positive carrier nor refuted.

## 6. Route verdicts at this analytic rung

- Route A's fixed observed-propagator theorem is established for every fixed
  sufficiently thin Cao member once its charge obeys
  `C_mon(delta)g^2<c_HF delta`: the exact returned phase, real expanding
  multiplier, constrained DA/tag/Gauss packet and Units A--F give (36bg).
  Exact reciprocity at nonzero `g` is not claimed.
- Route A's frozen all-time uniform-Lipschitz open-ball refutation is blocked
  at Unit G.  The exact coadjoint/tag/Gauss path (36bh) has the right tangent,
  but the fixed-`P_tot` leaf is not typed on `mathcal X_in^s` and no explicit
  carrier-specific row-witness matrix (36bh.2) has been constructed.  A
  result on the enlarged leaf would not imply the frozen smaller-leaf claim.
- Route A for `C_mon(delta)g^2` comparable with or larger than `delta` remains
  active.  A uniform `g^2<=c_*delta` refutation is also open because no
  `delta`-uniform bound on `C_mon` has been proved.
- Route B constant-lambda Beltrami is refuted by (37)--(38).
- Route B on the actual compact Gavrilov swirl carrier is refuted by the
  reviewed accessible growth theorem (39).
- Route B Hill remains active at the explicit bulk/interface tangent (41),
  whose full constrained Hessian is the next construction.
- Route C has no LP2 verdict.  Choi's exact axisymmetric Hill theorem is a
  restricted all-time comparator; any finite-time or symmetry-restricted
  continuation remains below LP2.

The strongest exact Cao result is therefore linear and observed, not a
particle or nonlinear-persistence claim: on every fixed sufficiently thin
Cao member the compact-interior first-curvature parametric resonance and a
real expanding multiplier survive for `0<|g|<g_0(delta)`, producing the
fixed-`J_loc` essential-norm lower bound (36bg).  Its conversion to the frozen
nonlinear Lipschitz contradiction is blocked by the exact fixed-row leaf
construction.  The complementary charge/thinness regime is unresolved.  The best smooth
compact swirl alternative tested so far is likewise genuinely hyperbolic and
cannot be coercive.  Route B remains active because Hill survives as a full
bulk/interface three-dimensional problem and the broader materially different
carrier class has not been exhausted.

# Analytic trace and invariant-metric reduction

## Fixed returned cocycle

Fix the activated 0032 carrier, a compact regular Gavrilov pressure plateau,
and Baldi's chart `x=Phi(sigma,beta,I)`.  On the plateau

    F_t=[[1,0,t Omega_1'],[0,1,t Omega_2'],[0,0,1]],
    k_z(t)=F_t^(-T)(n,-m,0).

Choose `I_0,n,m` with `n Omega_1'(I_0)=m Omega_2'(I_0)`.  The covector is
then returned in action coordinates after each meridional period; the
remaining beta-shift is identified by the actual rigid axial rotation.

In a returned physical orthonormal frame `E(t)` inside `k(t)^perp`, the exact
Kelvin--Leray amplitude equation is

    dot y=C(t)y,
    C=E^*[ -L+2 k tensor(k^T L)/|k|^2 ]E-E^* dot E,
    C(t+T)=C(t).

The exact ordered-integral object is

    M=Pexp integral_0^T C(t) dt.

No chart Euclidean metric is used in this definition.

## Conserved identities and pressure-localizable sector

For `c=k cross A`, direct differentiation of the full pressure-retaining
system gives

    dot c=L c+(varpi dot k) A,
    A=(c cross k)/|k|^2,
    d(varpi dot k)/dt=0.

The second term is essential and survives exactly when a general transverse
polarization is dynamically accessible.  The passive Cauchy shortcut
`dot c=Lc` is therefore not an admissible trace derivation.

For the invariant pressure phase `k=grad P`, steady Euler and localizability
imply

    L u_*=-k,       u_* dot k=0,       u_* dot grad|u_*|^2=0.

Baldi's dense irrational tori and continuity give `|u_*|^2=S(P)`.  Hence

    A_1=u_*,
    A_2=(k cross u_*)/(|k|^2 |u_*|^2)

are two exact bounded amplitudes.  The normalized complement has defect

    dot A_2-C A_2
      =[A_2 dot grad|u_*|^2/(2|u_*|^2)] A_1=0.

Thus the pressure-normal principal cocycle has an explicit positive periodic
metric and identity return.  This is a genuine bounded sector, not a claim
about angle phases.

## Polarization area and metric classification

For any two full-pressure solutions, including moving `k(t)`,

    area=(A_1 cross A_2) dot k

is constant.  After returned-frame identification this gives `det M=1`.

Let `nabla_t H=dot H+C^T H+H C`.  A positive `T`-periodic solution of
`nabla_t H=0` exists exactly when the returned map is conjugate to an
orthogonal matrix:

* if `|tr M|<2`, `M` is elliptic and diagonalizable; averaging a positive
  form over its compact rotation closure produces a periodic positive metric;
* if `M=+I` or `M=-I`, every positive initial form returns;
* if `|tr M|>2`, positivity is impossible because an H-orthogonal map cannot
  have reciprocal eigenvalues off the unit circle;
* if `|tr M|=2` but `M != +/-I`, positivity is impossible because an
  H-orthogonal map is diagonalizable while a nontrivial Jordan map is not.

This proves the exact metric alternative and separately tests scalar identity
cases.  It does not determine the source-specific trace.

## What reducible algebra does and does not determine

The action-angle equations determine `F_t`, `k(t)`, the conserved `varpi dot k`,
and the returned period.  They do not determine `L(t)` in the polarization
equation from frequencies alone: `L(t)` also contains the physical chart
metric and its connection.  The full trace remains the convergent
finite-dimensional ordered series

    tr M = 2 + sum_{r>=1} integral_{0<t1<...<tr<T}
      tr(C(tr)...C(t1)) dt1...dtr.

Conserved area fixes only the determinant; it does not select hyperbolic,
parabolic, or elliptic trace.  The small-shell comparison
`M -> [[1,0],[-sqrt(2)pi,1]]` fixes a limiting transient but not the sign of
`Delta=(tr M)^2-4`.

An attempted inference “determinant one plus near-Jordan limit implies
growth” is invalid.  The exact next object is the source-defined coefficient
`C(t;I_0,n,m)` and its ordered trace, not a boxed matrix power.

## Domain, accessibility, and global pressure

The physical space is the kinetic-norm closure of smooth coadjoint tangents,
modulo Euclidean symmetry distance.  The accessible leading vorticity is

    delta varpi_pr=i[(k dot varpi)d-(k dot d)varpi],   k dot d=0.

Whenever `varpi dot k !=0`, any transverse amplitude can be injected by
`d=(k cross A)/(i varpi dot k)`.  The exact packet is scalar-normalized by its
Biot--Savart velocity norm; its raw velocity scale is `N^-1`.  The full-space
Leray multiplier, pressure tail, collar commutators, exterior smoothing blocks,
and quotient leakage enter a fixed-time remainder relatively as
`C(T,delta)N^-1`.  No annulus wall is introduced.

If a later exact trace result gives hyperbolic or nontrivial parabolic gains,
the 0025 bridge `for every j exists delta_j exists N_j` transfers them to exact
accessible packets.  If the trace is elliptic, the positive periodic metric
route must still be extended over profiles and nonlocal terms before any
carrier stability conclusion.

## Route verdict

The pressure-normal invariant sector and the `SL(2,R)` metric classification
are established.  The source-specific resonant trace is blocked, not
refuted: the missing construction is an analytic evaluation or sign enclosure
of the ordered integral for `C(t;I_0,n,m)`, including the physical chart
connection.  No numerical activation is claimed by the present receipt; if
analytic reductions leave only that scalar remainder, a separately activated
small-ratio design is required.

## Append-only continuation: explicit metric-pulled resonant coefficient

The action-angle derivative matrix and resonant covector are

    B=D(Omega)=[[0,0,a],[0,0,b],[0,0,0]],
    kappa=(n,-m,0),       n a(I_0)-m b(I_0)=0.

Write `J=D Phi` and `g=J^T J`.  Since
`J^{-T}kappa=J(g^{-1}kappa)`, the cross-product identity

    (Jx) cross (Jy)=det(J) J^{-T}(x cross y)

pulls the exact physical equation to

    c_z_dot=B c_z + mu R_g(c_z),
    R_g(x)=det(J) g^{-1}[x cross h]/(kappa^T h),
    h=g^{-1}kappa,
    mu=varpi dot (J^{-T}kappa).

This is the requested coefficient reduction.  It contains no unevaluated
physical cross product or hidden Euclidean chart metric.

For Baldi's source functions, put `theta=beta+eta` and use the moving
orthonormal cylindrical frame.  The three columns of `J` are

    v_sigma=(rho_sigma, rho eta_sigma, zeta_sigma),
    v_beta =(0,rho,0),
    v_I    =(rho_I, rho eta_I, zeta_I),

so the metric entries and density are the explicit scalars

    g_ss=rho_s^2+rho^2 eta_s^2+zeta_s^2,
    g_sb=rho^2 eta_s,
    g_sI=rho_s rho_I+rho^2 eta_s eta_I+zeta_s zeta_I,
    g_bb=rho^2,
    g_bI=rho^2 eta_I,
    g_II=rho_I^2+rho^2 eta_I^2+zeta_I^2,
    det J=rho(rho_sigma zeta_I-zeta_sigma rho_I).

All functions in these equations are the analytic Baldi functions evaluated
at `(sigma,I_0)`; beta drops out by axial symmetry.  Set

    e=(m,n,0),  f=(0,0,1),  q=n^2+m^2,
    h=(h_s,h_b,h_I)=g^{-1}(n,-m,0),
    d=kappa^T h,
    alpha=det(J)/d.

The constraint `kappa dot c_z=0` is exactly `c_z=x e+y f`.  Since resonance
also gives `B f=(a/m)e`, equation (14) is the explicit scalar periodic system

    [x_dot]   [ mu alpha E_ee,    a/m+mu alpha E_ef ] [x],
    [y_dot] = [ mu alpha E_fe,    mu alpha E_ff       ] [y],

where

    E_ee=e^T g^{-1}(e cross h)/q,
    E_ef=e^T g^{-1}(f cross h)/q,
    E_fe=f^T g^{-1}(e cross h),
    E_ff=f^T g^{-1}(f cross h).

Equations (14)--(16) are the requested source-function coefficient.  They
retain the full accessible factor `mu`; setting it to zero removes the term
that dynamic accessibility supplies.

The coefficient is periodic in `sigma(t)=sigma_0+Omega_1(I_0)t`, and its
returned-frame monodromy is the ordered solution of this explicit 2x2 system.
The first trace correction about the solvable small-shell coefficient `C_0`
is the concrete Duhamel scalar

    T_1(I_0,n,m)=integral_0^T tr[U_0(T,t) C_1(t) U_0(t,0)] dt,

where `C_1=lim_{I->0} I^(-1/2)(C-C_0)` is obtained by inserting the Baldi
Taylor coefficients into the displayed metric entries, and `U_0` is the
known 0019 small-shell propagator.  Analyticity on the compact sigma-circle
gives the controlled remainder

    tr M=tr M_0+sqrt(I_0) T_1+O(I_0),
    |O(I_0)| <= K I_0,

with `K` bounded by the second Taylor-remainder sup norms and the finite-time
Duhamel estimate.  Since `tr M_0=2`,

    Delta=(tr M)^2-4=4 sqrt(I_0) T_1+O(I_0).

Substitution of the actual Baldi/Gavrilov source expansion removes that residual
at first order.  With `epsilon=sqrt(2 I)`, the returned-frame coefficients are

    C0 = [[0,0],[-1/sqrt(2),0]],
    C1(sigma) = [[3 cos(sigma)/2, sin(sigma)/sqrt(2)],
                 [sin(sigma)/(2 sqrt(2)), -cos(sigma)]].

The resonant transverse covector component is O(epsilon^3), hence does not
enter C0 or C1.  Since U0(t)= [[1,0],[-t/sqrt(2),1]], exact integration gives

    integral_0^(2 pi) U0(2 pi-t) C1(t) U0(t) dt = diag(pi,-pi),
    trace = 0.

Therefore the previously suggested sqrt(I) coefficient is exactly T1=0; this
is an evaluated cancellation, not an assumed normal-phase Jordan limit.  The
first potentially nonzero trace term is O(I) and requires both the second
ordered product of C1 and the genuine C2 coefficient.  C2 needs the cubic
Taylor data of g_c, gamma_c, action inversion, and the second-order frame
gradient, which are not present in the frozen source expansion (known only
through O(epsilon^2)).  The strongest source-supported statement is thus

    tr M = 2 + O(I),   Delta = O(I),

with the missing cubic source Taylor rung explicitly identified; no numerical
or sign claim for its coefficient is made.

### Cubic recurrence and the order-(I) obstruction

The exact source definitions impose the recurrence
\(F_c(g_c(\sigma))=\sigma\), where
\(F_c(\theta)=\int_0^\theta\partial_c\gamma_c(\tau)d\tau\), followed by
\[
\rho=1+\sqrt{2\gamma_c(g_c)}\sin g_c,qquad
\zeta={\sqrt{2\gamma_c(g_c)}\cos g_c\over
1+\sqrt{2\gamma_c(g_c)}\sin g_c},
\]
and the Baldi action equation for \(\eta\).  Expanding these identities to
third order is algorithmic, but the presently pinned defining data specify
\(\gamma_c\) only through its \(c^{3/2}\) term (and \(h(I)\) only through the
displayed order).  Writing the next coefficient as
\(\gamma_c=c/4+g_{3/2}(\theta)c^{3/2}+g_2(\theta)c^2+\cdots\), the inverse
recurrence shows that \(g_2\) enters \(\rho_3,\zeta_3,\eta_3\), hence enters the
single \(C_2\) trace integral linearly.  It is therefore not legitimate to
replace this datum by a printed lower-order Taylor truncation.

The replayable exact calculation (``verify_c2_recurrence.py``) evaluates the
entire order-prescribed double-Dyson contribution and returns
\[
 \int_{0<u<t<2\pi}\!\!\operatorname{tr}
 [U_0(2\pi-t)C_1(t)U_0(t-u)C_1(u)U_0(u)]\,du\,dt=2\pi^2.
\]
The complete order-(I) coefficient is this fixed \(2\pi^2\) plus the single
ordered integral of \(C_2[g_2]\), together with the known period/action
conversion term.  Since \(g_2\) (equivalently the cubic action-inversion datum)
is absent from the frozen source receipt, the first nonzero discriminant
coefficient is underdetermined, not numerically irreducible.  Supplying that
next source equation is the precise remaining construction.

## Sol-High correction audit and completed source recurrence

This section supersedes the two preceding small-shell conclusions while
preserving them as audit history.  Baldi (4.29)--(4.33) determines every
`W_j` recursively.  In particular (4.31) gives

    W2 = sqrt(2)*(80 sin(sigma)^6-40 sin(sigma)^4
                  -3 sin(sigma)^2-11)/128,

so the claimed free `g_2` datum does not exist.  The full-`mu` repair in
(16) is algebraically correct, but it did not validate the later displayed
physical coefficient: that coefficient omitted rotation of the cylindrical
orthonormal frame.  The old `C1` and its matrix-valued first Dyson integral
are therefore withdrawn as physical statements.

Put `epsilon=sqrt(2I)` and assume the fixed regular shell lies in a flat
cutoff plateau, so `omega(K)=1`.  From (3.39), (3.48), (4.39), and (4.31),
coefficient comparison in `F_c(g_c(sigma))=sigma` gives

    g_c(sigma)=sigma+epsilon*g1(sigma)+epsilon^2*g2(sigma)+O(epsilon^3),
    g1=cos(sigma)^3-9 cos(sigma)/4+5/4,

with `g2` printed exactly by `verify_c2_source_recurrence.py`.  Substitution
in `sqrt(2 gamma)=sqrt(c) w`, `c=h(I)=2 epsilon^2+O(epsilon^6)`, and Baldi
(3.74) supplies the physical trajectory through cubic order.  Gavrilov's
defining ODE/PDE (4.2)--(4.4) supplies `alpha_2` through degree four and
`H(s)=4s-21s^2/2+39s^3/32+...`.  Differentiating the physical velocity,
including cylindrical-frame angular velocity, and using
`grad(sigma)||(z_I,-rho_I)` yields

    C0 = [[0,0],[-1/sqrt(2),0]],
    C1 = [[cos(sigma)-3 cos(3 sigma)/2, sqrt(2) sin(sigma)],
          [sin(sigma)/(4 sqrt(2))+3 sin(3 sigma)/(4 sqrt(2)),
           -cos(sigma)]],

and the source-derived `C2` printed by the verifier.  The resonant ratio is

    m/n = Omega1'(I)/Omega2'(I)
        = (3195/256) I^(3/2) [1+O(I)],

so its physical azimuthal covector correction is relative `O(epsilon^4)` and
does not alter `C0,C1,C2`.  Also
`Omega1=K'(I)=1+(3195/1024)I^2+O(I^3)` on the flat plateau; consequently the
return-period correction first enters beyond order `I`, not in `C2`.

For `U0(t)=I+t C0`, exact ordered integration now gives

    single-C2 trace  = 9*pi^2,
    double-C1 trace  = 2*pi^2,
    epsilon^2 trace  = 11*pi^2.

Analyticity of the source functions and the compact sigma-circle makes the
coefficient remainder uniform.  Hence, along the derivative-resonant small
shell sequence,

    tr M = 2 + 11*pi^2*epsilon^2 + O(epsilon^3)
         = 2 + 22*pi^2*I + O(I^(3/2)),
    det M = 1,
    Delta = (tr M)^2-4 = 88*pi^2*I + O(I^(3/2)) > 0

for sufficiently small positive resonant `I`.  Thus the actual returned
polarization cocycle is hyperbolic on that sequence.  This is a linear
geometric-optics/accessible-packet mechanism; it does not by itself prove
nonlinear instability or create a particle/quantum interpretation.

## Bounded 0038 correction: cubic spatial jet in `C2[1,0]`

The source-recurrence oracle originally applied `mp(z,3)` to the physical
velocity-gradient jet before substituting `X,Z=O(epsilon)`.  That truncation
retains spatial degree at most two.  In the `u_phi,r` entry, however,
`H'(c)/sqrt(H(c))=O(epsilon^-1)`, so spatial degree three contributes at order
`epsilon^2` and belongs to `C2`.  The corrected oracle applies `mp(z,4)` and
proves exactly that the restored increment is

    Delta C2[1,0]
      =sqrt(2)*(3*cos(sigma)^4/2
                 -13*cos(sigma)^2/8-9/16).                 (20)

No other displayed `C0`, `C1`, or `C2` entry changes.  In particular the
single-`C2` trace functional depends on `C2[0,0]`, `C2[0,1]`, and `C2[1,1]`,
not `C2[1,0]`.  Therefore its value remains `9*pi^2`; the double-`C1` term
remains `2*pi^2`, and the established conclusions

    tr M=2+22*pi^2*I+O(I^(3/2)),
    Delta=88*pi^2*I+O(I^(3/2))

are unchanged.  The first corrected execution is preserved in
`c2-corrected-spatial-jet.*`; the original successful receipt remains as
append-only provenance.

The same corrected oracle was then executed with the repository interpreter
`/home/dan/substrate-framework/.venv/bin/python`; its command, identical
stdout, empty stderr, and exit-zero files are named
`c2-repository-interpreter.*`.  All pre/post and execution hashes are pinned
in `c2-spatial-jet-correction-receipt.md`.

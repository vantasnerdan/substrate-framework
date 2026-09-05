# Actual affine axial second jet: radial control, whole-history projection and current

This continues the same stationary material-tag construction. It concerns
the exact affine Euler comparison, not yet the nonlinear periodic cell.
All macro-K statements below are axial; the generic common laboratory K
and nonlinear-cell transfer remain named subsequent constructions.

## 1. The measured inertia is independent of the radial mode profile

Replace the Gaussian vorticity radial transform by any sufficiently
decaying smooth F(k), vanishing to adequate order at k=0, retaining the
same fourth angular harmonic. The exact full pressure inverse gives

    phi_2=tau integral F(k)J2(kr)dk/k,
    phi_-2=tau² phi_2,
    phi_0=-tau² integral F(k)J0(kr)dk/k.

The Bessel recurrence proves

    partial_r phi_0=tau(partial_r phi_2+2phi_2/r).     (1)

For the SAME stationary positive tag define
L_chi[F]=integral chi(x)[phi_2+x partial_x phi_2]dx.
The exact physical rows are

    theta0=-i pi(1+tau²/3)L_chi[F]/(2Omega Q),
    S0=2pi rho tau L_chi[F],
    G0=-i S0/(4Omega)=Delta theta0.                   (2)

Thus every coherent radial superposition keeps the same positive
Delta=rho Q tau/(1+tau²/3). Signed field amplitudes do not change the
positive tag density or introduce signed inertia. Choose L_chi[F]!=0.

## 2. The actual axial Euler correction and its radial homogeneity

Remove only the explicitly retained axial advection phase exp(-iKW0 t)
when writing the coefficient equations; the original field contains it.
Actual material-tag evaluation along X=X0+W0t cancels this factor in the
registered material Fourier amplitude. The laboratory energy still keeps
the axial Noether connection derived in0205.

For a full Kelvin Fourier characteristic with normal wave q,

    q_dot=-A_core^T q,
    omega_vector=(zeta,K b)+O(K²),
    b_dot=A_core b-h zeta Jq/|q|², h=1+Omega².

The full vorticity equation is
omega_dot=Du_aff omega+hK (k cross omega)/|k|², k=(K,q).
This follows from the complete Euler velocity
w=i k cross omega/|k|² and the actual constant base vorticity -h e_X.
In particular the pressure contribution has not been suppressed.
Writing I(t)=integral_0^t |q(s)|^(-2)ds gives the exact first correction

    b(t)=exp(tA_core)b(0)-h zeta Jq(t)I(t).           (3)

Its contribution to zeta at second order includes
-h²K² zeta integral_0^t I(s)ds. This is an actual branch-memory term,
not a chosen dispersion coefficient. The free initial b(0) and the
divergence condition q.b(0)=-zeta are retained.

The normal Fourier radius k in strain coordinates is invariant along
this characteristic. Start with the nonresonant displacement preparation
of the companion proof. Its initially invisible resonant row may be zero,
or the homogeneity-preserving phase preparation
s_4(0)=i c(-Delta_R)^(-1)zeta_4 described below. All initial normal
velocity and displacement amplitudes then
have radial homogeneity F(k)/k. Exact P_K completion and the equations
above imply that their intrinsic normal second-K coefficients have
radial factor F(k)/k³ times an angular/time function. The first-K axial
coefficients have factor F(k)/k². This is a homogeneity statement about
the complete Fourier solution, including its resonant time polynomials,
not a truncation to one angular mode.

Every stationary-tag NORMAL first-moment observation pairs these fields
against a radial Fourier derivative chi_hat'(k) times an angular linear
factor. Hence all intrinsic second-K angle and G terms are proportional
to the single radial functional L_chi[F/k²]. The same is true of the
spin connection. Indeed, with Fchi_psi=chi,

    S-G_t=2iK rho integral Fchi xi_X,
    Fchi_hat=-Omega ell² chi_hat'(k)/k.

Thus its intrinsic second coefficient uses that same inverse-wave
functional. The literal current, rather than a frequency average, makes
this cancellation check possible.

## 3. Cancel the intrinsic second jet by actual radial preparation

Choose F with

    L_chi[F]=1, L_chi[F/k²]=0.                        (4)

This cancels the complete intrinsic observed second jet just identified,
including periodic and secular angular terms. It does not set the
unobserved Euler state to zero.

There is an explicit nondegenerate physical moment construction. For
chi(x)=exp(-b x), b>0, and F_s(k)=k^4 exp(-s k²/2), s>0, direct Hankel
integration gives

    L_chi[F_s]=16b tau/(1+2bs)^3,
    L_chi[F_s/k²]=2tau/(1+2bs)^2.                     (5)

Their ratio is (1+2bs)/(8b), strictly increasing in s. Two distinct
positive s therefore solve(4) by signed coherent amplitudes. The
response profile and its pressure field are actual functions, not two
oscillators assigned different frequencies: both retain the exact same
fourth-harmonic affine Euler dynamics.

If a compact tag is required, choose a sufficiently remote smooth radial
cutoff of chi. The actual two-by-two moment matrix remains nonsingular
by dominated convergence. Solve its ACTUAL integrals, not the limiting
Gaussian numbers, to impose(4) exactly. Equation(1) and the inertia in(2)
remain exact for that fixed compact tag.

## 4. A real first-K axial current supplies positive curvature

Add an independent first-K axial Euler velocity iK G, with radial
fourth-harmonic transform G=i g F, g>0. Real ±K phases are paired by
conjugation. This is actual initial velocity/circulation data. It is not
an imposed change of the generator's clock. Its Eulerian axial amplitude
is advected, and the complete Fourier equations give the additional
normal velocity coefficient

    w_perp,2=e^(i4Omega t)[h t J grad phi_G+grad phi_G]. (6)

The second term supplies the divergence required by the axial velocity.
It cannot be dropped as a longitudinal pressure term.

For ANY normal displacement define its actual first-moment matrix
M=integral chi r xi^T dA. Since the reference chi is stationary,

    M_t=A_core M+M A_core^T+integral chi r w^T dA.    (7)

This closes the observation exactly, not the full Euler field. In units
where I=L_chi[F], the defining full-pressure moment integrals give

    Kmat=integral chi r (grad phi_F)^T
       =[2pi I/(1+Omega)^2] [[1,i],[iOmega²,-Omega²]],
    Vmat=Kmat J^T.

Let R_omega be the exact matrix inverse of
i(4Omega)M-A_core M-M A_core^T. It is nonsingular on these first
moments (their homogeneous frequencies are0,±2Omega). Put

    M0=R_omega Vmat,
    M1=R_omega(Kmat-h M0).

The particular first moments in(6) are exactly
exp(i4Omega t)[h t M0+M1]. All homogeneous first-moment terms are
retained until the actual phase preparation below removes them.

## 5. Whole-history phase projection, not a clock replacement

Let P=2pi/Omega. Draw a real tau from the FOUR-fold convolution of
the centered uniform probability on[-P/2,P/2]. Prepare each actual
Euler/Lin initial column by its TRUE time-tau propagator, followed by
the fixed real two-phase amplitude rotation exp(-i4Omega tau). The
stationary reference tag chi is unchanged by the reference flow. This
is a positive declared family of actual prepared histories.

For every nonzero integer n, this probability satisfies

    E[tau^j exp(i n Omega tau)]=0, j=0,1,2,3.         (8)

Its characteristic function is sinc(n pi)^4, with zeros of order4;
differentiation proves(8). The affine Fourier/Lin second jet is a
polynomial of time degree at most3 with periodic coefficients of
period P. Consequently the actual phase-correlated average removes
all nonresonant periodic first-moment terms. It retains the polynomial
resonant coefficients, including any real branch-memory term; those
were canceled by(4), not wished away by phase averaging.

Average the ACTUAL phase actions on the common initial columns first.
The Euler/Lin propagator preserves their full phase, and every real
two-phase rotation preserves its two-by-two skew form. Thus the averaged
phase is that same form, not an average of squared frequencies. The
actual mean angle and rate then define the common initial physical chart.

After that initial normalization, (6),(7) give

    theta_K(t)=theta_initial(t)
                 with frequency Gamma(K)=4Omega+h g K²+O(K³),
    Gamma(K)^2=16Omega²+8Omega h g K²+O(K³).          (9)

The curvature is strictly positive. The finite ensemble has not changed
the actual field, its laboratory K or its measured tag angle.

## 6. The derived gradient inertia is part of the physical current

Define theta_i=(Mi_12+Mi_21)/Q and G_i=rho(Mi_12-Mi_21). Then
G0=Delta theta0, but the COMPLETE constant forcing row gives

    D=G1-Delta theta1
     =pi I rho (Omega-1)(2Omega^4+5Omega^3+4Omega²+5Omega+2)
       /[24Omega²(Omega+1)(Omega²+Omega+1)].           (10)

It is nonzero. This is not discarded; after the same initial angle
normalization it is the physical gradient inertia

    Delta_K=Delta+Delta2 K²,
    Delta2/Delta=g(2Omega^4+5Omega^3+4Omega²+5Omega+2)
                         /[12Omega(Omega²+Omega+1)]>0. (11)

The free axial displacement in(6) has fourth angular harmonic, so the
primitive-weight gyro-current integral is zero. The intrinsic current
was canceled by(4). Therefore the actual rows satisfy

    G_K=Delta_K theta_K+O(K³),
    S_K=Delta_K (theta_K)_t+O(K³).                    (12)

The positive equation curvature in(9) coexists with this gradient mass.
Both enter the action and energy. An equation-level order reduction may
use the leading optical equation; it may not silently erase(11) from
the inherited physical phase or current.

An actual on-site preparation already normalizes the leading phase.
Take s_4(0)=i c(-Delta_R)^(-1)zeta_4, where the inverse has radial
Fourier multiplier1/k². Its full displacement has radial factor F/k,
so the cancellation in(4) remains intact. Its initial and transported
homogeneous fourth harmonic is invisible to theta,G,S at K=0. Yet the
same integration by parts as the companion proof gives

    beta=rho c integral |zeta_hat_4|²/k² dA_hat.     (13)

This integral is finite and strictly positive for any nonzero selected
F. For F=k^4 exp(-k²/2), it is6pi zeta0² ell². Choosing c so that
beta=Delta(4Omega)C² supplies the actual leading positive physical
phase without a remote return or altered clock. The secular fourth-
harmonic state remains present. A further order-K² change of c has
invisible leading observation and can adjust the scalar even phase jet;
its coefficient is nonzero by(13). More generally the complete joint
cross forms, physical chart and energy can use0210 and0205's actual
returns at their own stated domain and error scope.

In particular the complete Jacobi energy is
rho/2 integral[|xi_t|²-|u_aff.grad xi|²+xi.Hess(p_aff)xi], with its
constraint and axial Noether connection, not an assigned oscillator
energy. It need not equal the energy inferred from the measured scalar
action merely because(13) matches phase. That remaining finite quadratic
form requires an actual energy return or an explicit evaluation. The
scalar action inferred after the actual phase normalization has mass
Delta_K and stiffness Delta_K Gamma(K)^2 through the shown spatial order;
this statement alone makes no equality claim for the full Jacobi energy.

At K=0 that difference can be evaluated exactly. Let N_F denote the
positive norm in(13), without rho c. The actual fourth-harmonic pressure
forcing is phi_4=-ell²(-Delta_R)^(-1)zeta_4, since the zero Fourier
coefficient of1/g is1. In the full stationary Euler/Lin phase space the
two real columns y1,y2 and their actual pure-label forcing z1,z2 obey

    L y1=-omega y2+z1, L y2=omega y1+z2,
    L z1=-omega z2,    L z2=omega z1, omega=4Omega.

Integration by parts gives Omega(y1,z1)=Omega(y2,z2)=-a,
a=rho ell² N_F/2>0; the crossed y,z rows vanish and Omega(z1,z2)=0.
The full Hamiltonian bilinear form is -Omega(.,L.), so its restriction
to the initial y columns is

    H_yy=(omega beta+a) Id.                          (14)

The factor1/2 in the actual quadratic energy then gives an excess
a|input|²/2 over the physical scalar oscillator. This is a genuine
positive contribution of the unobserved resonant state. It is not removed
by the stationary tag or the whole-history phase filter, both of which
preserve the complete conserved energy. Actual negative-energy passive
configuration returns are the available method repair; their own action,
cross rows and finite-K costs must accompany any implementation. Formula
(14) makes the full-energy obligation quantitative rather than ceremonial.

There is an exact K=0 repair on this same affine field. Choose two
disjoint smooth annular supports outside the fixed compact tag. On each
put one actual axial passive column g_i=G_i(r)cos(Ntheta),
h_i=T^(-1)g_i=-G_i(r)sin(Ntheta)/(N Omega), with
w_i=e^(-tT)g_i e_X and xi_i=e^(-tT)(h_i+t g_i)e_X. Their initial
cotangent is rho g_i e_X and T h_i=g_i. Therefore their full energy
matrix is diagonal with entry -rho||g_i||². Their phase matrix vanishes
by sine/cosine parity and disjoint support. All phase and energy cross
rows with the horizontal optical columns vanish by physical component
orthogonality at K=0. These are exact full-field statements.

Set ||g_i||²=a/rho for each column. Their added energy matrix is -a Id,
canceling(14)'s excess without changing beta. The annuli are invariant
under the actual affine flow and avoid the tag, so every tagged angle,
spin, G and shape row remains unchanged for all finite times. Whole-
history propagation and real column rotation preserve this phase-null
isotropic energy correction. Thus leading positive physical angle,
spin, phase AND complete energy can be joined on the affine comparison.
No measured normalization has been changed. The corresponding finite-K
energy and cross rows remain retained: K=0 orthogonality alone does not
set their second jet to zero. The finite homogeneous cross controls and
own-scale return estimates of0205/0210 identify the next repair, rather
than asserting that repair has already been applied here.

## 7. Scope still to complete

Equations(4)–(12) are a constructive observed axial second-jet interface
for actual affine Euler histories, with a positive finite phase ensemble
and literal currents. Full joint action requires applying the actual
phase/energy returns at their stated error scope, not substituting them
as exact identities. The nonlinear fixed-cell transfer, generic K
geometry, and optical density as ell changes are still separate work.
In particular no nonzero j limit is hidden: Delta scales as rho ell^4
on this fixed cell. Neither compact Euclidean EPS geometry nor parent
completion follows from this affine axial interface.

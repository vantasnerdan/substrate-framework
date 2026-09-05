# Complete correlated whole-law current on the actual two-wave field

## 1. Actual Euler field, common laboratory input and preparation

Use0179's exact stationary periodic field, in its fixed body coordinates:

    u=(psi,sin Z,-sin Y), psi=cos Y+cos Z,
    alpha=cos Y-cos Z, A=sin Z partial_Y-sin Y partial_Z,
    curl u=-u, p=-|u|²/2.

The same unit laboratory kappa and unit transverse displacement D are
used in every whole-field realization. In body coordinates their joint
law is the Haar law of an orthonormal pair. Let

    d=kappa_Y D_Y-kappa_Z D_Z,
    e=kappa_Y D_Y+kappa_Z D_Z=-kappa_X D_X.

The actual translation-leading Kelvin/velocity expansion from0179 is
w_D=T_D+ik(q_D+z)+..., with q_D fixed and z_t=Lz+dF0. No stationary
z is required. Allow arbitrary measurable orientation-correlated smooth
initial z0, with finite averaged energy/enstrophy; it need not factor as
d times one field. Each such choice has its actual initial circulation.
The material preparation chi0=0 and mean-zero microscopic velocity retain
the exact common initial rho J form as in0179. The current of q_D is
constant in time and retained below as part of the initial constant.

Any x-dependent microscopic correction can first be averaged in the
translation-invariant x direction: Euler's linear operator, forcing and
the physical current commute with that average. Other axial Fourier
sectors do not feed the measured mean. Thus the following planar fields
describe the COMPLETE observed part of that three-dimensional correction.

Write z=(b,-phi_Z,phi_Y), eta=Delta phi, H=-Delta, B=H-1,
and r=b+eta=b-Hphi. The exact equations are

    H phi_t=-A Bphi+d Aalpha,
    r_t+A r=-(d/2)Aalpha.                              (1)

The second source is computed from F0=P(u_Zwave cross u_Ywave):
F0_x=-sin Y sin Z, curl_x F0=2 sin Y sin Z. Their sum equals
sin Y sin Z=-Aalpha/2. It is not silently declared passive. The
homogeneous difference r_h is actual passively transported axial data.

## 2. Combine the literal axial and planar current before averaging

The scalar physical stress row in the direction D is <c.z>, where
c=(kappa.u)D+(D.u)kappa. Its axial coefficient is

    c_x=2kappa_X D_X psi
       +(kappa_Y D_X+kappa_X D_Y)sin Z
       -(kappa_Z D_X+kappa_X D_Z)sin Y.

The planar contribution is <l_c phi>, with

    l_c=partial_Z c_Y-partial_Y c_Z
       =2kappa_Y D_Y cos Z+2kappa_Z D_Z cos Y
        -(kappa_X D_Y+D_X kappa_Y)sin Z
        +(kappa_X D_Z+D_X kappa_Z)sin Y.

Since b=r+Hphi and Hc_x=c_x, the COMPLETE current is

    <c.z>=<c_x r>+<(Hc_x+l_c)phi>
          =<c_x r>-e<psi phi>-d<alpha phi>.              (2)

All mixed sine rows cancel exactly. The last equality uses only the
physical transverse constraint kappa.D=0. Dropping the axial-vorticity
connection b=r+Hphi would miss this cancellation and falsely leave several
independent planar output sectors. This is the representation repair that
extends the separable class of0179 to arbitrary correlated data.

The psi moment is constant for EVERY solution(1), because Hpsi=psi,
Apsi=0 and A is skew. Correlating it with e changes only a constant.
No per-orientation stress stationarity has been imposed.

## 3. Whole-law correlations reduce to one actual forced range state

The actual orthonormal-pair fourth moment gives

    E[d²]=1/5, E[d c_x]=0.                             (3)

Put phibar=E[d phi], which is a well-defined finite-enstrophy field by
Cauchy--Schwarz. Since the same body-frame operator acts in every
realization, (1) gives

    H phibar_t=-A Bphibar+(1/5)Aalpha.                  (4)

This is an actual weighted average of prepared Euler fields, not a chosen
scalar cancellation function. Its initial value can contain ALL the
correlated planar controls; no first-shell restriction remains.

The forced part of r is d times one fixed solution of(1). Its contribution
to E<c_x r> vanishes by(3). The remaining correlated passive current

    P_h(t)=E<c_x,exp(-tA)r_h0>

is bounded for all real t by unitarity and finite averaged second moments.
It is real analytic: c_x has finite Fourier support and is an analytic
vector of A; move the group to this test row, not to the arbitrary data.

Consequently the complete whole-law physical D stress is

    R_D,iso(t)=C0-<alpha,phibar(t)>+P_h(t),             (5)

where C0 includes q_D and the conserved psi sector. This identity covers
general correlated planar and passive axial preparations in the stated
translation-leading class, not merely the separable forced ansatz.

## 4. Its exact range-energy implication

On the complement of B's first-shell kernel define

    Zbar=(HB)^(1/2)phibar, Q=(B/H)^(1/2),
    G=-QAQ, f=(1/5)Q Aalpha.

The actual full-pressure range equation and observed-current derivative are

    Zbar_t=G Zbar+f,
    d[-<alpha,phibar>]/dt=-5<f,Zbar>
                           =-5 d(||Zbar||²/2)/dt.      (6)

Thus

    R_D,iso(t)=C1-(5/2)||Zbar(t)||²+P_h(t).             (7)

The first-shell phi components are retained in the current; they were
not set to zero because Zbar excludes them. Equation(6) is what relates
that physical current to the range norm. The unitary group/domain and
analytic-vector estimate for the finite Fourier forcing f are exactly
those proved in0179, with its changed positive factor retained.

If the actual scalar stress were bounded on all real times, (7) would
make the affine-unitary range trajectory bounded. Its Cesaro averages
then yield GZ_*+f=0, as in0179's closed-graph argument. The even-even /
odd-odd parity split gives

    A[QZ_*,even-even-alpha/5]=0.

Every L2 first integral on the connected regular streamlines is a
function of psi. Pairing with alpha gives zero for that function and
QZ_*, but -<alpha²>/5=-1/5 for the displayed expression: contradiction.
This is an actual whole-law observed-current mechanism. No sign has been
assigned to the indefinite0180 Jacobi form, nor has every orientation
been required stationary or even to have constant stress.

Time reversal and whole-field reflection do not evade the result. Each
time-reversed range generator remains skew and carries the same positive
energy-current sign; finite positive weights give a sum of nonnegative
range energies. The full reflected field and polar inputs give the same
scalar contraction by Euler covariance. No independent reflection of just
the angle, tag or velocity is used.

## 5. Finite positive physical closure and honest scope

The physical stress output and the passive remainder are analytic in time
for arbitrary finite-energy prepared range data, by the analytic-vector
argument. Hence equality on an open interval to a bounded finite sum of
constant-frequency stable optical responses extends to all real time,
contradicting(7). At the second spatial jet, a positive constant-coefficient
finite Cosserat system has exactly such a bounded displacement-column
acceleration: acoustic displacement contributes its constant restoring
row, and any nonzero optical gaps contribute sines/cosines. Eliminating
Phi is allowed and creates this bounded memory; no isolated-mean autonomy
is imposed as an extra condition.

For example in a transverse curl-helicity sector write U=D+k²U2+...,
Phi=k phi1+... and omega²=4alpha/j. The actual positive canonical
Cosserat pencil gives

    phi1_tt+omega² phi1=(omega²/2)D,
    rho U2_tt=-(mu+alpha)D+2alpha phi1
              =-mu D+2alpha[(phi1(0)-D/2)cos(omega t)
                                  +phi1_t(0)sin(omega t)/omega].

This includes arbitrary actual initial optical displacement/rate at this
order. Adding finitely many positive optical pairs or the retained
time-independent derivative mass map adds only finitely many bounded
oscillatory terms to this D acceleration row. It is that precise necessary
joint-response property, not a claim of zero mean memory, which conflicts
with(7) for the present microscopic preparation.

Thus this translation-leading two-wave preparation class cannot produce
that EXACT finite positive physical stress response by arbitrary correlated
first-cell velocity controls and passive returns. This is not a theorem
about different leading microstates, different Euler fields, an arbitrary
moving physical observation, or controlled finite-window approximation.
It does not refute the same-EPS parent: the two-wave field is the explicit
comparison that generated the present correction route, not C-CST-011's
full inserted geometry. The finite physical angle/spin/current construction
remains the goal, and new field or leading-state routes remain active.

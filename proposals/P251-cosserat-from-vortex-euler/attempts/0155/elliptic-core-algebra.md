# Full elliptic strain: principal symbol, physical clocks and candidate repair

This is an analytic construction receipt, not a completed0155 theorem.
The parent validated263/12 before0151 was opened. Coordinates are
normal (a,b)=(y-pi/2,z), longitudinal s=x, positively oriented.

## 1. Exact strain and biorthogonal variables

For0<d<1/10 let Omega=sqrt(d), h=1+d and

    A=[[0,1],[-d,0]], S=diag(d^(1/4),d^(-1/4)),
    X=S(a,b), G=S²=diag(Omega,1/Omega), J=[[0,-1],[1,0]].

Then SAS^-1=-Omega J, but A is NOT antisymmetric in the physical
Euclidean metric. The exact linear core has A²=-dI and isotropic
physical pressure Hessian dI. Its axial speed expansion is

    W=1+d-(Omega/2)|X|²+(X1^4+d X2^4)/24+...,

where the quartic formula follows from d a^4+b^4 after the S map.
The full nonlinear field remains the exact0151 periodic Euler field,
not its linear jet.

For a velocity perturbation put V=S v_perp. Its principal normal
Euler operator is Omega partial_theta V+Omega JV. The pressure
gradient is G grad_X P, and incompressibility is
div_X V+i k v_s=0. It is NOT legitimate to turn G into the identity.
On V=e_+ F_l, e_+=(1,i)/sqrt(2), F_l=e^(il theta)f(r), the
laboratory frequency relative to kW(0) is

    omega_0=(1-l)Omega=(2-m)Omega, m=l+1.                    (1)

The opposite polarization and its different angular sectors must
be kept when solving subsequent orders. The leading Kelvin map is

    eta=S xi_perp=i G V/h,
    beta_fiber=-(rho/(2h))integral |F_l|² dX <0               (2)

in the Re/Im convention. This is the actual negative-curl KKS,
not a sign selected from an eigenfrequency label.

## 2. Pressure-resolved first correction

For carrier k=p>0, the transverse projection of the leading pressure
is obtained from the full longitudinal equation, not a local pressure
ansatz. Its resonant e_+ component has coefficient
(trG)/4 times Delta_X. Balance with k(W-W0) gives

    ell^4=(trG)/p³, delta_e=sqrt((trG)/p), cD=Omega delta_e,
    [-Delta_R+R²]F_l=2(2n+m)F_l,
    F_l=R^(m-1)e^(-R²/2)L_n^(m-1)(R²),
    omega=(2-m)Omega-(2n+m)cD+O_d(delta_e²).                 (3)

This leading equation does not alone establish an all-order packet.
Its helicity/angular couplings and remainder are the next construction
step. In particular frequency(1) must not be replaced by particle
Omega or by the circular-core value2Omega without an actual observer.

## 3. A genuine simple-marker normalization obstruction

For m=2, choose a stationary elliptic material marker with normal
density chi(|X|), centered covariance proportional to S^-2. Its
physical Euclidean covariance angle is well defined for d!=1.
For radial profile f, set

    C0=integral chi r² dX,
    AF=integral chi r² f dX,
    N=integral_whole r² f² dX.

Direct variation of position and velocity gives, up to the shared
mode amplitude convention,

    angle row c=-AF/[sqrt(2)(1-d)C0],
    physical spin row i rho(1-d)AF/[sqrt(2)h],
    beta=-rho N/(2h),
    S/Pi=AF²/(C0 N)<=1.                                    (4)

The inequality follows from Cauchy--Schwarz and0<=chi<=1.
Equality requires f constant on the observed material and no
unobserved f² action. Thus eta=1 saturation is not available for a
nonconstant localized mode in THIS stationary radial-marker class.
An eta<1 target has a different existence question; it is not obtained
by claiming equality in (4).

Also, a fast painted-ellipse angle cannot replace this slow branch:
at the exact linear core pressure Hessian dI is radial, so leading
total material spin is constant. An angle winding at2Omega from
background paint alone does not supply matching canonical momentum.

## 4. Registered odd-lobed physical moments: a new candidate

The preceding failure motivates an actual m-lobed physical moment
with odd m>=3. The unperturbed isotropic elliptic marker has no odd
physical moment. Add a bounded, radially controlled cos(m theta)
lobe marking, keeping the material density nonnegative and centered.
Let z_phys=a+i b=A_+ Z+A_- conjugate(Z),

    A_+=(d^(-1/4)+d^(1/4))/2,
    A_-=(d^(-1/4)-d^(1/4))/2, r_m=(A_-/A_+)^m<1.

The actual Euclidean mth reference moment on the linear core is

    Q_m(t)=Q0[A_+^m e^(-imOmega t)+A_-^m e^(imOmega t)].       (5)

It is nonzero, with a strictly positive physical chart margin.
The leading velocity/displacement polarization in (2), integrated
against the isotropic part of the tag, gives the actual angle-row
complex factor, apart from a nonzero real normalization,

    F(t)=e^(i omega_0 t)
      [(A_+^(2m)+A_-^(2m))e^(imOmega t)
           +2A_+^m A_-^m e^(-imOmega t)]/|Q_m(t)/Q0|².        (6)

The physical angle is (1/m) arg Q_m after varying the ACTUAL
moment; (6) is that linear observation row, not a chosen Floquet
logarithm. Its phase-rate minimum is exactly

    min gamma=Omega[2-4m r_m/(1+r_m)²].                     (7)

When2m r_m/(1+r_m)²<1 the whole linear-core observation has gamma>0,
so (2) gives positive physical scalar mass -beta/(gamma |c|²).
For any fixed d>0 sufficiently large odd m satisfies this inequality.
This establishes the principal physical sign, not its nonlinear-
background or finite-packet error budget.

The leading total material spin for this high-carrier sector is
constant, as above. With a pure m-lobe marker it vanishes by angular
selection: its principal spin density has harmonic m-2, whereas
the pressure torque first supplies the m harmonic. The actual
nonzero oscillating spin therefore starts at the pressure order
delta_e. Radial reference-moment cancellation, analogous to0147,
is the candidate to match that SMALL physical spin to the same
action; no leading-order kinematic spin is invented.

For the isotropic part of a fixed radial marker chi proportional
to r², the Laguerre Laplace derivative gives

    (mean x)/2=(2n+m)+m/(2n+m),
    gamma=2Omega+[m/(2n+m)]cD
       +counter-rotating physical terms+O_d(delta_e²).        (8)

Its leading pressure contribution has positive carrier curvature.
The counter-rotating terms are fixed by (5)--(6), not erased. Their
size, the large-m moment conditioning, finite carrier derivatives,
the small reference moment and the global Euler energy constant
must be compared before choosing a joint d,m,p hierarchy. Merely
letting m grow until(7) holds is not that comparison.

## 5. Current constructive boundary

Equations(1)--(8) register the actual physical candidate and the
mechanisms the next proof must retain. The simple stationary m2
eta=1 route is refuted at its leading normalization; the odd-lobed
route has a positive principal physical action with explicit strain
and clock. Its pressure-order spin and full finite-action/parameter
hierarchy are active constructions. No parent completion, exact
isolated spectrum, or absolute-director identification is claimed.

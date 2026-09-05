# Fixed non-small eigenvalue: an actual robust Beltrami torus

This is a new construction, not an extension of EPS §7's small-eigenvalue
asymptotics. The registered axisymmetric candidate is established below. It
produces an **unknotted** invariant tube, an elliptic core, and a physical
gyro-dominant Kelvin amplitude on that same core. It does not prove the
fixed-eigenvalue normal form for every prescribed knotted centerline, or
identify scalar canonical inertia with a material parcel's spin inertia.

All limits are parameter choices in an existence proof. First fix a nonzero
physical eigenvalue lambda and U>0. By reversing orientation it suffices to
write lambda>0. Choose a positive minor radius a with delta=lambda*a small.
Next choose a sufficiently large but finite major radius R, with a<R/2.
Finally fix every geometric parameter before global approximation and the
high-carrier finite-time Euler construction. No eigenvalue tends to zero.

## 1. Exact elliptic construction

In cylindrical coordinates (r,theta,z), set x=r-R and let D_a be the disk
x^2+z^2<a^2. Write kappa=1/R. Solve the scalar Dirichlet problem

    [Delta_(x,z) - kappa/(1+kappa*x) d_x + lambda^2] phi = -lambda*U,
    phi|boundary(D_a)=0.                                      (1)

The Dirichlet Laplacian is invertible. On scaling the disk to unit size,
its perturbations have operator norms controlled by |kappa|a and delta^2
from C^(m+2,alpha)_0 to C^(m,alpha). A Neumann series gives a unique solution
for sufficiently small parameters, with analytic dependence on kappa in
these spaces. Elliptic bootstrapping gives all finite derivative norms.
The coefficients and boundary are analytic; analytic boundary regularity
gives an analytic extension across the circular boundary, as required by
the neighborhood version of the KAM and approximation theorems.

Define psi=R*phi, C=R*U and

    u_r     = -psi_z/r,
    u_theta = (lambda*psi+C)/r,
    u_z     = psi_r/r.                                      (2)

Direct cylindrical differentiation gives div u=0 and

    (curl u)_r=lambda*u_r,
    (curl u)_z=lambda*u_z,
    (curl u)_theta=-(psi_rr-psi_r/r+psi_zz)/r=lambda*u_theta.

Thus (2) is an exact smooth stationary Euler field, with
p=-rho*|u|^2/2 plus a constant. Also u·grad psi=0 exactly. The boundary
psi=0 is an exact invariant torus, rather than an approximately tangent
surface. The pressure and velocity are the full three-dimensional fields;
there is no filament cutoff or neglected mutual velocity in this assertion.

At kappa=0 the normalized solution is explicit:

    phi_0(s)= U/lambda [J0(lambda*s)/J0(delta)-1],
    u_ax(s)= U J0(lambda*s)/J0(delta),
    u_pol(s)=U J1(lambda*s)/J0(delta),  s=sqrt(x^2+z^2).       (3)

The quotient at s=0 is defined by its analytic limit. For small delta,
phi_0=U*lambda*(a^2-s^2)/4+O(U*lambda^3*a^4). Its Hessian is uniformly
negative definite on the closed disk. Small kappa preserves strict
concavity, a unique interior critical point q_*(kappa), positive u_theta,
and a foliation by nested analytic level circles outside q_*. Revolving
q_* gives an exact circular periodic core, not merely an approximate
streamline. The transverse linear particle matrix is elliptic: in the
oriented normal plane it is -J Hess(phi)/(1+kappa*x_*) and has determinant
det Hess(phi)/(1+kappa*x_*)^2>0 and trace zero.

## 2. Actual return angle, flux action, and twist

An axisymmetric orbit completes a poloidal circuit on each regular level
curve C_c={phi=c}. The toroidal angle advanced during that circuit is

    Delta_theta(c) = integral_Cc
         kappa*(U+lambda*phi)/[(1+kappa*x)*|grad phi|] dl.

Consequently the rotation angle of the toroidal-section return map is

    Theta(c)=4*pi^2/Delta_theta(c).                         (4)

This includes the true nonuniform transit time. A convenient area action is

    I(c)=(1/(2*pi)) integral_inside(Cc) u_theta dx dz,

oriented to increase outward. The return map preserves u_theta dx dz;
the canonical angle obtained from this measure is conjugate to rigid
rotation by Theta(I). For the straight limit define the reduced rotation
F_kappa(I)=Theta(I)/(2*pi*R). Equation (4), divided by kappa before the
limit, is a regular analytic contour integral. On a closed annulus about
the boundary, phi, its level curves, F_kappa, and their action derivatives
depend analytically on kappa. Their limits are

    F_0(s)=J1(lambda*s)/[s*J0(lambda*s)]
          =lambda/2 +lambda^3*s^2/16+lambda^5*s^4/96+O(s^6),
    dI_0/ds = U J0(lambda*s)*s/J0(delta),
    dF_0/dI = [lambda^3/8+O(lambda^5*s^2)] J0(delta)/[U J0(lambda*s)]. (5)

At the boundary this is lambda^3/(8U)+O(lambda^5*a^2), strictly positive
for small delta. Choose delta to make the error less than half the leading
term, then kappa sufficiently small to retain half that strict margin.
This is an analytic continuity argument on a fixed positive-radius
annulus, not a claim that an unspecified O(epsilon) error is dominated by
epsilon^2. Thus dTheta/dI=2*pi*R*dF_kappa/dI is nonzero. The normal torsion
in EPS Theorem 7.6 is a nonzero positive-coordinate multiple of this
action derivative: in action-angle variables the map is
(I,vartheta)->(I,vartheta+Theta(I)), whose off-diagonal derivative is
Theta'(I). Nonvanishing is coordinate invariant.

The core reduced rotation F_core(kappa) extends analytically to kappa=0
by the elliptic linearization, with F_core(0)=lambda/2. The boundary value
F_boundary(0)=J1(delta)/(a J0(delta)) is also positive. Since

    Theta_boundary(R)=2*pi*R*F_boundary(1/R),
    d_R Theta_boundary=2*pi*[F_boundary(kappa)-kappa F_boundary'(kappa)],

the boundary angle is strictly monotone for sufficiently large R. The
same assertion holds for the core angle. Exclude the discrete set where
Theta_core belongs to pi*Z; the remaining set contains open intervals.
Diophantine angles have full Lebesgue measure, and the boundary-angle map
is a local diffeomorphism. Therefore within any such interval there is a
finite R for which the boundary angle is Diophantine and the core is
strictly elliptic with multipliers different from both +1 and -1.
This is a geometric parameter selection, not a choice of Floquet logarithm.

## 3. Robust global realization

For this now fixed torus and field, EPS Theorem 7.6 applies directly:
analytic boundary conjugacy follows from the axisymmetric action-angle
coordinates, the selected angle is Diophantine, and (5) supplies its
nonzero twist. Nearby divergence-free return maps preserve nearby positive
flux measures. The explicit near-identity Moser measure identification in
EPS equations (7.26)-(7.32) conjugates these to the same measure before
Theorem 7.6 is applied. Merely claiming that the two flux measures are
identical would omit this step.

The thickened closed solid torus has connected complement. EPS Theorem
8.3, with this fixed nonzero lambda, approximates (2) to arbitrary C^k
accuracy by a global curl u=lambda*u field whose derivatives decay O(1/r).
Take the approximation error below the preceding KAM threshold and below
the strict elliptic-core implicit-function threshold. The resulting
global stationary Euler field has a nearby invariant torus and an
elliptic periodic core. All constants are finite at the fixed R,a,lambda;
no uniform-in-R approximation estimate is being assumed.

This uses the general KAM theorem, measure adjustment and global Runge
theorem in the source. It does not invoke Theorem 7.8's epsilon^3 twist
formula, nor Theorem 6.8 as if it already supplied KAM persistence.
The torus is unknotted. A source knotted local seed can additionally be
rescaled to the same lambda and placed disjointly before global
approximation; that gives a second, knotted tube, not a transfer of this
rotor to the knot. The latter identification is explicitly not made.

## 4. Physical, rather than torsion-only, Kelvin rotation

On the exact circular core use the Frenet frame with its actual geometric
lift fixed in physical space. Its torsion is zero. Axisymmetry makes the
core velocity-gradient and frame-connection matrices constant in that
frame. As kappa tends to zero at fixed a,lambda they tend to the straight
core values. Write U_c=U/J0(delta)>0 and Omega=lambda*U_c/2.
The particle transverse matrix tends to +Omega*J, while the stationary
periodic covector normalized by k·u=1 tends to t/U_c. The latter exists
and remains bounded because the transverse particle eigenvalues are
separated from zero in the constant-frame algebraic equation. Its tilt
is O(kappa/|lambda|), with a finite constant depending on fixed delta.

For true Euler geometric optics the material amplitude equation is

    a_dot = (Du - frame_connection) a
            -lambda/|k|^2 * k cross a,    k·a=0.          (6)

It is the equation derived in 0112 from the Euler amplitude, not a
particle variational equation substituted for Euler. At kappa=0 its
restriction is exactly

    a_dot = -Omega*J*a.                                  (7)

For finite sufficiently small kappa it is a real trace-zero 2x2 matrix
arbitrarily close to -Omega*J, hence elliptic with a definite signed
Hamiltonian of the same sign. Curvature contributions, including the
tilt of k, are retained as continuous matrix perturbations. The small
global approximation and the finite-time continuation of the frame and
covector retain this property over a fixed number of optical periods.
No circulation-frequency multiple or frame winding is added to (7).

The actual packet KKS sign is sign(lambda). In the physical registration
of 0114/0115 the comoving coefficient r is close to -Omega, so
-beta*r>0 for beta with sign(lambda). The moving-frame carrier terms
cancel by the exact pullback identity before this small energy is bounded.
The high-order divergence-free WKB and exact Kelvin initialization in
0112 then give genuine linearized Euler packets in this same tube, with
full Euler and material reconstruction error O(N^(-m-1)) on a fixed
finite interval. Select the global approximation and kappa first, then
the finite N needed for that error to be below the positive action and
physical moment margins. The near-identity rectifier of 0114 produces the
positive autonomous two-dimensional *packet* action on that interval;
physical observation rows retain their time-dependent frame factors.

Actual parcel moments use the transported tag and the full formula
S=rho integral_D [r cross xi_dot+2 xi cross u]. The local straight-limit
shear packet in 0109 therefore gives a nonzero, nonconstant moment, and
its strict endpoint change survives sufficiently small curvature/global
approximation and the high-order packet error. This is a transported
parcel inside the invariant tube; it is not the integrated spin of the
whole solid torus. Pressure traction on that parcel is part of the Euler
balance, not an omitted external torque. Full compact-packet returns can
cancel its spin outside the observed parcel.

## 5. Executed mechanical-inertia countercheck

Even at the exact straight-core leading limit, the physical shear
xi=-J*Phi*z is not a rigid parcel rotation. With u=Omega*t cross r,
the Kelvin equation gives Phi_dot=-Omega*J*Phi. For an isotropic small
parcel let j_D=rho integral_D z^2. Direct evaluation including the
moving boundary gives

    S_D=j_D*Phi_dot.

But the KKS density of the two shear-angle columns integrated over that
same parcel is 2*Omega*j_D. A scalar elimination in this two-angle
first-order system therefore gives I_canonical=2*j_D, not j_D. The
parcel-restricted density is not itself a closed global symplectic leaf;
the full packet also has its return/collar contributions. Neither this
factor of two nor those ambient terms may be dropped in identifying a
mechanical rotor. Making lambda dominate curvature does not remove it.

The failure-derived rigid-rotation candidate can also be evaluated
exactly: xi=Phi cross r has
v=(Phi_dot-Omega*t cross Phi) cross r. Its fixed-Kelvin curl equation is
2(Phi_dot-Omega*t cross Phi)=-2*Omega*t cross Phi, hence Phi_dot=0.
That candidate is the static rotational symmetry, not the optical mode.

Thus the new positive theorem is a robust same-tube gyro-dominant Euler
packet with actual nonconstant material-parcel spin. A whole-tube
mechanical/canonical matching theorem and a fixed-lambda arbitrary-knot
twist theorem remain different constructions, not consequences of this
one. Parent continuum completion is not claimed.

## Evidence and route verdicts

`verify.py` is an exact differentiation/series/sign oracle; it does not
numerically approximate the elliptic solution or its KAM threshold. The
existence, parameter-continuity, measure matching and source applicability
proofs are the analytic evidence above. Numerical error-budget,
eigenpair-floor and mesh prescriptions have no unresolved numerical
remainder here: all smallness selections are against fixed strict analytic
margins. The first execution is captured separately.

- Axisymmetric robust-torus/physical-gyro route: `route_verdict:
  established`; `evidence_scope: EXACT_LOCAL_EXISTENCE_PLUS_GLOBAL_SOURCE_
  CONTINUATION_AND_CONTROLLED_FINITE_TIME_EULER_PACKET`.
- Dominant-lambda-alone mechanical canonical matching: `route_verdict:
  refuted`; mechanism is the two-angle/shear factor and ambient current,
  with the rigid-rotation repair yielding a zero-frequency symmetry.
- Arbitrary-knot fixed-lambda route: `route_verdict: blocked` at the named
  missing fixed-eigenvalue return/twist normal form. It is not refuted or
  exhausted by the circular construction.

These are route-specific conclusions. The stronger parent objective and
its remaining physical coupling/moment identification stay active.

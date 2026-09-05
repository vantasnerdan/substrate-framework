# Actual helical Euler optical modes and their positive action

Fix c,C>0, d=c^2+r^2, h=r e_theta+c e_z, and u=f h with f=C/d.
The axial period is Lz=2pi c; the domain initially is R^2 times this
periodic axis, with no radial wall. The background energy diverges
logarithmically at radial infinity; the perturbation action and every
tagged observable constructed below are finite. Neither background
finite energy nor a Euclidean closed knotted core is asserted here.

The exact stationary Euler field and invariant scalar are the0183 input.
This proof independently derives the operator, modes, action and observer.

## 1. Full pressure and the helical momentum

Use s=theta-z/c. For helical-equivariant vector fields [h,v]=0, the
Cartesian derivative h.grad v=Jv, J=e_z cross. Therefore the exact linear
Euler equation is

    v_t+2f Jv+f_r v_r h=-grad pi, div v=0.                (1)

Both the basis rotation and the radial term are present. Write
tau=v.h=r v_theta+c v_z. Since h.grad pi=0, equation(1) gives

    tau_t=-[2rf+d f_r]v_r=-(d f)_r v_r.                  (2)

For the stated profile d f=C, tau_t=0. Every nonzero-frequency mode
therefore has tau=0. This is not the assertion that its mechanical axial
spin vanishes: tau is angular momentum plus c times axial momentum.

On this sector write the exactly solenoidal velocity as

    v_r=phi_s/r, v_theta=-c^2 phi_r/d, v_z=cr phi_r/d.     (3)

Its physical kinetic metric is

    integral |v|^2 r dr ds
       =integral [phi_s^2/r^2+c^2 phi_r^2/d] r dr ds.

The two independent pressure equations from(1) are

    pi_r=-phi_st/r-2f c^2 phi_r/d,
    pi_s=rc^2 phi_rt/d-2f c^2 phi_s/d.

Equating their mixed derivatives gives

    H phi_t=beta phi_s,
    H=-r^-1 partial_r(rc^2/d partial_r)-r^-2 partial_s^2,
    beta=8Cc^2/d^3>0.                                   (4)

Thus the pressure has been eliminated through its compatibility equation,
not omitted or replaced by a radial boundary condition.

There is an exact nonlinear geometric check. For any equivariant velocity
with total helical momentum C, the identity
h cross curl v=grad(v.h) makes its horizontal vorticity parallel to h.
Its total vorticity is

    omega_total=(H phi+b)/c * h, b=2Cc^2/d^2.

Euler vorticity transport becomes actual advection of H phi+b by(3).
The vertical velocity C h/d contributes the fixed energy integral C^2/d;
the variable reduced kinetic energy is precisely the positive metric above.
Equation(4) is the linearization of this actual Euler reduction.
The finite relative energy is defined by its pointwise difference: the
background is vertical and v is perpendicular, so that difference is
|v|^2/2. This is not an adjustable infrared subtraction. Helical material
maps commute with h and preserve the constant transported helical momentum;
their higher-order horizontal corrections cannot introduce a hidden
background cross term into this quadratic energy.

## 2. A finite-energy spectral problem, including both endpoints

Take integer m>=1 and a mode phi=phi_m(r) exp(i m s+i sigma t).
Its physical laboratory angular frequency is omega=-sigma in the
exp(-i omega t) convention. Equation(4) is

    sigma H_m phi_m=m beta phi_m, sigma>0.               (5)

The energy space is the completion of regular compact radial profiles in

    E_m(phi)=integral_0^infinity
           [rc^2/d |phi'|^2+m^2/r |phi|^2]dr.            (6)

At the origin finite energy selects the regular r^m branch; the r^-m
branch is excluded. On each bounded disk this is the ordinary m-sector
H1 energy with comparable coefficients, so Rellich compactness applies
to its weighted L2 beta norm. It can also be proved directly by the
one-dimensional radial energy estimate away from0 and the m^2/r term
near0.

At infinity the radial coefficient is c^2/r^2, not a uniform elliptic
constant. Compactness still follows from the angular energy:

    integral_R^infinity beta |phi|^2 r dr
      <=sup_(r>R)(r^2 beta/m^2) E_m(phi)
      <=C_(c,C,m) R^-4 E_m(phi).                          (7)

Thus the embedding of(6) into L2(beta r dr) is compact. The positive
operator H_m^-1/2 beta H_m^-1/2 is compact and has infinite rank. The
variational maximum gives a positive largest eigenvalue; subsequent
orthogonal variations give positive eigenvalues tending to0. Their
corresponding frequencies sigma_j are real and positive. The ground
radial mode is simple and strictly positive by the one-dimensional
Sturm variational argument. No numerical eigenvalue or boxed wall enters.

At radial infinity equation(5), after multiplication by r^2/c^2, is

    -phi''+phi'/r+(m^2/c^2)phi+O(r^-2)phi'
                                      +O(r^-2)phi=0.

Comparison with exponential supersolutions shows that its finite-energy
branch and every fixed derivative decay as exp(-a r), for any a<m/c.
The regular mode is analytic at finite r. Its finite positive frequency
is specified by the variational eigenproblem(5), not supplied as a fit.

## 3. Exact Kelvin reconstruction and KKS sign

Use the real phase convention

    phi(r,s,t)=phi_m(r)[a(t)cos(ms)+b(t)sin(ms)],
    a_t=sigma b, b_t=-sigma a.                            (8)

For the complex mode at nonzero frequency, the actual Lin displacement is

    xi=v/(i sigma)+f_r v_r h/(i sigma)^2.                 (9)

Indeed [u,xi]=-f_r xi_r h, so v=xi_t+[u,xi]. The second term is parallel
to h and has zero divergence; the first term is solenoidal. This is a
material displacement, not merely a velocity polarization.

Here omega0=(2cC/d^2)h. Direct substitution of(1) yields

    xi cross omega0=v+grad pi/(i sigma).

Hence the full Leray projection gives P(xi cross omega0)=v exactly.
The mode is on the actual fixed-Kelvin orbit. No independent hidden
circulation oscillator is attached to it.

Let

    I_m=E_m(phi_m), h_m=rho*pi*Lz*I_m.

The full phase Hamiltonian is H_phase=h_m I_2. The vertical part of(9)
does not contribute to omega0.(xi_a cross xi_b). For the two real phase
velocities one finds

    h.(v_a cross v_b)=m c phi_m phi_m'/r.

Integration by parts, with both endpoint terms zero, and(5) give the
ACTUAL KKS form

    Omega_ab=rho integral omega0.(xi_a cross xi_b)
       =rho*pi*Lz*m/sigma^2 integral beta phi_m^2 r dr
       =h_m/sigma>0.                                    (10)

With Omega=[[0,h_m/sigma],[-h_m/sigma,0]], the convention
L=-x^T Omega x_t/2-x^T H_phase x/2 gives exactly(8).
This computes the Euler action sign in the physical clock; it is not a
Floquet-branch choice or a magnetic stability analogy.

## 4. A stationary material angle and literal mechanical spin

Choose m=2 for an actual transverse helical quadrupole. The scalar
F_m=(x+i y)^m exp(-i m z/c)=r^m exp(i m s) is exactly stationary along u.
Take a compact radial material fraction

    w0=chi(r)[1+epsilon b_tag(r)cos(ms)],
    0<=chi<=1/2, |epsilon b_tag|<1.

It is nonnegative and transported without reference dephasing. Let

    B_tag=integral chi b_tag r^(m+1)dr !=0,
    Q0=rho*pi*Lz*epsilon B_tag.

Its actual material moment Q=integral w F_m has this nonzero stationary
reference. The first variation satisfies delta Q_t=integral w0 v.grad F_m,
including moving labels. The helical perturbation with matched initial
Lin displacement(9) gives

    delta Q=(rho*pi*Lz*m I_obs/sigma)(a+i b),
    I_obs=integral chi r^(m-1)(m phi_m+r phi_m')dr,
    theta=delta(arg Q)/m=c_obs b,
    c_obs=I_obs/(sigma epsilon B_tag).                    (11)

Choose chi in a core interval where I_obs!=0; the positive regular ground
mode has m phi_m+r phi_m'>0 near0. An arbitrary initial label perturbation
could add a constant offset to(11). We choose the actual particular
displacement(9), so that offset is explicitly fixed, not discarded.

Eliminating a from the actual phase action gives

    L_theta=M/2(theta_t^2-sigma^2 theta^2),
    M=h_m/(sigma^2 c_obs^2)>0,
    Pi=M theta_t=-h_m a/(sigma c_obs).                    (12)

The angular momentum here must be computed mechanically. Its first
variation about the material tag's transverse centroid is

    S_z=rho integral w0[r v_theta+(r^2 f)'xi_r].          (13)

For m=2 the transverse reference centroid and momentum vanish by angular
orthogonality, so no origin/centroid correction is hidden in(13). The
tag's axial mechanical momentum obeys the exact companion identity

    S_z=-c delta P_z,
    delta P_z=rho integral w0[v_z+c f_r xi_r].            (14)

Zero helical charge is not zero mechanical spin. Substitution of(3),(9)
in(13) gives

    S_z=rho*pi*Lz*epsilon a J_tag,
    J_tag=integral chi b_tag
       [-c^2 r^2 phi_m'/d+2mCc^2 r phi_m/(sigma d^2)]dr.  (15)

This is a literal finite material moment. A finite material axial parcel
may also be used on the covering space: all scalar integrands above are
helically invariant, and integration over complete label phases gives
the same factor Lz. The mode's action itself is per periodic cell until
the finite-packet construction in axial-continuation.md is applied.

Fix eta>0 in advance, for example eta=1. The exact physical row condition
S_z=eta Pi is the linear radial equation

    rho*pi*Lz*J_tag+eta h_m B_tag/I_obs=0.                (16)

It can be solved with B_tag!=0 using two smooth radial controls. To see
their independence, divide the integrand of J_tag by r^(m+1). This analytic
function is not constant on any core interval. If it were constant, decay
at infinity would make it0 everywhere. Then J's pointwise density would
vanish and the mode would solve

    phi_m'/phi_m=2mC/(sigma r d).

Its solution tends to a nonzero constant at infinity, contradicting(6).
Therefore choose two interior radial bumps whose J/B ratios differ, solve
the exact2x2 moment system with a prescribed nonzero B_tag, and decrease
epsilon so the physical density remains nonnegative. No mode frequency
or inertia was fitted: sigma comes from(5), and every coefficient in(16)
is the actual Euler action or material moment. Equality holds throughout
the full linear mode history, not just at initial time.

## Scope

This proves actual neutral optical Euler modes, positive fixed-Kelvin
action, a stationary physical material-quadrupole angle, and exact tagged
mechanical-spin normalization on the stated periodic helical column.
It explicitly retains the linked axial momentum. It is not an independent
translation/rotation continuum, a Euclidean closed EPS tube, or a carrier-
curvature theorem. Smooth radial localization and actual finite-action
axial packets are separate constructions in the next file.

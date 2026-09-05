# An actual stationary-domain tilt mode with positive action and exact current jets

## 1. A different physical observation and the fixed target

This implements0181's smooth finite-core spectral candidate. The source
0137 supplies the full radial Euler system and thin-annulus transfer
method, not an assumed pole for an arbitrary EPS ambient. We extend
that method to an explicitly identified RADIAL-NODAL m=1 branch.

Start with the ordinary Rankine column of radius a and angular speed
Omega>0, with zero axial base velocity. Select a fixed nonzero axial
carrier p and the phase exp(i theta+i p z-i omega t). The actual
material tag is axisymmetric and stationary, not painted with rotating
m2 lobes. Its centered covariance has

    I_perp=rho integral_tag x²=rho integral_tag y²,
    I_z=rho integral_tag z², Delta=I_z-I_perp>0.

Its measured tilt is the ordinary physical orientation of this
distinguished covariance axis:

    theta_x=-delta Q_yz/Delta,
    theta_y=delta Q_xz/Delta.                            (1)

The response to a rigid physical tilt is exactly one. A smooth positive
mass fraction chi(r,z)=chi_r(r)chi_z(z), even in z and supported inside
the rigid core, is an invariant material density because u.grad chi=0.
Its finite reference shape is therefore stationary. This removes the
painted-marker clock by an actual symmetry of the observed material,
not by renaming an intrinsic phase.

The positive result is an actual smooth-column Euler eigenmode with
positive physical tilt action and squared-frequency curvature, plus
a ONE fixed positive stationary tag with the full mechanical and
displacement current identities through the second carrier jet. Its
measured spin overlap may be any resulting eta>0; eta=1 is not required.
The whole-Euclidean stationary EPS and coherent band-edge continuation
are separate, explicitly retained below.

## 2. A genuine nodal Euler branch and its transverse zero-frequency crossing

Let x=p a, tau=-sigma/Omega in(0,2), sigma=omega-Omega, and

    l=x sqrt(4-tau²)/tau,
    J(l)=l J_1'(l)/J_1(l), K(x)=x K_1'(x)/K_1(x).

The COMPLETE Rankine pressure/normal-displacement determinant of0135
is equivalently

    D(l,x)=J(l)+(l²/x²)K(x)+sqrt(1+l²/x²)=0.           (2)

Its exterior perturbation is the actual decaying irrotational Euler
response, with its radially varying Doppler pressure. No wall or
discarded exterior pressure is involved.

Here is a direct monotonicity proof, valid between Bessel poles. Put
y_in(r)=J_1(lr)/J_1(l) and y_out(r)=K_1(xr)/K_1(x).
Differentiating their radial Sturm equations in l² and x², multiplying
by y and integrating the Wronskian gives

    J'(l)=-2l integral_0^1 r y_in² dr<0,
    K'(x)=-2x integral_1^infinity r y_out² dr<0.        (3)

The parameter derivatives of y at r=1 vanish because its Dirichlet
normalization is fixed; its radial derivatives need not vanish.
Also K=-1-x K_0/K_1<-1. Consequently

    D_l=J'(l)+(l/x²)[2K+1/sqrt(1+l²/x²)]<0.          (4)

In every interval(j_1,n,j_1,n+1), J decreases from+infinity to
-infinity, so(2) has one simple real radial root l_n(x). The pressure
P(r)=J_1(l_n r/a) has n interior zeros. Select n>=2.

At a LABORATORY zero frequency, tau=1 and l=sqrt(3)x. The exact
crossing function is

    D0(x)=J(sqrt(3)x)+3K(x)+2.                        (5)

It decreases strictly from+infinity to-infinity on
(j_1,n/sqrt(3),j_1,n+1/sqrt(3)). Thus it has one zero x_* there.
Equations(3),(4) show that it is transverse, without any assumption
about the arithmetic or numerical value of Bessel zeros. In fact

    l_n'(x_*)-sqrt(3)=-D0'(x_*)/D_l<0,
    tau'(x_*)=sqrt(3)[sqrt(3)-l_n'(x_*)]/(4x_*)>0.

Hence omega'(x_*)=-Omega tau'(x_*)<0. Choose a fixed x0>x_* close
enough that omega(x0)<0 but tau<2. At the crossing

    (omega²)''=2 omega'²>0,

so this squared-frequency curvature remains strictly positive on a
one-sided finite neighborhood. The frequency itself is nonzero at x0.
All choices and margins are fixed before smoothing. This is an actual
simple Euler mode, not a quasimode promoted to a spectral pole.

## 3. Same full KKS and positive laboratory action

In the rigid core the actual displacement components before the common
Fourier phase are

    xi_r=A=(2Omega P/r-sigma P')/[sigma(4Omega²-sigma²)],
    xi_theta=iB=i(2Omega P'-sigma P/r)
                                      /[sigma(4Omega²-sigma²)],
    xi_z=i p P/sigma².                                (6)

They give v=-i sigma xi in cylindrical components and satisfy full
Lin reconstruction, including the rotating Cartesian basis. For the
two real columns the COMPLETE orbit pairing is

    beta=4pi rho Omega L_z integral_0^a A B r dr.       (7)

There is no exterior vorticity volume contribution; the exterior
velocity/pressure response is already in(2). Put tau=-sigma/Omega.
The integral in(7), multiplied by the positive denominator
tau² Omega⁴(4-tau²)², is

    2tau integral_0^a [r(P')²+P²/r]dr
                     +(tau²+4)P(a)²/2>0.             (8)

This proves beta>0 even for the nodal pressure. The nodal signs were
not mistaken for a negative norm. With the conventions of0135, the
actual LABORATORY mode Hamiltonian is -beta omega I, and for a measured
scalar tilt column of nonzero length c_theta,

    M=-beta/[omega c_theta²]>0,
    L=M(theta_dot²-omega² theta²)/2.                  (9)

This uses the stationary physical reference(1), so its clock is omega,
not sigma or a rotated material paint clock. The spatial kinetic jet
is retained: normalizing the mass leaves the exact propagation
combination M0 partial_p²(omega²)/2>0. Neither M_p nor M_pp is set
to zero. A positive lab Hamiltonian holds on this selected branch.

## 4. Complete stationary-domain shape, spin and displacement moments

Write Z(p)=integral chi_z(z)exp(i p z)dz, real for the even axial tag,
and define actual radial integrals

    R1=integral chi_r r(A+B)dr,
    R2=integral chi_r r²P dr.

Let

    Dq=rho pi[-Z' R1+p Z R2/sigma²],
    Dg=rho pi[ Z' R1+p Z R2/sigma²].                  (10)

Angular integration of the complete displacement gives, in complex
notation with e=(1,i),

    delta Q_hz=i Dq e exp(-i omega t),
    theta=(Dq/Delta)e exp(-i omega t),
    G_h=rho integral_tag(r cross xi)_h
                              =-Dg e exp(-i omega t). (11)

These are centered moments; reference centroid is zero, and the
centroid-variation terms vanish at first order because the reference
mean velocity and first position moment vanish. The tag is invariant,
so the full mechanical spin can also be written exactly as

    S=rho integral_tag[r cross xi_t+2xi cross u0].

Equivalently it includes BOTH the position and material-velocity
variations. Using(6), its transverse row is

    S_h=i rho pi[omega Z'R1+(sigma-Omega)p Z R2/sigma²]
                                           e exp(-i omega t). (12)

Thus the sole radial condition

    R2=0                                                     (13)

implies the exact all-time identities on the actual eigenmode

    G_h=Delta theta, S_h=Delta theta_dot,
    G_h(0)+integral_0^t S_h=Delta theta(t).              (14)

No initial-G adjustment is being supplied as an extra coordinate.
It is the ACTUAL displacement moment. For a single carrier a positive
tag on opposite sides of one pressure zero solves(13); the R1 row
has one sign in a sufficiently small such neighborhood and is nonzero.

The measured spin/action overlap is eta=Delta/M>0. It can be small,
and scales with the actual positive tag fraction. It is not imposed
equal to1. All coefficients in(9),(14) come from the same Euler mode,
tag and full-cell action. The optical translational coupling after the
physical current map is eta*kappa/2, as derived in0176, not zero merely
because eta differs from1.

## 5. One fixed POSITIVE tag controls the current through carrier two

For the actual mode pressure normalization write
P(r,p)=A0(p)J_1(l(p)r/a), A0(p)!=0. Define

    F(l)=integral chi_r(r) r² J_1(lr/a)dr.

If F,F_l,F_ll vanish at the selected l0, then R2 and both of its p
derivatives vanish there, including the pressure normalization factors.
No material support is moved when differentiating the carrier.

The Bessel equation gives the exact identity

    l² F_ll+l F_l-F=-(l²/a²) F4,
    F4=integral chi_r r⁴J_1(lr/a)dr.                  (15)

It is enough to impose F=F_l=F4=0. Choose two consecutive pressure
zeros r1,r2 inside the core; n>=2 ensured their existence. Begin with
positive point-annulus weights w1,w2 at those radii. F and F4 vanish
individually. Since J_1' alternates sign at successive simple zeros,
there is a positive ratio with

    w1 r1³ J_1'(l0r1/a)+w2 r2³ J_1'(l0r2/a)=0.        (16)

This solves F_l=0. Vary the two physical annulus centers and one
relative positive weight. The Jacobian of(F,F4,F_l) is nonsingular:
the two center columns in the first two rows are nonzero multiples
of(1,r1²) and(1,r2²), while the weight column has zero entries there
and nonzero F_l entry. Their determinant contains r2²-r1²!=0.

Replace each point annulus by a smooth narrow NONNEGATIVE bump of
fixed physical width. Its moment map is a smooth small perturbation;
the IFT adjusts the two centers and positive ratio to solve the three
EXACT equations. A common scale makes0<=chi_r<=1. A small positive
stationary connecting/background density may be added before the same
IFT if a connected material support is desired. It does not change
the moment count, positivity margins or stationary property.

The tilt is not annihilated by these current constraints. At the
point-annulus seed, R1 is a nonzero factor times

    w1 r1 J_1'(l0r1/a)+w2 r2 J_1'(l0r2/a).

Using(16) makes this proportional to1/r1²-1/r2², hence nonzero.
It stays nonzero for the smooth tag. Choose an even finite axial
profile with Delta>0 and Z'(p0)!=0. Such profiles exist: take a
sufficiently elongated positive profile and avoid the discrete zeros
of its nontrivial analytic Fourier derivative by a small width change.
The axial period can be chosen to contain this fixed finite tag.
It is an actual stationary material region because W=0; there is no
unobserved axial label drift.

Equations(14) now hold through the COMPLETE second carrier jet with
the SAME physical tag, for ALL times of the eigenmode. More precisely
their discrepancy is O((p-p0)³), uniformly on any fixed time interval
with its needed time derivatives. The exact mode clock, positive
inertia, curvature and full spin/G current have no painted dephasing
at that order. Higher carrier coefficients are not claimed zero.

## 6. Transfer to an actual smooth compact-vorticity column

Choose a smooth cutoff and smooth the Rankine vorticity only in
[a,a+epsilon], keeping it exactly2Omega on r<=a and exactly zero
outside the annulus. The actual angular speed is
O_e(r)=r^-2 integral_0^r s Z_e(s)ds. This is smooth stationary Euler.

The exact first-order radial pressure/displacement system from0137 is

    f'=-[1/r+2O_e/(r s)]f+(1/r²+p²)P/s²,
    P'=(s²-2O_e Z_e)f+2O_e P/(r s),
    s=omega-O_e(r).                                   (17)

Near the selected mode, |s| is bounded below in the smoothing annulus;
no Z_e' occurs. Its transfer matrix differs from identity by O(epsilon)
in the required omega,p parameter derivatives. Match it to the exact
decaying irrotational exterior pressure/normal-displacement vector.
The simple determinant root of(2),(4) then continues by the ordinary
IFT, with its frequency and first two carrier derivatives converging.
This extends the0137 method to this specified nodal branch; it does
not cite its small-carrier surface conclusion as this different pole.

Because the selected laboratory omega is NEGATIVE, there is not even
an exterior particle critical radius: omega-O_e(r)<0 for all r.
The invariant perturbation class has vorticity supported in the active
core/annulus; exterior base vorticity is zero, so this class is preserved
by full linear Euler. Its velocity and pressure tails are retained.
The actual material generator may be extended through the vorticity-
free exterior without altering the KKS volume form. No general gap
against arbitrary exterior-vorticity initial data is asserted.

The full KKS core integral and annular correction converge, preserving
beta>0. The laboratory gap, curvature and two interior pressure zeros
persist. Construct the positive fixed tag directly for the smooth
mode's actual l(p0), using(15),(16), entirely inside the EXACT rigid
core. Thus(10)--(14) and all their second carrier jets remain exact
for this smooth-column Euler mode. Smoothness of the axial and radial
tag can include rounded stationary support; all finite moment conditions
are retained by the same positive IFT.

## 7. Precise next transfer, not an inferred EPS pole or band edge

This establishes the positive stationary-domain optical mode and
literal current through carrier two on the actual smooth ordinary
column. It is stronger than an initial marker match and does not
require eta=1. The mode is a full Euler eigenmode on its stated active-
vorticity class, not a fitted oscillator or a stationary trial tangent.

The selected branch has omega'(p0)!=0. Its coherent whole-orientation
average therefore still has0172's group-variance term; constant
INDIVIDUAL clock/current does not make that ensemble term vanish.
A genuine nonzero physical band edge or another actual shared-mode
preparation is a next construction, not a conclusion of this source.

Adding a force-free axial W(r), embedding in a closed stationary tube,
or replacing the exterior by analytic constant-curl EPS ambient changes
both the material-domain invariance and the allowed spectral class.
In particular an EPS exterior has nonzero vorticity; pressure/velocity
tails can then generate exterior vorticity and recouple critical-layer
or other continuous-spectrum channels. Local approximation alone is
not an exact-pole continuation theorem. Its full operator/observation
transfer must be supplied or kept as a finite-time history statement.
The live next routes are a stationary invariant closed-domain mode
construction with its actual exterior response, or a profile/array
spectral band-edge construction on the admitted Euler class. These are
the positive continuation targets, not a scientific exhaustion claim.

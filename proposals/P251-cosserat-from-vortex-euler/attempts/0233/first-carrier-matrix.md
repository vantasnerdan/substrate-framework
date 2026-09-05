# First actual carrier matrix and physical G/spin correction

This body computes the first carrier coefficient of the actual Kelvin
preparation, not a prescribed oscillator. The ring identifies k=n/R,
n=-1. It is useful to keep k and the geometric metric coefficient
kappa=1/R separate until their respective first coefficients have been
computed. The continuous k calculation here is a local Fourier-fiber
derivative; its eventual common-laboratory-K realization is a separate
consumer, not supplied by an integer harmonic.

## 1. Actual full-pressure first coefficient

Use the 0222 right-handed local coordinates. In the literal-curl core
u=(r Omega e_theta,W e_zeta), W'=-lambda r Omega and
Z=2Omega+rOmega'=lambda W. Let S=A(r)exp(-i theta),
A=Delta_1 f, with f compact and regular. The exactly divergence-free
Kelvin preparation is

    Xi_k=(J grad S-ik lambda grad f, lambda S)exp(ik zeta).

Its actual velocity is P_k(Xi_k cross omega). At k=0,
w0=b e_zeta, b=-lambda Omega S_theta. Writing w1 for partial_k w,
zeta1=curl_perp w1,perp and d1=div_perp w1,perp gives EXACTLY

    d1=-i b,
    (partial_t+Omega partial_theta)zeta1+Z' w1,r=i Z b,
    w1,perp=J grad Delta_1^-1 zeta1+grad Delta_1^-1 d1.     (1)

These equations include the pressure, not a transverse incompressibility
constraint imposed separately. The initial vorticity coefficient is

    zeta1(0)=i lambda^2[W A+W'f'],                        (2)

and the actual initial displacement coefficient is -i lambda grad f.
The radial pressure inverses use the compact moments of 0226; otherwise
their actual tails remain in (1).

Now expand only the small inner-core parameter epsilon=(lambda a)^2,
on the fixed scaled preparation/tag support. Put Omega=Omega0 at leading
order, but retain the axial shear

    W(r)=Wc-lambda Omega0 r^2/2, Wc=2Omega0/lambda.        (3)

Its gradient matters to the first carrier row even though W-Wc is a
small relative change. It cancels a divergence term in Lin reconstruction.
Equations (1)-(2) yield at this order, suppressing exp(iOmega0 t),

    zeta1=(2i lambda Omega0-2lambda Omega0^2 t) S,
    d1=lambda Omega0 S,
    Xi1,perp=lambda(-i+Omega0 t)grad f
         +(2i lambda Omega0 t-lambda Omega0^2 t^2)Jgrad f
         -i W(r)t Jgrad S.                              (4)

The coefficient in the final term is the actual W(r). Replacing it
by Wc before differentiating would leave a false divergence residual.
The full equations (1) specify the O(epsilon) corrections to (4), with
fixed-time bounds obtained from their actual pressure/Lin equations.
Thus (4) is an exact coefficient of the stated joint inner-core
expansion, not an exact finite-epsilon solution asserted by omission.

## 2. The measured row and the new secular coefficient

Define actual radial tag functionals

    B(q)=integral r chi'(r)q(r)dr,
    D(q)=integral r^2 chi'(r)q'(r)dr,
    B_W(q)=B(Wq).

The leading density-quadrupole angle row is proportional to iB(A).
Direct integration of (4) gives its first carrier row

    exp(iOmega0 t){-i lambda D(f)
       +t[B_W(A)+lambda Omega0(D(f)-2B(f))]
       -i lambda Omega0^2 t^2 B(f)}.                    (5)

The compact return identity is

    integral chi r^2 A=B(f)-D(f).                       (6)

Consequently the direct-current return of 0226, B(f)=D(f), does not
remove the new t^2 coefficient in (5). This is the concrete new response
row exposed by the full pressure calculation. B(f)=D(f)=0 removes both
unwanted rows while preserving B(A), when those three actual functionals
are independent. The remaining linear-time B_W(A) is an actual carrier
frequency coefficient, not a fitted clock.

For a tag equal to exp(-r) on a fixed correction annulus, after dividing
each point-evaluation column by its positive exp(-r) factor, the three
weights for B(f), D(f), B(Delta_1 f) are

    -r,  2r-r^2,  -r+1+1/r.                             (7)

At r=1,2,3 their determinant is 17/3. Narrow smooth bumps at those fixed
points therefore give an actual invertible moment matrix. The tag is
constant and radial-smooth near the axis and smoothly joins the fixed
superflat positive cutoff of 0226 outside this annulus. The point values
are an exact rank calculation; sufficiently narrow fixed bumps inherit
that nonzero minor by continuity. A preselected inner f=r^3/8 near the
axis can be retained and the three bumps solve the inhomogeneous rows,
so A=r and the physical core displacement are not erased to obtain rank.

## 3. Full mechanical spin requires an additional actual lift

Deleting the direct rows in (6) does not by itself imply S_mech=G_t.
The complete identity remains S_mech-G_t=2rho integral chi Xi cross u.
For the general actual Kelvin initial lift

    Xi=(Jgrad S,h), h=Delta_1 g,
    Xi_k=(Jgrad S-ik grad g,h)exp(ik zeta),               (8)

the leading transverse cross-current at the reference time is

    integral chi (Xi cross u)_x
       =-pi[B_W(A)+Omega0(N(h)-lambda N(A))],
    N(q)=integral chi r^2 q=B(Delta_1^-1 q)
                                      -D(Delta_1^-1 q). (9)

The large W-weighted row is real physical transport. It cannot be
removed as an angular-impulse convention. The new lift (8) leaves the
initial tag density unchanged in its axial component and supplies an
actual field-changing Kelvin return, with its pressure retained.

Put F=lambda f-g and let H=Delta_1^-1 F and G=Delta_1^-1 g be compact
potentials. At leading inner-core order the exact triangular Euler/Lin
system gives

    w_perp=i lambda Omega0 Jgrad F exp(iOmega0 t),
    w_z=[i lambda Omega0 A-lambda^2 Omega0^2 t F]
                                                   exp(iOmega0 t),
    Xi_perp=[Jgrad A+i lambda Omega0 t Jgrad F]
                                                   exp(iOmega0 t),
    Xi_z=h exp(iOmega0 t).                              (10)

The axial shear contributions cancel in the final Lin equation; the
axial component is not discarded. The zero-carrier observed secular
coefficient is B(F). Its cross-current has a constant and a linear-time
coefficient proportional to

    C0=B_W(A)+Omega0[N(h)-lambda N(A)],
    C1=B_W(F)-lambda Omega0 N(F).                       (11)

For the first carrier row the coefficients, after the common exponential,
are

    t^0: -iD(g),
    t^1: A1=B_W(A)-lambda B_W(g)
                    -lambda^2 Omega0 D(G)+lambda Omega0 D(f),
    t^2: i lambda Omega0[B_W(F)+lambda Omega0 D(H)
                                                   -Omega0 B(f)],
    t^3: lambda^2 Omega0^3 B(H)/3.                       (12)

The derivation uses the actual divergence row and the exact compact
identities Delta_1^-1(r F')=rH'-2H and
Delta_1^-1(W Delta_1 H)=WH+2lambda Omega0(rG_H'-G_H),
where G_H=Delta_1^-1 H. It retains the pressure-generated polynomial
time terms before projecting them onto the physical tag.

## 4. Explicit finite physical-row reduction

Equations (10)-(12), together with (6), reduce simultaneous zero-carrier
clock purity, C0=C1=0, and removal of the first-carrier t^2,t^3 rows to

    B(f)=0,
    B(H)=0, B(Delta_1 H)=0,
    E(H):=D(H)-B(r^2 Delta_1 H)/2=0,
    D(Delta_1 H)=-B_W(A)/Omega0.                        (13)

Here g=lambda f-Delta_1 H and h=lambda A-Delta_1^2 H. These are actual
preparation rows. In particular C0 fixes D(g), rather than silently
deleting its contribution to the initial physical angle.

For chi=exp(-r) on the correction annulus the four H-weights, divided
by exp(-r), are

    B(H):             -r,
    B(Delta H):       -r+1+1/r,
    E(H):             (r^3-7r^2+7r)/2,
    D(Delta H):       -r^2+5r-2-2/r.                    (14)

At r=1,2,3,4 their exact four-point minor is 37/2. This was computed
before a bump selection was made. The matrix, not the dimension of a function space,
licenses solving (13). Compact f and H preserve all required radial
inverse returns. Further finite pressure moments use separate off-tag
bumps and their explicitly checked independent rows.

The surviving first carrier time row is A1. Its frequency and initial
amplitude conversion are the actual physical observation (12); they are
not set equal to Wc or to a desired Cosserat coefficient. The correction
can have inverse-small-parameter amplitude. Its cost, the full phase
and energy, and the finite-time remainder must therefore be retained
before any density/common-K limit. The next body supplies ring-native
phase/energy controls, rather than importing a different cell's action.

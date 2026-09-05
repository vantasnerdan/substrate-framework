# Exact rotating-Euler sidebands and the averaging distinction

Use the stationary comparison u=Omega ez cross x+U ez, Omega!=0.
Its physical pressure is rho Omega²(x²+y²)/2. This is not an EPS
field or a finite-energy global vortex tube. A bounded local transfer
needs0120's separate construction. All formulas below are exact for
this comparison. No numerical spectrum or curvature tolerance is used.

## Derive the pressure-projected amplitude and name the clock

Let A=Omega[ez cross], R(t)=exp(At). An exact Kelvin plane wave has

    k(t)=R(t)k0,
    phi(t,x)=k(t).x-U k0,z t,
    v=R(t)b(t) exp(i phi),
    k0.b=0.

The phase solves (partial_t+u.grad)phi=0. Substitution in the COMPLETE
linear Euler equation, with its pressure amplitude eliminated using
incompressibility, gives in this physical co-rotating frame

    bdot=-2 P_k0 A b,
    P_k=1-kk^T/|k|².                                    (1)

On k0-perp, P_k[ez cross] is (ez.k/|k|)[k/|k| cross].
Thus the exact squared frequency in this specified frame is

    f(k)=omega_rot²=4 Omega² (ez.k)²/|k|².               (2)

The factor two arises by differentiating the physical frame as well as
the velocity amplitude; dropping either A term is an exposed wrong
convention. The lab observation is R(t)b(t)exp(i phi), not b(t).
Its fixed-probe lines include the rotation matrix and transported
carrier. Equation (2) is not silently renamed the lab optical gap of
the knotted medium or the frequency of a material core-angle observation.

## Scalar sideband jet with co-rotated geometry

Fix |p|=N and c=ez.p/N with 0<|c|<=1, then add a slow sideband K.
Rotate BOTH ez and p by R in SO(3), preserving c. Write

    f_R(p_R+epsilon K)=f0+epsilon a_R+epsilon² b_R+O(epsilon³),
    f0=4 Omega² c².

Exact differentiation followed by the invariant second moment gives

    <a_R>=0,
    <b_R>=4 Omega²(1-3c²)|K|²/(3N²),
    <a_R²>=64 Omega^4 c²(1-c²)|K|²/(3N²).              (3)

In particular the AVERAGED scalar curvature is positive for c=1/2,
and negative for axial c=1. These are exact curvature signs, not
evidence that the averaged physical dynamics has a single normal mode.

## The false-positive closure exposed by the actual evolution

For each orientation let x_R solve x_R,tt+f_R x_R=0 with x_R(0)=1,
x_R,t(0)=0. Its exact initial derivatives are x_R,tt(0)=-f_R and
x_R,tttt(0)=f_R². Therefore an oscillator built from the averaged
frequency reproduces the second derivative but misses the fourth by

    <x_R,tttt(0)> - <f_R>²
      =epsilon² <a_R²>+O(epsilon³).                   (4)

For the positive-curvature oblique candidate c=1/2 this is strictly
positive at the SAME second spatial order as the desired stiffness.
No simulation is needed: the initial Taylor coefficient already exposes
the missing branch dispersion/memory. The positive averaged curvature
is a valid scalar fact but not the requested exact dynamical closure.
It cannot be used as a replacement for the true Euler response oracle.

The axial candidate has a_R=0 for every orientation, so this specific
obstruction is absent to second spatial order. For fixed finite time,
analytic dependence of the oscillator propagator on f shows that its
scalar average is governed through that order by <b_R>; variance first
enters at fourth spatial order. Its scalar curvature is negative,
not a source of positive curvature hidden by a convention change.
The c=0 candidate also has a_R=0 but loses the optical gap and the
leading fixed-Kelvin displacement/velocity invertibility, so it cannot
substitute for this positive-gap construction.

## A vector-polarization screening identity, not a mode transfer

If a proposed axial preparation were represented by P_t Phi, t=R ez,
the exact fourth spherical moment would give

    <P_t>=2 I/3,
    <[|K|²-(t.K)²]P_t>=2|K|² I/5+2KK^T/15,
    <P_t>^-1 <[|K|²-(t.K)²]P_t>
       =3|K|² I/5+KK^T/5.                              (5)

This is the transverse/longitudinal weighting that such a preparation
would carry. It is NOT asserted that its actual Kelvin velocity,
vorticity-angle, moving-tag spin and Fourier pressure projections all
have that unaltered polarization away from K=0. Their changing projectors
and the physical frame remain additional terms in an actual vector join.
The scalar and tensor screens are preserved rather than promoted as
Cosserat coefficients.

## Route result and continuation

Established: the exact pressure-projected comparison dynamics, named
frame frequency, sideband curvature and variance; the oblique averaged-
frequency closure is refuted at its proposed second-jet scope by (4).
The scalar axial finite-time closure survives that particular test but
has negative curvature. Neither result is a no-go for the parent.

Next actual routes are the fixed-lambda curved-core dispersion0120 and
the full periodic Bloch response0116, or retaining the oblique branch
distribution as an explicit extra state rather than fitting one oscillator
to its averaged frequency. A physical-current/field map is required before
any of these identities can enter the autonomous Cosserat equations.

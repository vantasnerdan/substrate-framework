# Neutral vortex street: physical impulse and exact acoustic symbol

This failure-derived candidate was registered in README before selection.
It is a planar point-vortex theorem first; the smooth-core continuation
below has its own explicitly smaller dynamical scope. No external walls,
magnetic force, or prescribed sound speed is introduced.

## Physical registration and the zero-wavenumber mass

Let a>0, Gamma>0. Put +Gamma vortices at (ma,b/2), and -Gamma vortices
at ((m+1/2)a,-b/2). The laboratory far-field velocity is zero. The street
travels with V=Gamma tanh(pi b/a)/(2a) in the +x direction. The fixed
street frame subtracts this constant V; this is not an adjustable
frequency branch. Work per period and per unit axial length, density rho.

Use coordinates

    x_+=X-d/2, x_-=X+d/2,
    y_+=Y+b/2, y_-=Y-b/2.

The exact vortex symplectic form becomes

    rho Gamma (dx_+ wedge dy_+ - dx_- wedge dy_-)
      =rho Gamma(dX wedge db+dY wedge dd).

Thus P_X=rho Gamma b, P_Y=rho Gamma d. The first is also **actual fluid
momentum**, not merely a canonical label. The x-average of the velocity
is Gamma/a between the two rows and zero outside. Consequently

    rho integral_[0,a]xR u_x dxdy=rho Gamma b.

More generally, integrate omega=partial_x u_y-partial_y u_x against y:
for any neutral smooth periodic street decaying at y infinity,
rho integral u_x=rho integral y omega. The same formula holds exactly
for smooth cores with circulations +/-Gamma and centroid separation b.
In a translating frame it is momentum relative to the fixed uniform
far-field background; no infinite boost mass is counted.

The row interaction and self regular part of periodic Biot--Savart give

    H(a,b,d)=gamma log[a sqrt(cosh(2pi b/a)-cos(2pi d/a))]+constant,
    gamma=rho Gamma^2/(2pi).

The core constant is independent of a,b,d. At d=a/2, write t=pi b/a:
H=gamma log[a cosh t]+constant. Its b derivative is rho Gamma V.
H_bb=gamma pi^2 sech^2(t)/a^2>0, whereas H_dd=-H_bb.
Eliminating the separation momentum in the longitudinal sector gives

    M_X=(rho Gamma)^2/H_bb
       =2rho a^2 cosh^2(t)/pi >0.

This equals dP_fluid/dV. The transverse common-coordinate sector has
opposite mass sign; no positive full four-coordinate action is claimed.

## Full lattice symbol, not a local-energy assumption

Set a=1 temporarily; physical wavenumber K corresponds to k=Ka. Write
gamma=rho Gamma^2/(2pi), g=rho Gamma and h=b/a=t/pi. Use a Bloch phase
at each vortex's actual x location, thereby including the half-period
row offset. For 0<k<2pi define

    D=pi k-k^2/2-pi^2 sech^2(t),
    F=pi/cosh(t) [(pi-k)cosh((pi-k)h)
                     -pi tanh(t)sinh((pi-k)h)],
    G=pi/cosh(t) [(pi-k)sinh((pi-k)h)
                     -pi tanh(t)cosh((pi-k)h)].

The exact point-vortex energy Hessian in the longitudinal pair (X,b) is

    H_L(k)=gamma [[2(D+F),-iG],[iG,(F-D)/2]],
    L_L=g Re(b* X_t)- (1/2)(X*,b*) H_L (X,b).

The transverse pair (Y,d) gives the complementary sign-reversed block.
These formulae retain all lattice interactions. To derive them, the pair
Hessian is gamma s_i s_j [[Re z^-2,-Im z^-2],[-Im z^-2,-Re z^-2]].
The same-row sum is pi k-k^2/2. The cross-row sum follows by differentiating

    sum_m exp(ikm)/(m+z)=pi exp(i(pi-k)z)/sin(pi z), 0<k<2pi,

at z=1/2+/-ih. Symmetric summation defines the first identity; the squared
denominator sum is absolutely convergent. Multiplication by exp(ik/2)
implements the physical row offset. Its real- and imaginary-kernel
combinations give F and iG. Negative k is the conjugate symbol.

Let A=-1+2t tanh(t)+t^2 sech^2(t). Exact expansions give

    F+D=(A/2) k^2+O(k^3),
    F-D=2pi^2 sech^2(t)-2pi k+O(k^2),
    G=-pi[tanh(t)+t sech^2(t)]k+O(k^2).

Therefore the two longitudinal frequencies are

    omega=-(Gamma/(2pi a^2)) [G +/- sqrt(F^2-D^2)],

and, if A>0, their physical long-wave speeds are

    omega/K = Gamma/(2a)[tanh(t)+t sech^2(t)
                             +/- sqrt(sech^2(t) A)] + O(|K|a).

At the classical spacing cosh^2(t)=2, A=-1+sqrt(2)t+t^2/2>0: indeed
t=asinh(1)>1/sqrt(2), directly from its positive integral. Both roots
are real and distinct for all sufficiently small nonzero |K|. This is
an exact point-vortex spectral statement with a controlled analytic
small-k remainder, not a whole-spectrum or smooth-Euler stability import.
The positive H_bb(k) supplies a finite dynamic mass approaching M_X;
the physical longitudinal impulse variation is g b at k=0.

The registration is important. The longitudinal action after separation
elimination is, through its leading slow derivatives,

    L= M_X/2 (X_t-v_d X_s)^2 - C_X/2 X_s^2,
    v_d= -Gamma/(2a)[tanh(t)+t sech^2(t)],
    C_X=gamma A,   C_X/M_X=(Gamma/(2a))^2 sech^2(t) A,

where s labels the reference street positions. Thus the propagation
speeds are -v_d +/- sqrt(C_X/M_X), as obtained from the full symbol.
The locally varying period calculation H(a(1+X_s),b) gives the same
coefficients, but it is not used instead of that symbol derivation.

Crucially H_aa H_bb-H_ab^2=-gamma^2 pi^2/a^4. Equivalently

    v_d^2-C_X/M_X=(Gamma/(2a))^2>0.

The fixed street-frame Hamiltonian is indefinite; the two real branches
have a convective bias larger than their relative acoustic speed.
A coordinate change can center them, but it does not create positive
energy in the preselected physical street frame. The result is a real,
finite-impulse-mass **convected longitudinal acoustic pair**, not a
positive isotropic elastic Euler continuum. Its phase coordinate tracks
a vortex-row displacement; it is not a uniform whole-space boost.

## Consequences for the candidate ladder

The point-street positive-mass acoustic pair is established in its stated
row-pattern scope. Its nonzero k=0 physical impulse is stronger than a
pure canonical label, but it does not by itself establish a nonzero-k
bulk mean chart. For a longitudinal Bloch perturbation v(y)exp(iKx)
decaying at y infinity, exact incompressibility gives

    iK integral_R v_x dy + [v_y]_-infinity^infinity=0.

Thus the **complete-strip longitudinal fluid momentum vanishes for every
K!=0**. The k=0 impulse limit is nonuniform: ambient return flow spreads
over distances of order 1/|K| and cannot be deleted from the observation.
The row displacement and a finite strip/tag observation remain nonzero,
with that ambient exchange retained. In particular this result is not a
bulk longitudinal mean-acoustic mode and does not meet the positive
fixed-frame, isotropic smooth-Euler target.
The negative transverse mass and fixed-frame energy mechanism are explicit.
Finite-core stability/action transfer requires its own dynamics proof;
mere persistence of steady cores does not transfer the lattice spectrum.
The next executable variants are an independently tested smooth-core
isospectral/adiabatic street family, and a neutral multi-row geometry
whose full physical-momentum Schur matrix removes this sign defect.
Neither has been silently excluded from the parent candidate universe.

# Positive fourth-harmonic material response on the actual nonlinear cell

This establishes a fixed-tag response, not yet a changed vortex-structure
mode or the common acoustic/optical second spatial jet. Every field below
is on the actual fixed C016 cell. Its nonlinearity is essential to the
positive observed branch; an affine replacement deletes its angle.

## Exact homogeneous Euler/Lin and material observations

Use coordinates (X,a,b), psi=cos b+Omega^2 cos a, Omega=1/10,
u=(psi,sin b,-Omega^2 sin a), and T=u_perp.grad. Choose a fixed smooth
nonnegative stationary chi(psi) supported inside the elliptic island,
with Q=integral chi(a^2-b^2)dA>0. It can be chosen decreasing in the
small energy E=1+Omega^2-psi, hence chi_psi>0 on an inner regular annulus.
The tag occupies the complete axial circle, with all expressions below
per axial length; densities divide by the actual fixed cell area.

For an actual compactly supported stream s0 in that annulus, set

    xi_perp(0)=J grad s0, xi_X(0)=0, w(0)=0.

The full linear Euler solution is exactly w=0. The full Lin solution is

    s(t)=exp(-tT)s0, xi_perp(t)=J grad s(t),
    xi_X(t)=-t T s(t).                                  (1)

Indeed [u_perp,J grad s]=J grad(Ts), and
xi_perp.grad psi=-Ts. Thus xi_t+[u,xi]=0 in all three components.
The displacement remains solenoidal, smooth and supported over the same
invariant annulus. Its actual coordinate velocity is not zero. The
stationary tag variation is delta_chi=chi_psi Ts and is passively
transported. These are field-preserving material-label variations, not
changes of Eulerian vorticity; a physical vortex-mode interpretation
requires an additional construction.

The literal linear angle, spin and displacement moment are

    theta = integral chi_psi a b Ts / Q,
    S = rho integral chi_psi (Dpsi) Ts,
    G = rho integral chi Ds
      = -rho integral (2chi+chi_psi Dpsi)s,                (2)

where D=a partial_a+b partial_b. Both initial centroid components vanish
for the even angular sectors used below. The term involving 2chi has zero
integral for each nonzero streamline harmonic. Differentiating (2) using
s_t=-Ts gives EXACTLY G_t=S, without a missing initial G or velocity spin:
the Eulerian velocity variation term is zero because w=0. This also follows
from the complete stationary-tag identity in0209. The axial component in
(1) remains in the full energy, phase and centroid response.

## Actual time-angle orbit, not an affine ellipse

The elliptic one-degree-of-freedom analytic Hamiltonian admits regular
time-angle coordinates on each sufficiently small positive-energy annulus.
Choose the increasing flow angle theta_o so T=omega partial_theta_o and
the orbit phase fixed at b=0,a>0. Analytic periodic-orbit continuation,
or the energy integral and analytic period, gives convergent small-
amplitude expansions on such an annulus. Use epsilon as the coefficient
of the first cosine harmonic of a (not its turning-point value). Directly
substituting in omega a'=sin b, omega b'=-Omega^2 sin a gives

    a = epsilon cos theta_o
        +epsilon^3(3Omega^2-1)cos(3theta_o)/192+O(epsilon^5),
    b = -Omega epsilon sin theta_o
        +epsilon^3[Omega(1-Omega^2)sin theta_o/16
                   +Omega(3-Omega^2)sin(3theta_o)/192]
        +O(epsilon^5),
    omega = Omega-epsilon^2 Omega(1+Omega^2)/16
                 +O(epsilon^4).                         (3)

Remainders and their required finite energy/angle derivatives are analytic
on a fixed sufficiently small regular annulus. The zero-energy center is
not used as a selected orbit. The orbit-energy check verifies that (3)
uses a time angle, rather than imposing a circular parametrization.

Write a b = sum A_l(E) sin(l theta_o), and
Dpsi = B_0(E)+sum B_l(E) cos(l theta_o). Reflection/time-reversal
symmetries make these sine/cosine forms exact. From (3),

    A_2=-Omega epsilon^2/2+O(epsilon^4),
    B_2=Omega^2(1-Omega^2)epsilon^4/24+O(epsilon^6),
    A_4=Omega(1-Omega^2)epsilon^4/96+O(epsilon^6),
    B_4=Omega^2(1+Omega^2)epsilon^4/96+O(epsilon^6).         (4)

In particular both A_4 and B_4 are strictly positive on a sufficiently
small but FIXED positive-energy annulus. Omitting orbital distortion
changes B_4 by a factor of two and makes A_4 zero. The full actual orbit
is therefore essential, not a minor error in an affine clock.

## Smooth-band physical clock with fixed positive inertia

The exact coarea measure is dA=dE dtheta_o/omega(E), with the positive
orientation chosen for area. Put s0=f(E)cos(l theta_o), or its sine
quadrature, where f is smooth and supported in the selected annulus.
Since T preserves E, equations(1)--(2) become exact one-dimensional
smooth-band integrals. For the cosine column,

    theta(t)=-pi/Q integral chi_psi l A_l f cos(l omega t)dE,
    G(t)=-pi rho integral chi_psi B_l f/omega cos(l omega t)dE,
    S(t)= pi rho integral chi_psi l B_l f sin(l omega t)dE. (5)

The sine column has the corresponding opposite sine angle and cosine
spin. Select E_* strictly inside the fixed annulus, with chi_psi(E_*)>0.
Choose smooth nonnegative band profiles of width h tending to zero and
normalize their signed common amplitude so theta(0)=1. Formula(5) and
smoothness give, uniformly on every fixed finite time window and through
any fixed number of time derivatives,

    theta_1 -> cos(nu_* t), theta_2 -> -sin(nu_* t),
    G_i -> j_* theta_i, S_i -> j_* theta_i,t,
    nu_*=l omega(E_*),
    j_*=rho Q B_l(E_*)/[l omega(E_*) A_l(E_*)].             (6)

Symmetric bands give the usual O(h^2) error after normalization; an O(h)
bound suffices without imposing that symmetry. All actual finite
preparation norms are retained: amplitude normalization grows inversely
with band width, and Sobolev norms may diverge. Neither density nor a
mass coefficient is inferred from a vanishing norm. The tag chi and Q,
the selected E_*, the cell and nu_* are FIXED before h is taken small.

The l=2 route has negative leading j and is unsuitable for the requested
positive spin branch near this center. The l=4 route instead has

    j_* = rho Q (1+Omega^2)/[4(1-Omega^2)]
                          +O(rho Q epsilon_*^2)>0.        (7)

This is a derived physical spin/angle-rate ratio on a fixed positive tag.
It does not vanish as the observation error h tends to zero. Frequency
selection narrows preparation support, not the physical material tag.
The two real observed columns have nonzero limiting angle/rate determinant.

## Full action and spatial continuation remain explicit

The actual initial cotangent is pi=rho Du xi, not zero despite w=0.
Full phase and Jacobi energy are evaluated from

    Omega_12=integral (xi_1.pi_2-xi_2.pi_1),
    H=rho/2 integral [|xi_t|^2-|T xi|^2+xi.Hess(p)xi].    (8)

The actual axial displacement (1), full pressure Hessian, ambient support
and initial G stay in these forms. Equality of (8)'s physical-chart mass
with (7), the conserved-energy match and acoustic cross forms are further
calculations, not asserted by the observed ratio. Reviewed0210 supplies
actual phase controls only after those finite forms and cross constraints
are resolved. A positive observed spin is not itself an inherited inertia.

At purely axial Bloch wave number K_X, an initial curl-generated horizontal
displacement stays solenoidal without a transverse projection correction.
Its passive material scalar has the exact phase frequency
l omega(E)+K_X psi(E). Thus symmetric whole-field averaging alone leaves
a genuine K_X^2 t^2 clock-spread term. The common spatial second-jet
construction must cancel or derive that term; it is not erased by taking
h small. The full field-changing/angular-structure and EPS parent
obligations also remain active. This result selects a positive fixed-tag
clock route and exposes its next concrete action/current requirements.

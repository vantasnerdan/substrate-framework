# Nonlocal finite-time transfer of the physical moment match

This theorem transfers 0123 to the **actual** globally realized torus of
0120. The conclusion is an arbitrarily accurate physical spin/action
normalization on a fixed finite optical interval, with an exact initial
scalar moment match if the profile root is retuned. It is not exact
all-time closure of a finite Euler sector.

## 1. Parameter order and the exact affine intermediate

Fix an optical interval 0<=t<=T with Omega*T bounded. The torus construction
permits a local core jet arbitrarily close, on this interval, to rotation
with nonzero Omega in its **fixed physical frame**, after following the
core translation. Choose that geometric/approximation error eta_g first.
Fix one actual smooth global stationary field u and its periodic core
trajectory X(t). This field and every finite derivative bound used below
are now fixed, before choosing a packet size. The decaying EPS field has
finite such bounds. No estimate below assumes a uniform bound on the
whole family of global approximants.

Use the translated, nonrotating position coordinate x=world_position-X(t).
The actual background is

    w(x,t)=u(X(t)+x)-X_dot(t),  w(0,t)=0,
    A(t)=D u(X(t)),  tr A=0.

The acceleration of the translated origin is absorbed in a linear
pressure. Differentiating the actual Euler equation on the trajectory
gives

    A_dot+A^2=-D^2 p(X(t))/rho.                         (1)

Thus w_aff(x,t)=A(t)x is an EXACT, generally time-dependent Euler solution,
with quadratic pressure. It is not an arbitrary linearization declared to
be a new stationary solution. Its vorticity omega(t)=curl A is exactly
the actual core vorticity, and omega_dot=A omega. Let M_dot=A M, M(0)=Id;
det M=1. The affine and actual backgrounds have the same center, first
jet and core vorticity, and on a fixed ball |x|<L,

    |w-w_aff|<=K2 |x|^2,  |D(w-w_aff)|<=K2 |x|.         (2)

All constants L,K2, sup|u| and sup|D^j u| are fixed properties of this
chosen field. In contrast, directly comparing the global EPS field to
uniform rotation would leave a fixed curvature error multiplied by a
global growth constant chosen only afterward. Equation (1) removes that
circularity.

## 2. The exact affine Fourier propagator includes pressure

For linearized Euler about A(t)x, a Fourier covector follows
k(t)=M(t)^(-T) k0. The velocity amplitude solves

    b_dot=-A b+2 k(k·A b)/|k|^2,   k·b=0.               (3)

The material amplitude solves a_dot=A a+b with
b=P_k(a cross omega). These are exact equations for the affine flow,
not a high-frequency asymptotic system. In material coordinates y=M^-1 x,
write a=M q and G=M^-1 M^-T. The displacement generator has symbol

    q_dot=L_t(k0)q,
    L_t(k)=[G-(Gk)(Gk)^T/(k·Gk)] [q -> q cross omega(0)]. (4)

The matrices in (3)-(4) are smooth homogeneous degree zero functions of
nonzero k. On the finite interval their sphere derivatives have finite
bounds depending only on the local affine history and Omega*T. Their
propagators have the same property. The inverse Fourier kernel of a
degree-zero smooth multiplier, away from its delta distribution at zero,
has |D^j K(x)|<=C_j |x|^(-3-j). This follows by radial integration and
integration by parts on annuli; the constant is bounded by finitely many
sphere derivatives of its symbol. Applying it to a curl adds one power
of decay. Pressure is inside these symbols throughout; no local projector
is substituted for Leray.

Take 0123's compact potential and displacement, with radial scale a,
core scale ell, tail scale R, ell<<R<<a:

    P_h(q)=-(Jq,0) g(x_perp/a) H(z),  H'=h,
    Xi_h(q)=curl P_h(q),  ||P_h||_1<=C |q| a^2 B,
    B=integral z h(z) dz ~ ell^3.                       (5)

The matched small root satisfies integral h^2~ell^3 and, for each fixed
j, ||h^(j)||_2<=C_j ell^(3/2-j). These estimates include its disjoint
low-amplitude tail: its amplitude is O(ell^3/R^3). They hold uniformly
for the root in a fixed relative neighborhood of the simple root.

The exact affine Kelvin solution initialized by Xi_h has, outside the
advected support ball |x|>C_T a,

    |D^j xi_aff|<=C_T a^2 B |x|^(-4-j),
    |D^j v_aff|<=C_T |omega(0)| a^2 B |x|^(-4-j).        (6)

The curl potential in (5), rather than an estimate of Xi's unsigned
first moment, is what gives the useful a^2 B coefficient. It preserves
the cancelling return circulation and the entire exterior velocity.

The affine extension of 0123's separated-collar estimate is also exact.
In y coordinates the planar profile q(t)h(z) obeys (4) at k=e_z, for
EVERY axial shape h. Its compact divergence completion is Xi_h(q(t)).
The difference between its evolution and (4) is a degree-one symbol
which vanishes on k parallel to e_z. A smooth partition of the frequency
sphere factors it as sum_(j=1,2) k_j T_j(t,k), where each T_j is smooth
degree zero. Near either axial pole this is the integral form of the
Hadamard factorization; away from the poles divide by |k_perp|^2.
Therefore the residual is a sum of degree-zero multipliers applied to
(partial_j g)H, supported before projection in the radial collar.
Composing with the exact propagator retains the degree-zero kernel bound.
Its source L1 norm is O(a B). The material tag and collar remain separated
by c_T a after the same affine map M. Hence on the transported tag

    |xi_aff-xi_planar_completed|<=C_T B/a^2,
    |D^j difference|<=C_T B/a^(2+j).                   (7)

The time derivative version has the corresponding local frequency
factor. Integrating the actual tag spin, including the boundary and
shape contributions, gives relative O(R/a), exactly as in 0123. The
core-angle error is O(B/a^3). No global EPS norm appears in (7).

The local affine matrices and physical tag maps depend continuously on
the actual core jet. Comparing these material-coordinate matrices with
the uniform rotating matrices costs O(eta_g) on this interval, without
an a/ell dephasing loss: the profile argument is the transported coordinate
k(t)·x, and the tag is transported by the same M. Holding h(z) fixed in
the laboratory while tilting k would not have this property. Thus the
exact affine solution inherits the matched physical action/spin rows
with error O(eta_g+R/a+B/a^3).

## 3. Actual EPS evolution: an explicit global remainder bound

Initialize the actual solution by the same compact Xi_h and its EXACT
Kelvin velocity

    v(0)=P[Xi_h cross curl u].                         (8)

The affine initial velocity is P[Xi_h cross omega(0)]. Their difference
has L2 norm at most C K2 a^2 ell^(3/2), by Leray contraction and (2).
This includes the initial noncompact pressure tail.

Put d=w-w_aff. Substitution of the exact affine velocity in the actual
Euler equation gives residual

    P[d·grad v_aff+(D d)v_aff].                        (9)

Choose a so small that C_T a<L/2. Inside C_T a, the global H1 bound from
(3) is ||v_aff||_(H1)<=C_T a ell^(1/2), with the fixed physical frequency
included in C_T. Equation (2) bounds (9) there by C K2 a^3 ell^(1/2).
Between C_T a and L, (6) bounds its L2 norm by C K2 a^(1/2) B: explicitly
the residual is bounded pointwise by C K2 a^2 B |x|^-3 and its squared
radial integral is O(a^-3). This is smaller than the preceding term.
Outside L, use the actual GLOBAL bounds

    |d|<=2 sup|u|+sup|A| |x|,
    |D d|<=sup|Du|+sup|A|.

Together with (6), this gives the finite contribution

    C_T a^2 B [sup|u| L^(-7/2)
                     +(sup|Du|+sup|A|)L^(-5/2)].       (10)

Thus the whole-space residual, not a locally truncated pressure, obeys

    ||residual||_2<=C_u,T [a^3 ell^(1/2)+a^2 ell^3].    (11)

The standard divergence-free Euler energy estimate has growth at most
exp(T sup|Du|), since the translated transport is skew in L2. The same
estimate for the Lin reconstruction, with (9)'s analogous displacement
residual, gives

    sup_t (||v-v_aff||_2+||xi-xi_aff||_2)
        <=C_u,T a^3 ell^(1/2).                        (12)

Frequency/unit conversion factors in the sum are fixed, not silently
set equal to one. The smaller initial and far-field terms have been
absorbed for small a. The possibly large constant C_u,T is now fixed;
the powers of a below really can dominate it.

Higher regularity is used only as an interpolation bound, not as a false
uniform local-to-nonlocal inference. Smooth Euler transport and the
Leray H8 bound give for each of the two exact solutions

    ||xi||_(H8)+||v||_(H8)<=C_u,T a ell^(-13/2).

Interpolating (12) with this bound at s=11/4 (theta=11/32), followed by
H^(11/4)(R3) embedding into C1, yields

    ||xi-xi_aff||_(C1)+||v-v_aff||_(C1)
       <=C_u,T a^(37/16) ell^(-61/32).                 (13)

The physical core angle and its first time derivative are controlled by
this estimate. At the followed center w=0; the linearized vorticity
equation there uses curl v, v and the background first/second jets, not
a second spatial derivative of the displacement error. This is why (13)
suffices for the phase-action coefficient and angular-velocity maps.

## 4. Actual transported tags, not a fixed integration window

Let D(t) be the image of the chosen initial cylinder under the actual
unperturbed material flow; D_aff(t)=M(t)D(0). These flows preserve volume.
Equation (2) and finite-time ODE comparison give a position difference
O(C_u,T a^2) on this tag. Evaluating v_aff or xi_aff at the two positions
costs at most C a^2 times their local H1 norm in L2, of the same order as
(12). Centers and mean velocities are subtracted separately.

For ANY transported parcel, the useful full observation formula is

    delta S=rho integral_D [
       (xi-delta X) cross (u-U_D)
       +r cross (v+(Du)xi-delta U_D)],                 (14)

where r=x-X_D, delta X is the parcel mean of xi, and delta U_D the
variation of its mean velocity. This follows by differentiating the
material integral before any boundary integration by parts. In a
stationary invariant reference tag it reduces to the 0123 formula;
for a generally transported tag one should use (14), not discard the
Eulerian flux of its moving boundary.

On a cylinder of transverse scale a and axial scale R,
||r||_(L2(D))<=C a^2 sqrt(R). Apply (12), the flow comparison and (14),
and retain the smaller background Taylor terms. With the normalization
S_*=rho*Omega*a^2*ell^3 supplied by the matched profile,

    sup_t |S-S_aff|/S_*
       <=C_u,T a^3 sqrt(R)/ell^(5/2).                  (15)

The same comparison controls the first translation and symmetric shape
moments at their explicit natural scales. They are not identified with
spin or removed from the centroid/ambient momentum current. The full
compact packet may still have cancelling spin outside D(t).

## 5. An executed hierarchy, exact initial root, and finite-time claim

Choose a fixed reference length L0, put delta=a/L0, and set

    ell=a*delta^(1/8),   R=a*delta^(1/16).              (16)

As a shrinks, ell/R->0 and R/a->0. Equations (13),(15) become, up to
fixed dimensionful constants,

    C1 core error = O(delta^(43/256)),
    relative material spin transfer = O(delta^(23/32)).

The separated-collar error is O(delta^(1/16)); its core-angle error is
O(delta^(3/8)). The direct Lundquist dephasing parameter lambda*a^2/ell
also tends to zero, but the proof uses the stronger exact affine
reference rather than relying on that one parameter. Thus, for any
chosen tolerance eta, first take eta_g sufficiently small in the actual
torus construction, fix the resulting global field, and then take a
sufficiently small that every C_u,T times the displayed power is <eta.

For the physical scalar angle selected in 0123, its initial core row is
unchanged by the tail: h'(0)=1 and Xi_h(0)=0. Its true initial spin is a
linear functional of Xi_h and its exact Kelvin velocity (8); the exact
global KKS is a quadratic functional of Xi_h. Therefore the actual
initial condition det([angle row; spin row])=beta is a quadratic equation
in the profile parameter, with all coefficients defined by the actual
field, full P, and material tag. Unlike the uniform problem, overlap of
primitive H terms can add a cross coefficient; it is retained. Rescale
the tail parameter by ell^3/R^3 and the equation by a^2 ell^3. The
uniform simple root has a fixed nonzero derivative in these variables.
The preceding estimates and their parameter derivatives give continuity
to the actual coefficients. The implicit function theorem, or signs on
two fixed neighboring root values, supplies an exact actual initial
root. No frequency is fitted and no fixed Kelvin rotor mass is added.

Over the entire optical interval, let Phi be the ACTUAL core angle and
S its ACTUAL transported-tag spin. After this initial normalization,

    |S-j Phi_dot| <= C [eta_g+delta^(1/16)
                         +C_u,T delta^(43/256)] S_*,  (17)

with the natural unit-amplitude angular-velocity normalization understood;
at zeros of Phi_dot this is an absolute uniform bound, not division by a
vanishing observable. The exact KKS of the two genuine Euler solutions
is conserved. The physical two-angle map and its derivative are close
to the fixed rotating reference. Pulling back the action with the exact
0115 formula therefore gives a positive nondegenerate phase action close
to that reference, and its scalar elimination has positive inertia.
Canonical momentum and measured S agree initially and differ afterward
by the controlled bound, not by an asserted identity. A near-identity
time-dependent canonical normalization retains this observable remainder
and the physical clock. There is no winding/frequency branch selection.
Only the first observation derivative is controlled here. Positivity is
asserted for the TWO-angle phase Hamiltonian and the scalar kinetic
coefficient. Differentiating a time-dependent scalar momentum shift can
involve a second observation derivative; no separate scalar potential
coefficient or second-order force remainder is inferred without that
additional estimate. The retained first-order phase/action remainder is
the direct downstream object.

This proves a controlled finite-time spin/action normalization for actual
Euler solutions on the same invariant tube. Exact all-time mechanical
matching, an invariant finite-dimensional ansatz and a finite-k continuum
action do not follow and are not used as premises.

## 6. Periodic variant and scope

The estimate from (1) onward needs a fixed smooth final field with finite
global norms, local core jet/persistence margins, and the relevant Leray
operator. A periodic realization supplied at the actual scope of 0116
has these properties too. In a small coordinate ball the periodic Leray
kernel is the whole-space singular kernel plus a smooth bounded kernel;
the latter contributes an additional O(a^2 B) term, absorbed in (10)-(12).
More explicitly, cut off the known affine reference at a fixed radius
inside a cell and apply periodic Leray. On this cutoff collar all its
derivatives are O(a^2 B) by (6), so the cutoff, divergence and equation
commutators are O(a^2 B) in each fixed norm. The smooth difference of the
two Green kernels has the same bound after integrating the compact curl
potential by parts. These add to (10)-(12). One does not pretend the
unbounded affine velocity is a periodic background. A spatial
translation law of this fixed periodic stationary field gives identical
cell density, but it does not make the evolving perturbations stationary
or remove their physical phase error. Rescaling the source eigenvalue is
a declared geometric choice made before the final field is fixed.

## 7. Genuine small-amplitude smooth Euler flows

The preceding action and moments are quadratic/first-variation objects,
but their finite-time trajectories are realized as variations of actual
smooth Euler flows. After fixing the background, packet geometry and
optical interval, take the volume-preserving map exp(e Xi_h) and push
forward the background vorticity and the initial material tag by this
same map. More precisely, pull the initial velocity one-form back by
the inverse map and apply the complete Leray projection. This fixes the
Kelvin periods as well as the curl. On a periodic domain its harmonic
mean component is the one determined by this pullback; it is retained,
not independently reset to its background value. On R3 the projected
compact one-form difference fixes the decaying perturbation. Since Xi_h is compact
and smooth, the velocity difference is smooth and H^s for every finite
s, with expansion e*v(0)+O(e^2). This is an exactly isovortical initial
curve, not a Kelvin constraint imposed only to first order.

For the perturbation w_e of the fixed steady solution the actual equation
is w_e,t+P[u·grad w_e+w_e·grad u+w_e·grad w_e]=0. Integer H^s energy,
s>5/2, gives d||w_e||_(Hs)/dt<=C_u||w_e||_(Hs)+C||w_e||_(Hs)^2.
Choose the disturbance amplitude e last, so this solution exists on the
fixed interval. Applying the corresponding difference estimate at one
lower derivative, with the smooth initial Taylor remainder, gives
||w_e-e*v||<=C_u,T,packet e^2. Higher initial smoothness supplies the
C1 velocity/material-flow control used by the observations. The same
transport equation advances the tag, so its angle, spin and shape
variations differ from the constructed first variations by O(e).
Take e sufficiently small relative to their nonzero normalization scales.
This realizes (17), with an additional arbitrarily small nonlinear error,
by genuine finite-time smooth Euler flows. It does not assert global
regularity in time or an exact finite-dimensional nonlinear ansatz.

`route_verdict: established` for the conditional-on-0120/0123 construction
and quantitative finite-time estimate above.
`evidence_scope: ACTUAL_EULER_NONLOCAL_FINITE_TIME_PHYSICAL_NORMALIZATION`.
The analytic proof supplies the operator bounds; the exact oracle checks
the affine Euler/Kelvin algebra, interpolation powers and hierarchy. It
does not purport to numerically solve the global Euler PDE or its spectrum.

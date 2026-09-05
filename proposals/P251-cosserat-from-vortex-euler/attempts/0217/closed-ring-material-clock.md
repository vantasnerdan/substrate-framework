# Positive actual material tilt on the same closed Beltrami ring

The selected route is A: actual homogeneous Euler/Lin transport. This
constructs an observed material-tag clock with unchanged Eulerian velocity
and vorticity. It is not advertised as a changing vortex-structure mode.
The geometric input is the independently reviewed 0211 global steady ring.
All frames below are its actual steady laboratory frame, including its
uniform far velocity. No separate boost or Floquet winding is used.

## 1. Exact material and stationary-tag equations

Let u be that smooth steady Euler field, omega=curl u, T=u.grad, and
D=Du. For any smooth compactly supported divergence-free initial xi,
transport it by the actual volume-preserving flow:

    xi_t+[u,xi]=0,       w=xi_t+Txi-Dxi=0.                 (1)

The Eulerian linear velocity and pressure variations are identically zero.
The material map varies, and its Kelvin label data vary accordingly. This
does not impose a fixed-circulation leaf on an arbitrary initial xi.
It is a genuine homogeneous Lin solution, not a projected oscillator.

Choose a smooth nonnegative tag chi=chi(phi), compactly supported in the
open constant-curl region. It is stationary because Tphi=0. Its perturbed
density is q=-xi.grad chi and obeys exactly q_t+Tq=0. Thus for every
fixed Euclidean test F,

    delta integral chi F = integral q F,
    partial_t delta integral chi F = integral q T F.      (2)

The actual canonical cotangent is pi=rho Dxi. The complete homogeneous
phase and conserved material Jacobi energy are

    Omega_12=-rho integral omega.(xi_1 cross xi_2),
    H[xi]=rho/2 integral omega.(xi cross xi_t).            (3)

These are the direct homogeneous restrictions of the full Jacobi action,
not the energy of w alone. To obtain the second identity, expand
rho[|xi_t|^2-|Txi|^2+xi.Hess(p)xi]/2, insert (1), use
Hess(p)=-D^2-TD, and integrate T by parts. The actual compact support
removes the boundary term. Root's independent 0216 derivation supplies
the same general identity; the column calculation below independently
retains every component in its application to this geometry.

## 2. Physical global tilt, spin and centroid

The reference tag is axisymmetric and z-even. Its mass M_tag is positive,
its centroid is zero, and its covariance has I_xx=I_yy. Put

    D_tag=I_xx-I_zz>0,
    Q=delta integral rho chi (x+i y)z,
    theta_+=Q/(i D_tag).                                 (4)

For a rigid physical infinitesimal rotation a about a horizontal axis,
Q=i D_tag(a_x+i a_y). Hence (4) has EXACT unit rigid-rotation response.
It is an actual global Euclidean tilt, not a freely assigned poloidal
angle or an angle whose stationary reference quadrupole vanishes.

The complete first-order mechanical spin is

    S=integral rho q x cross u.                           (5)

Indeed the reference mean momentum vanishes:
integral chi u=integral T(chi x)=0. Therefore the moving-centroid
correction to spin is zero at first order, but the centroid itself and
its derivative are retained as

    delta X=M_tag^-1 integral rho q x,
    delta P=partial_t(M_tag delta X)=integral rho q u.     (6)

The symmetric shape row is equally the actual integral of q x_i x_j.
Since w=0, the complementary ambient tag has the opposite density,
momentum, spin and shape variations. No uncancelled complete-fluid
Eulerian current is inferred from (5) or (6).

## 3. Full homogeneous column phase and energy

First fix the small-core size and pass to the large-R straight limit of
0211. Use right-handed local cylinder coordinates (s,theta,zeta), with
zeta along the ring tangent; locally the global vertical coordinate is
-s sin(theta). The exact inner column is

    u=V(s)e_theta+W(s)e_zeta,
    Omega_p=V/s, Z=2Omega_p+s Omega_p',
    Z=lambda W,       W'=-lambda s Omega_p,
    Omega_p(0)=Omega_0=A lambda^2/2>0.                   (7)

For the meridional scalar generator A_0(s) exp(i m theta+i k zeta),
the displacement amplitude is

    xi=(-i m A_0/s,A_0',0),       q_f=m Omega_p+kW.

It is solenoidal without a radial wall. A_0 is smooth, regular as an
m=1 scalar at the origin, and compactly supported in the inner region.
For its real cosine/sine pair the exact homogeneous rate is

    xi_t=(-i q_f xi_s,
          -i q_f xi_theta+s Omega_p' xi_s,
          W' xi_s).                                     (8)

The last component matters. With alpha=m A_0/s and b=A_0', the angular
average of omega.(xi_cos cross xi_cos,t) is

    Z q_f alpha b + (W'^2+Z s Omega_p')alpha^2/2.

Integration by parts, retaining (7), gives the complete forms per axial
length L:

    beta=Omega_cos,sin=rho pi L m integral Z' A_0^2 ds,
    H_cos=H_sin=-rho pi L m integral Z' q_f A_0^2 ds.      (9)

The off-diagonal energy vanishes by the angular integral. In particular
the energy MATRIX is twice the displayed unit-column energy. This is
not the matrix obtained by silently assigning generator -q_f to a
closed two-dimensional phase plane: the toroidal/shear components of
the true displacement need not lie in that plane.

Choose local m=-1 and global ring harmonic n=-1, so k=-1/R and
q_f=-Omega_p-W/R<0 on the chosen inner disk. Since Z'<0 there,
beta>0 and the full H is strictly positive. All fields are real after
taking the cosine/sine pair; negative harmonic number is a physical
polarization, not negative density or a selected logarithm branch.

## 4. Literal moment calculation and positive physical spin

Choose fixed smooth dimensionless profiles c(t)>=0 and f(t)>=0 supported
in t<1, with c nonincreasing, positive near zero, and f nonzero where
c'<0. Take c radial-smooth at the origin and flat at its outer edge.
Set chi_a(s)=c(s/a), A_0(s)=s f(s/a), with lambda a small but positive.
The fixed positive radius a is selected before R. Define

    I_c=integral s chi_a ds>0,
    I_A=integral s A_0^2 ds>0,
    B=integral s A_0 chi_a' ds<0.                        (10)

In the actual ring use the exactly solenoidal meridional preparation

    S=-R A_0(s) exp(-i theta-i varphi),
    xi=(-S_z/r,0,S_r/r),                                (11)

with the smooth Cartesian expression A_0(s) exp(-i theta)
=(x-i y_local) f(s/a). The minus sign accounts for y_local=-z in the
global meridional curl convention. The initial generator is compact in
the physical core. The
factor R/r is part of the field, not discarded before taking moments.
On the actual ring choose chi as the corresponding smooth function of
phi, converging to chi_a; this makes the reference tag EXACTLY stationary.

In the straight-core limit, retaining the toroidal volume Jacobian and
the R/r factor in (11), the complex mode gives

    Q= C (a_1-i a_2),       C=rho pi^2 R^2 B,
    S_+=rho pi^2 R^2 integral s A_0 chi_a'
                          (Omega_p-W/R) ds (a_1-i a_2),
    M_tag delta X_+=-i C/R (a_1-i a_2),
    D_tag=2rho pi^2 R^3 I_c+O(rho R a^4).                (12)

Transport multiplies the scalar density at each s by
exp(i[Omega_p(s)+W(s)/R]t). Thus the measured frequency is the positive
LAB quantity nu=Omega_p+W/R, not merely a poloidal comoving frequency.
The initial scalar angle row is theta_x=-c_theta a_2 and its rate is
theta_x,t=c_theta nu a_1 at leading order, where c_theta=C/D_tag.
Consequently the actual scalar phase mass is

    M_phase=beta/(nu c_theta^2)>0,
    S_+=j_tag theta_+,t,
    j_tag=D_tag (Omega_p-W/R)/(Omega_p+W/R)>0             (13)

at the leading narrow-core order. The centroid in (12) is nonzero and
linked to the same amplitudes; it is not independent translation data.
Its momentum/spin ratio is O(1/R), not an omitted zero current.

For finite a and R all rows mean their exact integrals (2), not an
assumption that different radial clocks coincide. Axisymmetry makes
the two-angle observation a single complex function of time multiplying
(a_1-i a_2). Reflection across z=0 makes Q(0) real, Q_t(0) imaginary,
and S_+(0) real in this registration. Thus the initial spin/rate ratio
is a real scalar, and is strictly positive by the displayed limit.

## 5. Positive tag fraction and controlled finite-time normalization

For the reference tag at full fraction, the preceding derived moments
give

    M_phase/j_tag
      =4lambda^2 I_c I_A/B^2 [1+O((lambda a)^2)+o_R(1)]
      =O((lambda a)^2).                                 (14)

The last equality follows from I_c=a^2 c_0, I_A=a^4 a_0,
B=a^2 b_0 with fixed c_0,a_0>0 and b_0<0. It is not a fitted
coefficient. For sufficiently small fixed a and then large R this
ratio lies strictly between zero and one.

Multiplying the ACTUAL stationary material fraction chi by the exact
number alpha=M_phase(0)/j_tag(0) leaves theta and its time derivatives
unchanged, multiplies its mechanical spin, mass, centroid moment and
shape by alpha, and does not change the complete-fluid phase. Therefore
the same fixed positive tag satisfies EXACT initial spin/phase-inertia
normalization. This is a fractional material tag with density rho chi,
bounded by the fluid density, not a redefinition of the microscopic rho.

Here is the actual finite-time license. For any fixed T and finite
number r of time derivatives, Taylor's Bessel equation in the chosen
core gives Omega_p=Omega_0[1+O((lambda a)^2)] and the same controlled
derivatives. Formula (2) then bounds the normalized moment histories
and their time derivatives by those of the common rotation at Omega_0,
with error O_{T,r}((lambda a)^2)+o_R(1). For the actual ring, the
0211 local C^k convergence in axisymmetric coefficients, including
the exact R/r factors, controls the ordinary flow and its finite
derivatives on this invariant core. Duhamel for this LOCAL transport
equation gives o_R(1) after a is fixed. There is no omitted nonlocal
pressure comparison: w and its pressure are exactly zero on both
backgrounds, while the full material displacement is transported.

Choose R after all the fixed-a phase and tag denominators so the
phase, energy and normalized observation errors are smaller than their
strict limits. In particular the small phase beta is compared with its
own value proportional to R lambda^2 Omega_0 I_A, not with an order-one
unscaled error. The required R can depend strongly on a and T; none of
these limits is interchanged.

The actual two real scalar angle/rate columns therefore have a positive
Wronskian on [-T,T]. Pulling the conserved phase through THESE rows
gives positive mass and stiffness, uniformly close to
M_phase(0) and M_phase(0) Omega_0^2. The same physical spin is close to
M_phase(0) theta_t with controlled relative error, and agrees exactly
at the reference time by alpha. All time-dependent mass, amplitude,
frequency and shape connections remain in the exact moving chart.
This is an actual finite-time positive angular action, not an assertion
of an exact autonomous oscillator at finite a and R.

The separate conserved energy (9) is positive and has leading MATRIX
2beta Omega_0 I. Its factor of two relative to the simple moving-angle
oscillator energy is retained. A same-ring invisible energy control
would need its own construction to change that conserved form. The
0210 control lives on another geometry and is not imported here by name.

## 6. Scope, alternative route and next exact construction

This earns a positive stationary-tag tilt clock, actual mechanical
spin/centroid/shape, nondegenerate phase, and positive full energy on
the SAME globally smooth ring with a literal nonzero Beltrami core.
The normalization uses one fixed positive material tag and actual full
homogeneous Lin histories. Its moving-angle action and conserved energy
are both derived and distinguished. Accuracy has been ordered at fixed
finite time; no nonzero-density continuum or all-time isochrony is inferred
from shrinking a or increasing R.

The Eulerian vorticity is unchanged. Thus route A establishes a material
observation response, not by itself the parent's field-changing optical
structure. Candidate B remains a distinct active construction: the exact
axisymmetric angular-momentum perturbation delta l=r w_varphi obeys

    (partial_t+u_p.grad)delta l+w_p.grad F(psi)=0.

Unless the meridional response is derived too, a proposed passive delta l
changes the centrifugal forcing and is not an exact independent Euler
sector. The actual meridional pressure/vorticity equation must be retained.
The parallel 0213 curved-mode construction addresses a genuine field-
changing branch on this same ring; it is not replaced by the label result.

Failure-derived spatial candidate, registered here before execution:
a compact homogeneous generator on a contractible closed ring can be
periodized with a common lab Bloch phase, because its exact transported
support stays within the invariant tube. A prescribed change of its
actual flux support with K can change the observed transport clock;
that requires new derivative-level moment and action estimates and is
not supplied by integer n alone. The root's spatial/interface continuation
can consume the present finite-time rows without assuming those estimates.

Route verdict A: established as the stated positive material-tag response.
Evidence scope: actual Euler/Lin finite-time construction with controlled
physical normalization, not a changed vorticity mode. Candidate B and
the parent same-core coupled continuum remain active.

# A nonempty actual Euler spin-to-centroid response window

This is the explicit first-spatial-jet join of0124/0126 and0117.
It supplies actual coupled finite-time motions on a smooth stationary
Euler background, including a genuine small nonlinear realization.
It does not establish the remaining autonomous second-gradient pencil.

## 1. A physical preparation on actual identical cells

Use the periodic variant of0124 on a fixed smooth stationary Beltrami
field with one robust invariant torus and its material tag per chosen
cell. Its finite radius, cell volume, profile and field are fixed
geometric inputs. Rotate and translate the COMPLETE background, tags
and preparation to form the stationary isotropic law. Time reversal
pairs the complete backgrounds and independently prepared histories;
it does not impose a common microscopic trajectory. The original
knotted-tube construction remains available separately;0120's rotor
torus is unknotted and is not renamed a prescribed knot.

For the reference matched optical profiles, prepare the transverse
angles q_n(0)=P_n Phi0 and q_n_dot(0)=0, |Phi0|=1. The physical macro
angle is reconstructed by (E P_n)^-1 E q_n, as in0126. Let the derived
positive reference inertia density be j=(2/3)nu M and choose
T=pi/(2|Omega|). The actual0124 flows can be chosen to satisfy, uniformly
on this fixed interval,

    |S(t)-j Phi_ref_dot(t)| <= eta j |Omega|,
    Phi_ref(t)=cos(|Omega|t) Phi0,

for any fixed eta>0. Hence

    Delta S=-j |Omega| Phi0+e_S, |e_S|<=2eta j|Omega|.       (1)

This is an absolute bound at the natural nonzero spin scale, valid also
at instants with zero angular velocity. The actual angle, action and
their measured corrections remain attached; Phi_ref is a comparison,
not a replacement for those observations. Parameter order is that of
0124: select geometric error, fix the final global field, then select
the finite packet size, then the nonlinear disturbance amplitude.

An angular preparation may initially have a nonzero full-fluid mean
velocity. Remove it with its ACTUAL Galilean solution if a zero-mean
Bloch initialization is desired. A boost V has

    xi_G=t V, v_G=V-t(V.grad)u,

so its mean is V, its actual centered material spin is zero, and its
material vorticity-direction variation is zero. Choosing V opposite the
initial mean therefore preserves (1). It changes a declared harmonic/
translation datum, not the compact fixed-Kelvin perturbation's periods
by stealth. At nonzero slow wave number the complete Bloch projection
and initialization from0116/0125 are used, not a boost asserted to solve
a spatially modulated equation. This keeps the ray-limit mean projector
from introducing an unprepared direction-dependent zero mode.

## 2. Arbitrary slow directions have a genuine Euler realization

There is no common finite supercell for every Haar rotation and a fixed
laboratory k. That arithmetic issue does not obstruct the finite-time
flow. Write a rotated periodic field as u0(Bx), where B is an invertible
3x3 reciprocal-cell matrix. Add the slow phase phi=k.dot x and set

    theta=(Bx, k.dot x) in T^4,
    D_i=sum_(a=1)^3 B_ai partial_theta_a+k_i partial_phi.

The lifted velocity U(theta,t) has THREE physical components. For a
Fourier label (m,n) in Z^3 x Z, its physical wave vector is
p=B^T m+n k. Define P_D on this Fourier coefficient by

    P_D(p)=identity-p p^T/|p|^2 when p!=0;
    P_D(0)=identity.                                      (2)

This includes all resonant zero PHYSICAL modes, not just the zero label.
The symbol is orthogonal with norm at most1 and commutes with every
theta derivative. Thus it is bounded on every H^s(T^4), uniformly even
when nonzero physical p approach zero arbitrarily quickly. No lower
Diophantine bound on |p| and no inverse pressure estimate are used.

Solve the lifted, fully projected equation

    U_t+P_D[(U.dot D)U]=0, D.dot U=0.                      (3)

The four-component transport velocity is (B U,k.dot U), whose ordinary
four-dimensional divergence is D.dot U=0. Since P_D is orthogonal and
commutes with derivatives, the standard differentiated energy argument
reduces to the transport commutator. On the4-torus, integer s>=5 gives

    d||U||_Hs/dt <= C_(s,B,k) ||U||_Hs^2.

For a perturbation W of the fixed smooth lifted stationary background,
the corresponding bound is C_u||W||_Hs+C||W||_Hs^2. Smooth Fourier
Galerkin truncation commutes with P_D and preserves the constraint.
These uniform bounds, compactness at one lower derivative, and the same
energy estimate for differences give local existence and uniqueness;
applying the estimate at successive s propagates smooth initial data.
This uses the usual torus Sobolev product/embedding inequalities, not
ellipticity of D as an operator on four variables. The restricted field
u(x,t)=U(Bx,k.dot x,t) is smooth and solves physical3D incompressible
Euler, by the chain rule.

The omitted force is a gradient in physical space: each nonzero Fourier
coefficient is parallel to p, hence its physical curl vanishes; its
zero-physical-mode coefficient is zero by(2). Its smooth restriction on
simply connected R3 has a smooth pressure potential. That potential
itself need not be a bounded quasiperiodic function on T^4; only its
actual gradient enters Euler. Dividing by arbitrarily small |p| to claim
a smooth pressure on the4-torus would be an unjustified extra step.

Initial compact-cell Kelvin fields and their Bloch preparation are
smooth periodic functions of theta with finitely many initial phi
harmonics. Real +/- slow harmonics are paired. Choose the exact
volume-preserving physical initial map generated by that smooth
quasiperiodic displacement, apply the Kelvin one-form pullback and
P_D, and push the material tags by the same map. The corresponding
lifted transport has zero4D divergence, so this map is also constructed
by a smooth flow on T^4. Its expansion supplies the linear Kelvin
initial field, plus the separately declared Galilean/harmonic data.
Equation(3) retains every nonlinearly generated slow harmonic.

At fixed background, packet and finite interval, take the amplitude e
last. The perturbation estimate ensures existence through T, and the
difference from e times the linearized solution is O(e^2) in H^(s-1).
Take s>=6 to control the C1 material observations after restriction.
Constants are finite uniformly over proper rotations (a compact set)
and k in a fixed bounded set because the norms of(2) stay bounded.
This realizes the full preparation as genuine smooth Euler flows,
not merely a rational-wavevector approximation or a prescribed field.

## 3. The actual centroid signal dominates all retained errors

The complete isotropic linear preparation has the finite first Bloch
jet supplied by0116/0125, with its mean data prepared as above.0117
then gives, for k perpendicular to Phi0,

    Delta J_H=-(i/2)k cross Delta S+R_k,
    |R_k|<=C_k |k|^2,                                   (4)

where J_H contains all actual tag centroids AND the entire ambient.
C_k includes full stress response, shape-rate, tag radii and time
integration. Isotropy cancels only the coherent k=0 symmetric tensor
rows, not microscopic shape fluctuations or their finite-k response.
The complex expression is the Fourier amplitude of a real paired wave.

Fix eta<=1/16. All geometric parameters in j and C_k are now fixed,
with j>0 finite. Choose positive upper bound C_k and a regularity radius
k0 for the actual first-jet family. The interval

    0<|k|<=min(k0, j|Omega|/(16 C_k))                    (5)

is nonempty. Equations(1),(4) give

    |Delta J_H| >= |k|j|Omega|(1/2-eta)-C_k |k|^2
                 >= (3/8)|k|j|Omega|.                 (6)

Let C_nl>0 bound the difference between the amplitude-normalized actual
nonlinear current and its first variation at these now fixed inputs.
Choosing additionally

    0<e<=|k|j|Omega|/(16 C_nl)

preserves a lower bound (5/16)|k|j|Omega| for the normalized actual
current change. This is a physical nonzero translational response,
derived from the same material spin rather than an unrelated stress seed.
The corresponding centroid-plus-ambient velocity is J_H/rho, with the
TOTAL density rho. The initial Galilean datum can separately fix its
constant velocity part; it is not an additional rotor mass.

All intervals are finite but genuine. As packet size tends to zero,
j may also tend to zero and the admissible k/amplitude window may shrink.
No nonzero limiting microinertia, uniform-in-geometry modulus or
autonomous second-gradient action follows from(5). For each selected
finite geometry, however, the same actual stationary Euler realization
has nonconstant mechanical spin and a nonzero actual coarse centroid
response on the specified optical/coarse window. The remaining positive
spatial stiffness and full constitutive identification are the next
parent obligations, not hidden premises of this result.

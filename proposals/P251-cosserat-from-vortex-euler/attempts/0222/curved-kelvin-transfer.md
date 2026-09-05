# Actual field-changing Kelvin clock on the same closed Beltrami core

The new route is the compact Kelvin lift with an actual whole-space
pressure estimate. The background is precisely the reviewed 0211 ring,
in its steady frame. This is not a new globally constant-curl field.
The displacement support lies inside its literal-lambda region; its
Euler velocity and pressure may extend throughout space.

## 1. Exact column lift and the complete cotangent

In the straight reference column use the right-handed coordinates of
0217 and write u=(V e_theta,W e_zeta), T_p=Omega_p partial_theta.
On the preparation support, Z=lambda W and W'=-lambda V. For a smooth
compact scalar S independent of zeta, take

    xi=(J grad_perp S,lambda S),
    xi cross u=grad(WS)-(T_p S)e_zeta,
    w=P(xi cross omega)=-lambda(T_p S)e_zeta.             (1)

These are literal vector identities. The final projection does not
discard pressure: the displayed remainder is already solenoidal, and
the preceding gradient is compact. With S_t=-T_p S, the full Euler and
Lin equations hold exactly. The meridional displacement and the density
perturbation of a radial stationary tag are those of 0217. The Eulerian
vorticity variation curl w is nonzero.

The full canonical cotangent is pi=rho(w+Du xi), not the homogeneous
formula rho Du xi. Its axial component vanishes because
w_zeta=-lambda T_p S and (Du xi)_zeta=lambda T_p S.
For the real m=-1 cosine/sine pair, the complete phase is consequently
the same positive beta as 0217:

    beta=rho pi L m integral Z' A_0^2 ds>0,
    H_cos=H_sin=rho pi L lambda^2 m^2/2
                         integral s Omega_p^2 A_0^2 ds.  (2)

The energy is the complete material Jacobi energy, equivalently here
rho||w||^2/2. The reference velocity has zero self-helicity. Direct
material calculation gives the same result: adding xi_zeta=lambda S
subtracts its equal transport-rate square from the old homogeneous
energy. Thus the conserved energy MATRIX tends to beta Omega_0 I,
not 2beta Omega_0 I. It naturally matches the leading physical clock
energy. The vorticity-cross part of phase alone would have the wrong
answer; the xi.w cross terms restore the displayed full cotangent.

## 2. Compact radial return and exact divergence completion

Keep the fixed positive observed tag of 0217 inside s<a. Choose the
radial m=-1 amplitude A_0(s) there as before, and add a smooth signed
return on 2a<s<3a such that

    integral_0^infinity s^2 A_0(s) ds=0.                 (3)

This is one explicit nonzero moment row: a fixed positive bump on that
annulus has a strictly positive s^2 integral, so its coefficient is
minus the inner moment divided by that integral. The return does not
alter the tag's inner moment B. All profiles remain smooth, and 3a is
chosen inside the literal-curl region. No negative tag density is used.

The regular radial inverse of the m=1 Laplacian gives a smooth compact
f=f_0(s) exp(-i theta) with Delta_perp f=S. Indeed its only decaying
exterior term is proportional to s^-1 integral t^2 A_0(t)dt, and (3)
sets it to zero. Interior regularity and flat source matching make f
smooth across both the axis and outer boundary. This is a compact
potential construction, not an imposed wall or an inverse spectral gap.

Set zeta=R varphi, k=-1/R, and include the actual factor exp(ik zeta).
In straight tube coordinates define

    Xi_k=(J grad_perp S-i k lambda grad_perp f,lambda S).

Its divergence is exactly -ik lambda Delta_perp f+ik lambda S=0.
The map

    F_R(X,Y,zeta)=((R+X)cos(zeta/R),
                  (R+X)sin(zeta/R),-Y),
    J_R=1+X/R

has positive Jacobian on the support. The Piola field
xi_R(F_R y)=DF_R Xi_k/J_R is exactly divergence-free and compact in
the actual core. It is smooth and single-valued because kR=-1 is an
integer harmonic. The extra meridional component is O_a(1/R); no
nonlocal projection of xi_R introduces support outside the Beltrami
region. Use real and imaginary parts as the actual two preparations.

The actual initial velocity is now, globally,

    w_R(0)=P_R3(xi_R cross omega_R).                     (4)

On the support omega_R=lambda u_R, so the literal-curl identity applies
there without being assumed in the exterior. This is an actual Kelvin
initial condition. Evolve (4) with full linearized Euler, and xi_R with
Lin reconstruction. Then delta omega=curl(xi_R cross omega_R) holds
for the entire finite interval, and the material tag is advected by
the same displacement. In particular its vortex tangency is varied
consistently, unlike the w=0 route.

The exact toroidal cross identity explains the geometric remainder.
For the uncompleted field eta=(-S_z/r,-lambda S/r,S_r/r) and
F(psi)=lambda psi+constant,

    eta cross u=-grad(FS/r^2)+(TS/r)e_varphi
                                      -(2FS/r^3)e_r,
    div eta=-lambda S_varphi/r^2.                       (5)

The last two curvature/divergence facts are nonzero. The Piola
construction repairs the divergence exactly, and the pressure estimate
below controls, rather than deletes, the remaining terms.

## 3. Uniform actual-field and full pressure estimates

Here and below a, lambda, the fixed outer taper, the profile/return,
the finite time T and all required derivative orders are selected FIRST.
Constants may depend on them, but not on the subsequently large R.

The global fields u_R have uniform C^j bounds for each fixed j. To see
this from the actual 0211 geometry, omega_R has uniform local C^j bounds
on a tube of fixed cross radius about a circle of radius R. The volume
of that tube in any ball of radius d above the cross scale is at most
C d until d~R. The far part of the Biot-Savart kernel and its derivatives
is therefore uniformly integrable. On the near part subtract the local
Taylor value in the singular derivative kernels and use the uniform
local smooth norms. Higher derivatives use the same subtraction with
derivatives of omega. The actual uniform far velocity U_R=O(log R/R)
is bounded as well. This proves uniform global C^j bounds, not just
local convergence near the core. The full Euler H^s growth constants
on a fixed interval are consequently uniform in R.

Let S_t be the exact reference transported scalar from (1), with the
single ring factor exp(-i varphi), and let W_R^app be its compact
tangential passive velocity on the curved tube. It is not declared
solenoidal. Its divergence and every fixed Cartesian derivative have
size O_{a,T}(1/R), supported on a tube of length O(R) and fixed width.
For such a source d, the actual Newton energy obeys

    ||d||_(dot H^-1)^2
      <=C integral integral |d(x)d(y)|/|x-y| dxdy
      <=C_{a,T}(1+log R)/R.                             (6)

For each x the integral of |x-y|^-1 over the fixed-width circular tube
is O(1+log R); multiplying its volume O(R) by the two O(1/R) amplitudes
proves the final bound. The same estimate applies to each fixed
derivative of d. Thus the FULL whole-space Leray correction satisfies

    ||(I-P_R3)W_R^app||_(H^s)
                      <=C_{a,T} sqrt((1+log R)/R).       (7)

It includes the ambient tail and the lowest global azimuthal harmonic.
It uses no high-n multipole suppression, no Poincare constant, and no
radial wall. Relative to the norm O(sqrt R) of the prepared column,
the pressure error is O(sqrt(log R)/R).

Put V_R=P_R3 W_R^app. The compact local part of its Euler residual is
O_{a,T}(log R/R) per unit length: it consists of the actual 0211
coefficient difference, the tube metric, and the n/R derivative terms.
The remaining part is exactly

    P[u_R.grad(P-I)W_R^app+((P-I)W_R^app).grad u_R].

By (7) and the uniform global bounds this has H^s norm at most
C_{a,T} sqrt R log R/R. The initial velocity (4) differs from V_R(0)
by the same order: expand the compact identity (1) through the Piola
map, retain its exact gradient and apply P. All metric, return-completion
and actual coefficient errors are compact O_a(log R/R) terms before
projection. The projection has norm one on every H^s.

Full Euler Duhamel, followed by the actual Lin transport equation, now
gives for every selected finite s and number r of time derivatives

    ||w_R-V_R||_(C^r_t H^s)
       +||xi_R-Xi_R^app||_(C^r_t H^s)
                  <=C_{a,T,s,r} sqrt R log R/R.          (8)

Xi_R^app is the mapped exact reference displacement (1), with its
toroidal lambda S_t component; its compact transport residual is
estimated in the same way. A compact inverse for f_t at later times
is NOT assumed: f was needed only for the exact initial completion.
All later reconstruction is the actual Euler/Lin evolution. This avoids
requiring the moment constraint (3) to survive radial dephasing.

## 4. Same actual tag, currents and action normalization

Use the exactly stationary chi_R(phi_R) from 0217, with its positive
global covariance difference D_tag. Its first variations are computed
from the actual xi_R,w_R. In particular

    delta X=M_tag^-1 integral rho chi_R xi_R,
    G=integral rho chi_R x cross xi_R,
    S=integral rho chi_R[x cross xi_R,t+2xi_R cross u_R],
    delta I_ij=integral rho chi_R(x_i xi_j+xi_i x_j).       (9)

These integrated material identities retain the moving centroid and
the ambient exchange; they are equivalent to the full moving-domain
formula, not fixed-domain angular impulse substitutions. The reference
centroid and momentum vanish. For estimating the centroid use the FIRST
identity, rather than bounding x times the small density divergence;
the latter loses an unnecessary power of R.

The toroidal displacement/velocity in (1) has no leading effect on the
tag density or global tilt. Its extra transverse mechanical-spin and
G moments have relative size O_a(1/R). The leading reference rows are

    theta_+=Q/(iD_tag),
    G_+=D_tag theta_+ [1+O((lambda a)^2)+o_R(1)],
    S_+=D_tag theta_+,t [1+O_{T}((lambda a)^2)+o_R(1)].    (10)

For example, with the 0217 real-phase registration their extra direct
column rows are G_toroidal=-i rho pi^2 lambda R integral chi_a A_0 s^2 ds
and S_velocity=rho pi^2 lambda R integral chi_a Omega_p A_0 s^2 ds.
The leading rows have size R^2 B; these explicit integrals display the
suppression by one power of R. The compact meridional completion has
the same or smaller fixed-a relative order, and is included in (8).

All statements about errors in (10) mean bounded linear observation
rows on the two initial amplitudes, not division by a sinusoid at its
zeros. The exact remainder G_t-S=-2 integral rho chi_R xi_R cross u_R
is retained. The polar centroid and symmetric shape remain separate
rows; no complete-fluid current is inferred from the tag alone.

Estimate (8), Cauchy-Schwarz on the actual tube, and (9) control all
these normalized rows at their own fixed-a scales. The Q test has size
O(R a), not R^2; its signal is O(R^2 B), and its denominator D_tag is
O(R^3 I_c). The resulting relative errors tend to zero as R grows.
The complete phase and energy errors are O_{a,T}(R log R/R). Compare
them with beta~R lambda^2 Omega_0 I_A and H~beta Omega_0, after a and
the actual return are fixed. This is an explicit nonzero own-scale
ordering; no absolute order-one estimate decides the small phase sign.

The signed outer return increases I_A but leaves B and I_c unchanged.
It therefore preserves the strictly positive beta and the scaling

    M_phase/j_tag=4lambda^2 I_c I_A/B^2 [1+o(1)]
                 =O((lambda a)^2).                     (11)

Choose a small enough and then R large enough that the exact positive
ratio is below one. Its value alpha scales the fixed stationary material
fraction as in 0217, leaving the physical angle unchanged. To fix exact
initial real phase conventions, use the ring's reversing half-turn
(x,y,z)->(x,-y,-z): the background transforms to its time reversal,
and chi is invariant. Project the first compact initial generator onto
its prescribed parity and obtain the second by a quarter azimuthal
rotation. This preserves compactness, solenoidality, Kelvin preparation
and the leading pair, and fixes Q(0) real, Q_t(0) imaginary and S_+(0)
real. Hence the initial scalar spin/rate ratio is a real positive number.

The exact alpha=M_phase(0)/j_tag(0) then supplies exact initial spin/
phase-inertia matching. Throughout any fixed chosen interval the actual
spin, G, positive scalar mass/stiffness and conserved energy are as close
to the same positive physical oscillator normalization as prescribed,
by first reducing lambda a and then choosing finite R using (8).
Time-dependent phase, frequency, shape and current connections remain
in the exact finite-preparation action. Unlike the homogeneous-label
route, its conserved-energy limit already has the correct factor one.

## 5. Earned scope and remaining parent interface

This constructs a genuine field-changing Kelvin/Euler response and a
fixed positive material-tag optical clock on the SAME globally smooth
Euler ring with a literal nonzero closed Beltrami core. The vorticity
and the tag are transported by the same displacement. Its complete
phase, energy and physical normalization are positive with a controlled
fixed-time error; the initial spin normalization is exact. The full
ambient pressure and the actual lowest azimuthal harmonic are included.

The limits are ordered, and every selected preparation has finite
radius, finite action and a nonzero tag fraction. This does not assert
a uniform positive action density when accuracy, radius and cell packing
are taken to their limits together. It does not give an isolated global
spectral pole or all-time autonomous oscillator, nor a common-K second
jet merely from the integer n=-1. Actual spatial derivatives, physical
hybrid/shape currents and the acoustic joining remain the next parent
constructions. No label-only or foreign-geometry phase control was used.

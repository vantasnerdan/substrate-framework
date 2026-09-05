# Full homogeneous material action and the field-changing continuation

This closes the full initial phase/energy evaluation for the fixed-tag
response, while keeping its distinction from a vortex-structure mode.
Density rho is constant; all integrations include the complete fluid.
Compact normal support and periodic X remove boundary terms.

## Exact full three-component initial forms

For any steady Euler field, write B=Du, T=u.grad, omega=curl u.
Homogeneous Lin data obey xi_t=Bxi-Txi and pi=rho Bxi. The full phase is

    Omega(xi1,xi2)=-rho integral omega.(xi1 cross xi2).     (1)

Differentiate the steady Euler equation to obtain
Hess(p)=-B^2-TB. Integrate its TB term by parts in the exact Jacobi
Hamiltonian. The resulting identity is

    H=rho/2 integral omega.(xi cross xi_t).               (2)

Equation(2) includes all components of xi_t; it is not the kinetic energy
of w, which is zero for this family. It is the full second-variation
Hamiltonian in the stationary Euler/material chart. Its invariance under
time evolution does not turn the selected two initial stream columns
into an invariant two-dimensional subspace.

For the actual C016 field curl u=-u and xi(0)=(0,J grad s), formulas
(1)--(2), with the ACTUAL xi_X,t=-Ts, reduce to

    Omega_12=-rho integral s1 T s2,
    H= rho integral (Ts)^2.                             (3)

For example omega.(xi cross xi_t)=psi {s,Ts}+(Ts)^2,
and integral psi {s,Ts}=integral(Ts)^2. Removing the axial coordinate
velocity removes half the energy. The exposing verifier computes the
full pressure Hessian, Du and vector rates first and checks (3) on three
independent actual periodic stream profiles. It also checks a nonzero
full cotangent phase pairing. The proof uses integration by parts for
arbitrary admissible streams; the examples are exact sign/omission anchors.

For s1=f(E)cos(l theta_o), s2=f(E)sin(l theta_o),

    Omega_12=-pi rho l integral f^2 dE,
    H1=H2=pi rho l^2 integral omega(E)f(E)^2 dE,
    H12=0.                                             (4)

Here Hi is the scalar energy of one unit column, not the diagonal entry
of its Hessian, which is twice Hi. The observed angle columns have the
limiting form (C cos(nu t),-C sin(nu t)); their angle/rate determinant is
-nu C^2. The inherited limiting phase mass is therefore

    M=pi rho l integral f^2/[nu C^2]>0.                 (5)

It generally differs from the finite literal spin coefficient j_*.
For bands narrowed with fixed C, M grows inversely with width. The full
energy in (4) tends to M nu^2 C^2, twice the mechanical oscillator's
one-column energy M nu^2 C^2/2. This is a derived actual excess, not a
claim that a positive observed action already equals full-fluid energy.
Reviewed0210 phase returns and0205 energy returns are possible finite-
norm repairs at their scopes. This attempt does not silently apply them
without the common acoustic cross-form and spatial-current construction.

## Why the next route changes the actual Euler field

The exact planar vorticity equation and material density equations give

    z=Ts, r=Hphi-z, r_t+Tr=0,
    Hphi_t=-T(H-1)phi.                                 (6)

The present positive clock uses phi=0 and nonzero r. It changes the
material tag while leaving the Eulerian vortex field unchanged. On the
fixed Kelvin leaf, by contrast, delta omega_X=-xi.grad omega_X=-Ts,
while delta omega_X=-Hphi, so Hphi=Ts and r=0. A nontrivial homogeneous
tag displacement with phi=0 does not satisfy that vortex-transport
condition. This is a precise route distinction, not a refutation of the
positive material response or a claim that every Euler initial variation
must stay on one Kelvin leaf.

The field-changing continuation now executes (6) with r=0, or a separately
derived nonzero Euler variation whose moved vortex geometry matches its
actual material tags. The canonical cell is particularly useful because
H=-Delta on the mean-zero normal torus has H>=1. Its actual field-changing
operator can be written using S=(1-H^-1)^(1/2), with the finite first-shell
nullspace treated explicitly. This supplies a positive energy-space
representation for the next actual spectral/output calculation; it is not
by itself an observed-mode existence or cyclicity theorem.

Thus the attempt establishes a positive fixed-density physical material
clock with its complete initial forms. The user's field-changing vortex
mechanism and common spatial action continue on their original scope.

# Same-orbit common and relative angle action

## Exact finite-domain stationary field

Use physical cylindrical coordinates (r,theta,z), with z of period L.
Let lambda>0 and choose R with J_2(lambda R)=0. Such positive simple roots
exist; the Bessel equation and positive-zero facts are recorded in NIST
DLMF 10.2.1 and 10.21. Define

    psi=A J_0(lambda r)+B J_2(lambda r) cos(2 theta),
    u=(psi_theta/r)e_r-psi_r e_theta+lambda psi e_z,
    omega=lambda u, p=pstar-rho |u|²/2.

The Bessel equation gives Delta_perp psi=-lambda² psi, hence curl u=lambda u
and stationary incompressible Euler exactly. The apparently polar functions
are analytic at the axis. On r=R, psi is constant, so u_r=0. Thus the finite
cylinder is an invariant material fluid domain, with no sidewall mass flux.
The axial direction is periodic. The domain itself has an exact SO(2) action.

For A>0 and 0<|B|<2A, the quadratic core streamfunction is

    psi=A+lambda²[(-A/4+B/8)x²+(-A/4-B/8)y²]+O(r⁴).

It has an elliptic, noncircular core. Nearby psi-level ovals crossed with
the axial circle are actual invariant smooth vortex tubes. Rotation of this
quadratic core is a physical orientation, not an axisymmetric relabeling.
No singular point vortex or imposed core-polarization modulus is involved.

The cylinder is a declared Euler material-domain candidate, not yet a
boundary/gluing theorem for an arbitrary prescribed EPS knot. The ambient
pressure and momentum transfer remain explicit when it is used as an RVE.

## The Euler orbit and boundary terms

Generators are divergence-free and tangent to the cylindrical boundary.
Keep the boundary circulation periods fixed. For these generators and
omega.n=0, delta u=P(xi cross omega) has tangential boundary trace equal to
a surface gradient. The background tangential velocity one-form is closed
because (curl u).n=0. Consequently the boundary pairing of u with delta u
in curl integration by parts vanishes: the integral of a closed one-form
wedged with an exact one-form on the boundary torus is zero.

The Leray projection includes the harmonic component determined by this
coadjoint leaf. It is not replaced by an arbitrary zero-mean inversion.
The exact orbit Hessian and symplectic form therefore are

    H(xi,eta)=rho integral_D[delta u_xi.delta u_eta
                            -delta u_xi.curl(delta u_eta)/lambda],
    Omega(xi,eta)=rho integral_D omega.(xi cross eta).

Both come from the same complete Euler action. The identity E'=0 on this
orbit follows pointwise from u parallel omega. Boundary periods are part
of the leaf, not a source of a fitted angular stiffness.

## Explicit common-angle momentum partner

Let K=e_z cross x, the global rotation generator, and

    v_K=P(K cross omega)=e_z cross u-(K.grad)u.

This is the physical tangent of the rotated stationary field. It is nonzero
because B is nonzero. Rotation invariance of energy on the same cylinder
and E'=0 imply H(K,eta)=0 for every admissible eta, including K itself.

Choose a nonnegative smooth cross-sectional cutoff chi supported where
curl v_K is nonzero and away from the core-jet neighborhood, and put

    eta_0=curl(chi curl v_K).

It is smooth, divergence-free, independent of z and supported away from
the sidewall. The common angular-momentum derivative is explicitly

    B_body=Omega(K,eta_0)
          =-rho integral_D v_K.eta_0
          =-rho integral_D chi |curl v_K|² != 0.

Thus the common angle has a genuine conjugate shape on the Euler orbit;
it has not acquired a supplied rigid-body inertia. The sign fixes momentum
orientation only. Its nonzero magnitude is an actual fluid integral.
The analytic nonaxisymmetric tangent is nonzero on open sets away from the
axis, so the stated support choice exists. Both momentum partners can then
be supported away from the core, preserving the physical section-angle map
when their canonical momenta are eliminated.

The energy H(eta_0,eta_0) need not be positive. Add a compact negative-helicity
cage generator eta_2 at a large nonzero axial harmonic k_2. Its transverse
support lies in a region with omega_z of one sign. The full Helmholtz/Leray
calculation of 0045 applies with the axial Fourier Neumann resolvent on this
cylinder: its L² bound is at most k_2^{-2}, and compact source support removes
sidewall integration terms. It gives H(eta_2,eta_2)>0 growing linearly with
|k_2| for a fixed envelope. All cross terms with eta_0 or K vanish exactly
by axial Fourier orthogonality. Set eta=eta_0+eta_2. For a sufficiently large
finite integer harmonic,

    h_body=H(eta,eta)>0,
    Omega(K,eta)=B_body.

The coefficients are the exact finite-harmonic integrals above. The large-k
bound proves the existence of a usable finite choice; it does not replace
its coefficient by a fitted asymptotic constant.

## A separate internal core/cage pair

Choose a distinct nonzero axial harmonic k_1. The internal pair is a compact
negative-helicity cage pair at k_1, with the first generator also carrying
a disjoint compact transverse core rotation chi_core(e_z cross x)cos(k_1 z).
The latter is divergence-free for a radial chi_core and rotates the actual
noncircular quadratic core at z=0, with the opposite rotation at a half-period
section. Its leading self-helicity vanishes; its full energy is bounded at
large k_1 as in the invariant-tube construction 0042. The cage gives the
positive leading energy and nondegenerate symplectic form. Keep the complete
positive matrix H_internal and nonzero B_internal, not merely one diagonal.

Every common/internal H and Omega cross entry is exactly zero because the
background is z-independent and k_1 differs from zero and |k_2|. The Leray
operator commutes with axial translations, so this includes its nonlocal
kinetic term, rather than just local energy products.

The resulting four-dimensional constrained Euler action is

    L=B_body s_body Bdot - h_body s_body²/2
      +B_internal s qdot
      -(H_qq q²+2H_qs q s+H_ss s²)/2.

Eliminating both canonical momenta gives, up to a total time derivative,

    L=I_body Bdot²/2+I_q qdot²/2-K_q q²/2,
    I_body=B_body²/h_body>0,
    I_q=B_internal²/H_ss>0,
    K_q=H_qq-H_qs²/H_ss>0.

This is the sought intrinsic common/relative construction on one Euler
orbit: the two inertias and locking coefficient are all fluid integrals.
It does not mix an isovortical action with a separately added rigid rotation
velocity, which would change Kelvin data.

## Physical affine-cage field map

The two section orientations are theta_plus=B+q and theta_minus=B-q.
Put beta=theta_minus, a=1+I_q/I_body, and Psi=beta+a q. Then Psi transforms
as a physical angle under a common frame rotation. Direct substitution yields

    T=J_Psi Psidot²/2+J_beta betadot²/2,
    V=K_Psi(Psi-beta)²/2,
    J_Psi=I_body²/(I_body+I_q),
    J_beta=I_body I_q/(I_body+I_q), K_Psi=K_q/a².

The field agrees with a specified linear combination of the actual section
angles; it is not an angular rate renamed as an angle. Under a declared
affine material cage beta=n.curl(U)/2, it has the standard relative-angle
spring and absolute microrotation inertia plus the explicitly retained
macroscopic gradient-inertia term. Spatial momentum gradients must also be
retained before identifying optical k² coefficients, as established in 0041.

## Current verification boundary

The positive common-rotor construction and its coefficient map are exact
conditional statements once the stated finite-domain high-harmonic bound
and the internal pair are carried over with their boundary conditions.
The Bessel PDE, boundary circulation argument, all symmetry/momentum
identities and four-coordinate Routh algebra are checked independently next.
Neither EPS domain transfer nor the full continuum gradient action is
inferred solely from this free-rotor construction.

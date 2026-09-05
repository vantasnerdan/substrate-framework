# Shaped stationary force-free cores and an exact local closed tube

## 1. A freely shaped Euler core with derived axial flow

Let zeta(r)>=0 be a smooth radial vorticity profile, smooth as a function
of r² at the axis, with zeta(0)=2 Omega>0 and compact support r<=c.
It may equal 2 Omega on an inner disk and taper smoothly to zero. Define

    V(r)=r^-1 integral_0^r s zeta(s) ds,
    D(r)=2 integral_0^r zeta(s)V(s) ds,
    W(r)=sqrt(U0²-D(r)), U0²>D(infinity),
    u=V(r)e_theta+W(r)e_z, f(r)=zeta(r)/W(r).          (1)

Every quantity is smooth at r=0; V=Omega r+O(r³), D=O(r²), W>0.
Direct differentiation gives

    V'+V/r=zeta, W'=-zeta V/W,
    curl u=(0,-W',zeta)=f u, div u=0.                 (2)

This is an exact smooth stationary constant-density Euler field with
p=-rho|u|²/2+constant. In particular p'=rho V²/r. Its pressure is the
same radial pressure as the ordinary pure-swirl column: the added axial
flow creates no extra radial force. Outside c, V=Gamma/r, W=W_infinity,
and vorticity vanishes identically. There is no sheet or wall force.

An integration by parts gives the alternative formula

    D(r)=V(r)²+2 integral_0^r V(s)²/s ds.              (3)

Thus the axial flow is derived from the chosen swirl, not fitted to an
optical frequency. It is a generalized Beltrami field; f is not silently
substituted for a constant eigenvalue in any accepted constant-curl API.
The infinite column is not finite total energy or a closed knotted tube.

A globally smooth explicit noncompact-vorticity example is

    V=Omega r/(1+r²/a²),
    zeta=2 Omega/(1+r²/a²)²,
    W²=U0²-Omega²a²+Omega²a²/(1+r²/a²)²,
    U0>Omega a.                                     (4)

This supplies an exact shaped-core comparison; compact vorticity is
provided by (1), not falsely attributed to (4).

## 2. Actual finite-time Euler response transfers from the pure column

Apply the genuine Galilean axial boost -U0. The stationary velocity is
V e_theta+h e_z, h=W-U0. At fixed profile,

    ||h||_Cm <= C_m/U0,                              (5)

for all finite m and sufficiently large U0, including the exterior
constant h. This follows directly from the square root and bounded
derivatives of D. For a fixed axial Fourier wave number q, the velocity
linearized Euler generators on divergence-free L² cross-sectional
fields differ by the bounded operator

    E_q v=-P_q[i q h v+h' v_r e_z],
    ||E_q|| <= |q| ||h||_infinity+||h'||_infinity.     (6)

P_q is the COMPLETE pressure projector, not a local pressure model.
The coefficients are independent of z, so the comparison uses the same
Fourier sector and the same projector. Smooth bounded background
gradients give the usual L² Euler energy/semigroup bound. Duhamel then
gives a C_T/U0 difference between the two actual linearized solutions
on each fixed finite interval. Higher Sobolev versions follow by the
same differentiated energy estimate, with constants for the fixed
profile and q. Lin transport and transported material tags are compared
in those same Sobolev norms; genuine angle/spin moments therefore also
transfer when their defining reference chart is nonsingular.

For a simple ISOLATED column eigenvalue with a supplied resolvent contour,
the bounded perturbation (6) also preserves its Riesz projection when
||E_q|| max_contour||(z-L_column)^-1||<1. This is a conditional
perturbation implication, not an assertion that a desired mode is isolated.
It does not move an embedded critical-layer mode out of the continuum
by fiat, or transfer a stability theorem for zero axial shear to all U0.

## 3. Exact local toroidal continuation of a compact-vorticity profile

Choose a disk radius b>c. Let phi0(s)=integral_s^b V(t)dt, and define
F(phi0(s))=zeta(s). This defines a smooth F on the attained interval:
near the center phi0(0)-phi0(s)=Omega s²/2 if the core is flat, and
near the boundary F=0 identically. Extend F smoothly outside that
interval. Define C(phi)>0 by

    C(phi)²=C(0)²+2 integral_0^phi F(t)dt.             (7)

Then C(phi0(s))=W(s) if C(0)=W_infinity, and
Delta phi0+F(phi0)=0 on the disk. The linearized Dirichlet operator is
Delta+F'(phi0). Its potential Q(r)=-zeta'(r)/V(r) is compactly supported
and smooth, including the flat central region.

The freedom to choose b makes this operator invertible. For angular
number m>=1, V>0 solves

    V''+V'/r+(Q-1/r²)V=0.

The ground-state identity for radial Dirichlet h is

    integral_0^b [r h'²+(1/r-rQ)h²]dr
      = integral_0^b r V²[(h/V)']² dr.               (8)

There is no boundary term since h(b)=0, V(b)>0 and regular h vanishes
at the axis in these sectors. Adding (m²-1) integral h²/r proves that
no m>=1 sector has a Dirichlet zero mode. For m=0, the regular homogeneous
radial solution is unique up to scale; outside c it is A+B log r, with
A,B not both zero by ODE uniqueness. It has at most one zero radius
b>c. Choose any other finite b. The self-adjoint elliptic Fredholm
operator is then invertible on the usual Dirichlet Schauder spaces.

Set x=r_toroidal-R, kappa=1/R. The implicit-function theorem applied to

    [Delta-kappa/(1+kappa x) partial_x]phi+F(phi)=0,
    phi|disk_boundary=0,                             (9)

gives a smooth exact solution for sufficiently small finite kappa,
close in every fixed C^m norm to phi0. Define the actual 3D velocity

    u_r=-phi_z/(1+kappa x),
    u_theta=C(phi)/(1+kappa x),
    u_z=phi_x/(1+kappa x).                           (10)

Full cylindrical differentiation yields div u=0 and
curl u=C'(phi)u; hence this is an exact stationary Euler field in the
closed solid torus, with its full pressure and no filament approximation.
The level phi=0 is an invariant boundary. The central nondegenerate
maximum persists, giving a circular elliptic core. The radial reference
has no other critical point; compactness away from the center preserves
this nested level-set structure for sufficiently small kappa.

The straight-limit reduced section rotation is V(s)/(s W(s)). In the
outer annulus this equals Gamma/(s² W_infinity), with strictly nonzero
radial derivative. The true flux action has I'=s W_infinity>0 there,
so its twist is nonzero. As in0120, actual contour transit integrals
and the flux measure give continuous twist and allow a finite major
radius with Diophantine boundary return and nonresonant elliptic core.
This is a local stationary smooth torus construction, not a global
extension theorem or an arbitrary-knot claim.

## 4. Exact source licenses and the remaining global construction

Gallay–Smets, [Spectral stability of inviscid columnar vortices](https://www-fourier.univ-grenoble-alpes.fr/~gallay/SpecStab.pdf),
Theorem1.3 treats pure swirl under strict profile assumptions H1/H2,
in fixed nonzero axial Fourier sectors. Remark1.4 retains essential
spectrum and Kelvin eigenvalues. These statements do not assert that
our compact-core taper or axial-flow extension satisfies that theorem.

Enciso–Poyato–Soler, [arXiv1605.06626v2](https://arxiv.org/abs/1605.06626), Theorem4.2 supplies generalized
Beltrami fields with robust knotted tubes on an EXTERIOR domain. Its
factor perturbation is supported in a boundary-to-boundary stream tube,
not arbitrarily prescribed inside the invariant knots. Theorem5.3 is
a small-ball local perturbation result, not a global Runge extension.
PDF SHA256: `929342806b7ac20ba209ac6d3f8f5adbbad910694c9bd45a1b21514f13813399`.

Consequently neither theorem alone embeds (10) in the accepted global
stationary ensemble. A global extension, or a direct global stationary
construction with the same physical core response, remains an active
route. The constant-factor EPS Runge theorem cannot approximate an
arbitrary nonconstant-factor field to arbitrarily small C¹ error.
The positive results here are (1), the actual finite-time comparison
(6), and the exact local closed tube (10).0135 investigates the actual
optical action/observation that this construction can preserve.

## 5. A global Bernoulli lift, without a single-valued vorticity law

There is a different exact global construction. Let v=(v_x,v_y) be a
smooth stationary planar incompressible Euler field, with actual global
physical pressure p and constant rho>0. Suppose B=p/rho+|v|²/2 is bounded
above on its domain. For any constant C>sup B define

    W=sqrt(2(C-B)), u=(v_x,v_y,W),
    zeta=partial_x v_y-partial_y v_x.                (11)

The planar Euler identity gives

    grad B=(zeta v_y,-zeta v_x), v.dot(grad B)=0.

Consequently W_x=-zeta v_y/W and W_y=zeta v_x/W, and direct Cartesian
curl gives

    curl u=(W_y,-W_x,zeta)=(zeta/W)u.                (12)

The divergence remains zero, the full three-dimensional Euler pressure
is the original p, and p/rho+|u|²/2=C exactly. These conclusions hold
globally wherever the given planar fields and the strict inequality
hold. Periodic planar fields with global periodic pressure have bounded
B and yield a genuine periodic three-dimensional field. A stationary
street with bounded B yields a global nondecaying columnar street.

Crucially, no globally single-valued relation zeta=F(psi) is required.
Different streamline components may carry different vorticity laws;
the actual globally defined pressure supplies B and its smooth joining.
On a connected irrotational exterior grad B=0, so W is constant there.
This repairs0134's global-streamfunction-law obstruction without adding
spurious vorticity in the far field. It does not impose a pressure where
none is known to exist. A pressure gauge p->p+rho g is accompanied by
C->C+g and changes no velocity or curl factor.

## 6. The planar dynamics, Kelvin form and action are preserved exactly

The z-independent subspace of three-dimensional Euler is invariant:

    v_t+(v.dot(grad_xy))v=-grad_xy p/rho, div_xy v=0,
    w_t+v.dot(grad_xy)w=0.                            (13)

Thus the full planar Euler evolution is unchanged, while the axial
velocity is a genuinely advected scalar. Equation(11) specifies the
stationary background; a perturbed trajectory is NOT required to stay
force-free. Its exact evolution follows(13).

Restrict to area-preserving planar Kelvin maps with z fixed. For planar
generators xi,eta the full three-dimensional KKS form reduces to

    Omega_3(xi,eta)=rho integral zeta (xi_x eta_y-xi_y eta_x) dx dy dz,

since xi cross eta is axial. The planar induced velocity is the usual
two-dimensional Euler coadjoint variation; its axial component is
-xi.dot(grad w). There is no omitted axial vorticity contribution to
this restricted two-form. On this orbit w=w0 composed with the inverse
area-preserving map, hence integral w² is constant. The kinetic
Hamiltonian restricted to the orbit differs from the planar Euler
Hamiltonian by that constant only. The local phase action therefore has
the same two-form and Hamiltonian variations, up to a boundary term in
its potential. No additional axial inertial term is attached to a planar
mode. The invariant subsystem(13) independently verifies this action
restriction rather than inferring dynamics from the restriction alone.

For periodic cells all integrals are per axial length/period. In an
unbounded street, use compactly supported planar variations and the
background-subtracted axial energy, with vanishing boundary flux; the
claim is not finiteness of its total kinetic energy. Physical in-plane
velocity, coarse in-plane momentum and axial material angular momentum
are exactly the corresponding planar quantities. Axial displacement
and nonplanar generators are outside this action restriction and are
retained separately when constructing optical modes.

The reusable `planar_bernoulli_lift` API checks stationary planar Euler
and exposes the remaining strict domain inequality. Cartesian residual,
pressure-gauge, wrong-sign/domain and passive-energy tests exercise it.
This supplies a common force-free background for actual planar acoustic
and three-dimensional core-mode work. It does not itself supply a
globally knotted EPS ensemble, an isotropic same-field array, or a closed
autonomous Cosserat field equation. Those are the next constructions,
not consequences of the word force-free.

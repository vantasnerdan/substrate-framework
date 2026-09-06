# Exact nonlinear material balances and the ambient-pressure closure test

This is exact calculus for classical Euler on its smooth existence interval, not a small-angle approximation or an existence/stability theorem for a particle. Physical pressure is p, constant density is rho>0, and D_t u=-grad p/rho. A material tag is a transported bounded volume, not a vortex boundary imposed as a wall. All volume and surface integrals below are at the same time.

## 1. Full material state and exact balance

Write dm=rho dx, M=integral dm, X=M^-1 integral x dm, V=M^-1 integral u dm, r=x-X, c=u-V. Use unnormalized tensors I=integral rr^T dm, B=integral cr^T dm, T=integral cc^T dm and spin L=integral r cross c dm. With F=-integral grad p dx and P=integral grad p r^T dx, transport and incompressibility give

    Mdot=0, Xdot=V, M Vdot=F,
    Idot=B+B^T,
    Bdot=T-P,
    Iddot=2T-P-P^T,                                      (1)
    Ldot=-integral r cross grad p dx
        =-integral_boundary p r cross n dA,               (2)
    Tdot=-integral (grad p c^T+c grad p^T) dx,
    E_int=tr(T)/2, E_int_dot=-integral c dot grad p dx
                           =-integral_boundary p c dot n dA. (3)

For the total tag energy E_tag=M|V|^2/2+E_int,

    E_tag_dot=-integral_boundary p u dot n dA.             (4)

The centrifugal or restoring force is not supplied by assigning I or L: actual pressure work and the full tensor T enter. The pressure moment identity is

    P_ij=integral_boundary p n_i r_j dA-delta_ij integral p dx. (5)

It makes (1) insensitive to an additive constant pressure and gives the full scalar virial I_trace_ddot=2 tr(T)-2 integral_boundary p r dot n+6 integral p. Equation(2) has the pressure torque sign appropriate to outward n. Shape eigenvalues may evolve; Idot is not assumed to be a rigid commutator. Full-space conservation requires the corresponding decay/boundary conditions and is distinct from tag conservation.

These identities follow by D_t r=c, D_t c=-grad p/rho-F/M, integral r dm=integral c dm=0. Differentiating the centered quantities before dropping these mean terms is essential. Equation(3) also specifies precisely the next unresolved pressure/velocity correlation. Existing euler_observation.material_tag_moments supplies the simultaneous I/B/L algebra and is reused by the additive balance API; it does not already evolve these moments through physical pressure.

## 2. Smooth finite-energy ambient construction

Let f(s)=g(s^2) be a nonzero smooth radial bump with support s<a, flat at its outer edge. On R^3 choose a center D and a unit axis n and set

    y=x-D, s=|y|, u_0(x)=f(s) n cross y,
    J=integral_0^a s^4 f(s)^2 ds>0.                       (6)

This is C-infinity compactly supported, divergence free, and lies in every H^s. It is admissible Euler initial data; it is not asserted steady. Its unique decaying physical pressure at t=0 solves

    -Delta p=rho tr[(Du_0)^2].                            (7)

For n=e_z, direct differentiation gives tr[(Du_0)^2]=-2f^2-2ff' s sin^2(theta). Its only angular harmonics are

    S_0=-2f^2-(4/3)s f f',
    S_2=(4/3)s f f' multiplying P_2(cos theta).           (8)

The decaying radial Poisson Green formula outside support is p_l(s)=rho s^(-l-1)/(2l+1) integral_0^a t^(l+2) S_l(t)dt. Integration by parts gives integral t^2 S_0=0 and integral t^4 S_2=-(10/3)J. Consequently the exterior field is EXACTLY

    p(x)=-(2rho J/3) P_2(n dot y/|y|)/|y|^3
        =-(rho J/3)[3(n dot y)^2-|y|^2]/|y|^5, s>a.      (9)

There are no omitted far-field multipoles for this radial swirl. No cutoff is applied to pressure. A polynomial radial test can check the integration signs, while C-infinity bump admissibility is the analytic construction in(6).

Smooth local Euler existence for these data is the only PDE import needed. See Terence Tao, [254A Notes3, local well-posedness](https://terrytao.wordpress.com/2018/10/09/254a-notes-3-local-well-posedness-for-the-euler-equations/), pressure formula(2) and the local existence section for H-infinity data. The pressure computation and moment consequences here are rederived, not imported particle claims. No global 3D smoothness assertion is needed for this initial-acceleration test.

## 3. Identical local velocity state, different exact future

Compare the zero initial field with(6) at D=d e_z, n=e_z, d>a+b, and use the identical initial material tag B_b(0). Both have exactly the same velocity (zero) on an open neighborhood of the tag, hence identical position/velocity distributions, I, B=0, T=0, L=0, and every local spatial velocity derivative. In the zero solution its centroid acceleration vanishes. In the swirl solution p is harmonic throughout the tag and the harmonic mean-value property applies to each component of grad p. Thus

    Xddot(0)=-grad p(0)/rho=(2J/d^4)e_z !=0.              (10)

This distinguishes whole nonlinear Euler solutions with identical local observed initial states. It rules out an autonomous law of only those tag position/velocity/spin/shape moments on the unrestricted smooth finite-energy Euler class. It does not rule out a law that includes ambient pressure, a proved invariant family, or an exact nonlocal/memory closure.

The shape consequence is equally exact. For harmonic h in a neighborhood of B_b,

    average_B [x_j partial_i h]=b^2 partial_ij h(0)/5.     (11)

To prove it, expand harmonic homogeneous polynomials; the spherical integral of x_j times a harmonic of degree l vanishes unless l=1. For(9) on the symmetry axis,

    Hess p(0)=(rho J/d^5) diag(4,4,-8),
    Iddot(0)=-(2 M b^2/(5rho)) Hess p(0)
             =(M b^2 J/(5d^5)) diag(-8,-8,16).           (12)

The tag begins with zero spin and shape rate but nonzero trace-free shape acceleration. This removes the unbounded affine-flow limitation of the earlier P252 counterexample. It is not an argument that a restricted ring manifold cannot close; that would require a separate construction and stability test.

## 4. Constructive continuation

The useful state is enlarged by F, P and the pressure-velocity work in(1)–(3), with pressure computed from the same full velocity using(7). This yields exact nonlinear balances immediately. Eliminating the ambient velocity requires deriving its conditional evolution and initial unresolved state; replacing it by a fitted torque loses(10). The next R3 route is an exact projected dynamics with its unresolved initial-data term retained, followed by a test on an actual R1 carrier family. There is no irreducible numerical remainder in(1)–(12); a force fit or long integration would weaken this exposing oracle.

The exact pressure tail is a derived interaction channel, but(10) is an initial acceleration of a resting material tag near a nonsteady swirl. It is not an all-time two-particle force, Coulomb law, electron charge, restoring stability or a neutrino mechanism.

## 5. Two disjoint swirls: exact initial interaction and retained ambient momentum

The same construction provides an exact initial two-excitation interaction, without a prescribed future path. Choose two smooth radial swirls with disjoint velocity supports. At t=0, the product u_1 tensor u_2 vanishes pointwise, so the pressure source for u_1+u_2 is exactly the sum of their sources. Decay fixes p=p_1+p_2. This additivity is an INITIAL disjoint-support statement; the velocity develops pressure-driven tails immediately and later cross terms must be retained.

For spherical material tags centered on the two supports and not intersecting the other support, each self-pressure force averages to zero by inversion symmetry. The external pressure is harmonic in the other tag, so its exact centroid acceleration is minus its center gradient divided by rho. Writing r=X_1-X_2, h=n_2 dot r, s=|r|^2 gives

    X_1ddot(0)=J_2 [2h n_2 s+(s-5h^2)r]/s^(7/2).        (13)

For collinear centers and parallel axes, r=-d n,

    X_1ddot=+2J_2 n/d^4, X_2ddot=-2J_1 n/d^4.            (14)

For an equatorial placement h=0, (13) is outward J_2 r/|r|^5. Thus the interaction is anisotropic and its initial acceleration scales d^-4. Reversing a swirl's circulation does not reverse its quadratic pressure source. These are derived features of THIS initial-data family, not the charge/sign/range properties of electromagnetism.

The tag forces M_1 X_1ddot and M_2 X_2ddot need not sum to zero: their remainder is balanced by the untagged ambient fluid and pressure flux. Assigning arbitrary tag mass to the excitation and discarding that ambient momentum would create a false Newtonian pair law. Here the complete Euler field provides the balance. The profile energy is E_j=(4pi rho/3)J_j, so (9)/(13) also give an exact energy-normalized initial pressure channel. Neither fixed labels nor this identity selects J_j or a particle mass.

This is positive P3 input at exact initial-response scope. Its next achievement is persistence and a controlled interaction law on an actual stable R1 family; the swirls here are admissible initial data, not asserted enduring objects.

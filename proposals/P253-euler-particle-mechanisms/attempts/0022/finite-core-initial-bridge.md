# Exact finite-core Euler initial bridge on a straight vortex background

This construction executes the initial-state part of the background route.
It does not prescribe the later Euler curve from a binormal solution.
Use Cartesian coordinates (x,y,z), fluid mass density rho_m>0, a smooth
radial axial vorticity w(r^2)e_z with compact transverse support, and
circulation Gamma=2pi int_0^infinity r w(r^2)dr. The exact stationary column is

    U=V(r)e_theta, V(r)=r^-1 int_0^r s w(s^2)ds,
    p_0'(r)=rho_m V(r)^2/r, p_0(infinity)=0.             (1)

The velocity is smooth and bounded with bounded derivatives, and has the
Gamma/(2pi r) exterior swirl. Its total energy on R3 is infinite; the
localized perturbation is measured relative to this specified state.

## Volume-preserving bend and actual vorticity

For smooth rapidly decaying a(z),b(z), the map

    Phi(X,Y,Z)=(X+a(Z),Y+b(Z),Z), det DPhi=1

is a global volume-preserving diffeomorphism. Its pushforward gives

    omega_a=w((x-a(z))^2+(y-b(z))^2) (a'(z),b'(z),1).   (2)

The derivatives of the shifted profile cancel the divergence of the
transverse components, so div omega_a=0 exactly. This is vorticity
pushforward; pushing U pointwise would not be the Euler velocity.

For the embedded soliton of filament-construction.md, take z=gamma_1(s,0)
and (a(z),b(z))=(gamma_2(s,0),gamma_3(s,0)). When xi^2>eta^2, z_s is
bounded below by a positive constant, so this defines a smooth rapidly
decaying graph without an implicit multi-valued center. A circular transverse
support gives an actual finite-thickness tube for every positive core size.
Its sections are horizontal; a Frenet-normal circular section is not assumed.

The difference delta_omega=omega_a-omega_0 is smooth, bounded in transverse
support, and rapidly decaying in z with all derivatives. Moreover

    int_R2 delta_omega dxdy=Gamma(a',b',0),
    int_R3 delta_omega dV=0.                           (3)

Define the actual velocity by

    v=curl(-Delta)^-1 delta_omega, u_a=U+v.             (4)

Then div u_a=0 and curl u_a=omega_a. The Newton kernel expansion using (3)
and all finite moments gives v=O(|x|^-3) at infinity, with the differentiated
bounds and v in every H^m. The tail estimate also holds in the axial ends of
the support cylinder because the source there decays rapidly; splitting the
convolution into |source|<|x|/2 and its complement makes this explicit.

## Finite actual excess energy and nonzero relative impulse

The initial excess kinetic energy is the convergent, physically normalized
integral

    E_rel=(rho_m/2) int (2U dot v+|v|^2)dV.             (5)

For the only potentially problematic cross term, |U| is bounded by
C min(r,1/r) and |v| by C(1+sqrt(r^2+z^2))^-3. At large r the cylindrical
measure cancels 1/r, and integrating z leaves an integrable O(r^-2) dr;
at bounded r the z integral is finite. Thus the cross term is absolutely
integrable, not a subtraction of two infinite energies. No positivity or
Euler conservation of the separate NLS curvature deficit is inferred.

Transverse integration of the actual impulse difference gives

    (1/2)int_R2 x cross delta_omega dxdy
      =(Gamma/2)(b-zb', za'-a, ab'-ba').               (6)

Consequently

    I_rel=(Gamma int b dz, -Gamma int a dz,
           (Gamma/2)int(ab'-ba')dz).                 (7)

Physical impulse is rho_m I_rel. The exact Gaussian check integrates (6)
from the three-dimensional field rather than assigning a filament moment.
For a=e^-z^2,b=z e^-z^2, it gives
I_rel=(0,-Gamma sqrt(pi),Gamma sqrt(pi/2)/2).

## The bend is dynamically accessible, but does not create relative helicity

The interpolation Phi_tau=(X+tau a(Z),Y+tau b(Z),Z) is a smooth orbit path.
Although its literal generator (a(z),b(z),0) is not transverse-compact, a
planar Hamiltonian chi=a(z)y-b(z)x times a radial cutoff has velocity
(partial_y chi,-partial_x chi,0) agreeing with it on a neighborhood of the
union of all transported supports. Choose the cutoff identically one on
that fixed compact transverse disk. Its z decay is rapid. It generates the
same vorticity path and supplies a legitimate volume-preserving lift.

With u_tau=U+curl(-Delta)^-1(omega_tau-omega_0), define

    H_rel(tau)=int (u_tau-U) dot (omega_tau+omega_0)dV.  (8)

All terms converge: the background vorticity is confined transversely and
the perturbation velocity decays axially. Integration by parts has vanishing
boundary terms (U=O(1/r), perturbations dipolar/rapid), so differentiating
(8) gives 2 int u_tau dot partial_tau omega_tau. As
partial_tau omega_tau=curl(v_tau cross omega_tau), this is exactly
2 int omega_tau dot(v_tau cross omega_tau)=0. Thus H_rel=0 along this pure
bend. This is an exact relative-Casimir calculation, not a statement that
all swirling or twisted backgrounds have zero helicity. A geometrically
chiral tube by itself does not earn an extra helicity sector.

## Local Euler time evolution is an initial-value problem

For an integer m>=3, write the actual Euler evolution as U+v(t). Leray
projection gives the full equation

    v_t+P[U dot grad v+v dot grad U+v dot grad v]=0.    (9)

U is stationary, so there is no omitted background forcing. In the H^m
energy estimate, the leading transport by divergence-free U+v integrates
to zero. Leibniz commutators and H^m embedding into W^(1,infinity) give

    d||v||_Hm/dt <= C_m (||U||_W^(m+1,infinity)+||v||_Hm)||v||_Hm. (10)

The highest U dot grad derivative cancels before estimating; multiplication
by bounded derivatives of U is the relevant background assumption, not
U in L2. Friedrichs regularization gives a uniform short-time bound from
(10), difference estimates in H^(m-1) give convergence and uniqueness, and
smooth initial data propagate higher regularity on this local interval.
Therefore (2)-(4) has an actual smooth local-time Euler lift on this specified
background. Its existence time can shrink with core size; (10) provides no
thin-core-uniform or all-time concentration estimate.

The initial finite excess energy, relative impulse, and orbit path are exact.
A later-time relative-energy conservation assertion additionally requires
propagated decay sufficient for its boundary flux; local H^m regularity alone
is not used to smuggle in that stronger assertion. The local kinetic
perturbation balance that follows without those extra weighted assumptions is

    (1/2)d||v||_L2^2/dt=-int v_i v_j partial_j U_i.     (11)

This retains the physical exchange with the background. It is not a positive
restoring deficit.

## Route verdict and immediate next construction

The smooth finite-core initial geometry, Hodge velocity, finite excess energy,
relative impulse and local Euler evolution are established as stated. The
persistent particle supplier remains active. The next positive task is to
construct a relative traveling wave or a controlled modulation on this same
column and prove persistence under a neighborhood of initial perturbations,
retaining its full pressure, exterior velocity and relative invariants.
A reduced-model soliton cannot close that task by changing the energy's name.

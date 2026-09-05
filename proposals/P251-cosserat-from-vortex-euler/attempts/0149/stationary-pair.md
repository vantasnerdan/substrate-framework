# A stationary periodic vorticity pair with a positive physical director mode

## 1. Exact background and Green Hessian

Let Lambda=L Z^2, A=L^2, and normalize the real zero-mean Green function
by Delta G=delta-1/A. Use the same Euler/Green convention as0139:
u=J grad psi and H=-rho Gamma^2 G(r) for two equal vortices with
relative separation r; all position-independent self terms are omitted.
The actual total vorticity is

    omega=Gamma(delta_z1+delta_z2)-2Gamma/A, Gamma>0.          (1)

The constant part is transported Euler vorticity, not a rotating container.
Take z1=-r0/2,z2=r0/2 with r0=(L/2,L/2). Since -r0=r0 modulo Lambda,
evenness of G gives grad G(r0)=0. Thus both vortices are exactly stationary.
The quarter-turn also fixes r0 modulo Lambda. Consequently Hess G(r0)
is a scalar matrix; its trace is the defining smooth Poisson value -1/A.
This determines the full Hessian exactly, without summing a truncated
Fourier series or assigning a numerical lattice constant:

    Hess G(r0)=-I2/(2A), Hess H(r0)=h I2,
    h=rho Gamma^2/(2A)>0.                                   (2)

The unreduced four-position Hessian is h[[I,-I],[-I,I]]. Its only
kernel consists of simultaneous translations; both relative directions
are strictly positive. Square symmetry and the Poisson trace together
are essential: trace alone would not determine the sign of both entries.

This is also a single finer square vortex lattice viewed in a doubled
cell. The optical mode below is a folded finite-wave-number lattice mode,
not a newly invented independent particle rotation. Selecting the two
neighboring vorticity cores gives an actual pair director in a continuous
local lift of their separation; exchanging the labels changes its angle
by pi. The angle is not an arbitrary orientation of passive mass labels.

## 2. Full point-Euler relative action and geometric angle

Put Z=(z1+z2)/2 and r=z2-z1. The actual point-vortex KKS form is

    Omega=2rho Gamma dZx wedge dZy
                   +(rho Gamma/2) drx wedge dry.             (3)

The relative mode has fixed Z and no harmonic impulse change. Let
d=|r0|, r=(d+a)(cos theta n+sin theta Jn), n=r0/d. At the equilibrium
a=theta=0, its quadratic action, in the conventions of0139, is

    L2=-B a theta_dot-h(a^2+d^2 theta^2)/2,
    B=rho Gamma d/2.

The conjugate momentum is p=-B a. Elimination gives

    Ltheta=I theta_dot^2/2-kappa theta^2/2,
    I=B^2/h=rho A d^2/2,
    kappa=h d^2=rho Gamma^2 d^2/(2A),
    omega_opt^2=kappa/I=Gamma^2/A^2.                          (4)

Both the physical angle and the conjugate radial reaction are actual
vorticity-centroid motions. There is no pinned external core or prescribed
strain. The relative point system is genuinely one degree of freedom,
not a guessed invariant subspace: H depends only on r and Z is constant.
The nondegenerate minimum gives nearby closed energy curves and periodic
point trajectories by the local Morse lemma and the nonzero symplectic
form. Global fluid rigid rotation rotates the lattice as well; (4) is
relative to that lattice, not an absolute-space locking energy.

The simpler dilute signed-center polygon is not automatically a substitute.
For six positive ring vortices, stationarity against an isolated central
vortex at leading order requires Gamma_c=-5Gamma/2. On the alternating-
triangle fixed-total-radius mode of0036, its added radial Hessian is
3rho Gamma Gamma_c/(2pi), while the original radial Hessian is
3rho Gamma^2/(8pi). Their sum is -27rho Gamma^2/(8pi), whereas the
relative-angle Hessian remains positive. This candidate has a hyperbolic
mode. The periodic positive pair avoids that specific sign mechanism.
The argument does not exhaust central polygons with other surroundings.

## 3. Exact smooth stationary Euler and actual finite-time histories

There is a direct smooth-core construction using the bordered inverse
already proved in0036/0139. On a fixed rescaled disk near z1 the extra
regular Green kernel is

    K_epsilon(y,z)=R(epsilon(y-z))+G(r0+epsilon(y-z)),
    R=G-log|x|/(2pi).

Subtract K_0 into the constant border and solve

    U-mu+N*F(U)+integral[K_epsilon-K_0]F(U)=0,
    integral F(U)=Gamma.                                    (5)

The same profile is copied to z2. Both kernels are quarter-turn invariant
and their linear jets at zero vanish. Restricting U to C4-invariant
functions excludes the translation kernel exactly as C6 did in0139;
the radial mass border is unchanged. The perturbing operator is O(epsilon^2)
on the fixed disk. Its inverse and the negative collar persist by IFT.
The fixed smooth flat F and elliptic bootstrap yield globally C-infinity
vorticity

    omega_epsilon=epsilon^-2 sum_i F(U((x-zi)/epsilon))
                                              -2Gamma/A.    (6)

Equation(5) makes each excess core a function of its actual streamfunction;
outside, vorticity is constant. Hence (6) is exactly stationary Euler,
not a smooth approximation inserted into the steady equation.
This derivation uses smooth vorticity. The separate Sakajo--Zou torus
patch theorem corroborates the equilibrium setting but only supplies
smooth patch boundaries, not the smooth-vorticity assertion in(6).

Opposite compact Hamiltonian translations of the two cores preserve
all Casimirs, total harmonic impulse and the compensating vorticity.
Their KKS is exactly(3), independent of core radius. Translation self
energies cancel and every remaining Green kernel is smooth at the
separated centers; its C2 interaction error is O(epsilon^2).
Thus the exact restricted finite-core Hessian retains both positive
relative directions for sufficiently small nonzero epsilon.

Actual evolution is supplied separately, not inferred from this restricted
action. The positive-blob confinement proof transferred to the torus in
0139 applies verbatim to these two separated paths on each fixed finite
interval. Choose the desired number of optical periods first, perturbation
amplitude delta next, and epsilon last. The exact transported excess-core
centroids obey

    sup_t |Z_i,epsilon(t)-z_i,linear(t)|
                          <= C_T delta^2+C_T epsilon^alpha,
    0<alpha<1/3.                                            (7)

Since d>0, the actual geometric director angle of their separation has
the same controlled error divided by d. Choosing both errors small
relative to delta realizes nonzero smooth Euler optical director histories.
The exact finite-core phase action and these actual histories are both
present. This does not claim an exact fixed-core spectral pole or an
invariant two-dimensional smooth Euler manifold.

The exact0136 Bernoulli lift extends this smooth planar stationary base
to a stationary three-dimensional generalized force-free column, and its
z-independent dynamical subsystem is actual Euler. Variable curl factor
is not silently promoted to constant-lambda EPS compatibility.

## 4. Physical spin and the remaining join

The canonical momentum p in(4) is an Euler orbit momentum. Its equality
to an actual selected material tag's mass angular momentum is not asserted.
For two very small core mass tags of equal mass M with initially stationary
centroids, the centroid spin alone is M d^2 theta_dot/2. Matching I from(4)
would require M=rho A, larger than the mass available to each disjoint
half-cell tag. Surrounding fluid and internal/boundary spin therefore
matter; assigning core mass to I is not a valid normalization.

There is a useful exact Eulerian identity, also exposing the distinction.
On the centered square fundamental domain, periodic velocities satisfy

    integral r cross u = -1/2 integral |r|^2 omega,

because opposite-edge |r|^2-weighted tangential boundary integrals cancel.
For the relative pair perturbation the right side varies by -Gamma d a/2,
so rho times this Eulerian cell spin equals p initially. A material cell
subsequently has a moving boundary; that boundary and the unperturbed
transported reference must be included before using this fixed-domain
identity as a physical spin/current law. This is the next executable
observation construction, not a missing arbitrary coefficient.

Established route: exact stationary positive pair, full relative director
action, smooth stationary continuation and controlled actual optical
director histories. Remaining: actual material-spin/centroid observation,
long-wave coupled three-dimensional dynamics, and constant-curl knotted
stationary realization. The stronger parent objective remains active.

Primary contextual source: Sakajo and Zou, Steady Vortex Patches on Flat
Torus with a Constant Background Vorticity, Theorem2.6,
https://link.springer.com/article/10.1007/s00332-025-10236-6.
Its rank2N-2 equilibrium condition agrees with(2); it is not used as
a smooth-vorticity or dynamical-mode theorem.

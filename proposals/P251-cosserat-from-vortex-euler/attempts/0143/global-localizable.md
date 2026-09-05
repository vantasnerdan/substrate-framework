# Exact stationary global cells and the interaction they do not supply

## 1. Source license and the pressure localization identity

Constantin--La--Vicol, [2019 primary manuscript](https://cims.nyu.edu/~vicol/CLV1.pdf),
Theorems1/2 and the construction following(108), gives nonzero C-infinity
compactly supported stationary Euler velocities in toroidal shells.
The local pre-cutoff velocity at the central circle is only Hölder;
its swirl is proportional to sqrt(psi). The smooth construction places
the cutoff in an annulus excluding that circle. Thus this import does
not supply a smooth nonzero central Rankine core by assertion.

The localization mechanism is elementary and exact. Let smooth v,p solve
steady Euler, div v=0, rho(v.dot(grad))v=-grad p, and v.dot(grad p)=0.
For smooth chi define

    u=chi(p)v, P(p)=integral^p chi(s)^2 ds.

Then div u=0 and rho(u.dot(grad))u=-grad P. The missing term in an
arbitrary cutoff would be rho chi chi' (v.dot(grad p))v; it vanishes
here by the actual pressure first integral. Choosing chi supported
between compact regular pressure levels inside the smooth annulus and
flat at its endpoints yields a global C-infinity zero extension. The
pressure is constant in the exterior, and its gauge can be fixed there.
No surface force, wall, or vorticity sheet is inserted.

The exact local torus in0136 only has constant Bernoulli, not a proved
pressure first integral. Therefore this operation does not directly
globalize that torus. In particular a generalized Beltrami field with
finite total energy cannot be inserted as the nonzero compact template.

## 2. One global smooth isotropic stationary Euler ensemble exists

Fix one such nonzero compact template v_*,p_*, with support in a ball
of radius a, and p_*=0 outside the ball. Take a locally finite set of
centres with separation greater than2a and independent rotations R_j.
The single actual field

    u(x)=sum_j R_j v_*(R_j^T(x-X_j)),
    p(x)=sum_j p_*(R_j^T(x-X_j))

is globally C-infinity and solves stationary Euler. At every point at
most one velocity is nonzero; hence every mixed convective product
vanishes identically. Smooth flat support boundaries justify the same
statement there, and the displayed pressure supplies the complete
global Euler equation, not a local pressure replacement.

A stationary isotropic hard-core process of positive intensity nu
provides centres. One explicit choice is a homogeneous Poisson process
retained only at points having no other point within distance d>2a;
its intensity is nu=lambda exp(-lambda volume(B_d))>0. This rule is
translation/rotation covariant. Independent Haar rotations of each
template therefore give a stationary isotropic field law. This is a
single Euler field per sample, not an average of incompatible local
equations. Fixed bounded template derivatives give uniform local
smoothness, and its finite positive energy density is

    E_density=(rho nu/2) integral |v_*|².

The pressure-localized toroidal shells are genuine invariant regions.
Their existence does not supply the arbitrary knotted tube topology or
constant-factor Gaussian support used by the accepted EPS construction.
The 2023 Enciso--Peralta-Salas survey lists arbitrary-topology compact
stationary tubes as an open problem there. The newer primary preprint
[Peralta-Salas--Slobodeanu2606.13462](https://arxiv.org/html/2606.13462v1),
Theorem1.1, proves axisymmetry for analytic localizable fields tangent
to a smooth bounded boundary with constant boundary Bernoulli and
nonzero boundary Bernoulli gradient. It concerns that precise analytic
class, not all C-infinity stationary Euler fields or all global joins.
No parent no-go is inferred from it.

## 3. Rigid-pose stiffness vanishes exactly in this construction

Independent small translations/rotations of each template, keeping the
supports disjoint, give another EXACT stationary Euler field. They are
also genuine Kelvin preparations: extend the rigid Euclidean map near
each entire support to a volume-preserving map in its disjoint ball.
On the support its derivative is orthogonal, so coadjoint transport of
the velocity one-form gives the rotated/translated template; off the
support the one-form is zero. The resulting velocity is already
divergence free, and the full pressure projection changes nothing.

Consequently the kinetic Hamiltonian is constant on the whole
finite-dimensional rigid-pose family:

    H(X_1,R_1,...,X_N,R_N)=N H_*.

All its first and second pose derivatives, including every intercell
cross derivative, vanish. Differentiating the actual stationary Euler
equation also proves each pose derivative v_a is a zero-frequency
solution of the complete linearized Euler equation:

    P[(u.dot(grad))v_a+(v_a.dot(grad))u]=0.

This is an exact mechanism, not a numerical near-zero eigenvalue. It
refutes the route that would derive positive centre/frame elastic
stiffness solely from independent rigid poses of these disjoint cells.
It does not refute internal shape modes, dynamical nonlocal response,
or interacting noncompact stationary backgrounds.

## 4. Stress and physical observations are retained

For a compact stationary template, integration of
partial_j(rho u_i u_j+p delta_ij)=0 against x_k gives

    rho integral u_i u_k=-delta_ik integral p.

Thus its integrated velocity stress is already isotropic, even before
Haar averaging. This is not an elastic modulus: the pose family above
has exactly constant energy. An affine field constraint can have a
nonzero Hessian while free rigid repositioning does not; the actual
admissible family decides which proposition is established.

Also integral u=0, since u_i=partial_j(x_i u_j) for a compact solenoidal
field. Pose derivatives therefore have zero total linear momentum.
They are not a positive-mass bulk displacement mode merely because
their support centre moves in a family of equilibria. Actual transient
perturbations can induce noncompact velocities and pressure; that
dynamical problem is not replaced by this zero-frequency family.

## 5. Route verdicts and executed continuation

Global smooth stationary isotropic field assembly: established as stated.
Rigid-pose elasticity from disjoint compact supports: refuted by the
exact stationary Kelvin family and constant Hamiltonian. Neither
verdict settles the full Euler--Cosserat objective.

Method repair: retain the true pressure first integral and exclude the
nonsmooth central circle in the source construction. Representation
change: use the exact full Kelvin pose family to distinguish stationarity
from elasticity. Materially different continuation is already executing:
0141 uses an interacting noncompact periodic Bernoulli-lifted array;
0142 uses a constant-factor field with direct EPS approximation;0144
constructs the general three-dimensional mean/ambient response. These
routes preserve interactions instead of treating disjoint support as
the desired constitutive mechanism.

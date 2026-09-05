# Material spin produces an actual centroid/ambient first-gradient response

This is a balance/observation theorem, not a constitutive assumption.
The actual smooth incompressible Euler fluid has density rho>0. Select
disjoint material tags D_a(t), keep their entire continuous ambient
complement, and use Fourier convention exp(-i k.x). No artificial wall
is applied to a tag. All boundary pressure and advective exchanges remain.

## Exact moving-tag identity and finite-radius error

For one tag let M=int_D rho, X=M^-1 int_D rho x, r=x-X,
P=int_D rho u=M Xdot, I_ij=int_D rho r_i r_j, and
S=int_D rho r cross u. Reynolds transport gives

    C_ij=int_D rho u_i r_j,
    Idot=C+C^T,
    C_ij=Idot_ij/2-epsilon_ijl S_l/2.

Here I is the central second mass moment, not the rotational inertia
tensor tr(I)1-I. In particular non-rigid tags generally have Idot!=0.
Direct Taylor expansion of the defining Fourier integral gives

    int_D rho u exp(-i k.x)
      =exp(-i k.X)[P+i k cross S/2-i Idot k/2+R_2],
    |R_2| <= |k|^2 int_D rho |u| |r|^2/2.             (1)

This bound follows from |exp(-iy)-1+iy|<=y^2/2, not from a numerical fit.
The centroid phase is retained. The exact formula applies before
linearization, so variation also differentiates the domain, centroid,
phase and ambient indicator. The shape rate and spin are actual material
observations, not the fixed-domain integral of an Eulerian perturbation.

Define J_E as the fully resolved point momentum and J_H as the momentum
of actual tag centroids plus the complete ambient point momentum. Summing
(1), with S and Idot understood as their tagged Fourier densities, yields

    J_E=J_H+i k cross S/2-i Idot k/2+O(|k|^2).        (2)

The ambient cancels from this comparison, but is present in BOTH currents.
Consequently a vanishing complete-packet dipole is not a proof that
the actual centroid current has no spin-coupled response.

## Euler balance, isotropy and an actual changing spin

Freeze a stationary ensemble of actual backgrounds/tags and a covariant
linear preparation driven by an axial vector a. Suppose its full finite-
time response has a controlled first jet along k=epsilon n. This is an
application hypothesis, not supplied by a formal moment identity. In a
periodic construction0116 supplies the Bloch family and finite-time
energy estimate; its mean projection at k=0 is prepared explicitly.
SO(3) averaging is over the COMPLETE response and its actual preparation.

The exact point momentum equation is

    d_t delta J_E=-i k_j delta T_ij,
    T=rho u tensor u+p 1,       T=T^T.              (3)

At k=0 an isotropic linear map from an axial vector to a symmetric
rank-two tensor vanishes: its only SO(3)-invariant rank-three tensor is
epsilon_ijl, whose symmetric part is zero. This applies separately to
the full Euler stress response and the central shape-rate response.
It neither removes their fluctuations nor sets the microscopic tensors
to zero. Equations (2)-(3) therefore imply, on the fixed finite interval,

    delta J_H(t)-delta J_H(t0)
      =-i k cross [delta S(t)-delta S(t0)]/2
         +O(|k|^2 C_T).                              (4)

The remainder includes the actual Taylor radii, response regularity and
time-integrated first stress jet. A reflection-paired smooth response can
improve the axial-to-polar parity remainder, but (4) claims only the
first jet and does not need that improvement. Longitudinal components
are handled by the actual incompressibility/mean projection; the displayed
spin term is transverse. The total leading mass is rho, not a filling
fraction. If U_H is defined by delta J_H=rho U_H,dot, then (4) is a
physical centroid-plus-ambient velocity response to changing spin.

The relation also holds after any temporal linear observation annihilating
constants, with its own finite norm multiplying the error. For a genuine
single nonzero frequency omega and an independently derived physical law
S=j Phi_dot, it yields

    U_H=-i j k cross Phi/(2 rho)+O(|k|^2).

The latter constitutive hypothesis is NOT inferred from (4). In0112/0114
the actual spin has a time-dependent geometric observation row. That row,
not a renamed canonical momentum, is inserted in (4).

## A nonzero isotropic response can be prepared without fitting a modulus

Suppose one actual tag has a nonzero transverse endpoint spin row D and
an invertible transverse core-angle observation A. Choose a unit normal
axis n for which n^T D is nonzero; n^T A is then nonzero. There is a phase
vector z0 with n^T A z0=1 and n^T D z0!=0: an affine hyperplane on which
the latter functional vanished identically would force that functional
to be zero (or proportional to the former with zero proportionality).
Use the fixed actual preparation z0(n.a) and rotate the COMPLETE marked
background, tag and preparation by R in SO(3). Haar averaging gives

    <R(Az0)n^T R^T>a=a/3,
    <R(Dz0)n^T R^T>a=(n^T Dz0)a/3 !=0.              (5)

Rescale the declared coherent input by three if unit mean angle is the
chosen convention. The spin coefficient is a derived moment of that
preparation, not a supplied inertia.0114 constructs D!=0 for an actual
periodic-covariance tag through its elliptic Euler return. The coherent
first jet (4) is then nonzero for wave directions not parallel to a.
Uniform orientation marginals alone would not prove this: the preparation
is correlated covariantly with each entire marked state.

## What the construction supplies next

The exact physical spin/shape map and (4) bridge actual same-core angular
dynamics to a genuine coarse centroid response. The surrounding pressure
reaction is essential, not omitted because the packet is compact initially.
Positive two-dimensional phase energy and a nonconstant observed spin do
not yet fix the full vector inertia/constitutive pencil of C-CST-009.
That final coefficient and spatial-dynamics bridge remains distinct.

The new euler_observation API reuses the existing discrete_mass_moments
definition and tests (1) against independently differentiated Fourier
sums. A deforming tag exposes the missing shape term; reversing the spin
sign exposes the orientation convention. These exact finite fixtures
verify the algebra, not the analytic Euler/Bloch application hypotheses.

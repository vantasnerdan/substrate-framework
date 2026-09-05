# Nineteen exact physical tube/ambient response constraints

## Stationary physical tags and a complete one-fluid partition

Let u0 be the smooth stationary Beltrami field, curl u0=lambda u0, and D
one of its actual invariant solid-torus EPS domains, with a nontrivial
trefoil core. Set omega=lambda u0. The material tag chi_D obeys

    u0 dot grad chi_D=0.

Because u0 dot n=0 on its boundary,

    integral_D u0 = integral_boundaryD x (u0 dot n) = 0.

The mass centroid X_D is stationary. The domain and every mass moment
integral_D rho f(x-X_D) are stationary in the background. This is an
actual material parcel, not a fixed Eulerian box substituted for one.
For disjoint selected domains, chi_A=1-sum chi_D is the invariant
ambient tag. Its material action uses its normalized volume measure;
no finite centroid is assigned to this generally unbounded phase.
Masses of all tubes plus the ambient mass give exactly rho per unit
total volume. The shared pressure and velocity are those of the same
fluid, and their interface fluxes remain in the action.

## Actual moment functionals

Put r=x-X_D. For a compact curl-generated isovortical direction xi with
support away from the tube boundary, set

    F_xi=xi cross omega, v_xi=P_R3 F_xi.

Its physical tag variation is -xi dot grad chi_D=0. Its material centroid
position variation is zero if every interior piece is a compact curl.
At this reference tangent, there is no moving-boundary contribution to
the following velocity moments:

    M_D(a,h;xi)=rho integral_D (a+h r) dot v_xi,
    tr h=0.

The eleven components are three actual tube centroid momenta, five
tracefree symmetric first moments, and three actual intrinsic angular
momenta. In particular h r=b cross r gives b dot S_D, not the whole-space
vortex angular impulse. The trace part of the first moment is retained
as a pressure/shape scalar; it does not pair with incompressible affine U.

By Leray self-adjointness these are compact test functionals of xi:

    M_D(a,h;xi)=integral xi dot f_D(a,h),
    f_D(a,h)=rho omega cross w_D(a,h),
    w_D(a,h)=P_R3[1_D (a+h r)].

Also retain the eight whole-space translation and STF affine KKS rows

    M_G(a,E;xi)=rho integral (a+E r) dot F_xi,
    E=E^T, tr E=0.

These expressions are finite on compact F before any decay condition
is imposed. After their moments are set to zero, the induced velocities
have the 0069 compact-curl plus double-divergence representation. They
then have the integrable angular tails needed for the equivalent
velocity-integral interpretation. No global angular row is constrained.

## Exact independence of the nineteen rows

Use the same permitted prototype as 0059: an EPS Beltrami field whose
smooth Herglotz density gives O(1/r) field and derivative bounds, plus
arbitrarily small same-helicity ABC atoms in three linearly independent
wave directions k_i=lambda e_i. The perturbation is small enough to
preserve the chosen KAM invariant trefoil tube. The proof below concerns
this prototype; a strict local Gram margin transfers to the bounded
positive-probability good-patch event by continuity, not by assuming
Gaussian samples have its far-field expansion.

Suppose a linear combination of the nineteen response kernels has zero
curl on an exterior open ball. They are analytic throughout the
connected exterior of D. Analytic continuation therefore gives

    curl[omega cross (a+E r+w_D(a_D,h_D))]=0

on that entire exterior. Both vector fields are divergence-free, so
this is a vanishing commutator (its overall sign is immaterial).

Along a generic ray x=r n, the three positive and three negative ABC
frequencies are distinct. The E r term times grad omega has leading
order r. All derivatives of E r, and the EPS contribution, have lower
order. Dividing by r and taking the elementary long-ray Fourier/Cesaro
mean for each distinct frequency forces (E n) dot k_i=0 for every i.
The set of generic n is open and dense; consequently E=0. The next,
order-one, coefficients similarly force a dot k_i=0, hence a=0.

The remaining w_D is harmonic and curl-free in the exterior, and is
O(r^-3). Unless it vanishes there identically, its convergent exterior
harmonic multipole series has a first nonzero term r^-m W(n), m>=3.
In the commutator, w_D dot grad omega_atom is order r^-m; the derivative
of w_D and all EPS terms are O(r^-m-1). Isolating each of the three
distinct frequencies gives W(n) dot k_i=0. Their independence forces
W=0, contradicting its selection as the first nonzero multipole.
Thus w_D vanishes throughout the exterior.

Write w_D=1_D(a_D+h_D r)-grad p. In the exterior p is a constant, chosen
zero. The Newtonian single-layer potential p is continuous across the
smooth boundary, and Delta p=tr h_D=0 in D. The maximum principle gives
p=0 inside too. Distributional incompressibility now implies

    (a_D+h_D r) dot n=0 on boundary D.

Hence D is invariant under this affine one-parameter flow. A bounded
domain with nonempty interior admits no unbounded affine group: its
action on an interior ball bounds both exp(t h_D) and the translation,
for all positive and negative t. The closure is a compact affine group.
Haar averaging gives a fixed point and an invariant positive inner
product. Relative to that point, h_D is skew in this inner product. In
three dimensions every nontrivial such one-parameter group is conjugate
by an invertible linear map to ordinary axial rotations.

A solid torus invariant under ordinary axial rotations is unknotted.
Here is the elementary geometry needed for this particular assertion,
not a claim about arbitrary circle actions on S^3. Its boundary torus
cannot meet the axis: a connected compact surface of revolution meeting
the axis has interval orbit space and collapsed endpoint circles, and
is a sphere rather than a torus. Away from the axis the domain is
Sigma cross S^1 in cylindrical coordinates. Since its fundamental
group is Z, the planar cross-section Sigma is a disk; a round rotation
orbit is its core. Invertible affine maps do not turn that core into a
trefoil. The actual trefoil tube D therefore has no such nonzero affine
generator. It follows that a_D=h_D=0.

All nineteen coefficients are zero. Thus, for every nonnegative smooth
chi supported in an exterior ball and positive on a nonempty open set,

    G_AB = integral chi curl f_A dot curl f_B

is strictly positive definite. Its smallest eigenvalue is a
geometry-dependent constant. The good-patch event records a positive
lower bound for this actual matrix and upper derivative norms; there is
no universal formula depending only on circulation and core radius.

The compact-group averaging step is the usual explicit invariant-inner-
product construction; see Brian Conrad's primary lecture notes
[Compact subgroups](https://virtualmath1.stanford.edu/~conrad/210CPage/handouts/maxcompact.pdf).
The multipole and affine/topological arguments above are given here in
full rather than imported as an unstated symmetry theorem.

## Disjoint dual responses and exact physical spin prescription

Choose nineteen pairwise disjoint finite exterior response balls,
also disjoint from every raw core and carrier support. On ball j use
its positive Gram matrix to form one compact curl field eta_j satisfying

    M_A(eta_j)=delta_Aj.

Explicitly, eta_j is the jth normalized column of the fields
curl(chi_j curl f_A). Each eta_j is supported in its own ball, so

    Omega(eta_i,eta_j)=0 for all i,j,
    Omega(raw Q,eta_j)=Omega(raw S,eta_j)=0.

Define Pi=I-sum_j eta_j M_j. For the fixed physical core-angle generator
Q0 and its negative-helicity reaction cage S0, put B=Omega(Q0,S0)!=0 and

    Q=Pi Q0,
    S=Pi S0+B sum_i n_i eta_(tube-spin,i).

The exact conclusions are

    Omega(Q,S)=B,
    M_G(Q)=M_G(S)=0,
    M_D(Q)=0,
    tube-centroid(S)=tube-STF(S)=0,
    S_D(S)=B n.

No response changes either raw physical core jet. All Gram corrections
are fixed finite fields once the actual good patch and geometry are
chosen. Their complete cross energies and ambient Leray velocities
remain in H. The high-carrier positivity argument transfers because
the fixed responses and their mixed terms are bounded while the local
negative-helicity contribution grows linearly in the carrier. The
parent 0083 supplies the full kinetic-Schur version; no isolated
reaction inverse or pressure-free interior energy is introduced here.

Further gradient-only cages can be Pi-projected with zero prescribed
nineteen rows. Their KKS pairings with these Q/S fields still vanish
exactly by the same disjoint-response identity. Thus the full gradient
transfer uses the actual tube-spin construction, not coefficients from
an unrelated background.

## Exact reference affine/cotangent pairing

For a divergence-free affine U(x)=U(X_D)+h r write
h=E+[beta cross], E symmetric tracefree. Define the selected affine
configuration lift by

    Gamma_U=U-beta Q

for this angle direction, with the analogous sum over the marked
directions in the ensemble. It has the actual centroid displacement
delta X_D=U(X_D), since the compact interior Q pieces integrate to zero.
Near the boundary it is U, so the tube's tag has the ordinary material
affine variation -U dot grad chi_D. At the observed physical core, the
rotation of Gamma_U is subtracted by Q; an independent absolute core
angle Phi multiplies Q. Common rigid rotation satisfies
Gamma_K+Q=K exactly. Slow envelopes are formed by the complete curl
potential construction, including all gradient returns, not by
multiplying a divergence-free field and discarding its new divergence.

For every retained reaction v_s the moment prescription proves

    Omega(Gamma_U,Ss)
      =rho integral U dot v_s-beta dot D*s
      =rho integral_ambient U dot v_s
          +sum_D U(X_D) dot integral_D rho v_s.

The STF contribution is zero by construction and the tube spin equals
D*s. The same equation holds for Q because all its tube moments vanish.
This is precisely the hybrid physical tube-centroid plus ambient
momentum functional, restricted to the selected tangent space. It is
not the whole-fluid point mean from 0072. The equality identifies an
explicit volume-preserving affine configuration lift and its actual
moment pairing; it is stronger than simply renaming the momentum field.

The physical tube spin on a general material tangent also contains

    rho integral_boundaryD r cross u0 (delta g dot n),

and any centroid/background terms. These vanish on Q/S here because
those configuration generators vanish near the boundary, but not on
arbitrary boundary/shape directions. They must be retained in the full
material reconstruction and its subsequent momentum elimination.

## Remaining material joining, kept explicit

These results establish the stationary tagged geometry, actual tube
moment functionals, nineteen-row response construction, and exact
reference affine/current pairing. They do not independently prove that
an arbitrary path confined to Q/S reconstructs the correct material
tag evolution. Its induced Leray velocity may have a nonzero normal
trace, even though the Q/S configuration generators vanish near the
boundary. A complete material pullback retains the boundary/shape and
ambient coordinates, or supplies the specified Cauchy--Born material
lift and its pressure response before any elimination. It cannot use
the free-space P unchanged while silently imposing a new no-normal-flow
boundary condition. The independent 0082 joining work and parent full
kinetic block address this distinct obligation. No autonomous physical
Cosserat claim is promoted solely by the moment matching above.

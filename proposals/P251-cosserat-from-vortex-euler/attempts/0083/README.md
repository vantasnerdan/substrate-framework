# 0083 — stationary phase mass and kinetic-centering coercivity

Owner `/root`; parent P251 / issue #198 remains active. This bounded
analytic continuation supplies two inputs to 0080/0082, not their missing
material-action identification. Accepted base v0.171.0 and the original
conditional slow-affine objective are unchanged.

Frozen route: (i) transport actual invariant EPS domains, retaining their
continuous ambient complement; (ii) bound an exact kinetic-Gram Schur
subtraction by the full Leray kinetic Gram before choosing carriers. The
alternative of arbitrary advected cubic parcels is preserved in 0075;
their reference-time identities do not establish a stationary shape law.
No empirical comparator, numerical eigenvalue, fitted coefficient, or new
all-wave-number Euler closure is involved. The oracle is transport calculus,
orthogonal projection, and exact noncommuting block algebra.

## 1. Stationary tagged mass, without an ambient centroid

Let a smooth stationary incompressible Euler field have disjoint bounded
invariant domains D_a: u.n=0 on every boundary. Write chi_a=1_Da and
chi_A=1-sum_a chi_a. Distributionally u.grad chi_a=0 and likewise for
chi_A. The reference tags solve the material transport equation exactly;
the domains, volumes, centroids and mass-moment tensors are time independent.
For each bounded domain, the divergence theorem gives

    integral_D u_i = integral_boundary x_i (u.n) = 0,
    integral_D (x_i u_j+x_j u_i)
      = integral_boundary x_i x_j (u.n) = 0.

Thus stationarity here controls the actual material shape moments, not
merely the Eulerian law of the velocity. Internal spin need not vanish.
For the continuous ambient phase there is no finite centroid assertion.
Its mass, action and momentum are kept per unit volume.

Use the stationary marked selection of 0071 with retained intensity nu,
bounded domains inside disjoint R-balls, and the actual domain volume V
as a mark. Campbell's identity gives tube volume fraction

    f = nu E_Palm[V],   rho_tube=rho f,
    rho_ambient=rho(1-f),   rho_total=rho.

Disjointness gives 0<=f<=1 pointwise after spatial averaging, and the
positive-probability nonempty domains give f>0. Sparse Poisson isolation
can ensure f<1, for example by choosing tau Vol(B_R)<1 and using
nu<=tau, V<=Vol(B_R). All retained phase labels and geometry are transported
under variations, not selected afresh to optimize energy. These identities
require the domain reconstruction in 0080 to be measurable and equivariant;
they do not construct that reconstruction themselves.

A common uniform translation adds exactly rho |Udot|^2/2 to the normalized
quadratic kinetic action: the base mean velocity is zero by the centered
stationary law and the phase fractions add to one. Finite-radius internal
inertia and gradient-dependent terms are not replaced by this zeroth-order
identity. Nor does the identity eliminate independent phase-relative
translation; its retention or the declared affine constraint belongs to
the explicit action in 0082.

## 2. The full kinetic Gram survives exact centering

Let V be the Hilbert space of all ambient-inclusive velocity variations,
with mass inner product, and let B:Z->V and A:Y->V be respectively the
retained macroscopic momentum-velocity and reaction-velocity maps in ONE
material kinetic action. Suppose M=B*B is coercive on its retained/gauge-fixed
space. Write G=A*A and C=B*A. Then

    G-C* M^-1 C = A*(1-Pi_B)A >= 0,
    Pi_B=B(B*B)^-1 B*.

This statement retains all interactions and is not an isolated-patch inverse.
The same conclusion uses the Moore-Penrose inverse on a closed retained
range when a declared gauge has a kernel. It is exact completion of the
kinetic square, not a license to add arbitrary momentum directions to a
fixed-Kelvin orbit. Identifying B,A,C with actual Euler phase variables is
the separate 0080/0082 construction.

If the reaction Hessian is P=G+L, with the local helicity contribution
L=-rho integral F.curl(F)/lambda, the correctly centered block is

    P_c=P-C* M^-1 C=L+A*(1-Pi_B)A.

Consequently a uniform lower bound L>=ell I implies P_c>=ell I,
independently of the size or nonlocality of the kinetic centering correction.
The same argument applies to the COMPLETE retained angle/reaction block,
using its full force map A and all mixed kinetic blocks. A positive local
helicity matrix therefore remains positive after exact kinetic centering.
This strengthens the small-subtraction bound in 0075: no large parcel radius
is required by this coercivity argument. It does not prove the physical
angle/centroid correspondence or positivity of background affine shear.

For the negative-helicity carriers and fixed dual responses of 0065/0080,
the lower bound has form ell(k0)>=c |k0|/|lambda|-C0, c>0. Once the entire
finite response geometry is fixed, a finite carrier makes ell>0. The
centering subtraction cannot remove this growing local term. The independent
gradient carrier uses the same reasoning for its complete kinetic Gram;
its fixed-base mixed remainders still require the operator-jet estimates
in 0065, rather than an assumption that they vanish.

## Scope and continuation

These transport and Gram propositions are established by the displayed
identities. `verify.py` corroborates exact phase normalization, noncommuting
Schur completion, a rank-deficient velocity example, and the mutation that
over-subtracting the kinetic cross destroys positivity. It is not a numerical
Euler or stationarity simulation. The next unsatisfied parent construction
is the actual hybrid material-action and field map in 0080/0082, followed
by full second-gradient closure and the individual parent claim review.

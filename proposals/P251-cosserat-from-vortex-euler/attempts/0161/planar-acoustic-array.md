# An actual smooth triangular Euler array on the acoustic time scale

This proves the planar candidate in the registered 0161 contract. It does
not replace the compensated compact-core array of 0156. All averages below
are normalized by the area of this same periodic cell. Pressure is divided
by the prescribed positive density rho. There is no numerical limit or
spectral sampling in the proof.

## 1. The field and its entire invariant-cell partition

Let b1,b2 have length lambda>0 and angle 120 degrees, b3=-b1-b2. On the
corresponding flat two-torus set

    psi=Psi [cos(b1.x)+cos(b2.x)+cos(b3.x)],   Psi>0,
    v=J grad psi,   zeta=Delta psi=-lambda^2 psi,
    J(x,y)=(-y,x).

The full stationary Euler pressure is

    p=-(|v|^2+lambda^2 psi^2)/2.

This is a smooth vortex array, not a compact-vorticity core and a separate
ambient fluid. Its velocity covariance, computed over the complete cell, is

    C_v=<v tensor v>=(3 Psi^2 lambda^2/4) I.                 (1)

Put a=b1.x, b=b2.x. The exact identity

    psi+Psi=4 Psi cos(a/2) cos(b/2) cos((a+b)/2)             (2)

exhibits a connected polygonal separatrix network at psi_s=-Psi. Its lines
are a=pi, b=pi, a+b=pi modulo 2pi. Every complementary polygon in the
universal cover is bounded. Its edges are invariant because v is tangent
to a level of psi; vertices are stationary saddles. Assign its geometric
centroid c_D and put r=x-c_D there. This defines a bounded, periodic,
piecewise smooth vector field with

    A r=v,    A=v.grad.                                    (3)

The jumps of r are across invariant edges, so their normal transport flux
is zero. Equation (3) holds in the weak transport domain, not merely on
regular contours. No positive lower bound on streamline frequency is
used. The centroid registration is equivariant under the sixfold rotation.

An important exact cell-boundary integral is

    <r_j partial_i psi>=psi_s delta_ij.                    (4)

Indeed on each polygon integrate partial_i(r_j psi), using the SAME
boundary value psi_s, and sum the resulting
delta_ij(psi_s area(D)-integral_D psi). The integral of psi over the torus
is zero. This retains the separatrix terms: setting them to zero would
destroy the nonzero quantity on the right of (4).

## 2. A bounded full horizontal Euler group in the excited sector

For mean-zero horizontal solenoidal perturbations a_h, write f=curl a_h.
The exact horizontal Euler equation at zero axial wave number is

    f_t=L_c f=-A H f,
    H=I-lambda^2 G,    G=(-Delta)^(-1).                    (5)

All inverses are those of the full periodic cell, with zero harmonic
scalar mode. Fourier wave numbers are n1 b1+n2 b2 and

    |n1 b1+n2 b2|^2/lambda^2=n1^2+n2^2-n1 n2.

The smallest nonzero value is 1, its shell has six elements, and the next
value is at least 3. Thus H is nonnegative globally and is at least 2/3
on the Fourier complement of the first shell. For the last assertion,
the integer quadratic form cannot equal 2: modulo 3 it is a square.

Use the real vector isotypic sector for rotation through pi/3. Common
horizontal velocity data transform in this sector, as do all fields
generated from them by the actual Euler operator. On the first shell,
rotation by pi acts as -1 in this sector, eliminating the cosine modes.
On the three sine modes the pi/3 rotation has characteristic polynomial
(z+1)(z^2-z+1). The z^2-z+1 space is exactly the two derivatives of psi;
the remaining sine sum has eigenvalue -1 and is a different representation.
Consequently the kernel of H IN THE EXCITED SECTOR is exactly

    t_x=lambda partial_x psi,  t_y=lambda partial_y psi.   (6)

There is no omitted stationary first-shell degree of freedom in that
sector. Rotations, rather than reflection invariance of a signed vortex
field, are used here.

Semidefinite energy alone would not control the coefficients of (6).
The actual bounded adjoint functions supplying the missing control are

    h_x=lambda r_y,   h_y=-lambda r_x,   A h_i=t_i.         (7)

Skew-adjoint transport, H t_i=0, and (5) give the two conserved rows

    d/dt <h_i,f>=<t_i,H f>=0.                             (8)

Their matrix on the translation kernel is, by (4),

    (<h_i,t_j>) = lambda^2 psi_s [[0,1],[-1,0]],           (9)

which is invertible. If f=f_perp+c_i t_i, energy conservation bounds
||f_perp|| by sqrt(3/2)<f(0),H f(0)>^(1/2), and (8)-(9) bound c by
the initial rows and ||h|| ||f_perp||. Hence

    sup_{t in R} ||exp(t L_c) f0||_2 <= C ||f0||_2        (10)

on this actual complete-cell sector. The constant is finite and is
computed from the gap, h, t and the nonzero matrix (9). It is not a
whole-group assumption imported from a compact core or an axisymmetric
whole-space theorem.

For a functional-analytic realization, -A generates a unitary transport
group on L2. The perturbation lambda^2 A G is bounded there. Thus (5)
already generates a C0 group; the energy/row identities follow first on
smooth data and then by density. The weak-domain fact (3) makes h a valid
adjoint-domain vector even at the separatrix. The preceding estimates
upgrade the group to (10).

Let Y be the solenoidal velocity space with norm ||curl a_h||_2. The
corresponding horizontal velocity group is bounded by (10). The remaining
mapping hypotheses of 0146 are also actual estimates here. With
d=grad G and F=P0[-A-(Dv)]d, curl F contains only products of smooth
zeta, grad zeta, d f and f; in particular F:L2->Y is bounded. The pressure
map Q_k:Y->L2 is bounded uniformly for small k by the periodic elliptic
estimate, and Q_k-Q_0=O(k^2). These statements keep the complete pressure
operator; there is no isolated-cell reaction or wall.

## 3. The physical acoustic history now follows for the actual array

Consider the FULL three-dimensional linear Euler system about u=(v,0),
with axial factor exp(ikz), k nonzero, and actual initial velocity

    w_h(0)=V0,    w_z(0)=0.

The real Fourier convention can use sqrt(2) cos(kz) for horizontal
components and the corresponding sine components vertically; its axial
mean-square normalization is one. Let m=<w_h>, X_t=m, X(0)=0.
Use exactly the full Hodge equations, pressure remainder pi_r, Z and
compensated current p_c of 0146, equations (1)-(4). Their bounded-group
hypothesis has now been PROVED in (10), including the response sector and
the forcing/pressure maps just identified. Therefore for every fixed
slow time T and fixed background Psi,lambda,

    sup_{0<=t<=T/|k|}
      |m(t)-cos(|k| c_b t)V0| <= C_T |k| |V0|,
    c_b=sqrt(3) Psi lambda/2,                             (11)

with the analogous bound for |k|X and the matrix-sine solution. Constants
may carry the fixed length/time units of the background; equivalently
rescale the cell to lambda=1 and time to Psi lambda^2=1, so the small
dimensionless parameter is |k|/lambda.

This is a uniform actual-Euler acoustic-window assertion. It is stronger
than a fixed-time k^2 Taylor coefficient or a compact low-frequency
resolvent limit. The proof uses no minimum transport frequency, no
contour deformation through an unexamined spectrum, and no estimate of
the form exp(C/|k|). The fast groups are controlled BEFORE Gronwall on
the slow interval. The current p_c and the point mean m remain distinct:

    p_c=m+i k <r Z>,    p_c-m=O_T(k),
    (|k|X)_tau=p_c+O_T(k),
    (p_c)_tau=-C_v (|k|X)+O_T(k).                         (12)

These are computed full-fluid current corrections. Rapid derivatives
of m itself are not inferred small from its uniform approximation.

The same argument permits a second, well-prepared displacement phase.
Given X0, take an initial material displacement eta0=(X0,0) exp(ikz)
and its EXACT Kelvin-prepared initial velocity

    w0=P_k(eta0 cross omega0)
       =T X0-i k P_k[(v.X0)e_z],   T X0=-(X0.grad)v.      (13)

Writing D_k=(-Delta_h+k^2)^(-1), its horizontal return is
k^2 grad D_k(v.X0) and its vertical component is
-i k(v.X0)+i k^3 D_k(v.X0). Both are retained. Relative to X(0)=X0,
the normal-form remainders are y0=O(k^2 X0), Z0=O(k^3 X0).
Thus (11)-(12) extend to |k|X0 and V0 bounded, with leading

    m(t)=cos(|k|c_b t)V0-|k| c_b sin(|k|c_b t)X0
          +O_T(|k|)(|V0|+|k X0|).                       (14)

No assertion that the two retained phase columns form an invariant
finite-dimensional Euler subspace is needed: they are actual solutions
and (11)-(14) are observation estimates on the specified time interval.

## 4. Full material action and its mass normalization

Use the exact Euler Jacobi action and cotangent momentum

    L2=rho/2 <|D_t eta|^2>-1/2 <eta.Hess(p_phys).eta>,
    pi=rho D_t eta,  p_phys=rho p.

For the common-velocity phase eta(0)=0, pi(0)=rho(V0,0).
For the prepared displacement phase (13), D_t eta(0)=
-i k P_k[(v.X0)e_z], because A eta0=0 and the Lin reconstruction
subtracts T X0. Therefore the actual phase symplectic pairing is
exactly rho X0.V0, with all displacement-displacement pairings zero:
the return momentum has zero horizontal cell average. This determines
the physical mass rho from the SAME action, not by appending a rotor or
an ambient kinetic term.

The common-velocity phase has its prescribed initial circulation-momentum
variation rho V0 exp(ikz); at nonzero k it is not a zero-circulation
perturbation. The displacement phase (13) is fixed-Kelvin prepared. Their
combination is the full material-Euler phase family with these explicit
initial circulation data, not a claim that both independent phases belong
to one fixed coadjoint leaf. Euler transports each phase's specified
circulation data exactly.

The Hamiltonian at this initial phase has no mixed X0,V0 term. The
spatial mean Hessian of the periodic pressure is zero, and A eta0=0.
Since v.X0 lies entirely in the first Fourier shell,

    H_phase=rho/2 |V0|^2
      +rho/2 k^2 [lambda^2/(lambda^2+k^2)] X0.C_v X0.     (15)

Indeed ||P_k[(v.X0)e_z]||^2=
lambda^2/(lambda^2+k^2)<(v.X0)^2>. This positive complete-fluid phase
energy and the exact symplectic pairing are independently calculable
from the actual initial fields. They retain the axial divergence return.

At later times the exact solution-column pullback retains the usual
moving-frame action and its connection. The observation chart
(|k|X,p_c) is uniformly invertible on this four-real-dimensional phase
family for small k: (12) and its two preparations make it C1-close on
slow time to the invertible oscillator fundamental matrix. Pulling the
exact symplectic form and connection through this chart consequently
gives the leading action with mass rho and stiffness rho k^2 C_v,
with O(k/lambda) corrections in the normalized slow chart. This is not
an assertion that the connection vanishes identically, nor a replacement
of p_c by the point mean before differentiating. The canonical
moving-frame identities in 0115 remain the exact finite-k rule.

Actual nonlinear Euler histories on each finite interval T/|k| follow by
taking perturbation amplitude sufficiently small LAST. This ordering does
not claim a uniform-in-k nonlinear amplitude or infinite-time nonlinear
stability.

## 5. Scope and the next same-field transfer

The planar actual-array acoustic candidate is established at its stated
scope. The full field (v,0) is smooth, periodic and has finite cell energy;
the result concerns its entire fluid, including the separatrix cells.
No compact-vorticity ambient decomposition is needed for this candidate.

The constant-curl lift W=zeta/lambda is not covered by (11). Its exact
horizontal mean includes -i k <W w_h+v w_z>, and its vertical mean and
pressure acquire corresponding terms. An O(k) microscopic change can
alter the leading slow equation. That registered next route must compute
these terms and its action before any same-field EPS/optical join.
The original compensated-array route, generic-wavevector isotropy and
the full coupled continuum remain active rather than being inherited
from this planar theorem.
